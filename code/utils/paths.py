"""Output directories and the provenance block stamped on shareable outputs."""
from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

MANIFEST = "manifest.json"

# Written into intermediate_phi/ at runtime so it survives `git clean -fdx`.
PHI_LABEL = """# intermediate_phi

Every patient-level artifact this pipeline produces: the cohort build's analytic
tables and the model fits later scripts read. One row per encounter block, or
per encounter-block-window.

SUBDIVIDED BY THE SCRIPT THAT WROTE THE FILE, the same way output/final_no_phi/
is -- `01_cohort/`, `05_titration/`, `07_outcome_model/`, `tabled/`. A reader
therefore says where a table came from, and the handoffs between scripts are
legible from the path. `_cache/` (the waterfall cache) and this label stay at
the root: neither is a script's output.

**This directory never leaves the site.** It is not part of any export bundle
and must not be committed, copied to a shared drive, or read into an analysis
transcript. Only `output/final_no_phi/` is shareable.
"""


def site_dirs(repo: Path) -> dict[str, Path]:
    """Create and return the directories this pipeline may write to.

    Per-script subfolders come from phase_dir() (shareable) and phi_dir()
    (patient-level), not from here.
    """
    d = {
        "out_phi": repo / "output" / "intermediate_phi",
        "out_final": repo / "output" / "final_no_phi",
        "logs": repo / "logs",
    }
    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)

    label = d["out_phi"] / "README.md"
    if not label.exists():
        label.write_text(PHI_LABEL)

    return d


def phase_dir(dirs: dict[str, Path], name: str) -> Path:
    """A per-script subfolder of the shareable output tree.

    The folder names the script that wrote the file. The old phaseN_ file
    prefixes were removed on 2026-09-23 -- the folder already says which script
    produced a file, so the prefix said it twice. `phase` in this function's
    name and in dirs["phase"] now means "this script's output folder" and
    nothing more.
    """
    d = dirs["out_final"] / name
    d.mkdir(parents=True, exist_ok=True)
    return d


def phi_dir(dirs: dict[str, Path], name: str) -> Path:
    """A per-script subfolder of the PATIENT-LEVEL tree. Mirrors phase_dir().

    Same rule as the shareable side: the folder names the script that WROTE the
    file, never the one that reads it. So a reader says where a table came from
    --- phi_dir(dirs, "01_cohort") / "trajectory_long.parquet" --- and the
    handoffs between scripts are legible from the path instead of from a flat
    directory of thirty files.

    `_cache/` and the README label stay at the root: the first is not a script's
    output and the second is about the whole tree.
    """
    d = dirs["out_phi"] / name
    d.mkdir(parents=True, exist_ok=True)
    return d


# The code the runner executes. code/tabled/ is excluded because the runner never
# calls it; tests/ because it produces no output.
CODE_GLOBS = ("code/*.py", "code/*.R", "code/utils/*.py", "code/utils/*.R")
# Protocol. MUST be byte-identical across sites or a pooled result is not
# meaningful. config.json is deliberately NOT here: it carries site_name and
# data_directory, so its digest is expected to differ between sites.
PROTOCOL_FILES = ("protocol.json", "outlier_config.json")
SITE_CONFIG = "config.json"


def _repo() -> Path:
    """This file is code/utils/paths.py, so the repo is two levels up."""
    return Path(__file__).resolve().parents[2]


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()[:16]


def code_digests(repo: Path | None = None) -> tuple[str, dict[str, str]]:
    """Content fingerprint of the code that produced a run: (overall, per file).

    Content-based rather than git-based: `git describe` returns "unknown" for a
    ZIP download and a bare "-dirty" for an uncommitted edit, while a digest is
    comparable across sites by equality. The overall digest hashes the
    (path, digest) PAIRS, so a rename moves it too.
    """
    repo = repo or _repo()
    files = sorted({p for g in CODE_GLOBS for p in repo.glob(g)},
                   key=lambda p: p.relative_to(repo).as_posix())
    per = {p.relative_to(repo).as_posix(): _sha(p.read_bytes()) for p in files}
    joined = "\n".join(f"{k} {v}" for k, v in per.items())
    return _sha(joined.encode()), per


def definition_versions(repo: Path | None = None) -> dict[str, str]:
    """The protocol version each config declares.

    One field per file, not one overall: they are bumped independently, so a
    single `definition_version` would be ambiguous about which it meant.
    """
    repo = repo or _repo()
    out = {}
    for name in PROTOCOL_FILES:
        f = repo / "config" / name
        if f.exists():
            out[name] = json.loads(f.read_text()).get("definition_version", "")
    return out


def provenance(config: dict) -> dict:
    """The block stamped onto every shareable output.

    `code_digest` is the authoritative answer to "did two sites run the same
    code?". `code_version` is the human-readable git label and may be "unknown".
    """
    repo = _repo()
    try:
        sha = subprocess.check_output(
            ["git", "describe", "--always", "--dirty"], text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        sha = "unknown"
    digest, per_file = code_digests(repo)
    cfg = config_digests(repo)
    return {
        "site_name": config["site_name"],
        "clif_version": config["clif_version"],          # CLIF SPEC version
        "definition_versions": definition_versions(repo),
        "code_version": sha,
        "code_digest": digest,
        "protocol_digests": {k: v for k, v in cfg.items() if k in PROTOCOL_FILES},
        "site_config_digest": cfg.get(SITE_CONFIG, ""),
        "file_digests": per_file,
        "generated": datetime.now(ZoneInfo(config["timezone"])).isoformat(),
    }


def config_digests(repo: Path) -> dict[str, str]:
    """SHA-256 of every file that governs what the outputs contain.

    require_manifest() uses these for WITHIN-site staleness; provenance() splits
    the same digests into protocol-vs-site for the CROSS-site check.
    """
    out = {}
    for name in (SITE_CONFIG,) + PROTOCOL_FILES:
        f = repo / "config" / name
        if f.exists():
            out[name] = _sha(f.read_bytes())
    return out


def clear_owned_outputs(dirs: dict[str, Path], owned: dict[str, list[str]],
                        retired: list[str] | None = None) -> int:
    """Delete this script's own outputs before it starts.

    A crash then leaves nothing rather than a MISMATCHED pair of files from two
    different code versions, which still looks valid.
    """
    n = 0
    for key, names in owned.items():
        for name in names:
            f = dirs[key] / name
            if f.exists():
                f.unlink()
                n += 1
    # Relocating an output leaves a stale twin at the old path, which the owned
    # list no longer names. Retired paths are cleared explicitly.
    for rel in retired or []:
        f = dirs["out_final"].parent.parent / rel
        if f.exists():
            f.unlink()
            n += 1
    # A renumbered SCRIPT leaves the whole directory behind, so remove a retired
    # directory once it is empty -- never one that still holds anything, and
    # never one this run writes to: a script may retire files in its OWN output
    # directory, leaving it legitimately empty.
    live = {Path(v).resolve() for v in dirs.values()}
    for d in {(dirs["out_final"].parent.parent / rel).parent for rel in retired or []}:
        if d.resolve() in live:
            continue
        if d.is_dir() and not any(d.iterdir()):
            d.rmdir()
    return n


def write_manifest(dirs: dict[str, Path], config: dict, repo: Path,
                   outputs: dict[str, int]) -> dict:
    """Written LAST. Its presence is what marks the outputs complete and current."""
    m = provenance(config)
    m["config_digests"] = config_digests(repo)
    m["outputs"] = outputs
    (dirs["out_final"] / MANIFEST).write_text(json.dumps(m, indent=2))
    return m


def require_manifest(dirs: dict[str, Path], config: dict, repo: Path) -> dict:
    """Fail loudly if the inputs on disk were not produced by this code and config."""
    f = dirs["out_final"] / MANIFEST
    if not f.exists():
        raise SystemExit(
            f"{MANIFEST} is absent, so the cohort build either never completed or "
            f"was cleared. Re-run code/01_build_cohort.py; do not read the parquet "
            f"files that may still be on disk."
        )
    m = json.loads(f.read_text())
    # REQUIRED, not .get(): reading an absent block as {} makes this guard pass
    # vacuously, a manifest without digests being uncheckable.
    stamped = m.get("config_digests")
    if not stamped:
        raise SystemExit(
            f"{MANIFEST} carries no config_digests, so the tables on disk cannot "
            f"be checked against the current config. Re-run code/01_build_cohort.py."
        )
    now = config_digests(repo)
    drift = {k: (v, now.get(k)) for k, v in stamped.items() if now.get(k) != v}
    if drift:
        raise SystemExit(
            f"config has changed since the cohort build ran: "
            + ", ".join(f"{k} {a} -> {b}" for k, (a, b) in drift.items())
            + ". Re-run code/01_build_cohort.py rather than analysing stale tables."
        )
    # Code drift is RECORDED, never refused: a site may have to edit _to_mcg_hr
    # to run at all. The digest in each provenance.json is what the coordinating
    # centre compares.
    was = m.get("code_digest")
    if was and was != code_digests(repo)[0]:
        print(f"  NOTE code changed since the cohort build ran ({was} -> "
              f"{code_digests(repo)[0]}); these outputs mix two code versions")
    return m
