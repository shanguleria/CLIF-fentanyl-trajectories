"""Output directories and the provenance block stamped on shareable outputs."""
from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

MANIFEST = "phase0_manifest.json"

# Written into intermediate_phi/ at runtime so it survives `git clean -fdx`.
PHI_LABEL = """# intermediate_phi

Every patient-level artifact this pipeline produces: the Phase 0 analytic tables
and the model fits that later phases read. One row per encounter block, or per
encounter-block-window.

**This directory never leaves the site.** It is not part of any export bundle
and must not be committed, copied to a shared drive, or read into an analysis
transcript. Only `output/final_no_phi/` is shareable.
"""


def site_dirs(repo: Path) -> dict[str, Path]:
    """Create and return the directories this pipeline may write to."""
    d = {
        "out_phi": repo / "output" / "intermediate_phi",
        "out_final": repo / "output" / "final_no_phi",
        "diagnostics": repo / "output" / "final_no_phi" / "diagnostics",
        "logs": repo / "logs",
    }
    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)

    label = d["out_phi"] / "README.md"
    if not label.exists():
        label.write_text(PHI_LABEL)

    return d


def provenance(config: dict) -> dict:
    """The block stamped onto every shareable output."""
    try:
        sha = subprocess.check_output(
            ["git", "describe", "--always", "--dirty"], text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        sha = "unknown"
    return {
        "site_name": config["site_name"],
        "clif_version": config["clif_version"],          # CLIF SPEC version
        "dataset_version": config.get("dataset_version", ""),  # conversion/ETL release
        "code_version": sha,
        "generated": datetime.now(ZoneInfo(config["timezone"])).isoformat(),
    }


def config_digests(repo: Path) -> dict[str, str]:
    """SHA-256 of every file that governs what the outputs contain."""
    out = {}
    for name in ("config.json", "covariates.json", "outlier_config.json"):
        f = repo / "config" / name
        if f.exists():
            out[name] = hashlib.sha256(f.read_bytes()).hexdigest()[:16]
    return out


def clear_owned_outputs(dirs: dict[str, Path], owned: dict[str, list[str]]) -> int:
    """Delete this script's own outputs before it starts.

    A crash then leaves nothing rather than a stale file that still looks valid.
    The dangerous case is a crash between two writes, which leaves a MISMATCHED
    pair from two different code versions.
    """
    n = 0
    for key, names in owned.items():
        for name in names:
            f = dirs[key] / name
            if f.exists():
                f.unlink()
                n += 1
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
            f"{MANIFEST} is absent, so Phase 0 either never completed or was "
            f"cleared. Re-run code/01_build_cohort.py; do not read the parquet "
            f"files that may still be on disk."
        )
    m = json.loads(f.read_text())
    now = config_digests(repo)
    drift = {k: (v, now.get(k)) for k, v in m.get("config_digests", {}).items()
             if now.get(k) != v}
    if drift:
        raise SystemExit(
            f"config has changed since Phase 0 ran: "
            + ", ".join(f"{k} {a} -> {b}" for k, (a, b) in drift.items())
            + ". Re-run code/01_build_cohort.py rather than analysing stale tables."
        )
    return m
