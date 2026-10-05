"""Where the pipeline is allowed to write, and that the PHI rule is by directory.

output/intermediate_phi/ is the single location for patient-level artifacts, per
CLIF convention. data/ and output/intermediate/ are retired; nothing may write to
them, and no source file may still name them.

Run standalone:  .venv/bin/python tests/test_paths.py
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "code"))
from utils.paths import (  # noqa: E402
    MANIFEST, clear_owned_outputs, config_digests, require_manifest, site_dirs,
)

SOURCE = list((REPO / "code").rglob("*.py")) + list((REPO / "code").rglob("*.R")) \
    + list((REPO / "validation").rglob("*.py")) + list((REPO / "validation").rglob("*.R"))
SOURCE = [f for f in SOURCE if "__pycache__" not in str(f)]


def test_site_dirs_returns_only_the_sanctioned_locations():
    d = site_dirs(REPO)
    assert set(d) == {"out_phi", "out_final", "logs"}, (
        f"site_dirs returns {sorted(d)}; the sanctioned set is out_phi, "
        f"out_final, logs. Per-script subfolders come from phase_dir()."
    )
    assert d["out_phi"] == REPO / "output" / "intermediate_phi"
    assert d["out_final"] == REPO / "output" / "final_no_phi"


def test_paths_r_and_paths_py_expose_the_same_surface():
    """paths.R's own header says the two must agree, and nothing checked it.

    They did not: R's site_dirs() returned three keys to Python's four, and R had
    no clear_owned_outputs() at all. This matters because Phase 0 (Python) writes
    what Phases 1-6 (R) read -- a function that exists on one side only is a
    guarantee that silently does not apply to half the pipeline.
    """
    r = (REPO / "code" / "utils" / "paths.R").read_text()
    r_fns = set(re.findall(r"^(\w+)\s*<-\s*function", r, re.M))
    py = (REPO / "code" / "utils" / "paths.py").read_text()
    py_fns = set(re.findall(r"^def (\w+)", py, re.M))
    # Private helpers are excluded on both sides: Python marks them with a
    # leading _ and R with a leading . (which the \w+ patterns above already
    # miss), so the two conventions differ and neither is part of the surface.
    py_fns = {f for f in py_fns if not f.startswith("_")}
    # write_manifest is Python-only by design: only Phase 0 writes the manifest.
    # Everything else must exist on both sides.
    py_only_by_design = {"write_manifest"}
    missing_in_r = py_fns - r_fns - py_only_by_design
    assert not missing_in_r, (
        f"paths.py has {sorted(missing_in_r)} and paths.R does not; either port "
        f"them or record them as Python-only in this test"
    )
    missing_in_py = r_fns - py_fns
    assert not missing_in_py, f"paths.R has {sorted(missing_in_py)} and paths.py does not"


def test_the_two_site_dirs_return_the_same_key_set():
    """A directory the Python half creates and the R half cannot name is a place
    Phase 0 writes and Phase 1 cannot read."""
    r = (REPO / "code" / "utils" / "paths.R").read_text()
    body = r[r.index("site_dirs <- function"):]
    body = body[:body.index("\n}")]
    r_keys = set(re.findall(r"^\s{4}(\w+)\s*=\s*file\.path", body, re.M))
    assert r_keys == set(site_dirs(REPO)), (
        f"paths.R site_dirs returns {sorted(r_keys)}; paths.py returns "
        f"{sorted(site_dirs(REPO))}"
    )


def test_the_two_languages_stamp_the_SAME_code_digest():
    """code_digest is the cross-site "did we run the same code?" check, so a
    Python/R disagreement makes every R-written provenance block incomparable to
    every Python-written one -- and the pipeline writes both in a single run.

    They disagreed on first implementation (2026-10-02), and the PER-FILE digests
    matched, which is what made it invisible: R's default sort() uses LOCALE
    collation and put code/utils/paths.py before paths.R, while Python's sorted()
    is byte order and gives paths.R first ('R' 0x52 < 'p' 0x70). Only the overall
    digest moved. Fixed with order(method = "radix").
    """
    from utils.paths import code_digests

    py_digest, py_files = code_digests(REPO)
    r = subprocess.run(
        ["Rscript", "-e",
         'suppressMessages(library(here)); source("code/utils/paths.R"); '
         'cd <- code_digests(); cat(cd$overall, length(cd$per_file))'],
        cwd=REPO, capture_output=True, text=True)
    assert r.returncode == 0, (
        "Rscript failed, so this check cannot run -- and it is the only thing "
        f"holding the two digests together:\n{r.stderr[-500:]}"
    )
    out = r.stdout.strip().split()
    assert out[0] == py_digest, (
        f"code_digest disagrees across languages: Python {py_digest}, R {out[0]}. "
        f"Compare CODE_GLOBS and the sort order in paths.py and paths.R."
    )
    assert int(out[1]) == len(py_files), (
        f"the two halves digest a different NUMBER of files: Python "
        f"{len(py_files)}, R {out[1]} -- CODE_GLOBS has drifted"
    )


def test_the_r_staleness_guard_does_not_cry_wolf_on_an_unchanged_tree():
    """require_manifest() in paths.R compares the manifest's code_digest to the
    tree's. .sha_string() returns a NAMED character, so identical() compared the
    names attribute too and was FALSE for equal digests: every R script printed
    `NOTE code changed since Phase 0 ran (424e4460bc02333e -> 424e4460bc02333e)`
    -- the same hash twice -- on every run from 2026-10-02 to 2026-10-05.

    A guard that can never say "no drift" cannot report drift either, so the
    bug disabled the check rather than just adding noise. The Python twin at
    paths.py:243 compares plain strings and was always right, which is how the
    two functions paths.R calls mirrors diverged unnoticed.

    This asserts the COMPARISON, not the digest value -- the existing
    cross-language test already covers the value and passed throughout.
    """
    r = subprocess.run(
        ["Rscript", "-e",
         'suppressMessages(library(here)); source("code/utils/paths.R"); '
         'cd <- code_digests(); '
         'cat(identical(as.character(cd$overall), unname(cd$overall)))'],
        cwd=REPO, capture_output=True, text=True)
    assert r.returncode == 0, (
        f"Rscript failed, so this check cannot run:\n{r.stderr[-500:]}")
    assert r.stdout.strip() == "TRUE", (
        "code_digests()$overall still carries attributes that make identical() "
        "false against its own string value, so require_manifest() reports code "
        "drift on an unchanged tree. Keep the unname() in paths.R.")

    src = (REPO / "code" / "utils" / "paths.R").read_text()
    assert "unname(code_digests(root)$overall)" in src, (
        "require_manifest() must unname() the digest before identical(), or the "
        "staleness NOTE fires on every run and stops meaning anything")


def test_the_code_digest_moves_when_any_covered_file_changes():
    """A fingerprint that does not move is worse than none: it would certify two
    different code bases as identical. Checked by perturbing a real covered file
    in a temp copy rather than by reasoning about the hash."""
    import shutil
    import tempfile

    from utils.paths import code_digests

    base, files = code_digests(REPO)
    assert files, "code_digests covered no files at all"
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / "repo"
        for rel in files:
            (root / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(REPO / rel, root / rel)
        assert code_digests(root)[0] == base, (
            "a faithful copy of the covered files must digest identically"
        )
        victim = root / sorted(files)[0]
        victim.write_text(victim.read_text() + "\n# one added comment\n")
        assert code_digests(root)[0] != base, (
            f"editing {sorted(files)[0]} did not move the code digest"
        )


def test_no_source_file_writes_to_a_retired_directory():
    retired = ("data/intermediate_phi", "output/intermediate/", 'dirs["data_phi"]',
               "dirs$data_phi")
    offenders = []
    for f in SOURCE:
        text = f.read_text()
        for token in retired:
            if token in text:
                offenders.append(f"{f.relative_to(REPO)}: {token!r}")
    assert not offenders, "retired paths still referenced:\n  " + "\n  ".join(offenders)


def test_phi_rule_is_by_directory_not_by_extension():
    """An extension rule misses the file type nobody thought of -- which is how a
    per-patient .rds became committable in this repo once already."""
    for name in ("x.parquet", "x.csv", "x.rds", "x.png", "x.txt", "x.json", "x"):
        path = f"output/intermediate_phi/{name}"
        rc = subprocess.run(["git", "check-ignore", "-q", path], cwd=REPO).returncode
        assert rc == 0, f"{path} is COMMITTABLE -- the PHI directory rule is not covering it"


def test_nothing_under_output_is_committable():
    """This repo ships to sites and every site generates its own outputs. Nothing
    under output/ is tracked: PHI artifacts never leave the site, and the PHI-free
    set reaches the coordinating centre by upload, not by git."""
    for path in ("output/intermediate_phi/trajectory_long.parquet",
                 "output/intermediate_phi/x.rds",
                 "output/final_no_phi/01_cohort/strobe.csv",
                 "output/final_no_phi/01_cohort/strobe.png",
                 "output/final_no_phi/01_cohort/diagnostics/missingness.csv",
                 "output/some_new_dir/anything.txt"):
        rc = subprocess.run(["git", "check-ignore", "-q", path], cwd=REPO).returncode
        assert rc == 0, f"{path} is COMMITTABLE; nothing under output/ may be tracked"


def test_clearing_removes_a_retired_path():
    """Moving an output leaves a stale twin at the old location that the owned
    list no longer names -- three survived a relocation before this existed."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "output" / "final_no_phi").mkdir(parents=True)
        stale = root / "output" / "final_no_phi" / "old_report.csv"
        stale.write_text("stale")
        d = {"out_phi": root, "out_final": root / "output" / "final_no_phi",
             "diagnostics": root, "logs": root}
        n = clear_owned_outputs(d, {}, retired=["output/final_no_phi/old_report.csv"])
        assert n == 1 and not stale.exists(), "a retired path must be cleared"


def test_no_output_file_is_tracked_in_the_index():
    """.gitignore does not affect files already added, so an output committed by
    mistake stays committed. Six were, once."""
    r = subprocess.run(["git", "ls-files", "output/"], cwd=REPO,
                       capture_output=True, text=True)
    tracked = [x for x in r.stdout.split("\n") if x.strip()]
    assert not tracked, (
        "output files are tracked in git despite the ignore rule:\n  "
        + "\n  ".join(tracked)
        + "\nUse: git rm --cached <path>"
    )


def test_the_phi_directory_carries_its_warning_label():
    d = site_dirs(REPO)
    label = d["out_phi"] / "README.md"
    assert label.exists(), "intermediate_phi must carry its PHI warning label"
    assert "never leaves the site" in label.read_text()


if __name__ == "__main__":
    import traceback

    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
        except Exception:
            failed += 1
            print(f"  FAIL  {t.__name__}")
            traceback.print_exc(limit=2)
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
