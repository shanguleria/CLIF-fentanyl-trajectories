"""Where the pipeline is allowed to write, and that the PHI rule is by directory.

output/intermediate_phi/ is the single location for patient-level artifacts, per
CLIF convention. data/ and output/intermediate/ are retired; nothing may write to
them, and no source file may still name them.

Run standalone:  .venv/bin/python tests/test_paths.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "code"))
from utils.paths import site_dirs  # noqa: E402

SOURCE = list((REPO / "code").rglob("*.py")) + list((REPO / "code").rglob("*.R")) \
    + list((REPO / "validation").rglob("*.py")) + list((REPO / "validation").rglob("*.R"))
SOURCE = [f for f in SOURCE if "__pycache__" not in str(f)]


def test_site_dirs_returns_only_the_sanctioned_locations():
    d = site_dirs(REPO)
    assert set(d) == {"out_phi", "out_final", "logs"}, (
        f"site_dirs returns {sorted(d)}; the sanctioned set is out_phi, out_final, logs"
    )
    assert d["out_phi"] == REPO / "output" / "intermediate_phi"
    assert d["out_final"] == REPO / "output" / "final_no_phi"


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


def test_shareable_outputs_are_not_ignored():
    for path in ("output/final_no_phi/phase0_strobe.csv",
                 "output/final_no_phi/validation/scaling_experiments.csv"):
        rc = subprocess.run(["git", "check-ignore", "-q", path], cwd=REPO).returncode
        assert rc != 0, f"{path} is ignored, but final_no_phi is the shareable set"


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
