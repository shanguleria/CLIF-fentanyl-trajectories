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
from utils.paths import (  # noqa: E402
    MANIFEST, clear_owned_outputs, config_digests, require_manifest, site_dirs,
)

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


def test_a_new_directory_under_output_is_ignored_by_default():
    """output/ denies by default and allows back explicitly, so a directory added
    later is PHI-safe without anyone remembering to add a rule. A bare `output/`
    would not work: git does not descend into an excluded directory, which kills
    every negation below it and takes final_no_phi with it."""
    for path in ("output/some_new_dir/x.csv", "output/some_new_dir/x.parquet",
                 "output/scratch/notes.txt"):
        rc = subprocess.run(["git", "check-ignore", "-q", path], cwd=REPO).returncode
        assert rc == 0, f"{path} is COMMITTABLE; output/ must deny by default"


def test_a_stray_parquet_in_the_shareable_set_is_still_blocked():
    rc = subprocess.run(["git", "check-ignore", "-q",
                         "output/final_no_phi/leak.parquet"], cwd=REPO).returncode
    assert rc == 0, "final_no_phi allows .csv and .json back, not patient-level formats"


def test_shareable_outputs_are_not_ignored():
    for path in ("output/final_no_phi/phase0_strobe.csv",
                 "output/final_no_phi/phase0_provenance.json",
                 "output/final_no_phi/validation/scaling_experiments.csv"):
        rc = subprocess.run(["git", "check-ignore", "-q", path], cwd=REPO).returncode
        assert rc != 0, f"{path} is ignored, but final_no_phi is the shareable set"


def test_absent_manifest_stops_a_downstream_phase():
    """A stale parquet on disk looks perfectly valid. The manifest is written last,
    so its absence means the run did not complete."""
    import json, tempfile

    with tempfile.TemporaryDirectory() as tmp:
        d = {"out_final": Path(tmp), "out_phi": Path(tmp), "logs": Path(tmp)}
        try:
            require_manifest(d, {}, REPO)
        except SystemExit as e:
            assert "never completed" in str(e)
        else:
            raise AssertionError("an absent manifest must stop the phase")


def test_a_changed_config_invalidates_the_outputs():
    import json, tempfile

    with tempfile.TemporaryDirectory() as tmp:
        d = {"out_final": Path(tmp), "out_phi": Path(tmp), "logs": Path(tmp)}
        stale = dict(config_digests(REPO))
        stale["covariates.json"] = "0" * 16
        (Path(tmp) / MANIFEST).write_text(json.dumps({"config_digests": stale}))
        try:
            require_manifest(d, {}, REPO)
        except SystemExit as e:
            assert "config has changed" in str(e) and "covariates.json" in str(e)
        else:
            raise AssertionError("a changed config must invalidate the outputs")


def test_a_matching_manifest_passes():
    import json, tempfile

    with tempfile.TemporaryDirectory() as tmp:
        d = {"out_final": Path(tmp), "out_phi": Path(tmp), "logs": Path(tmp)}
        (Path(tmp) / MANIFEST).write_text(
            json.dumps({"config_digests": config_digests(REPO)}))
        require_manifest(d, {}, REPO)


def test_clearing_removes_a_stale_output():
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        d = {"out_phi": Path(tmp), "out_final": Path(tmp)}
        (Path(tmp) / "trajectory_long.parquet").write_text("stale")
        n = clear_owned_outputs(d, {"out_phi": ["trajectory_long.parquet"]})
        assert n == 1
        assert not (Path(tmp) / "trajectory_long.parquet").exists(), (
            "a crash after clearing must leave nothing, not a stale file"
        )


def test_the_review_csvs_are_inside_the_phi_boundary():
    """They are one row per patient-window, so they are PHI regardless of format."""
    for name in ("trajectory_long.csv", "time_to_event.csv"):
        path = f"output/intermediate_phi/{name}"
        rc = subprocess.run(["git", "check-ignore", "-q", path], cwd=REPO).returncode
        assert rc == 0, f"{path} is COMMITTABLE -- a per-patient CSV is still PHI"


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
