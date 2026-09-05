"""
Tests for code/utils/outliers.py and config/outlier_config.json.

The behaviours worth testing here are the LOUD FAILURES, not the arithmetic.
Clipping a number to a range is not where these go wrong -- they go wrong by
silently doing nothing, which is indistinguishable from doing the job.

Run standalone:  .venv/bin/python tests/test_outliers.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "code"))
from utils.outliers import (  # noqa: E402
    OutlierConfigError, apply_long, apply_med_converted, apply_med_raw, apply_wide,
    fentanyl_sanity_ceiling, load_config, nee_sanity_ceiling,
)

CFG = load_config()
COV = json.loads((REPO / "config" / "covariates.json").read_text())


# ---------------------------------------------------------------- loud failures
def test_missing_config_file_raises_rather_than_skipping():
    try:
        load_config(REPO / "config" / "does_not_exist.json")
    except OutlierConfigError as e:
        assert "Refusing to continue" in str(e)
    else:
        raise AssertionError("a missing outlier config must be fatal, not a no-op")


def test_unknown_table_raises():
    df = pd.DataFrame({"cat": ["x"], "val": [1.0]})
    try:
        apply_long(df, "not_a_table", "cat", "val", config=CFG)
    except OutlierConfigError as e:
        assert "unbounded frame must never look bounded" in str(e)
    else:
        raise AssertionError("an unknown table must raise")


def test_a_category_with_no_bound_is_reported_not_silently_passed():
    """A category with no entry looks exactly like one checked and found clean."""
    df = pd.DataFrame({"lab_category": ["lactate"] * 3 + ["ferritin"] * 2,
                       "lab_value_numeric": [1.0, 2.0, 999.0, 5.0, 6.0]})
    out, rep = apply_long(df, "labs", "lab_category", "lab_value_numeric", config=CFG)
    assert "ferritin" in rep.unbounded, "an unbounded category must be named"
    assert "NO BOUND for ferritin" in str(rep)
    assert out.loc[3:, "lab_value_numeric"].notna().all(), "unbounded values pass through untouched"


def test_a_med_unit_with_no_bound_is_reported():
    df = pd.DataFrame({"med_category": ["fentanyl"] * 2,
                       "med_dose": [50.0, 60.0],
                       "med_dose_unit": ["mcg/kg/day"] * 2})
    _, rep = apply_med_raw(df, "med_category", "med_dose", "med_dose_unit", config=CFG)
    assert "fentanyl [mcg/kg/day]" in rep.unbounded


# ---------------------------------------------------------------- correctness
def test_long_bounds_are_applied_per_category():
    df = pd.DataFrame({
        "lab_category": ["lactate", "lactate", "bicarbonate", "bicarbonate"],
        "lab_value_numeric": [5.0, 999.0, 24.0, 999.0],
    })
    out, rep = apply_long(df, "labs", "lab_category", "lab_value_numeric", config=CFG)
    assert rep.n_nulled == 2
    assert out["lab_value_numeric"].tolist()[0] == 5.0
    assert pd.isna(out["lab_value_numeric"].iloc[1])
    assert pd.isna(out["lab_value_numeric"].iloc[3])


def test_wide_bounds_apply_only_to_known_columns():
    df = pd.DataFrame({"lactate": [5.0, 999.0], "some_other_col": [1e9, 2e9]})
    out, rep = apply_wide(df, "labs", config=CFG)
    assert rep.n_nulled == 1
    assert out["some_other_col"].tolist() == [1e9, 2e9], "unknown columns are untouched"


def test_raw_and_converted_med_bounds_are_different_functions():
    """Angiotensin charted in ng/kg/min: correct at 20, and 20x the converted
    ceiling. The raw layer must keep it; the converted layer would null it."""
    raw = pd.DataFrame({"med_category": ["angiotensin"], "med_dose": [20.0],
                        "med_dose_unit": ["ng/kg/min"]})
    out_raw, rep_raw = apply_med_raw(raw, "med_category", "med_dose", "med_dose_unit", config=CFG)
    assert rep_raw.n_nulled == 0, "a correct ng/kg/min dose must survive the raw layer"
    assert out_raw["med_dose"].iloc[0] == 20.0

    conv = pd.DataFrame({"med_category": ["angiotensin"], "dose": [20.0]})
    _, rep_conv = apply_med_converted(conv, "med_category", "dose", config=CFG)
    assert rep_conv.n_nulled == 1, (
        "the same number in the CONVERTED unit is 20x the ceiling -- which is why "
        "the two layers must not be interchanged"
    )


def test_angiotensin_is_bounded_here_since_clifpy_lacks_it():
    raw = {k: v for k, v in CFG["med_dose_raw"].items() if not k.startswith("_")}
    assert "angiotensin" in raw
    assert "ng/kg/min" in raw["angiotensin"], (
        "ng/kg/min is the unit some sites chart angiotensin in; without it the drug "
        "passes the raw layer unbounded"
    )


# ---------------------------------------------------------------- derived ceilings
def test_fentanyl_sanity_ceiling_is_derived_not_written_down():
    """Derived so it cannot drift when a bound changes."""
    expected = (CFG["med_dose_raw"]["fentanyl"]["mcg/hr"][1]
                / CFG["analysis_unit_bounds"]["vitals"]["weight_kg"][0])
    assert fentanyl_sanity_ceiling(CFG) == expected == 25.0


def test_nee_sanity_ceiling_matches_the_reference_implementation():
    """17.45 is the arithmetic ceiling CRRT-dose-lmtp reports for the same six
    drugs. Reproducing it exactly confirms our coefficient table and converted
    bounds both match theirs -- an independent check on two ported tables."""
    coef = COV["time_varying"]["nee"]["coefficients"]
    assert abs(nee_sanity_ceiling(coef, CFG) - 17.45) < 1e-9


# ---------------------------------------------------------------- config integrity
def test_every_nee_drug_has_both_a_raw_and_a_converted_bound():
    coef = COV["time_varying"]["nee"]["coefficients"]
    raw = {k for k in CFG["med_dose_raw"] if not k.startswith("_")}
    conv = {k for k in CFG["med_dose_converted"] if not k.startswith("_")}
    assert set(coef) <= raw, f"NEE drugs with no raw bound: {sorted(set(coef) - raw)}"
    assert set(coef) <= conv, f"NEE drugs with no converted bound: {sorted(set(coef) - conv)}"


def test_every_declared_lab_covariate_has_a_bound():
    labs = {k: v for k, v in COV["time_varying"].items()
            if not k.startswith("_") and (v.get("source") or {}).get("table") == "labs"}
    bounds = {k for k in CFG["analysis_unit_bounds"]["labs"] if not k.startswith("_")}
    missing = {v["source"]["category"] for v in labs.values()} - bounds
    assert not missing, f"declared lab covariates with no outlier bound: {sorted(missing)}"


def test_fio2_bound_matches_the_fraction_band_it_is_applied_after():
    band = COV["time_varying"]["oxygenation"]["fio2_scale"]["fraction_band"]
    assert CFG["analysis_unit_bounds"]["respiratory_support"]["fio2_set"] == band, (
        "the fio2 outlier bound and the fio2 fraction band are the same threshold "
        "written in two files; they must not drift"
    )


def test_declared_divergences_are_real():
    """The config claims three raw/converted divergences. If a bound is edited so
    one stops being true, the note becomes a false claim about the code."""
    raw, conv = CFG["med_dose_raw"], CFG["med_dose_converted"]
    assert raw["epinephrine"]["mcg/kg/min"][1] < conv["epinephrine"][1]
    assert raw["dopamine"]["mcg/kg/min"][1] > conv["dopamine"][1]
    assert raw["dobutamine"]["mcg/kg/min"][1] < conv["dobutamine"][1]


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
