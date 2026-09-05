"""Tests for the pure logic in code/01_build_cohort.py, on synthetic frames.

No CLIF data is read. The stages that need real tables are exercised in a
smoke run; what is tested here is the logic that is wrong-or-right independent
of the data: LOCF caps, absence-means-zero, the consumption assertion, unit
conversion, and the Severinghaus transform.
"""
from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "code"))
_spec = importlib.util.spec_from_file_location("build", REPO / "code" / "01_build_cohort.py")
B = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(B)


def _long(n_blocks=2, n_win=None):
    n_win = n_win or B.N_WINDOWS
    rows = []
    for b in range(n_blocks):
        for w in range(n_win):
            rows.append({"encounter_block": f"b{b}", "window_idx": w,
                         "window_start_hr": w * B.WINDOW_H, "at_risk": True})
    return pd.DataFrame(rows)


# --------------------------------------------------------------- LOCF caps
def test_locf_respects_the_declared_cap_in_windows_not_rows():
    """A 24h cap on a 4h grid is 6 windows, not 24."""
    df = _long(n_blocks=1)
    df["lactate"] = np.nan
    df.loc[0, "lactate"] = 4.2
    out, _, _ = B.apply_missingness(df)
    filled = out["lactate"].notna().sum()
    cap_h = B.COV["time_varying"]["lactate"]["locf"]["cap_hours"]
    assert filled == 1 + cap_h // B.WINDOW_H == 7, f"carried {filled} windows, expected 7"


def test_locf_flag_marks_only_filled_cells():
    df = _long(n_blocks=1)
    df["bun"] = np.nan
    df.loc[0, "bun"] = 30.0
    out, _, _ = B.apply_missingness(df)
    assert not out.loc[0, "bun_locf"], "an observed value must not be flagged as carried"
    assert out.loc[1, "bun_locf"], "a carried value must be flagged"
    assert out["bun_locf"].sum() == 6


def test_locf_does_not_cross_encounter_blocks():
    df = _long(n_blocks=2)
    df["lactate"] = np.nan
    df.loc[df["encounter_block"].eq("b0") & df["window_idx"].eq(B.N_WINDOWS - 1),
           "lactate"] = 9.9
    out, _, _ = B.apply_missingness(df)
    b1 = out[out["encounter_block"] == "b1"]
    assert b1["lactate"].isna().all(), "a value must never carry into the next block"


# ------------------------------------------------------- absence means zero
def test_absence_means_zero_variables_are_zeroed_not_carried():
    df = _long(n_blocks=1)
    df["nee"] = np.nan
    df.loc[0, "nee"] = 0.4
    out, _, _ = B.apply_missingness(df)
    assert out.loc[1, "nee"] == 0.0, "an absent vasopressor record means 0, not the prior dose"
    assert "nee_locf" not in out.columns, "absence-means-zero variables are not LOCF-eligible"


def test_imv_status_absence_means_not_ventilated():
    df = _long(n_blocks=1)
    df["imv_status"] = np.nan
    df.loc[0, "imv_status"] = 1.0
    out, _, _ = B.apply_missingness(df)
    assert out.loc[5, "imv_status"] == 0.0


def test_missingness_is_counted_before_the_fill():
    df = _long(n_blocks=1)
    df["lactate"] = np.nan
    df.loc[0, "lactate"] = 4.2
    _, rep, _ = B.apply_missingness(df)
    row = rep[rep["variable"] == "lactate"].iloc[0]
    assert row["n_missing_pre_locf"] == B.N_WINDOWS - 1, (
        "counting after the fill would hide how much was carried"
    )


def test_missingness_denominator_is_at_risk_rows_only():
    df = _long(n_blocks=1)
    df.loc[df["window_idx"] >= 10, "at_risk"] = False
    df["lactate"] = np.nan
    _, rep, _ = B.apply_missingness(df)
    assert rep[rep["variable"] == "lactate"].iloc[0]["n_at_risk"] == 10


def test_report_separates_pre_locf_from_final_missingness():
    """The point of the report: how much was carried, and what is still absent."""
    df = _long(n_blocks=1)
    df["lactate"] = np.nan
    df.loc[0, "lactate"] = 4.2
    _, rep, _ = B.apply_missingness(df)
    r = rep[rep["variable"] == "lactate"].iloc[0]
    assert r["n_observed"] == 1
    assert r["n_missing_pre_locf"] == B.N_WINDOWS - 1
    assert r["n_filled_by_locf"] == 6
    assert r["n_missing_final"] == B.N_WINDOWS - 7
    assert r["n_observed"] + r["n_filled_by_locf"] + r["n_missing_final"] == r["n_at_risk"]


def test_report_shows_zero_by_rule_separately_from_missing():
    """An absent vasopressor record is a zero, not a gap, and must not read as one."""
    df = _long(n_blocks=1)
    df["nee"] = np.nan
    df.loc[0, "nee"] = 0.4
    _, rep, _ = B.apply_missingness(df)
    r = rep[rep["variable"] == "nee"].iloc[0]
    assert r["n_zero_by_rule"] == B.N_WINDOWS - 1
    assert r["n_missing_final"] == 0
    assert r["pct_missing_final"] == 0.0


def test_report_percentages_use_the_at_risk_denominator():
    df = _long(n_blocks=1)
    df.loc[df["window_idx"] >= 10, "at_risk"] = False
    df["bun"] = np.nan
    _, rep, _ = B.apply_missingness(df)
    r = rep[rep["variable"] == "bun"].iloc[0]
    assert r["n_at_risk"] == 10 and r["pct_missing_final"] == 100.0


def test_pattern_table_pools_small_cells():
    min_cell = B.COV["missing_values"]["reporting"]["pattern_min_cell"]
    df = _long(n_blocks=1)
    for v in ("bun", "lactate"):
        df[v] = 1.0
    df.loc[0, "bun"] = np.nan          # a pattern seen once
    _, _, pat = B.apply_missingness(df)
    assert len(pat), "a pattern table must be produced"
    pooled = pat[pat["pattern"].astype(str).str.startswith("<other")]
    assert len(pooled) == 1, f"rare patterns must be pooled below {min_cell}"


# ------------------------------------------------- the consumption assertion
def test_consumption_assertion_fires_on_a_declared_but_absent_column():
    df = _long(n_blocks=1)
    try:
        B.assert_config_is_honoured(df)
    except SystemExit as e:
        assert "silent no-op" in str(e)
    else:
        raise AssertionError("a declared covariate with no column must fail the build")


def test_consumption_assertion_passes_when_every_declared_column_exists():
    df = _long(n_blocks=1)
    for k in list(B.COV["time_varying"]) + list(B.COV["time_invariant"]):
        if not k.startswith("_"):
            df[k] = 0.0
    for c in B.EXPOSURE["columns"]:
        df[c] = 0.0
    B.assert_config_is_honoured(df)


# ------------------------------------------------------- units and transforms
def test_infusion_unit_conversion_covers_the_charted_units():
    m = pd.DataFrame({
        "med_dose": [2.0, 100.0, 0.1, 0.05, 2.0],
        "med_dose_unit": ["mcg/kg/hr", "mcg/hr", "mg/hr", "mcg/kg/min", "mcg/min"],
        "weight_kg": [50.0] * 5,
    })
    out = B._to_mcg_kg_hr(m)
    assert out.iloc[0] == 2.0
    assert out.iloc[1] == 2.0            # 100 mcg/hr / 50 kg
    assert out.iloc[2] == 2.0            # 0.1 mg/hr -> 100 mcg/hr / 50 kg
    assert abs(out.iloc[3] - 3.0) < 1e-9  # 0.05 mcg/kg/min * 60
    assert abs(out.iloc[4] - 2.4) < 1e-9  # 2 mcg/min * 60 / 50


def test_an_unhandled_infusion_unit_raises_rather_than_dropping():
    m = pd.DataFrame({"med_dose": [1.0], "med_dose_unit": ["mcg/kg/day"],
                      "weight_kg": [70.0]})
    try:
        B._to_mcg_kg_hr(m)
    except SystemExit as e:
        assert "mcg/kg/day" in str(e)
    else:
        raise AssertionError("an unhandled dose unit must raise, not silently drop exposure")


def test_severinghaus_matches_the_reference_value():
    """Severinghaus(97) = 90.6 mmHg, per the reference implementation."""
    v = B._severinghaus(pd.Series([97.0]))
    assert abs(v.iloc[0] - 90.6) < 0.2, f"got {v.iloc[0]}"


def test_severinghaus_is_undefined_on_the_plateau():
    assert pd.isna(B._severinghaus(pd.Series([100.0])).iloc[0])


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
