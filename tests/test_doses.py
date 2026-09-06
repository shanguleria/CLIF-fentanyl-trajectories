"""Tests for code/utils/doses.py.

Every conversion is checked against a value computed by hand, because a factor
that is wrong by 60 or by 1000 produces numbers that still look like doses.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "code"))
from utils.doses import DoseUnitError, TARGET_OF, convert  # noqa: E402

W = 70.0


def one(drug, dose, unit, weight=W):
    df = pd.DataFrame({"med_category": [drug], "med_dose": [dose],
                       "med_dose_unit": [unit], "weight_kg": [weight]})
    out, rep = convert(df, "med_category", "med_dose", "med_dose_unit")
    return out.iloc[0], rep


def test_targets_are_declared_for_every_nee_drug():
    import json
    cov = json.loads((REPO / "config" / "covariates.json").read_text())
    for drug in cov["time_varying"]["nee"]["coefficients"]:
        assert drug in TARGET_OF, f"{drug} has no target unit"
    assert TARGET_OF["vasopressin"] == "u/min"
    assert TARGET_OF["norepinephrine"] == "mcg/kg/min"


def test_angiotensin_ng_per_kg_min_is_divided_by_1000():
    """100% of angiotensin at UCMC is ng/kg/min, and the NEE coefficient is 10,
    so an error here is a 1000x error on the sickest patients."""
    assert one("angiotensin", 20.0, "ng/kg/min")[0] == 0.02


def test_weight_based_targets_divide_by_weight():
    assert abs(one("phenylephrine", 50.0, "mcg/min")[0] - 50.0 / W) < 1e-12
    assert abs(one("norepinephrine", 10.0, "mcg/min")[0] - 10.0 / W) < 1e-12
    assert abs(one("norepinephrine", 0.7, "mg/min")[0] - 700.0 / W) < 1e-9


def test_per_hour_units_are_divided_by_60():
    assert abs(one("norepinephrine", 6.0, "mg/hr")[0] - 100.0 / W) < 1e-9
    assert abs(one("dopamine", 60.0, "mcg/hr")[0] - 1.0 / W) < 1e-12


def test_vasopressin_units_per_min_passes_through():
    assert one("vasopressin", 0.04, "units/min")[0] == 0.04


def test_vasopressin_unit_matching_is_case_insensitive():
    """UCMC charts 'Units/min'; the schema says 'units/min'. A case-sensitive
    lookup left 277,720 rows unconverted."""
    assert one("vasopressin", 0.04, "Units/min")[0] == 0.04
    assert one("vasopressin", 0.04, "  UNITS/MIN ")[0] == 0.04


def test_vasopressin_milli_units_are_multiplied_by_weight_and_divided_by_1000():
    got = one("vasopressin", 1.43, "milli-units/kg/min")[0]
    assert abs(got - 1.43 * W / 1000.0) < 1e-12
    assert abs(got - 0.1001) < 1e-4


def test_the_reference_dimensional_check_reproduces():
    """CRRT-dose-lmtp verifies its table with two conversions; ours must agree.
    angiotensin 10 ng/kg/min -> 0.01 mcg/kg/min, x10 = 0.1 NEE.
    vasopressin 0.04 U/min, x2.5 = 0.1 NEE."""
    import json
    coef = json.loads((REPO / "config" / "covariates.json").read_text()
                      )["time_varying"]["nee"]["coefficients"]
    assert abs(one("angiotensin", 10.0, "ng/kg/min")[0] * coef["angiotensin"] - 0.1) < 1e-12
    assert abs(one("vasopressin", 0.04, "units/min")[0] * coef["vasopressin"] - 0.1) < 1e-12


def test_an_uncovered_unit_raises_rather_than_dropping():
    try:
        one("norepinephrine", 1.0, "mcg/kg/fortnight")
    except DoseUnitError as e:
        assert "mcg/kg/fortnight" in str(e) and "rather than dropping" in str(e)
    else:
        raise AssertionError("an uncovered unit must raise")


def test_a_missing_weight_nulls_the_row_and_is_counted():
    val, rep = one("phenylephrine", 50.0, "mcg/min", weight=float("nan"))
    assert pd.isna(val)
    assert rep.n_no_weight == 1


def test_a_non_weight_unit_converts_without_a_weight():
    val, rep = one("norepinephrine", 0.1, "mcg/kg/min", weight=float("nan"))
    assert val == 0.1 and rep.n_no_weight == 0


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
