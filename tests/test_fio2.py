"""
Tests for code/utils/fio2.py.

The point of these is the distinction the module rests on: SCALE is decided per
column, BOUNDS are applied per value. A test suite that only checked "0.5 stays
0.5 and 50 becomes 0.5" would pass on an implementation that rescales individual
outliers, which is the behaviour this module exists to refuse.

Run standalone:  .venv/bin/python tests/test_fio2.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "code"))
from utils.fio2 import Fio2ScaleError, detect_fio2_scale, normalize_fio2  # noqa: E402


def test_fraction_column_is_detected_and_left_alone():
    v = pd.Series([0.21, 0.30, 0.40, 0.50, 1.0] * 20)
    out, rep = normalize_fio2(v)
    assert rep.scale == "fraction"
    assert not rep.rescaled
    assert rep.n_nulled == 0
    assert out.equals(v.astype(float))


def test_percent_column_is_rescaled_as_a_whole():
    v = pd.Series([21.0, 30.0, 40.0, 50.0, 100.0] * 20)
    out, rep = normalize_fio2(v)
    assert rep.scale == "percent"
    assert rep.rescaled
    assert rep.n_nulled == 0
    assert out.max() == 1.0 and out.min() == 0.21


def test_an_outlier_in_a_fraction_column_is_NULLED_not_rescaled():
    """The core rule. 88880 in a fractional column is a data entry error, not a
    unit. Rescaling it would manufacture a plausible-looking 888.8, or with a
    /100 heuristic an 888.8 that then looks like a real FiO2."""
    v = pd.Series([0.4] * 99 + [88880.0])
    out, rep = normalize_fio2(v)
    assert rep.scale == "fraction", "one wild value must not flip the column's scale"
    assert not rep.rescaled
    assert rep.n_nulled == 1
    assert pd.isna(out.iloc[-1]), "the outlier must be NULL, not a rescaled number"
    assert out.dropna().eq(0.4).all()


def test_a_percent_value_in_a_fraction_column_is_nulled_not_converted():
    """50 sitting in an otherwise-fractional column is ambiguous at the value
    level: it could be 50% or a typo. The column says fraction, so it is out of
    range and gets nulled. Converting it would be guessing."""
    v = pd.Series([0.4] * 99 + [50.0])
    out, rep = normalize_fio2(v)
    assert rep.scale == "fraction"
    assert rep.n_nulled == 1
    assert pd.isna(out.iloc[-1])


def test_mixed_units_raise_rather_than_guess():
    v = pd.Series([0.4] * 50 + [40.0] * 50)
    try:
        normalize_fio2(v)
    except Fio2ScaleError as e:
        assert "neither clearly fraction nor clearly percent" in str(e)
    else:
        raise AssertionError("a 50/50 mixed-unit column must raise, not pick a side")


def test_values_between_the_bands_push_toward_ambiguous():
    """(1, 21) is too high for a fraction and below room air as a percent. Such
    values must not be quietly absorbed into either band."""
    v = pd.Series([5.0] * 60 + [0.4] * 40)
    try:
        normalize_fio2(v)
    except Fio2ScaleError as e:
        assert "in neither" in str(e)
    else:
        raise AssertionError("a column dominated by (1, 21) values must raise")


def test_all_null_column_raises():
    try:
        detect_fio2_scale(pd.Series([None, None], dtype="float64"))
    except Fio2ScaleError as e:
        assert "no non-null" in str(e)
    else:
        raise AssertionError("an all-null fio2 column must raise, not report a scale")


def test_nulls_do_not_affect_the_scale_decision():
    v = pd.Series([0.4, None, 0.5, None, 0.21] * 20)
    out, rep = normalize_fio2(v)
    assert rep.scale == "fraction"
    assert rep.n_nonnull == 60 and rep.n_total == 100
    assert out.isna().sum() == 40


def test_non_numeric_values_are_coerced_to_null_not_crashed_on():
    v = pd.Series(["0.4", "0.5", "bad", "0.21"] * 25)
    out, rep = normalize_fio2(v)
    assert rep.scale == "fraction"
    assert out.isna().sum() == 25


def test_report_is_printable_and_names_what_it_nulled():
    v = pd.Series([0.4] * 99 + [88880.0])
    _, rep = normalize_fio2(v)
    text = str(rep)
    assert "scale detected: fraction" in text
    assert "nulled as out of range: 1" in text
    assert "88880" in text, "the report must show what it discarded, not just how many"


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
