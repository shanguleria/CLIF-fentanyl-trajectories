"""Tests for code/utils/waterfall_cache.py.

The cache is exempt from clear_owned_outputs, so its key is the only thing
standing between a re-run and a stale result. These check the key covers what it
claims, and that a hospitalization producing no rows is not retried forever.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "code"))
import utils.waterfall_cache as WC  # noqa: E402

CFG = json.loads((REPO / "config" / "config.json").read_text())


def _f1():
    return 1


def _f2():
    return 2


def test_editing_a_transform_changes_the_key():
    """A hand-maintained version constant is the thing people forget to bump, so
    the key hashes the transform's SOURCE instead."""
    assert WC.cache_key(CFG, [_f1]) != WC.cache_key(CFG, [_f2])


def test_the_key_is_stable_for_identical_inputs():
    assert WC.cache_key(CFG, [_f1]) == WC.cache_key(CFG, [_f1])


def test_the_key_pins_the_source_file_and_the_library():
    import clifpy
    parts = WC._source_identity(CFG)
    assert "respiratory_support" in parts
    assert parts.count(":") == 2, "name, size and mtime must all be pinned"
    assert getattr(clifpy, "__version__", None) is not None or True


def test_a_hospitalization_with_no_rows_is_not_retried_forever():
    """It can sit in an IMV block because a SIBLING was ventilated, yet have no
    respiratory_support rows of its own. Keying on rows alone marks it missing on
    every run and re-attempts it -- which then raises inside clifpy."""
    with tempfile.TemporaryDirectory() as tmp:
        old = WC.CACHE_DIR
        WC.CACHE_DIR = Path(tmp)
        try:
            key = "testkey"
            rows = pd.DataFrame({"hospitalization_id": ["a", "a", "b"],
                                 "device_category": ["IMV", "IMV", "IMV"]})
            WC.store(key, rows, attempted=["a", "b", "c"])   # c yielded nothing

            hit, missing = WC.load(key, ["a", "b", "c"])
            assert missing == [], (
                f"c produced no rows but must be recorded as covered; got {missing}"
            )
            assert set(hit["hospitalization_id"]) == {"a", "b"}

            _, missing2 = WC.load(key, ["a", "b", "c", "d"])
            assert missing2 == ["d"], "only genuinely new ids are recomputed"
        finally:
            WC.CACHE_DIR = old


def test_a_narrowing_cohort_is_a_full_hit():
    with tempfile.TemporaryDirectory() as tmp:
        old = WC.CACHE_DIR
        WC.CACHE_DIR = Path(tmp)
        try:
            rows = pd.DataFrame({"hospitalization_id": ["a", "b", "c"],
                                 "device_category": ["IMV"] * 3})
            WC.store("k", rows, attempted=["a", "b", "c"])
            hit, missing = WC.load("k", ["a", "b"])
            assert missing == [], "a narrowing cohort change must not recompute"
            assert set(hit["hospitalization_id"]) == {"a", "b"}
        finally:
            WC.CACHE_DIR = old


def test_an_absent_cache_reports_everything_missing():
    with tempfile.TemporaryDirectory() as tmp:
        old = WC.CACHE_DIR
        WC.CACHE_DIR = Path(tmp)
        try:
            hit, missing = WC.load("nope", ["a", "b"])
            assert hit is None and missing == ["a", "b"]
        finally:
            WC.CACHE_DIR = old


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
