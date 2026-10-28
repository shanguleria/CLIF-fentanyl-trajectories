"""The tracked configs must PARSE. Trivial, and it would have saved a bad tag.

On 2026-10-27 a syntactically invalid config/protocol.json -- one missing comma
-- was committed and TAGGED as v0.12.0. Nothing caught it: the 240-test suite
reads the WORKING TREE, which had already been repaired, while the commit did
not. Every script in the pipeline opens this file first, so the tag could not
execute a single line.

This runs in the preflight and reads the files from disk. It is deliberately
dumber than every other test here -- it asserts only that the bytes are valid
JSON and that the versions are present -- because the failure it guards is not
subtle, it is just unguarded.

Run standalone:  .venv/bin/python tests/test_config_parses.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
TRACKED = ("config/protocol.json", "config/outlier_config.json",
           "config/pooling_spec.json", "config/config_template.json")


@pytest.mark.parametrize("rel", TRACKED)
def test_the_tracked_config_is_valid_json(rel):
    f = REPO / rel
    if not f.is_file():
        pytest.skip(f"{rel} is not present in this checkout")
    try:
        json.loads(f.read_text())
    except json.JSONDecodeError as e:
        raise AssertionError(
            f"{rel} is not valid JSON: {e}. Every script reads this file before "
            f"it does anything else, so a committed syntax error means the whole "
            f"pipeline cannot run -- which is how v0.12.0 was tagged unrunnable."
        ) from e


def test_every_versioned_config_declares_its_version():
    """A protocol with no definition_version cannot be compared across sites, and
    check_provenance.py would report agreement by reading None == None."""
    for rel in ("config/protocol.json", "config/outlier_config.json",
                "config/pooling_spec.json"):
        f = REPO / rel
        if not f.is_file():
            continue
        d = json.loads(f.read_text())
        v = d.get("definition_version")
        assert isinstance(v, str) and v, f"{rel} must declare definition_version"
        assert v.count(".") == 2, (
            f"{rel} definition_version is {v!r}; the digest-comparison machinery "
            f"and every tag annotation assume major.minor.patch")


def test_the_site_config_is_valid_json_if_present():
    """config.json is gitignored, so this is about the LOCAL run, not the commit.
    A site editing it by hand is the most likely source of a syntax error here."""
    f = REPO / "config" / "config.json"
    if not f.is_file():
        pytest.skip("config/config.json absent -- copy config_template.json")
    try:
        json.loads(f.read_text())
    except json.JSONDecodeError as e:
        raise AssertionError(
            f"config/config.json is not valid JSON: {e}. Fix it before running; "
            f"every script reads it first."
        ) from e


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
