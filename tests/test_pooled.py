"""The pooled bundle must be arithmetically consistent with the site rows it came from.

output/pooled/ is the coordinating centre's product and had no test of its own
until 2026-10-31, when `pool_sum_counts()` was found summing a numerator over
the sites that published a cell while the denominator kept every site -- a
percentage over the wrong base, biased down, on a Table 1 race row. Nothing
would have caught it: the site suites validate site exports, and the centre
suite validates the inverse-variance maths.

Every check here re-derives the POOLED row from the per-site rows in the SAME
file, so it is independent of the engine that wrote it. Skips cleanly when
output/pooled/ is absent, which it is at a site.

Run standalone:  .venv/bin/python tests/test_pooled.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
POOL = REPO / "output" / "pooled"
SPEC = json.loads((REPO / "config" / "pooling_spec.json").read_text())


def _files():
    for rel, e in SPEC["files"].items():
        f = POOL / rel.replace("/", "__")
        if f.is_file():
            yield rel, e, pd.read_csv(f)


def _groups(d, keys):
    keys = [k for k in keys if k in d.columns]
    if not keys:
        return [("", d)]
    return list(d.groupby(keys, dropna=False))


def test_output_pooled_exists_or_this_suite_says_so():
    if not POOL.is_dir():
        print("  SKIP  output/pooled absent -- run code/pool_sites.R")
        return
    assert list(_files()), "output/pooled exists but claims no spec'd file"


def test_every_pooled_percentage_is_its_own_row_arithmetic():
    """pct must be 100 * numerator / denominator READ OFF THE SAME ROW.

    This is the check that fails when the numerator and the denominator are
    summed over different sets of sites: each is individually plausible and
    only their ratio is wrong.
    """
    checked = 0
    for rel, e, d in _files():
        den = e.get("denominator")
        if not den or den not in d.columns or e.get("denominator_within"):
            continue
        for pc, num in (e.get("percents") or {}).items():
            if pc not in d.columns or num not in d.columns:
                continue
            # The pooler rounds each percentage to the decimals the SITES
            # shipped for that column -- 1 in dose_summary, 3 in the
            # histograms -- so the tolerance is half a unit in the last place
            # of the value actually written, not a fixed epsilon.
            p = d[d.site == "POOLED"]
            for r in p.itertuples():
                n, dd, got = getattr(r, num), getattr(r, den), getattr(r, pc)
                if pd.isna(n) or pd.isna(dd) or pd.isna(got) or dd == 0:
                    continue
                txt = f"{got}"
                dec = len(txt.split(".")[1]) if "." in txt else 0
                want = 100 * n / dd
                assert abs(want - got) <= 0.5 * 10 ** -dec + 1e-9, (
                    f"{rel}: POOLED {pc} is {got} but {num}/{den} = "
                    f"{n}/{dd} = {want:.4f} -- the percentage is over a "
                    f"different base than the row reports"
                )
                checked += 1
    assert checked, "no pooled percentage was checked"
    print(f"  ..    {checked} pooled percentage(s) re-derived")


def test_a_pooled_count_is_the_sum_of_the_sites_behind_it():
    """And the sites behind it are the ones that did not suppress."""
    checked = 0
    for rel, e, d in _files():
        if e.get("method") not in ("sum_counts", "count_matrix"):
            continue
        flag = e.get("suppressed_flag")
        # THE DENOMINATOR IS RE-DERIVED TOO, and it is the one that matters: a
        # numerator summed over the sites that published, divided by a
        # denominator summed over all of them, is wrong while both halves look
        # right and the row's own pct/n/denominator stay self-consistent.
        nums = [n for n in list(e.get("numerators") or []) + [e.get("denominator")]
                if n and n in d.columns]
        for _, g in _groups(d, e.get("keys") or []):
            sites = g[g.site != "POOLED"]
            pooled = g[g.site == "POOLED"]
            if len(pooled) != 1 or sites.empty:
                continue
            use = sites
            if flag and flag in sites.columns:
                use = sites[sites[flag] != 1]
                if use.empty:
                    continue
            for nm in nums:
                want = use[nm].dropna()
                got = pooled[nm].iloc[0]
                if want.empty or pd.isna(got):
                    continue
                assert abs(want.sum() - got) < 1e-6, (
                    f"{rel}: POOLED {nm} is {got}, the contributing sites sum "
                    f"to {want.sum()}"
                )
                checked += 1
    assert checked, "no pooled count was checked"
    print(f"  ..    {checked} pooled count(s) re-derived")


def test_k_sites_counts_the_sites_that_actually_contributed():
    """A row pooled over fewer sites than the bundle has must say so, or the
    reader takes a two-site number for a three-site one."""
    checked = 0
    for rel, e, d in _files():
        if "k_sites" not in d.columns or e.get("method") not in (
                "sum_counts", "count_matrix"):
            continue
        flag = e.get("suppressed_flag")
        n_sites = d.loc[d.site != "POOLED", "site"].nunique()
        for _, g in _groups(d, e.get("keys") or []):
            pooled = g[g.site == "POOLED"]
            if len(pooled) != 1 or pd.isna(pooled["k_sites"].iloc[0]):
                continue
            sites = g[g.site != "POOLED"]
            if flag and flag in sites.columns:
                sites = sites[sites[flag] != 1]
            got = int(pooled["k_sites"].iloc[0])
            assert got == sites["site"].nunique(), (
                f"{rel}: k_sites is {got} but {sites['site'].nunique()} "
                f"site(s) contributed to that key"
            )
            assert got <= n_sites, f"{rel}: k_sites {got} exceeds {n_sites} sites"
            checked += 1
    assert checked, "no k_sites value was checked"
    print(f"  ..    {checked} k_sites value(s) re-derived")


def test_the_pooled_provenance_names_every_site():
    f = POOL / "pooled_provenance.json"
    if not f.is_file():
        print("  SKIP  pooled_provenance.json absent")
        return
    prov = json.loads(f.read_text())
    assert prov.get("sites"), "pooled_provenance.json names no sites"
    for _, _, d in _files():
        have = set(d.loc[d.site != "POOLED", "site"].dropna())
        assert have <= set(prov["sites"]), (
            f"a pooled file carries site(s) the provenance does not name: "
            f"{sorted(have - set(prov['sites']))}")


if __name__ == "__main__":
    if not POOL.is_dir():
        print("  SKIP  output/pooled absent -- run code/pool_sites.R")
        sys.exit(0)
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    bad = 0
    for fn in fns:
        try:
            fn()
            print(f"  PASS  {fn.__name__}")
        except AssertionError as exc:
            bad += 1
            print(f"  FAIL  {fn.__name__}\n        {exc}")
    print(f"\n{len(fns) - bad}/{len(fns)} passed")
    sys.exit(1 if bad else 0)
