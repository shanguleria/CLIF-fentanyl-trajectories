"""The federated-pooling exports must be poolable, exactly.

A median cannot be pooled across sites; a mean can, from n and the two sums. The
point of carrying sum and sum_sq rather than only mean and sd is that the pooled
figures are then exact rather than an approximation that assumes equal variances.
These checks run only when Phase 1 has produced the files.

Run standalone:  .venv/bin/python tests/test_pooling.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "output" / "final_no_phi"
CONT = OUT / "phase1_pooling_continuous.csv"
CAT = OUT / "phase1_pooling_categorical.csv"
TOL = 1e-5          # the exports are rounded to 6 decimals


def _skip_if_absent(f: Path) -> pd.DataFrame | None:
    if not f.exists():
        print(f"  SKIP  {f.name} absent -- run code/02_descriptive_trajectory.R")
        return None
    return pd.read_csv(f)


def test_mean_is_recoverable_from_sum_and_n():
    d = _skip_if_absent(CONT)
    if d is None:
        return
    d = d[d.n > 0].dropna(subset=["mean", "sum"])
    assert len(d), "no usable rows"
    err = (d["sum"] / d["n"] - d["mean"]).abs().max()
    assert err < TOL, f"sum/n disagrees with mean by {err}"


def test_sd_is_recoverable_from_the_two_sums():
    d = _skip_if_absent(CONT)
    if d is None:
        return
    d = d[d.n > 1].dropna(subset=["sd", "sum", "sum_sq"])
    var = (d["sum_sq"] - d["sum"] ** 2 / d["n"]) / (d["n"] - 1)
    err = (var.clip(lower=0) ** 0.5 - d["sd"]).abs().max()
    assert err < TOL, f"reconstructed sd disagrees by {err}"


def test_two_strata_pool_back_to_the_overall_row():
    """The whole point, demonstrated on real numbers: eligible and not-eligible
    are disjoint, so pooling them as if they were two sites must reproduce the
    overall row exactly."""
    d = _skip_if_absent(CONT)
    if d is None:
        return
    base = d[d.scope == "baseline"]
    checked = 0
    for var, g in base.groupby("variable"):
        s = {r.stratum: r for r in g.itertuples()}
        if not {"overall", "eligible", "not_eligible"} <= set(s):
            continue
        a, b, o = s["eligible"], s["not_eligible"], s["overall"]
        if any(math.isnan(x) for x in (a.sum, b.sum, o.mean)) or a.n + b.n == 0:
            continue
        n = a.n + b.n
        assert n == o.n, f"{var}: strata sum to {n}, overall says {o.n}"
        assert abs((a.sum + b.sum) / n - o.mean) < TOL, f"{var}: pooled mean differs"
        q = a.sum_sq + b.sum_sq
        pooled_sd = math.sqrt(max((q - (a.sum + b.sum) ** 2 / n) / (n - 1), 0.0))
        assert abs(pooled_sd - o.sd) < TOL, f"{var}: pooled sd differs"
        checked += 1
    assert checked >= 5, f"only {checked} variables had all three strata"


def test_small_cells_are_suppressed_not_published():
    """A mean over n = 1 is that patient's value. Anything below the site's
    small_cell_min_den must carry no statistics at all."""
    d = _skip_if_absent(CONT)
    if d is None:
        return
    cfg = json.loads((REPO / "config" / "config_template.json").read_text())
    lim = cfg["reporting"]["small_cell_min_den"]
    tiny = d[(d.n > 0) & (d.n < lim)]
    assert (tiny["n_suppressed_small_cell"] == 1).all(), (
        f"{int((tiny['n_suppressed_small_cell'] != 1).sum())} cells below n={lim} "
        f"are not flagged"
    )
    for col in ("mean", "sd", "sum", "sum_sq", "median", "min", "max"):
        assert tiny[col].isna().all(), f"{col} is published for a cell below n={lim}"


def test_categorical_counts_sum_to_their_denominator():
    d = _skip_if_absent(CAT)
    if d is None:
        return
    live = d[d.n_suppressed_small_cell == 0]
    for (var, st), g in live.groupby(["variable", "stratum"]):
        if (d[(d.variable == var) & (d.stratum == st)]
                .n_suppressed_small_cell == 1).any():
            continue                      # a suppressed level breaks the sum
        assert g["n"].sum() == g["denominator"].iloc[0], (
            f"{var}/{st}: levels sum to {g['n'].sum()} against a denominator of "
            f"{g['denominator'].iloc[0]}"
        )


def test_the_collapsed_race_matches_the_full_categories_level_by_level():
    """Collapsing must MERGE, not reshuffle.

    The first implementation used `unlist(map[v])`, which drops the NULLs for
    unmatched values -- so the replacement vector came back shorter than the
    input, ifelse recycled it, and every row was silently misaligned. The
    published table said 59.8% of the landmark cohort were Black against a true
    62.6%. Totals still matched, so only a level-by-level check catches it.
    """
    d = _skip_if_absent(CAT)
    if d is None:
        return
    full = d[d.variable == "Race"]
    coll = d[d.variable == "Race (collapsed)"]
    assert len(coll), "the collapsed race rows are missing"
    cmap = json.loads((REPO / "config" / "covariates.json").read_text())
    cmap = cmap["time_invariant"]["race"]["reporting_collapse"]
    named = {k: v for k, v in cmap.items() if not k.startswith("_")}
    for st in full.stratum.unique():
        f = full[full.stratum == st].set_index("level")["n"]
        c = coll[coll.stratum == st].set_index("level")["n"]
        expected: dict[str, int] = {}
        for lvl, k in f.items():
            tgt = named.get(lvl, cmap["_default"]) if lvl != "Missing" else "Missing"
            expected[tgt] = expected.get(tgt, 0) + int(k)
        for tgt, k in expected.items():
            assert int(c.get(tgt, 0)) == k, (
                f"{st}/{tgt}: collapsed says {c.get(tgt)}, the full categories "
                f"sum to {k}"
            )


def test_the_collapsed_race_preserves_the_total():
    d = _skip_if_absent(CAT)
    if d is None:
        return
    full = d[(d.variable == "Race") & (d.stratum == "overall")]
    coll = d[(d.variable == "Race (collapsed)") & (d.stratum == "overall")]
    assert len(coll), "the collapsed race rows are missing"
    assert full["n"].sum() == coll["n"].sum(), (
        "collapsing race changed the total, so a level was dropped rather than merged"
    )
    assert set(coll["level"]) <= {"Black", "White", "Other", "Unknown", "Missing"}


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
