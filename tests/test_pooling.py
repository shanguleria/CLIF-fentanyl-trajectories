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
# out_final is subdivided by the script that produced each file (2026-09-08).
CONT = OUT / "02_descriptive" / "pooling_continuous.csv"
CAT = OUT / "02_descriptive" / "pooling_categorical.csv"
# Phases 1 and 2 shared one pooling contract and both were held to it, the second
# through 06_landmark/ after the 2026-09-24 renumbers. This test caught the stale
# path the first time and is the reason it is worth keeping accurate.
#
# COVERAGE LOSS, 2026-09-30. landmark_cohort.R was tabled and is no longer run, so
# its pooling exports are not regenerated and cannot be held to the contract
# without re-running a tabled script by hand. Both entries are removed rather than
# left to fail, and the consequence is recorded here rather than silently
# absorbed: ONE producer is now checked, not two. The contract itself is
# unchanged. If a second live script ever emits pooling_*.csv -- 06_unit_variation.R
# deliberately does not, it suppresses directly -- add it to these lists in the
# same change that adds the writer.
ALL_CONT = [CONT]
ALL_CAT = [CAT]
TOL = 1e-5          # the exports are rounded to 6 decimals


def _skip_if_absent(f: Path) -> pd.DataFrame | None:
    if not f.exists():
        print(f"  SKIP  {f.name} absent -- run the phase that writes it")
        return None
    return pd.read_csv(f)


def test_the_pooling_exports_this_suite_checks_actually_exist():
    """Every other check here returns quietly when its input is missing, so a
    relocated output would make the whole suite pass vacuously. It did, once:
    out_final was subdivided by script on 2026-09-08 and these paths still
    pointed at the old flat locations. This is the check that cannot skip."""
    missing = [str(f.relative_to(REPO)) for f in ALL_CONT + ALL_CAT if not f.exists()]
    assert not missing, (
        "pooling exports are missing, so the rest of this suite is skipping "
        "rather than checking:\n  " + "\n  ".join(missing)
    )


def test_mean_is_recoverable_from_sum_and_n():
    for f in ALL_CONT:
        d = _skip_if_absent(f)
        if d is None:
            continue
        d = d[d.n > 0].dropna(subset=["mean", "sum"])
        assert len(d), f"{f.name}: no usable rows"
        err = (d["sum"] / d["n"] - d["mean"]).abs().max()
        assert err < TOL, f"{f.name}: sum/n disagrees with mean by {err}"


def test_sd_is_recoverable_from_the_two_sums():
    for f in ALL_CONT:
        d = _skip_if_absent(f)
        if d is None:
            continue
        d = d[d.n > 1].dropna(subset=["sd", "sum", "sum_sq"])
        var = (d["sum_sq"] - d["sum"] ** 2 / d["n"]) / (d["n"] - 1)
        err = (var.clip(lower=0) ** 0.5 - d["sd"]).abs().max()
        assert err < TOL, f"{f.name}: reconstructed sd disagrees by {err}"


def test_the_strata_pool_back_to_the_overall_row():
    """The whole point, demonstrated on real numbers: the Table 1 strata are
    disjoint and exhaustive, so pooling them as if they were separate sites must
    reproduce the overall row exactly.

    Deliberately agnostic to WHICH strata they are and HOW MANY. Table 1 was
    stratified on landmark eligibility (two strata) until 2026-09-24 and on
    predominant fentanyl intensity (four) after; a test that hardcodes the names
    fails on a change of stratification rather than on a change of arithmetic,
    which is not what it is for.
    """
    d = _skip_if_absent(CONT)
    if d is None:
        return
    base = d[d.scope == "baseline"]
    checked = 0
    for var, g in base.groupby("variable"):
        rows = {r.stratum: r for r in g.itertuples()}
        o = rows.pop("overall", None)
        parts = list(rows.values())
        if o is None or len(parts) < 2:
            continue
        if any(math.isnan(x.sum) or math.isnan(x.sum_sq) for x in parts) \
                or math.isnan(o.mean):
            continue
        n = sum(x.n for x in parts)
        if n == 0 or n != o.n:
            # a stratum suppressed for a small cell cannot be pooled back
            continue
        tot = sum(x.sum for x in parts)
        assert abs(tot / n - o.mean) < TOL, f"{var}: pooled mean differs"
        q = sum(x.sum_sq for x in parts)
        pooled_sd = math.sqrt(max((q - tot ** 2 / n) / (n - 1), 0.0))
        assert abs(pooled_sd - o.sd) < TOL, f"{var}: pooled sd differs"
        checked += 1
    assert checked >= 5, f"only {checked} variables were poolable"


def test_small_cells_are_suppressed_not_published():
    """A mean over n = 1 is that patient's value. Anything below the site's
    small_cell_min_den must carry no statistics at all."""
    for f in ALL_CONT:
        d = _skip_if_absent(f)
        if d is not None:
            _check_suppression(d, f.name)


def _check_suppression(d, name):
    cfg = json.loads((REPO / "config" / "config_template.json").read_text())
    lim = cfg["reporting"]["small_cell_min_den"]
    tiny = d[(d.n > 0) & (d.n < lim)]
    assert (tiny["n_suppressed_small_cell"] == 1).all(), (
        f"{name}: {int((tiny['n_suppressed_small_cell'] != 1).sum())} cells below "
        f"n={lim} are not flagged"
    )
    for col in ("mean", "sd", "sum", "sum_sq", "median", "min", "max"):
        assert tiny[col].isna().all(), f"{name}: {col} published below n={lim}"


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


def test_the_categorical_export_carries_every_stratum_the_continuous_one_does():
    """Both exports describe the same cohort under the same Table 1
    stratification, so their baseline stratum sets must agree.

    This fired on 2026-10-02 with 13 rows where 65 were expected. `pool_cat`
    built its stratum list as `c("overall", sort(unique(strata)))`, and `c()`
    dispatches on its first argument -- a character -- so a FACTOR `strata` was
    coerced to its INTEGER CODES. `strata == "1"` then matched nothing and every
    named stratum was silently dropped. The continuous twin was unaffected
    because its caller loops over `levels()`, which is already character, so the
    two files disagreed while each looked internally consistent and the counts
    that were present were all correct.

    Deliberately agnostic to which strata they are and how many, for the same
    reason as test_the_strata_pool_back_to_the_overall_row.
    """
    for fc, fk in zip(ALL_CONT, ALL_CAT):
        dc, dk = _skip_if_absent(fc), _skip_if_absent(fk)
        if dc is None or dk is None:
            continue
        want = set(dc.loc[dc["scope"] == "baseline", "stratum"].dropna())
        have = set(dk["stratum"].dropna())
        assert want <= have, (
            f"{fk.name} is missing stratum/strata {sorted(want - have)} and "
            f"carries only {sorted(have)}. A pooled categorical Table 1 cannot "
            f"be built without them."
        )


def test_no_suppressed_cell_is_recoverable_by_subtraction():
    """A suppressed cell that arithmetic recovers is not suppressed.

    Measured 2026-10-02: it recovered ALL FOUR. pooling_categorical.csv published
    a variable in two overlapping views (7-level Race, 4-level Race (collapsed))
    across two overlapping partitions (4 bands, overall), giving two systems of
    linear constraints on four unknowns. A brute-force solve returned exactly one
    consistent assignment, so every withheld count was determined.

    This reproduces the attack rather than trusting the fix: propagate the
    "a constraint with exactly one unknown determines it" rule to a fixpoint,
    which is precisely how the original leak unwound.
    """
    d = _skip_if_absent(CAT)
    if d is None:
        return
    cmap = json.loads((REPO / "config" / "covariates.json").read_text())
    collapse = cmap["time_invariant"]["race"]["reporting_collapse"]
    named = {k: v for k, v in collapse.items() if not k.startswith("_")}

    known: dict[tuple, float] = {}
    unknown: set[tuple] = set()
    for _, r in d.iterrows():
        key = (r["variable"], r["stratum"], r["level"])
        if r["n_suppressed_small_cell"] == 1 or pd.isna(r["n"]):
            unknown.add(key)
        else:
            known[key] = float(r["n"])

    # (total, [cells]) for every published identity the table exposes.
    cons: list[tuple[float, list[tuple]]] = []
    for (var, st), g in d.groupby(["variable", "stratum"]):
        # levels in a stratum sum to that stratum's published denominator
        cons.append((float(g["denominator"].iloc[0]),
                     [(var, st, lv) for lv in g["level"]]))
    strata = sorted(set(d["stratum"]) - {"overall"})
    for var, g in d.groupby("variable"):
        present = sorted(set(g["stratum"]) - {"overall"})
        if present != strata:
            continue          # not published per stratum, so nothing to sum
        for lv in sorted(set(g["level"])):
            if ("overall" in set(g["stratum"])):
                # named strata sum to overall, within a level
                cons.append(("overall_of", [(var, "overall", lv)]
                             + [(var, s, lv) for s in strata]))
    for st in sorted(set(d["stratum"])):
        src = d[(d.variable == "Race") & (d.stratum == st)]
        coll = d[(d.variable == "Race (collapsed)") & (d.stratum == st)]
        if src.empty or coll.empty:
            continue          # the overall-only rule has removed the overlap
        for _, c in coll.iterrows():
            members = [("Race", st, lv) for lv in src["level"]
                       if (named.get(lv, collapse["_default"])
                           if lv != "Missing" else "Missing") == c["level"]]
            if members and not pd.isna(c["n"]):
                cons.append((float(c["n"]), members))

    # Propagate to a fixpoint: a constraint with one unknown determines it.
    leaked = []
    changed = True
    while changed:
        changed = False
        for total, cells in cons:
            missing = [c for c in cells if c in unknown]
            if len(missing) != 1:
                continue
            if total == "overall_of":
                head, rest = cells[0], cells[1:]
                if head in unknown or any(c in unknown for c in rest):
                    if missing[0] == head:
                        val = sum(known[c] for c in rest)
                    else:
                        val = known[head] - sum(known[c] for c in rest
                                                if c != missing[0])
                else:
                    continue
            else:
                val = total - sum(known[c] for c in cells if c != missing[0])
            known[missing[0]] = val
            unknown.discard(missing[0])
            leaked.append((missing[0], val))
            changed = True

    assert not leaked, (
        "suppressed cells are recoverable by subtraction, so they are not "
        "suppressed:\n  " + "\n  ".join(
            f"{v} / {s} / {lv} = {val:g}" for (v, s, lv), val in leaked)
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
    # The two views no longer span the same strata: covariates.json
    # disclosure.full_levels_strata_scope ships the FULL levels at `overall`
    # only, so the collapse can be checked against its sources exactly where
    # both are published. Deriving the overlap rather than assuming it is what
    # keeps this a test of the invariant and not of today's disclosure policy.
    overlap = sorted(set(full["stratum"]) & set(coll["stratum"]))
    assert overlap, (
        "Race and Race (collapsed) share no stratum, so the collapse is never "
        "checked against its source levels anywhere"
    )
    n_checked = 0
    checked_strata = set()
    for st in overlap:
        g = full[full.stratum == st]
        c = coll[coll.stratum == st].set_index("level")["n"]
        expected: dict[str, int] = {}
        # A suppressed SOURCE level makes its target unverifiable from the
        # published file -- that is what suppression means. Only the stratum
        # `overall` is fully unsuppressed at UCMC, so without this the check
        # would crash on NaN at a site with any rare category inside a band.
        dirty: set[str] = set()
        for _, row in g.iterrows():
            lvl = row["level"]
            tgt = named.get(lvl, cmap["_default"]) if lvl != "Missing" else "Missing"
            if row["n_suppressed_small_cell"] == 1 or pd.isna(row["n"]):
                dirty.add(tgt)
                continue
            expected[tgt] = expected.get(tgt, 0) + int(row["n"])
        for tgt, k in expected.items():
            if tgt in dirty:
                continue
            assert int(c.get(tgt, 0)) == k, (
                f"{st}/{tgt}: collapsed says {c.get(tgt)}, the full categories "
                f"sum to {k}"
            )
            n_checked += 1
            checked_strata.add(st)
    # A level-by-level check that skipped every level would pass vacuously, which
    # is the failure this suite has hit before. The bound is the strata where
    # both views are published -- every one of them must contribute.
    missed = sorted(set(overlap) - checked_strata)
    assert not missed and n_checked, (
        f"verified {n_checked} collapsed level(s); strata published in both "
        f"views but never checked: {missed or 'none, but nothing was checked'} "
        f"-- the check is skipping, not passing"
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
