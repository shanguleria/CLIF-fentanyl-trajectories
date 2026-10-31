"""The federated-pooling exports must be poolable, exactly.

A median cannot be pooled across sites; a mean can, from n and the two sums. The
point of carrying sum and sum_sq rather than only mean and sd is that the pooled
figures are then exact rather than an approximation that assumes equal variances.
These checks run only when the scripts have produced the files.

Run standalone:  .venv/bin/python tests/test_pooling.py
"""
from __future__ import annotations

import json
import re
import math
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "output" / "final_no_phi"
# out_final is subdivided by the script that produced each file (2026-09-08).
CONT = OUT / "02_descriptive" / "pooling_continuous.csv"
# Table 1's categorical half MOVED to 07 on 2026-10-28 (protocol 0.17.0), when 07
# became the sole generator. 02 keeps pooling_continuous.csv for its by_window
# rows; the baseline rows and all of the categorical rows are 07's now.
CAT = OUT / "07_outcome_model" / "table1_categorical.csv"
T1_CONT = OUT / "07_outcome_model" / "table1_continuous.csv"
# The descriptive and landmark scripts shared one pooling contract and both were
# held to it, the second
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
        print(f"  SKIP  {f.name} absent -- run the script that writes it")
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


def test_the_table1_scopes_are_nested_not_a_partition():
    """Table 1's columns stopped being a partition on 2026-10-28.

    Until protocol 0.17.0 the strata were Landmark and its COMPLEMENT -- disjoint
    and exhaustive -- so pooling them as separate sites reproduced the overall row
    exactly, and this test checked that arithmetic. The scopes are now NESTED
    (overall > landmark_increase > model), so they do not sum to overall and any
    consumer that adds them is wrong. The invariant that replaces it: every scope
    is a subset of overall, so its n never exceeds overall's and its quantiles
    stay inside overall's observed range is NOT guaranteed -- a subset can be
    shifted -- but the count relation must hold.
    """
    f = OUT / "07_outcome_model" / "table1_continuous.csv"
    d = _skip_if_absent(f)
    if d is None:
        return
    checked = 0
    for var, g in d.groupby("variable"):
        rows = {r.stratum: r for r in g.itertuples()}
        o = rows.get("overall")
        if o is None:
            continue          # a model-only row has no overall by design
        for name, r in rows.items():
            if name == "overall":
                continue
            assert r.n <= o.n, (
                f"{var}/{name} has n={r.n} against overall n={o.n}; the scopes "
                f"are nested, so a subset cannot be larger than its superset"
            )
            checked += 1
    assert checked >= 3, f"only {checked} nested-scope comparisons were made"


def test_the_table1_scopes_do_not_sum_to_overall():
    """The converse, stated as a test so nobody re-adds the old assumption.

    If a future change made the scopes disjoint again this would fail loudly,
    which is the point: summing nested scopes double-counts every admission in
    the narrower one.
    """
    f = OUT / "07_outcome_model" / "table1_continuous.csv"
    d = _skip_if_absent(f)
    if d is None:
        return
    g = d[d.variable == "Age, years"]
    rows = {r.stratum: r for r in g.itertuples()}
    if "overall" not in rows or len(rows) < 3:
        return
    parts = sum(r.n for k, r in rows.items() if k != "overall")
    assert parts != rows["overall"].n, (
        "the scopes summed to overall, which means they are disjoint again -- "
        "if that is intended, restore the partition arithmetic test"
    )


def test_small_cells_are_suppressed_not_published():
    """A mean over n = 1 is that patient's value. Anything below the site's
    small_cell_min_den must carry no statistics at all."""
    for f in ALL_CONT:
        d = _skip_if_absent(f)
        if d is not None:
            _check_suppression(d, f.name)


def _check_suppression(d, name):
    cfg = json.loads((REPO / "config" / "protocol.json").read_text())
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


# ------------------------------------------------- the spec and the writers
SPEC = REPO / "config" / "pooling_spec.json"


def test_the_pooling_spec_claims_every_csv_the_pipeline_writes():
    """A path the spec does not claim is refused at the CENTRE, after the upload.

    pool_sites.R hard-stops on an unclaimed file -- deliberately, since a silent
    partial pool is worse -- but it does so on someone else's machine, days
    later, with a bundle that has to be rebuilt. This catches the same divergence
    here, on the site's own run. It earned its place on 2026-10-27, when four
    method checks moved into 05_titration/diagnostics/ and their spec keys had to
    move with them.

    Keyed on what is ON DISK rather than on the OWNED lists, because OWNED is R
    and parsing it from Python would test the parser. That makes this a check of
    a completed run; it skips when the tree is empty.
    """
    claimed = set(json.loads(SPEC.read_text())["files"])
    if not OUT.is_dir():
        print("  SKIP  no output/final_no_phi -- run the pipeline")
        return
    on_disk = {f.relative_to(OUT).as_posix() for f in OUT.rglob("*.csv")}
    # The cohort build's cohort_construction.csv and the tabled scripts' leftovers are not pooled.
    labels = set(json.loads(SPEC.read_text()).get("not_poolable_columns", []))
    unclaimed = sorted(on_disk - claimed)
    assert not unclaimed, (
        f"config/pooling_spec.json claims no entry for: {unclaimed}. Either add "
        f"one, or if the file is deliberately not pooled say so there -- "
        f"pool_sites.R will refuse the whole bundle otherwise.")
    # The converse -- a spec entry with no file -- is checked only inside phase
    # directories that produced SOMETHING. A phase mid-run has just cleared its
    # own outputs, which is indistinguishable from a relocation; asserting
    # globally made this test fail whenever it ran beside a live pipeline. Scoped
    # this way it still catches a file moved WITHIN a phase, which is the bug.
    produced = {f.parent.relative_to(OUT).as_posix() for f in OUT.rglob("*.csv")}
    missing = sorted(k for k in claimed
                     if str(Path(k).parent) in produced and not (OUT / k).is_file())
    assert not missing, (
        f"the spec claims files their own phase did not write: {missing}. A "
        f"relocated output leaves the spec pointing at the old path.")
    assert labels, "not_poolable_columns must be declared"


def test_every_sum_bins_file_declares_a_bin_key_the_pooler_can_group_on():
    """pool_bins() sums bins WITHIN the keys other than `bin`, and it used to
    hardcode `variable` for that grouping. 07's exposure_distribution.csv is keyed
    on `exposure` instead, so `tapply(out$n, out$variable, ...)` got a NULL second
    argument and pool_sites.R DIED -- after writing 36 of the files, which is the
    worst place to stop. Generalised to the declared keys; this is the check that
    the declaration and the pooler still agree."""
    spec = json.loads(SPEC.read_text())["files"]
    for rel, entry in spec.items():
        if entry.get("method") != "sum_bins":
            continue
        keys = entry.get("keys", [])
        assert "bin" in keys, (
            f"{rel} declares sum_bins but has no `bin` key; pool_bins() cannot "
            f"tell which key indexes the bins")
        assert len(set(keys) - {"bin"}) >= 1, (
            f"{rel} is keyed only on `bin`, so there is no group to total within")
        for col in entry.get("bin_edges", []):
            assert col, f"{rel} declares an empty bin_edges entry"
        f = OUT / rel
        if not f.is_file():
            continue
        cols = set(pd.read_csv(f, nrows=1).columns)
        missing = [k for k in keys if k not in cols]
        assert not missing, (
            f"{rel} is keyed on {missing} in the spec but does not ship those "
            f"column(s) -- the pooler would group on nothing")


def test_the_cohort_column_is_a_key_wherever_it_is_emitted():
    """`landmark` is a SUBSET of `all`. If cohort is emitted but not keyed, every
    sum_counts method adds the two rows and double-counts every landmark event --
    silently, in the plausible direction, with a pooled numerator that still
    looks like a count."""
    spec = json.loads(SPEC.read_text())["files"]
    for rel, entry in spec.items():
        f = OUT / rel
        if not f.is_file():
            continue
        cols = set(pd.read_csv(f, nrows=1).columns)
        if "cohort" not in cols:
            continue
        assert "cohort" in entry.get("keys", []), (
            f"{rel} ships a `cohort` column but the spec does not key on it; "
            f"pooling would add `all` and `landmark` together")
        d = pd.read_csv(f)
        # The value set is `all` plus any number of landmark strata. 05's four
        # files ship `landmark` (the protocol T); course_summary.csv sweeps T and
        # ships landmark_24h / landmark_48h / landmark_72h. Both shapes are a
        # NESTED FAMILY, which is the property that makes keying mandatory --
        # the count of strata is not.
        seen = set(d["cohort"])
        bad = {c for c in seen
               if c != "all" and not re.fullmatch(r"landmark(_\d+h)?", str(c))}
        assert not bad, (
            f"{rel} carries an unexpected cohort value: {sorted(bad)}. Allowed: "
            f"`all`, `landmark`, or `landmark_<N>h`.")
        assert "all" in seen and len(seen) > 1, (
            f"{rel} must carry `all` and at least one landmark stratum, or the "
            f"comparison the cohort column exists for is not available")


HISTS = [OUT / "02_descriptive" / "pooling_histograms.csv",
         OUT / "05_titration" / "pooling_histograms.csv"]


def test_every_declared_histogram_is_emitted_exactly_once():
    """A histogram declared in the config and written by nobody is the
    config-integrity failure CLAUDE.md names first; one written by two scripts
    would be double-counted the moment a centre concatenates the bundles."""
    cfg = json.loads((REPO / "config" / "protocol.json").read_text())
    declared = set(cfg["pooling"]["histograms"]["variables"])
    emitted: list[str] = []
    for f in HISTS:
        d = _skip_if_absent(f)
        if d is None:
            continue
        emitted += sorted(set(d["variable"]))
    assert emitted, "no histograms were emitted at all"
    assert set(emitted) == declared, (
        f"declared but not emitted: {sorted(declared - set(emitted))}; "
        f"emitted but not declared: {sorted(set(emitted) - declared)}"
    )
    dupes = [v for v in set(emitted) if emitted.count(v) > 1]
    assert not dupes, f"emitted by more than one script: {sorted(dupes)}"


def test_each_histogram_accounts_for_every_value():
    """The bins must partition the data: a value in no bin is silently dropped
    from the pooled distribution, and the pooled median moves with it.

    n_total is the count of NON-MISSING values, which is why it differs between
    variables measured at different rates (fentanyl is charted far more densely
    than lactate). A pooled histogram therefore describes the distribution
    among OBSERVED values -- pair it with missingness.csv, not with the window
    count.
    """
    for f in HISTS:
        d = _skip_if_absent(f)
        if d is None:
            continue
        for v, g in d.groupby("variable"):
            assert g["n"].sum() == g["n_total"].iloc[0], (
                f"{f.name}/{v}: bins sum to {g['n'].sum()} against n_total "
                f"{g['n_total'].iloc[0]} -- the bins do not partition the data"
            )
            assert g["n_total"].nunique() == 1, f"{v}: n_total is not constant"
            assert (g["n"] >= 0).all(), f"{v}: a negative bin count"
            # Contiguous and ordered, so a pooler can read a quantile off the
            # cumulative sum without re-sorting or checking for gaps.
            b = g.sort_values("bin")
            assert list(b["bin"]) == list(range(1, len(b) + 1)), f"{v}: bins not 1..n"
            lo, hi = b["bin_lo"].tolist(), b["bin_hi"].tolist()
            # The first bin may be a degenerate [0, 0] zero bin, whose hi equals
            # the next bin's lo -- that is the point of it.
            for i in range(len(b) - 1):
                assert hi[i] == lo[i + 1], (
                    f"{v}: gap or overlap between bin {i+1} and {i+2}: "
                    f"{hi[i]} then {lo[i+1]}"
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
    cmap = json.loads((REPO / "config" / "protocol.json").read_text())
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
    # No "strata sum to overall" identity: since protocol 0.17.0 the scopes are
    # NESTED (see test_the_table1_scopes_are_nested_not_a_partition), so a level
    # at overall is not the sum of the same level across scopes. Only the
    # within-stratum identities above hold.
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
    published table reported the wrong share of the landmark cohort for a race
    level. Totals still matched, so only a level-by-level check catches it.
    """
    d = _skip_if_absent(CAT)
    if d is None:
        return
    full = d[d.variable == "Race"]
    coll = d[d.variable == "Race (collapsed)"]
    assert len(coll), "the collapsed race rows are missing"
    cmap = json.loads((REPO / "config" / "protocol.json").read_text())
    cmap = cmap["time_invariant"]["race"]["reporting_collapse"]
    named = {k: v for k, v in cmap.items() if not k.startswith("_")}
    # The two views no longer span the same strata: protocol.json
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
        # `overall` is fully unsuppressed, so without this the check
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
            if tgt in dirty or pd.isna(c.get(tgt, 0)):
                continue          # suppressed in the collapsed view (primary or complementary)
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
    # Suppressed cells carry NaN, so published sums are lower bounds on the
    # stratum total and equality only holds when nothing was suppressed.
    den = coll["denominator"].iloc[0]
    assert full["denominator"].iloc[0] == den, "the two views disagree on the denominator"
    assert coll["n"].sum() <= den and full["n"].sum() <= den
    if not (full["n"].isna().any() or coll["n"].isna().any()):
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
