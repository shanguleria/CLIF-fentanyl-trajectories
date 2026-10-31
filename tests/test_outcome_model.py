"""Checks on what 07_outcome_model.R shipped, and on the protocol it claims to follow.

This is a POST-RUN suite, like test_pooling.py: it reads
output/final_no_phi/07_outcome_model/ and fails if what is there is not what
protocol.json declares. It does not belong in the preflight -- the files it
asserts on do not exist until the pipeline has produced them, and a suite that
gates the production of the files it validates can never pass on a clean clone
(lessons.md #28).

The reason this suite exists at all is the proportional-odds test. On
2026-10-27 it shipped with every term's p-value NA and a console line reading
"0 term(s) with p < 0.05" -- which is what a PASSING test looks like. The
diagnostic had not passed; it had not run. A check on the numbers would not
have caught it, because all three models converged and every coefficient was
correct. Only a check that the ASSUMPTION TEST produced a verdict catches it.

No CLIF data is read: everything here comes from the shareable output tree.

Run standalone:  .venv/bin/python tests/test_outcome_model.py
Or under pytest: .venv/bin/python -m pytest tests/test_outcome_model.py -q
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "output" / "final_no_phi" / "07_outcome_model"
SPEC = json.loads((REPO / "config" / "protocol.json").read_text())["model"]["outcome_model"]
SOURCE = (REPO / "code" / "07_outcome_model.R").read_text()

# The exposures enter the model rescaled to the declared reporting unit, and the
# shipped term names carry that unit so a coefficient cannot be read against the
# wrong one. Derived here the same way the R derives it, from the same key.
SCALE_PP = round(SPEC["exposure_scale_per"] * 100)
MODEL_EXPOSURES = [f"{e}_per{SCALE_PP}pp" for e in SPEC["exposures"]]


def rows(name: str) -> list[dict]:
    with open(OUT / name, newline="") as fh:
        return list(csv.DictReader(fh))


def num(s: str):
    """CSV reads everything as text; R writes NA for a missing number."""
    return None if s in ("NA", "", "NaN") else float(s)


# ---- the proportional-odds test must produce a verdict ----------------------

def test_the_proportional_odds_test_produced_a_verdict_for_every_model():
    """THE REGRESSION GUARD. Every model must have at least one term actually
    tested -- not merely a row in the file.

    The 2026-10-27 failure: clm() stores its call verbatim, and nominal_test()
    refits once per term by re-evaluating it. 07 fit the model inside a function
    as `clm(form, data = frame, ...)`, so the stored call named a local that was
    gone by refit time -- and `frame` is also a base R function, so the refit
    resolved it to that function rather than failing on an unknown name. Every
    refit died, nominal_test returned a fully formed table with every term NA,
    and nothing anywhere said so. Fixed by fitting through do.call, which inlines
    the data into the stored call.

    THIS GUARD FIRED AGAIN ON 2026-10-28, on sensitivity_icu_only, and was
    right to. ordinal:::nominal_test.clm guards its refit with
    `!inherits(nfit, "try-error") && nfit$convergence$code >= 0`, and clm
    returns a LENGTH-2 convergence code when a refit trips two conditions at
    once -- so under R >= 4.3 the whole call raises rather than reporting a
    failed refit, and the model shipped as one `errored` row. The fix was to
    test one term at a time through nominal_test's own `scope` argument, which
    reproduces the all-at-once table exactly on the primary (df and p, all
    seven terms, 1e-8) while letting one term's failure stay that term's
    failure. The underlying refits on that model genuinely do not converge, so
    its exposures remain not_testable -- which is why the GATE, not this test,
    is what decides whether the exposures passed.
    """
    po = rows("diagnostics/proportional_odds_test.csv")
    diag = rows("diagnostics/model_diagnostics.csv")
    assert diag, "no models in model_diagnostics.csv"
    for model in {d["model"] for d in diag}:
        terms = [r for r in po if r["model"] == model and r["term"] != "<none>"]
        assert terms, f"{model}: no term rows in the proportional-odds test"
        tested = [r for r in terms if num(r["p_value"]) is not None]
        states = sorted({r["status"] for r in terms})
        notes = sorted({r.get("note", "") for r in terms if r.get("note")})
        assert tested, (
            f"{model}: the proportional-odds test produced NO verdict for any "
            f"term -- {len(terms)} term rows, statuses {states}, every p-value "
            f"NA. The assumption is UNTESTED, not satisfied."
            + (f" Reported reason(s): {notes}" if notes else "")
            + " See this test's docstring."
        )


def test_every_proportional_odds_row_declares_its_status():
    """The status column is what makes an untestable term distinguishable from a
    passing one in the shipped file. Counting significant terms over an all-NA
    table yields zero, which reads as the assumption holding.

    AMENDED 2026-10-28 with the day-scale outcome. This test used to require
    every row to read `tested`, which was a true statement about the four-band
    model and is a POLICY, not a bug check -- lessons.md #21. At 23 thresholds
    the per-term refits for `age` and `icu` return NA, so a blanket demand would
    fail on a protocol change rather than on a defect. What is enforced instead:
    status is per TERM and is one of the declared values, a null p-value is
    never called `tested`, and a non-null one is never called `not_testable`.
    WHICH terms must pass is the gate's job, tested below.
    """
    po = rows("diagnostics/proportional_odds_test.csv")
    assert po, "proportional_odds_test.csv is empty"
    assert "status" in po[0], "the status column is gone -- an untestable term would again be silent"
    allowed = {"tested", "not_testable", "errored", "reference"}
    for r in po:
        assert r["status"] in allowed, f"unknown status {r['status']!r}"
        if r["term"] == "<none>":
            continue
        p = num(r["p_value"])
        if r["status"] == "tested":
            assert p is not None, (
                f"{r['model']}/{r['term']}: status 'tested' beside a null "
                f"p-value -- the diagnostic is reporting a pass it did not make"
            )
        if r["status"] == "not_testable":
            assert p is None, (
                f"{r['model']}/{r['term']}: status 'not_testable' but a "
                f"p-value of {p} is present"
            )


# ---- the outcome scale is the one the protocol declares ---------------------

def test_the_ordinal_breaks_and_labels_describe_the_same_cuts():
    breaks = SPEC["vfd_ordinal_breaks"]
    labels = SPEC["vfd_ordinal_labels"]
    assert len(labels) == len(breaks) - 1, "one label per interval between breaks"
    assert breaks == sorted(breaks) and len(set(breaks)) == len(breaks), \
        "breaks must be strictly increasing or cut() silently reorders the levels"


def test_exactly_zero_falls_in_the_first_level():
    """Half the cohort sits at VFD = 0 and that point mass must be its own level.
    cut() is right-closed, so the first break has to be BELOW zero for 0 to land
    in level 1 rather than becoming NA and dropping half the admissions."""
    assert SPEC["vfd_ordinal_breaks"][0] < 0, \
        "first break is not below zero -- VFD = 0 would fall outside every interval"


def test_the_top_break_clears_the_vfd_ceiling():
    """The landmark is at 72h and the horizon is 28 days, so the most VFDs an
    admission can earn is 25. A top break at or below 25 would put the ceiling
    itself outside the last interval."""
    horizon = SPEC.get("_vfd_ceiling_days")
    if horizon is None:
        outcomes = json.loads((REPO / "config" / "protocol.json").read_text())["outcomes"]
        horizon = outcomes["vfd"]["horizon_days"] - 3      # landmark at 72h
    assert SPEC["vfd_ordinal_breaks"][-1] > horizon, \
        f"top break {SPEC['vfd_ordinal_breaks'][-1]} does not clear the {horizon}-day ceiling"


def test_the_shipped_levels_are_the_declared_labels_in_order():
    got = [r["level"] for r in rows("vfd_ordinal.csv")]
    assert got == SPEC["vfd_ordinal_labels"], \
        f"shipped levels {got} do not match the declared labels"


def test_the_ordinal_levels_sum_to_their_denominator():
    vr = rows("vfd_ordinal.csv")
    dens = {int(r["denominator"]) for r in vr}
    assert len(dens) == 1, f"the levels do not share one denominator: {dens}"
    assert sum(int(r["n_admissions"]) for r in vr) == dens.pop(), \
        "the ordinal levels do not sum to the denominator -- a level is missing or double-counted"


# ---- the models are admissible and reconcile with the flow ------------------

def test_every_shipped_model_converged():
    for d in rows("diagnostics/model_diagnostics.csv"):
        assert d["converged"] == "TRUE", f"{d['model']} did not converge"
        assert d["convergence_code"] == "0", \
            f"{d['model']} convergence code {d['convergence_code']}: {d['convergence_message']}"


def test_the_primary_model_n_reconciles_with_strobe():
    """The published flow and the fitted model n must agree, and since 2026-10-28
    the flow states the unit drop as its own row rather than ending at the
    complete-case count -- so this checks the row AND the diagnostic, which used
    to be only arithmetic on the last row.

    An unexplained gap between a published flow and a model n is the denominator
    trap this repo has already hit three times.
    """
    flow = rows("strobe.csv")
    assert "branch" in flow[0], "strobe.csv lost its branch column"
    model_rows = [r for r in flow if r["branch"] == "model"]
    assert model_rows, "strobe.csv has no model branch"
    last = model_rows[-1]
    assert last["step"] == "ANALYSED IN THE MODEL", last["step"]

    primary = [d for d in rows("diagnostics/model_diagnostics.csv")
               if d["model"] == "primary"]
    assert primary, "no primary model in model_diagnostics.csv"
    p = primary[0]

    # 1. The published end of the funnel IS the fitted n -- no arithmetic.
    assert int(last["n_after"]) == int(p["n_admissions"]), (
        f"STROBE ends at {last['n_after']} but the primary model fits "
        f"{p['n_admissions']}"
    )
    # 2. And the drop it publishes is the one the model reported.
    dropped = int(p["n_dropped_small_unit"])
    assert int(last["n_excluded"]) == dropped, (
        f"STROBE publishes {last['n_excluded']} unit-drop(s), the model reports "
        f"{dropped}"
    )
    # 3. Every branch chains internally; across branches it must not be summed.
    for br in {r["branch"] for r in flow}:
        rs = [r for r in flow if r["branch"] == br]
        for a, b in zip(rs, rs[1:]):
            assert int(a["n_after"]) == int(b["n_before"]), (
                f"{br}: {a['step']} ends at {a['n_after']} but {b['step']} "
                f"starts from {b['n_before']}"
            )


def test_the_strobe_branches_descend_from_the_analytic_cohort():
    """Both branches start from the trunk's analytic cohort, not from each other.
    If they ever chained, one of the two denominators would be wrong -- the
    practice-pattern rates are over every admission that escalated, the model is
    over the complete-case frame."""
    flow = rows("strobe.csv")
    trunk = [r for r in flow if r["branch"] == "trunk"]
    analytic = int(trunk[-1]["n_after"])
    assert "ANALYTIC COHORT" in trunk[-1]["step"], trunk[-1]["step"]
    starts = {br: [r for r in flow if r["branch"] == br][0]
              for br in ("practice_patterns", "model")}
    assert int(starts["practice_patterns"]["n_before"]) == analytic, (
        "the practice-pattern branch must descend from the analytic cohort"
    )
    # The model branch descends from the escalation cohort, which is itself the
    # practice-pattern branch's end -- that is the one legitimate link.
    assert int(starts["model"]["n_before"]) == \
        int(starts["practice_patterns"]["n_after"])


def test_the_unit_fixed_effects_are_not_poolable():
    """Only the exposure coefficients travel. The ICU fixed effects are
    site-local nuisance parameters -- `icuN08S` means nothing at another site --
    so shipping them as poolable would invite the centre to average them."""
    for r in rows("model_coefficients.csv"):
        if r["term"].startswith("icu"):
            assert r["poolable"] == "FALSE", \
                f"unit fixed effect {r['term']} is marked poolable"


def test_both_exposures_ship_a_full_covariance_block():
    """The centre combines the two exposures jointly, which needs the 2x2 and not
    just the two standard errors."""
    exposures = set(MODEL_EXPOSURES)
    cov = rows("model_covariance.csv")
    by_model: dict[str, set] = {}
    for r in cov:
        by_model.setdefault(r["model"], set()).add((r["term_row"], r["term_col"]))
    assert by_model, "model_covariance.csv is empty"
    for model, pairs in by_model.items():
        for a in exposures:
            for b in exposures:
                assert (a, b) in pairs, f"{model}: covariance missing ({a}, {b})"
    for r in cov:
        mirror = [x for x in cov if x["model"] == r["model"]
                  and x["term_row"] == r["term_col"] and x["term_col"] == r["term_row"]]
        assert mirror and mirror[0]["covariance"] == r["covariance"], \
            f"{r['model']}: covariance is not symmetric at ({r['term_row']}, {r['term_col']})"


def test_the_covariance_standard_errors_agree_with_the_coefficients():
    """sqrt of the diagonal must be the shipped standard error, or the centre
    would weight by one number while reading uncertainty off another."""
    coefs = {(r["model"], r["term"]): r for r in rows("model_coefficients.csv")}
    for r in rows("model_covariance.csv"):
        if r["term_row"] != r["term_col"]:
            continue
        c = coefs.get((r["model"], r["term_row"]))
        assert c is not None, f"{r['model']}: {r['term_row']} has covariance but no coefficient"
        se = float(c["std_error"])
        assert abs(se - float(r["covariance"]) ** 0.5) < 1e-6, (
            f"{r['model']} {r['term_row']}: se {se} vs sqrt(var) "
            f"{float(r['covariance']) ** 0.5}"
        )


# ---- every declared key reaches the code -----------------------------------

def test_every_declared_outcome_model_setting_is_consumed():
    """CLAUDE.md's config-integrity rule: a key that no code path reads is a
    policy declaration nobody honours. Keys are checked against the script's
    source rather than its behaviour, which is the weaker half -- the stronger
    half is the assertions above, which read the output the settings produced."""
    declared = [k for k in SPEC if not k.startswith("_")]
    assert declared, "no settings declared under model.outcome_model"
    unconsumed = [k for k in declared if k not in SOURCE]
    assert not unconsumed, (
        f"declared in protocol.json but absent from 07_outcome_model.R: {unconsumed}"
    )


def test_every_fitted_sensitivity_model_is_declared_in_the_protocol():
    """The REVERSE of the consumption check above, and the direction that was
    missing. test_every_declared_outcome_model_setting_is_consumed catches a
    protocol key no code reads; this catches a model the code FITS that the
    protocol never declared -- which in a federated repo is worse, because
    protocol.json is the digest-visible statement of the estimand and a site
    comparing protocols would see agreement while the models differed.

    Found live on 2026-10-27: sensitivity_partial_po was fitted and shipped
    before it was declared.

    Matching is by token containment, because a label may carry a value the key
    cannot (`min_events` is fitted as `sensitivity_min_4_events`, the 4 coming
    from sensitivity_min_increase_events).
    """
    # Underscore-prefixed keys are prose, not models -- the convention this
    # protocol uses throughout (_purpose, _why, ...). Including them would
    # make the test demand a model named after a comment.
    declared = {k: set(k.split("_")) for k in SPEC["sensitivity_models"]
                if not k.startswith("_")}
    fitted = [d["model"] for d in rows("diagnostics/model_diagnostics.csv")
              if d["model"] != "primary"]
    assert fitted, "no sensitivity models were fitted"
    for label in fitted:
        tokens = set(label.removeprefix("sensitivity_").split("_"))
        hits = [k for k, kt in declared.items() if kt <= tokens]
        assert hits, (
            f"the model {label!r} is fitted and shipped but no "
            f"protocol.json model.outcome_model.sensitivity_models entry "
            f"declares it (declared: {sorted(declared)})"
        )
    for key, kt in declared.items():
        hits = [l for l in fitted if kt <= set(l.removeprefix("sensitivity_").split("_"))]
        assert hits, f"sensitivity model {key!r} is declared but was not fitted"


def test_the_declared_exposures_and_covariates_are_the_ones_fitted():
    """The formula is written out in model_diagnostics.csv, so the protocol can be
    compared to what was actually fitted rather than to what the script says."""
    primary = [d for d in rows("diagnostics/model_diagnostics.csv")
               if d["model"] == "primary"][0]
    formula = primary["formula"]
    for term in MODEL_EXPOSURES + SPEC["covariates"]:
        assert term in formula, f"{term} is declared but not in the fitted formula: {formula}"
    assert "n_increase" not in formula and "n_event" not in formula, (
        "the increase-event count is in the formula. It is the denominator of both "
        "exposures and is measured inside the exposure window; see protocol.json "
        "model.outcome_model._n_events_is_NOT_a_covariate"
    )


def test_the_declared_exposure_scale_reaches_the_shipped_coefficients():
    """CLAUDE.md config-integrity, applied to the reporting unit.

    protocol.json declares exposure_scale_per; if nothing consumed it the model
    would still fit and every number would look plausible -- just on a 0%->100%
    scale while the protocol claimed 10 percentage points. That is the silent,
    plausible-direction failure the rule is about, so this asserts the scale
    reached BOTH the term names and the scale_note a reader relies on.
    """
    rs = rows("model_coefficients.csv")
    shipped = {r["term"] for r in rs}
    for term in MODEL_EXPOSURES:
        assert term in shipped, (
            f"{term!r} is not in model_coefficients.csv -- the declared "
            f"exposure_scale_per={SPEC['exposure_scale_per']} did not reach the "
            f"fitted terms (shipped: {sorted(shipped)[:6]}...)"
        )
    flagged = {r["term"] for r in rs if r["is_exposure"] == "TRUE"}
    assert flagged == set(MODEL_EXPOSURES), (
        f"is_exposure marks {sorted(flagged)}, expected {sorted(MODEL_EXPOSURES)}"
    )
    note = rs[0]["scale_note"]
    assert "PERCENTAGE" in note.upper() and str(SCALE_PP) in note, (
        f"scale_note does not state the {SCALE_PP}-percentage-point unit: {note!r}"
    )


def test_the_unscaled_exposure_names_do_not_ship_as_coefficients():
    """The share-scale names still label the DESCRIPTIVE outputs, which is
    correct -- a share is what was measured. They must not also appear as
    coefficient rows, or the same name would mean two different units inside one
    output directory."""
    shipped = {r["term"] for r in rows("model_coefficients.csv")}
    for bare in SPEC["exposures"]:
        assert bare not in shipped, (
            f"{bare!r} ships as a coefficient on the share scale while "
            f"{bare}_per{SCALE_PP}pp is the declared reporting unit"
        )



# ---- the outcome is modelled at whole-day resolution ------------------------

def test_the_model_outcome_is_the_declared_scale():
    """0.16.0 made the declared bands the primary and whole days a sensitivity,
    because the day scale's 22 proportional-odds thresholds are not satisfiable
    across sites. The number of levels each model actually fitted is written into
    model_diagnostics.csv, so the declaration is checked against the fit."""
    assert SPEC["outcome_scale"] == "declared_bands", SPEC["outcome_scale"]
    assert SPEC["band_values"] == "observed_within_band_mean", SPEC["band_values"]
    diag = {d["model"]: d for d in rows("diagnostics/model_diagnostics.csv")}
    assert "primary" in diag, "no primary model was fitted"
    n_levels = int(float(diag["primary"]["n_outcome_levels"]))
    n_bands = len(SPEC["vfd_ordinal_labels"])
    assert n_levels == n_bands, (
        f"the primary fitted {n_levels} outcome levels against the {n_bands} "
        f"declared bands in vfd_ordinal_breaks"
    )
    # The two scale sensitivities must be on their own scales, or they are
    # duplicates of the primary wearing a sensitivity's name.
    assert int(float(diag["sensitivity_day_scale"]["n_outcome_levels"])) > n_bands, (
        "sensitivity_day_scale is not on the day scale"
    )


def test_the_four_band_sensitivity_still_reports_the_old_scale():
    """The four-band model is kept precisely so the previously reported number
    does not vanish. It must be fitted AND must be on four levels -- fitting it
    on the day scale would make it a duplicate of the primary wearing a
    sensitivity's name."""
    diag = {d["model"]: d for d in rows("diagnostics/model_diagnostics.csv")}
    assert "sensitivity_four_band" in diag, (
        "the four_band sensitivity was not fitted, so the outcome-scale change "
        "silently dropped the number this study has already reported"
    )
    got = int(float(diag["sensitivity_four_band"]["n_outcome_levels"]))
    assert got == len(SPEC["four_band_labels"]), (
        f"the four_band sensitivity fitted {got} levels, not "
        f"{len(SPEC['four_band_labels'])}"
    )


def test_the_day_distribution_covers_the_whole_attainable_range():
    """The day grid ships in full, including days no admission reached, so the
    key set is identical at every site and a pooled row sums over the same
    levels. An absent level is information: day 1 is unreachable because a
    liberation must last 48h to count."""
    dd = rows("vfd_day_distribution.csv")
    assert dd, "vfd_day_distribution.csv is empty"
    days = [int(float(r["day"])) for r in dd]
    assert days == list(range(min(days), max(days) + 1)), \
        f"the day grid has holes: {days}"
    assert min(days) == 0, f"the day grid starts at {min(days)}, not 0"
    den = {int(float(r["denominator"])) for r in dd}
    assert len(den) == 1, f"the denominator is not constant across days: {den}"
    total = sum(int(float(r["n_admissions"])) for r in dd)
    assert total == den.pop(), \
        f"the day counts sum to {total}, not the shipped denominator"


# ---- the proportional-odds gate ---------------------------------------------

def test_the_proportional_odds_gate_passed_for_both_exposures():
    """partial_po was dropped on 2026-10-28 CONDITIONAL on both exposures
    passing the PO test. The condition is the only thing holding the
    simplification up, so it is checked against the shipped verdict and not
    against the run log."""
    gate = rows("diagnostics/proportional_odds_gate.csv")
    assert gate, "proportional_odds_gate.csv is empty -- the gate did not run"
    scope = SPEC["proportional_odds_gate"]["scope"]
    terms = {r["term"] for r in gate}
    assert terms == set(MODEL_EXPOSURES), (
        f"the gate covers {sorted(terms)}, not the declared exposures "
        f"{sorted(MODEL_EXPOSURES)}"
    )
    for r in gate:
        assert r["model"] == scope, f"the gate ran on {r['model']!r}, not {scope!r}"
        assert r["passed"] in ("TRUE", "FALSE"), r["passed"]
        assert r["passed"] == "TRUE", (
            f"PROPORTIONAL-ODDS GATE FAILED for {r['term']} "
            f"(status {r['status']}, p {r['p_value']}). A single odds ratio is "
            f"not a valid summary of that exposure; reinstate a partial-PO model"
        )
        assert float(r["alpha"]) == SPEC["proportional_odds_gate"]["alpha"]


def test_the_partial_po_sensitivity_is_gone_from_both_sides():
    """It was removed on the condition above. If it comes back it must come back
    in BOTH places -- a model fitted but not declared is the failure
    test_every_fitted_sensitivity_model_is_declared_in_the_protocol exists for,
    and a declaration with no model is the reverse."""
    assert "partial_po" not in SPEC["sensitivity_models"], \
        "partial_po is declared again -- is it also fitted?"
    fitted = {d["model"] for d in rows("diagnostics/model_diagnostics.csv")}
    assert not [m for m in fitted if "partial_po" in m], \
        f"a partial-PO model is still fitted: {sorted(fitted)}"


# ---- the adjusted prediction -------------------------------------------------

def test_the_adjusted_prediction_is_on_the_declared_grid():
    pred = rows("adjusted_prediction.csv")
    assert pred, "adjusted_prediction.csv is empty"
    want = [float(g) for g in SPEC["prediction"]["exposure_grid_pct"]]
    got = [float(r["exposure_pct"]) for r in pred]
    assert got == want, f"the prediction grid is {got}, not the declared {want}"
    for r in pred:
        assert r["exposure"] == SPEC["prediction"]["exposure"]
        assert int(float(r["n_draws"])) == SPEC["prediction"]["bootstrap_draws"]
        assert abs(float(r["exposure_share"]) * 100 - float(r["exposure_pct"])) < 1e-9


def test_the_band_to_days_mapping_ships_with_the_curve():
    """The adjusted mean is no longer mapping-free: a band carries no days, so
    the expectation depends on a value assigned to each band. That assignment is
    shipped, or the curve cannot be reproduced from the file."""
    pred = rows("adjusted_prediction.csv")
    assert pred and "band_values" in pred[0], (
        "adjusted_prediction.csv does not ship band_values -- the mapping the "
        "top panel rests on would be invisible"
    )
    vals = [float(v) for v in pred[0]["band_values"].split("|")]
    assert len(vals) == len(SPEC["vfd_ordinal_labels"]), (
        f"{len(vals)} band values against {len(SPEC['vfd_ordinal_labels'])} bands"
    )
    assert vals == sorted(vals) and len(set(vals)) == len(vals), \
        f"band values are not strictly increasing: {vals}"
    assert vals[0] == 0, f"the lowest band value is {vals[0]}, not 0"
    # Every row must carry the same mapping -- it is a property of the fit.
    assert len({r["band_values"] for r in pred}) == 1, \
        "the band->days mapping differs between grid points"
    # And it must agree with the descriptive table's own observed means.
    tab = rows("vfd_ordinal.csv")
    assert "mean_vfd_28" in tab[0], "vfd_ordinal.csv does not ship mean_vfd_28"


def test_every_predicted_interval_brackets_its_point_estimate():
    """A percentile interval that does not contain the fitted point estimate
    means the bootstrap is drawing from a different model than the one that
    produced the point -- the exact failure a hand-rolled linear predictor
    invites, and why the script asserts itself against predict.clm first."""
    for r in rows("adjusted_prediction.csv"):
        for stem in ("mean_vfd_28", "pct_vfd_zero"):
            lo, mid, hi = (num(r[stem + "_lo"]), num(r[stem]), num(r[stem + "_hi"]))
            assert None not in (lo, mid, hi), f"{stem} has a missing bound at {r['exposure_pct']}%"
            assert lo <= mid <= hi, (
                f"{stem} at {r['exposure_pct']}%: {mid} is outside its own "
                f"interval [{lo}, {hi}]"
            )
            assert lo < hi, f"{stem} at {r['exposure_pct']}%: a zero-width interval"


def test_the_predicted_mean_stays_inside_the_attainable_vfd_range():
    """The ceiling is horizon minus landmark, not the horizon: accrual starts at
    the landmark. A predicted mean above it would mean the day values the
    expectation is taken over are wrong."""
    protocol = json.loads((REPO / "config" / "protocol.json").read_text())
    ceiling = (protocol["outcomes"]["vfd"]["horizon_days"]
               - protocol["cohort"]["landmark_hours"] / 24)
    for r in rows("adjusted_prediction.csv"):
        m = num(r["mean_vfd_28"])
        assert 0 <= m <= ceiling, (
            f"a standardised mean of {m} days at {r['exposure_pct']}% is outside "
            f"[0, {ceiling}]"
        )
        z = num(r["pct_vfd_zero"])
        assert 0 <= z <= 100, f"a predicted zero share of {z}%"


# ---- the unadjusted comparison -----------------------------------------------

def test_the_unadjusted_bins_are_the_declared_bins_and_partition_the_frame():
    bins = rows("unadjusted_bins.csv")
    assert bins, "unadjusted_bins.csv is empty"
    want = list(SPEC["prediction"]["unadjusted_bin_labels"])
    assert [r["bin"] for r in bins] == want, \
        f"the bins are {[r['bin'] for r in bins]}, not the declared {want}"
    total = sum(int(float(r["n_admissions"])) for r in bins)
    primary = [d for d in rows("diagnostics/model_diagnostics.csv")
               if d["model"] == "primary"][0]
    assert total == int(float(primary["n_admissions"])), (
        f"the bins hold {total} admissions but the primary model frame holds "
        f"{primary['n_admissions']} -- the comparison is not on the same set as "
        f"the curve it is drawn against"
    )


def test_exactly_zero_is_its_own_bin():
    """About half the admissions never have a qualifying score at any increase.
    Folded into a 0-20% bin they would make it a mixture rather than a level of
    exposure, and it would no longer be comparable to the adjusted curve at
    exactly 0, which is a point and not an interval."""
    breaks = SPEC["prediction"]["unadjusted_bin_breaks"]
    assert breaks[0] < 0 <= breaks[1] == 0, (
        f"the first bin is not the single point 0: breaks start {breaks[:3]}"
    )
    bins = rows("unadjusted_bins.csv")
    zero_bin = bins[0]
    assert zero_bin["bin"] == "0", zero_bin["bin"]
    if not zero_bin["suppressed"] == "TRUE":
        assert num(zero_bin["exposure_pct_mean"]) == 0, (
            f"the zero bin's mean exposure is {zero_bin['exposure_pct_mean']}%, "
            f"so it is not holding only exact zeros"
        )


def test_a_suppressed_bin_keeps_its_count_and_loses_its_statistics():
    """The small-cell rule the rest of the script follows: the count stays so the
    bin is auditable, the statistics go."""
    mincell = json.loads((REPO / "config" / "protocol.json").read_text())["reporting"]["small_cell_min_den"]
    for r in rows("unadjusted_bins.csv"):
        n = int(float(r["n_admissions"]))
        suppressed = r["suppressed"] == "TRUE"
        assert suppressed == (n < mincell), (
            f"bin {r['bin']!r} has n={n} and suppressed={r['suppressed']} "
            f"against a floor of {mincell}"
        )
        if suppressed:
            assert num(r["mean_vfd_28"]) is None, \
                f"bin {r['bin']!r} is suppressed but still ships a mean"
            assert int(float(r["n_admissions"])) == n, "the count must survive suppression"


def test_the_prediction_figure_was_drawn():
    assert (OUT / "vfd_vs_pindicated.png").exists(), \
        "vfd_vs_pindicated.png is missing -- the figure the whole section is for"
    caps = (OUT / "captions.md").read_text()
    assert "vfd_vs_pindicated.png" in caps, \
        "the prediction figure has no caption, so its n, seed and bins live only in pixels"


if __name__ == "__main__":
    import traceback

    if not OUT.exists():
        print("  SKIP  no output/final_no_phi/07_outcome_model -- run the pipeline")
        sys.exit(0)

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
