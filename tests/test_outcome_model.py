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
    """
    po = rows("diagnostics/proportional_odds_test.csv")
    diag = rows("diagnostics/model_diagnostics.csv")
    assert diag, "no models in model_diagnostics.csv"
    for model in {d["model"] for d in diag}:
        terms = [r for r in po if r["model"] == model and r["term"] != "<none>"]
        assert terms, f"{model}: no term rows in the proportional-odds test"
        tested = [r for r in terms if num(r["p_value"]) is not None]
        assert tested, (
            f"{model}: the proportional-odds test ran but tested NOTHING -- "
            f"{len(terms)} term rows, every p-value NA. The assumption is "
            f"UNTESTED, not satisfied. See this test's docstring."
        )


def test_every_proportional_odds_row_declares_its_status():
    """The status column is what makes an untestable model distinguishable from
    a passing one in the shipped file. Counting significant terms over an all-NA
    table yields zero, which reads as the assumption holding."""
    po = rows("diagnostics/proportional_odds_test.csv")
    assert po, "proportional_odds_test.csv is empty"
    assert "status" in po[0], "the status column is gone -- an untestable model would again be silent"
    allowed = {"tested", "not_testable", "errored"}
    for r in po:
        assert r["status"] in allowed, f"unknown status {r['status']!r}"
    assert not [r for r in po if r["status"] != "tested"], (
        "a model's proportional-odds assumption could not be tested: "
        + ", ".join(sorted({r["model"] for r in po if r["status"] != "tested"}))
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
    """STROBE ends at the complete-case count; the model then drops admissions in
    units below min_admissions_per_unit. The difference must be exactly that drop
    -- an unexplained gap between a published flow and a model n is the
    denominator trap this repo has already hit three times."""
    flow = rows("strobe.csv")
    primary = [d for d in rows("diagnostics/model_diagnostics.csv")
               if d["model"] == "primary"]
    assert primary, "no primary model in model_diagnostics.csv"
    p = primary[0]
    complete_case = int(flow[-1]["n"])
    dropped = int(p["n_dropped_small_unit"])
    assert complete_case - dropped == int(p["n_admissions"]), (
        f"STROBE ends at {complete_case} and the primary model fits "
        f"{p['n_admissions']}, which {dropped} unit-drop(s) do not explain"
    )


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
