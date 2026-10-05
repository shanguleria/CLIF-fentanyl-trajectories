"""
Integrity checks for config/covariates.json.

These run BEFORE Phase 0 exists. They cannot yet check that a declared covariate
reaches a column -- that assertion belongs in 01_build_cohort.py and is listed in
the config's _consumption_contract. What they check here is the half that can be
checked statically: that the declarations are internally consistent, that no rule
is a silent no-op, and that the two files declaring the window width agree.

Every failure mode below has been observed in a real CLIF repo.

Run standalone:  .venv/bin/python tests/test_covariates.py
Or under pytest: .venv/bin/python -m pytest tests/ -q
"""
from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
COV = json.loads((REPO / "config" / "covariates.json").read_text())
TEMPLATE = json.loads((REPO / "config" / "config_template.json").read_text())

TIME_VARYING = {k: v for k, v in COV["time_varying"].items() if not k.startswith("_")}
# A class is a block that declares a members list. "reporting", "imputation" and
# "no_variable_drop_threshold" live in the same object but are not classes.
CLASSES = {
    k: v
    for k, v in COV["missing_values"].items()
    if not k.startswith("_") and isinstance(v, dict) and "members" in v
}


def test_every_summary_rule_is_in_the_vocabulary():
    """A summary rule the dispatcher does not know must fail loudly, not default."""
    vocab = {k for k in COV["summary_rules"] if not k.startswith("_")}
    for name, spec in TIME_VARYING.items():
        how = spec.get("summary", "")
        if how.startswith("n/a"):
            continue
        assert how in vocab, f"{name}: summary {how!r} is not in summary_rules {sorted(vocab)}"


def test_every_time_varying_variable_has_exactly_one_missingness_class():
    """Two classes means whichever runs last silently wins; zero means undefined."""
    for name, spec in TIME_VARYING.items():
        assert "missing_class" in spec, f"{name}: no missing_class declared"
        cls = spec["missing_class"]
        assert cls in CLASSES, f"{name}: missing_class {cls!r} is not a declared class"

    seen: dict[str, list[str]] = {}
    for cls, spec in CLASSES.items():
        for member in spec.get("members", []):
            seen.setdefault(member, []).append(cls)
    dupes = {m: c for m, c in seen.items() if len(c) > 1}
    assert not dupes, f"variables in more than one missingness class: {dupes}"


def test_class_membership_lists_agree_with_the_per_variable_declarations():
    """The same fact written in two places is the failure mode that drifts."""
    for cls, spec in CLASSES.items():
        for member in spec.get("members", []):
            assert member in TIME_VARYING, (
                f"missing_values.{cls}.members lists {member!r}, "
                f"which is not a declared time_varying variable"
            )
            declared = TIME_VARYING[member]["missing_class"]
            assert declared == cls, (
                f"{member}: listed under missing_values.{cls} "
                f"but declares missing_class={declared!r}"
            )
    for name, spec in TIME_VARYING.items():
        cls = spec["missing_class"]
        members = CLASSES[cls].get("members", [])
        assert name in members, (
            f"{name} declares missing_class={cls!r} but is absent from that class's members list"
        )


def test_locf_eligibility_and_caps_are_consistent():
    """An eligible variable with no cap is the uncapped-ffill bug; a cap on an
    ineligible variable is a rule nothing will ever consume."""
    for name, spec in TIME_VARYING.items():
        locf = spec.get("locf")
        if locf is None:
            assert spec["missing_class"] == "not_applicable", (
                f"{name}: no locf block and not in the not_applicable class"
            )
            continue
        if locf["eligible"]:
            assert "cap_hours" in locf, (
                f"{name}: LOCF-eligible with no cap_hours. An uncapped lab ffill is "
                f"the bug that produced a 1,331h stale carry-in in a sibling repo."
            )
        else:
            assert "cap_hours" not in locf, (
                f"{name}: LOCF-ineligible but declares cap_hours -- a rule nothing consumes"
            )


def test_no_locf_cap_is_shorter_than_the_window_it_must_bridge():
    """A cap below one window width can never carry anything: it is a declared
    rule that is a guaranteed no-op, which reads as policy and is not."""
    width = COV["windows"]["granular"]["width_hours"]
    for name, spec in TIME_VARYING.items():
        locf = spec.get("locf") or {}
        cap = locf.get("cap_hours")
        if cap is None:
            continue
        assert cap >= width, (
            f"{name}: cap_hours={cap} is below the {width}h window width, so LOCF can never fire"
        )
        assert cap % width == 0, (
            f"{name}: cap_hours={cap} is not a whole number of {width}h windows"
        )


def test_locf_ineligible_variables_are_the_absence_means_something_ones():
    """The state/event/absence distinction is the one that matters. An
    absence-means-zero variable that was LOCF-eligible would invent exposure."""
    for name, spec in TIME_VARYING.items():
        cls = spec["missing_class"]
        locf = spec.get("locf") or {}
        if cls in ("absence_means_zero", "absence_means_not_ventilated"):
            assert locf.get("eligible") is False, (
                f"{name}: class {cls} must not be LOCF-eligible -- carrying it forward "
                f"invents exposure for a patient who had been weaned off"
            )


def test_nee_coefficients_match_its_declared_source_categories():
    """A drug in one list and not the other is silently dropped from, or silently
    absent from, the score."""
    nee = TIME_VARYING["nee"]
    cats = set(nee["source"]["categories"])
    coefs = set(nee["coefficients"])
    assert cats == coefs, (
        f"NEE source categories and coefficients disagree: "
        f"only in source={sorted(cats - coefs)}, only in coefficients={sorted(coefs - cats)}"
    )
    preferred = {k for k in nee["preferred_units"] if not k.startswith("_")}
    assert preferred <= coefs, (
        f"preferred_units names a drug that is not in the NEE table: {sorted(preferred - coefs)}"
    )


def test_repeat_block_decision_carries_its_required_diagnostics():
    """Keeping every encounter block is defensible only if the dependence it
    admits is measured. The decision and the diagnostics that bound it must travel
    together: a stated choice with no reporting obligation attached is how a known
    limitation becomes an unreported one."""
    eb = COV["encounter_blocks"]
    assert eb.get("blocks_per_patient") == "all", (
        "blocks_per_patient changed; the diagnostics below were written for 'all'"
    )
    diags = eb.get("required_dependence_diagnostics", [])
    assert len(diags) >= 3, (
        f"only {len(diags)} dependence diagnostics declared; keeping all blocks "
        f"requires at least the repeat-patient count, the same-patient-class check, "
        f"and the limitations sentence"
    )
    joined = " ".join(diags).lower()
    for token in ("more than one encounter block", "same patient", "limitations"):
        assert token in joined, f"required diagnostic missing: {token!r}"


def test_nee_carries_a_literature_citation():
    """These six factors are manuscript-bound constants. CRRT-dose-lmtp carries the
    identical table with no citation anywhere in that repo, which is how an
    unsourced number reaches a methods section. Keep the source attached to the
    values, not in a comment someone can delete."""
    nee = TIME_VARYING["nee"]
    cit = nee.get("citation", "")
    assert cit.strip(), "NEE coefficients have no citation"
    assert "doi:" in cit.lower(), f"NEE citation carries no DOI: {cit!r}"


def test_window_arithmetic_is_self_consistent():
    for grid in ("granular", "extended"):
        w = COV["windows"][grid]
        assert w["extent_hours"] % w["width_hours"] == 0, (
            f"{grid}: extent {w['extent_hours']}h is not a whole number of {w['width_hours']}h windows"
        )
        expected = w["extent_hours"] // w["width_hours"]
        assert w["n_windows"] == expected, (
            f"{grid}: n_windows={w['n_windows']} but extent/width={expected}"
        )


def test_window_width_agrees_with_the_site_config_template():
    """The same threshold written down twice, copies now disagreeing, is the
    second classic config-integrity failure. covariates.json and the site config
    both declare the grid; they must not drift."""
    cohort = TEMPLATE["cohort"]
    pairs = [
        ("granular.width_hours", COV["windows"]["granular"]["width_hours"], "window_hours", cohort["window_hours"]),
        ("granular.extent_hours", COV["windows"]["granular"]["extent_hours"], "granular_extent_hours", cohort["granular_extent_hours"]),
        ("extended.width_hours", COV["windows"]["extended"]["width_hours"], "extended_window_hours", cohort["extended_window_hours"]),
        ("extended.extent_hours", COV["windows"]["extended"]["extent_hours"], "extended_extent_hours", cohort["extended_extent_hours"]),
    ]
    for cov_key, cov_val, cfg_key, cfg_val in pairs:
        assert cov_val == cfg_val, (
            f"drift: covariates.json windows.{cov_key}={cov_val} but "
            f"config cohort.{cfg_key}={cfg_val}"
        )


def test_the_race_collapse_map_covers_every_permissible_value():
    """Race is collapsed for Table 1 display. A permissible CLIF value with no
    entry falls to _default silently, which is fine -- but a map that names a
    value NOT in the permissible list is a typo that will never fire, and the
    collapsed table would then show a category nobody ever lands in."""
    race = COV["time_invariant"]["race"]
    perm = set(race["permissible_values"])
    cmap = race["reporting_collapse"]
    named = {k for k in cmap if not k.startswith("_")}
    assert named <= perm, (
        f"reporting_collapse names values that are not permissible CLIF race "
        f"categories: {sorted(named - perm)}"
    )
    assert cmap.get("_default"), "no _default target for the uncollapsed values"
    # Unknown must not silently become Other: it means 'asked, not answered'.
    if not cmap.get("_default_includes_unknown", False):
        assert cmap.get("Unknown") == "Unknown", (
            "Unknown is not folded into _default unless "
            "_default_includes_unknown is set true"
        )


def test_every_declared_sedative_is_wired_end_to_end():
    """A sedative is declared in four places and must be live in all four.

    When a config governs N items, check all N. Adding dexmedetomidine touched
    config.json (the category to load), covariates.json exposure.sedatives
    (column and unit), covariates.json dose_units (how to convert it) and
    outlier_config.json (its bounds). Any one of those missing is a silent
    partial application: the run still exits 0 and the column is empty, wrong,
    or unbounded.
    """
    cfg = json.loads((REPO / "config" / "config.json").read_text()) \
        if (REPO / "config" / "config.json").exists() else TEMPLATE
    cats = cfg["medications"]["other_sedative_categories"]
    spec = COV["exposure"]["sedatives"]
    bounds = json.loads((REPO / "config" / "outlier_config.json").read_text())
    targets = {d for t, v in COV["dose_units"].items() if not t.startswith("_")
               for d in v.get("_target_for", [])}

    assert cats, "no sedative categories configured"
    for cat in cats:
        col = f"{cat}_dose"
        assert col in spec["columns"], (
            f"{cat} is loaded but {col} is not in exposure.sedatives.columns, so "
            f"Phase 0 raises rather than building it"
        )
        assert spec["units"].get(col), f"{col} has no declared unit"
        assert cat in targets, (
            f"{cat} has no dose_units target, so convert() raises on its first row"
        )
        assert cat in bounds["med_dose_raw"], (
            f"{cat} has no raw outlier bound; its rows would pass unbounded and "
            f"only appear as a NO BOUND line in the log"
        )
    # and nothing declared that is not loaded
    for col in spec["columns"]:
        assert col.removesuffix("_dose") in cats, (
            f"{col} is declared but its category is not in "
            f"config.json medications.other_sedative_categories"
        )


def test_the_sedative_target_unit_matches_its_declared_unit():
    """dose_units decides what the converter produces; exposure.sedatives.units
    is what the protocol, the figures and the docs promise."""
    spec = COV["exposure"]["sedatives"]
    target_of = {d: t for t, v in COV["dose_units"].items() if not t.startswith("_")
                 for d in v.get("_target_for", [])}
    for col, unit in spec["units"].items():
        drug = col.removesuffix("_dose")
        assert target_of.get(drug) == unit, (
            f"{drug}: dose_units target is {target_of.get(drug)!r} but "
            f"exposure.sedatives.units says {unit!r}"
        )


def test_the_extended_grid_is_labelled_while_nothing_consumes_it():
    """The declared-but-never-consumed failure mode, caught by construction.

    windows.extended and the two config extended_* keys are read by no code, and
    the mirror test above only proves the two DECLARATIONS agree with each other
    -- not that either reaches an operation. While that is true the block must
    carry a _STATUS saying so, or a reader takes it for a live setting.
    """
    src = []
    for d in ("code", "validation"):
        for pat in ("*.py", "*.R"):
            src += list((REPO / d).rglob(pat))
    src = [f for f in src if "__pycache__" not in str(f)]
    consumers = [
        f.relative_to(REPO)
        for f in src
        for key in ("extended_window_hours", "extended_extent_hours")
        if key in f.read_text()
    ]
    ext = COV["windows"]["extended"]
    if consumers:
        assert "_STATUS" not in ext, (
            f"windows.extended is consumed by {consumers}, so remove its _STATUS "
            f"rather than leaving the config claiming it is not built"
        )
    else:
        assert ext.get("_STATUS"), (
            "nothing reads the extended grid keys, so windows.extended must carry "
            "a _STATUS saying it is declared and not produced"
        )
        for f in (REPO / "config" / "config.json",
                  REPO / "config" / "config_template.json"):
            if not f.exists():
                continue
            cohort = json.loads(f.read_text())["cohort"]
            assert cohort.get("_comment_extended_grid"), (
                f"{f.name} declares extended_* with no note that nothing reads them"
            )


def test_every_declared_variable_is_documented():
    """A variable with no label and no source or algorithm cannot be implemented
    from this file, which is the only thing this file is for."""
    blocks = {"time_invariant": COV["time_invariant"], "time_varying": TIME_VARYING}
    for block, spec in blocks.items():
        for name, v in spec.items():
            if name.startswith("_"):
                continue
            assert v.get("label"), f"{block}.{name}: no label"
            assert v.get("source") or v.get("algorithm") or v.get("pipeline"), (
                f"{block}.{name}: neither a source, an algorithm, nor a pipeline"
            )


def test_small_cell_thresholds_do_not_contradict_the_site_config():
    pattern_min = COV["missing_values"]["reporting"]["pattern_min_cell"]
    site_min = TEMPLATE["reporting"]["small_cell_min_den"]
    assert pattern_min >= site_min, (
        f"pattern_min_cell={pattern_min} is below the site config's "
        f"small_cell_min_den={site_min}, so pattern reporting would be the laxer of the two"
    )


def _r_sources():
    return {f.name: f.read_text()
            for f in (REPO / "code").rglob("*.R")}


def test_the_dose_band_cuts_are_read_from_the_config_not_hardcoded():
    """Every consumer must reach the config. A literal 200 or 400 sitting in a
    cut() call is the drifted-duplicate failure: the consortium moves the band
    and one script keeps cutting at the old edge."""
    spec = COV["exposure"]["dose_states"]
    cuts = spec["cuts"]
    consumers = {n: t for n, t in _r_sources().items()
                 if "derive_dose_states(" in t and n != "states.R"}
    assert consumers, "nothing calls derive_dose_states; the config block is dead"
    for name, text in consumers.items():
        assert "DS_SPEC" in text and "dose_states" in text, (
            f"{name} calls derive_dose_states but never reads "
            f"exposure.dose_states from covariates.json")
        for c in cuts:
            assert f"c({c}" not in text and f", {c})" not in text, (
                f"{name} appears to hardcode the band edge {c}; read it from "
                f"the config instead")


def test_derive_dose_states_refuses_a_default_cut():
    """states.R must not default `cuts`: a default is what keeps a site cutting
    at the old edges after the consortium moves them."""
    src = (REPO / "code" / "utils" / "states.R").read_text()
    sig = src[src.index("derive_dose_states <- function("):]
    sig = sig[:sig.index(")")]
    assert "cuts," in sig or sig.rstrip().endswith("cuts"), sig
    assert "cuts =" not in sig, (
        "derive_dose_states must take `cuts` with NO default -- see "
        "covariates.json exposure.dose_states._cuts_MUST_BE_ABSOLUTE")


def test_the_dose_band_variable_is_produced_by_phase_0():
    """The band is defined on window_mcg, so Phase 0 must emit it and it must be
    the SUM over hourly cells -- not total_dose * window_hours, which overstates
    a window that straddles discharge."""
    spec = COV["exposure"]["dose_states"]
    build = (REPO / "code" / "01_build_cohort.py").read_text()
    assert spec["variable"] in COV["exposure"]["columns"], (
        f"{spec['variable']} is the dose-state variable but is not declared in "
        f"exposure.columns")
    assert 'inf_mcg="sum"' in build, (
        "window_mcg must come from a SUM over the hourly cells; reconstructing "
        "it as total_dose * window_hours overstates short windows")
    assert 'out["window_mcg"] = out["inf_mcg"] + out["bolus_mcg"]' in build


def test_the_dose_bands_and_labels_agree_in_length():
    spec = COV["exposure"]["dose_states"]
    assert len(spec["labels"]) == len(spec["cuts"]) + 2, (
        "labels = one zero band + one per interval the cuts create")
    assert spec["cuts"] == sorted(set(spec["cuts"])) and spec["cuts"][0] > 0


def test_the_exemplar_selection_rule_is_declared_and_consumed():
    """F1's selection rule is protocol, not a local preference.

    Baker et al., which F1 is modelled on, states no selection rule at all; the
    whole point of declaring ours is that a reader can check it. A rule that
    lives in config.json would be worse than none -- that file is gitignored and
    site-local, so two sites could silently draw exemplars by different rules and
    the figures would not be comparable. It must be here, and it must actually
    reach the code that applies it.
    """
    spec = COV.get("exemplar")
    assert spec is not None, "covariates.json must declare an `exemplar` block"

    build = (REPO / "code" / "01_build_cohort.py").read_text()
    figure = (REPO / "code" / "03_exemplar.R").read_text()

    # Every declared key must be claimed by someone. A threshold nobody reads is
    # a policy declaration that is not a policy -- the failure mode this repo
    # has already hit twice (small_cell_min_den, grid_resolution_minutes).
    live = {k for k in spec if not k.startswith("_")}
    unread = {k for k in live if k not in build and k not in figure}
    assert not unread, (
        f"exemplar keys declared but read by nothing: {sorted(unread)}. Either "
        f"consume them or delete them.")

    assert spec["draw"] == "uniform_random_from_eligible", (
        "the draw must be a seeded uniform draw from the eligible set; a "
        "weighted score needs weights nobody can defend, and hand-picking "
        "reintroduces the bias the rule exists to remove")
    for k in ("require_extubated_by_hours", "min_continuous_hours",
              "min_boluses", "min_rass_observations", "min_nvps_observations"):
        assert k in spec, f"exemplar rule is missing {k}"
        assert isinstance(spec[k], int) and spec[k] >= 0, (
            f"exemplar.{k} must be a non-negative integer, got {spec[k]!r}")

    assert spec["min_boluses"] >= 2, (
        "F1 must be able to show the near-simultaneous bolus staircase, which "
        "needs at least two administrations")
    assert spec["require_extubated_by_hours"] <= COV["windows"]["granular"]["extent_hours"], (
        "an exemplar cannot be required to extubate after the window ends")


def test_the_exemplar_figure_breaks_the_step_at_the_locf_cap():
    """The panel must not assert a score persisted longer than the analytic
    table lets it. Both must read the same key rather than restate the number."""
    figure = (REPO / "code" / "03_exemplar.R").read_text()
    assert "COV$time_varying$rass$locf$cap_hours" in figure, (
        "03_exemplar.R must read the step-breaking threshold from the same "
        "covariates.json key Phase 0 caps LOCF with, never restate it")
    assert COV["time_varying"]["rass"]["locf"]["cap_hours"] == \
           COV["time_varying"]["nvps"]["locf"]["cap_hours"], (
        "rass and nvps must share one cap, or the two panels break at different "
        "gap lengths for no stated reason")


def test_the_titration_rule_is_declared_and_consumed():
    """Pairing a bolus to a rate change is protocol: two sites using different
    windows or thresholds produce numbers that cannot be pooled. Every declared
    key must reach the code that applies it -- the same check that caught
    `require_extubated_by_hours` declared and never read."""
    spec = COV.get("titration")
    assert spec is not None, "covariates.json must declare a `titration` block"

    build = (REPO / "code" / "01_build_cohort.py").read_text()
    figure = (REPO / "code" / "05_titration.R").read_text()
    live = {k for k in spec if not k.startswith("_")}
    unread = {k for k in live if k not in build and k not in figure}
    assert not unread, (
        f"titration keys declared but read by nothing: {sorted(unread)}")

    assert spec["min_rate_change_mcg_hr"] > 0
    assert spec["bolus_window_minutes"] in spec["window_sensitivity_minutes"], (
        "the primary window must appear in the sensitivity list, or the curve "
        "cannot be read against the headline number")
    assert spec["adherence_landmark_hours"] < COV["windows"]["granular"]["extent_hours"], (
        "adherence must be measurable strictly inside the observation window, "
        "or it is not an exposure measured before outcome accrual")


def test_the_titration_analysis_never_reads_the_hourly_grid():
    """The grid bins to whole hours. A rate change moved up to an hour from the
    bolus that accompanied it makes a 30-minute pairing window meaningless, so
    this analysis must use raw charted timestamps -- and must say so loudly
    enough that a future edit cannot quietly reintroduce the grid."""
    build = (REPO / "code" / "01_build_cohort.py").read_text()
    figure = (REPO / "code" / "05_titration.R").read_text()
    assert "_infusion_events(" in build, (
        "titration_export must build events from _infusion_events(), the "
        "per-record frame, not from infusion_grid()")
    assert "titration_rate_events.parquet" in figure
    assert "infusion_grid" not in figure, (
        "05_titration.R must never touch the hourly grid")
    for src, name in ((build, "01_build_cohort.py"), (figure, "05_titration.R")):
        assert "%% 1" in src or "% 1 != 0" in src, (
            f"{name} must assert its timestamps are sub-hourly; whole-hour "
            f"timestamps would mean the grid leaked in")


def test_the_table1_stratification_is_declared_and_ordered():
    """Table 1 is stratified on predominant fentanyl intensity. Three things
    must hold or the strata mean something other than what the text will claim.
    """
    spec = COV.get("predominant_intensity")
    assert spec is not None, "covariates.json must declare predominant_intensity"

    # 1. it reuses the SAME bands the state figures use, so the two cannot drift
    assert spec["bands_from"] == "exposure.dose_states"
    assert spec["variable"] == COV["exposure"]["dose_states"]["variable"], (
        "the stratification must be cut on the same variable as the bands")

    # 2. the tie-break is declared, not left to sort order. 7.3% of episodes tie.
    assert spec["tie_break"] in ("higher", "lower"), (
        "a tie resolved by whatever order the levels happen to be in is not a rule")

    # 3. both artifacts stay visible as Table 1 rows
    assert spec["report_modal_share"] is True, (
        "modal share says how decisive each label is; a coin-flip label is not a stratum")
    assert spec["report_at_risk_hours"] is True, (
        "at-risk hours expose the duration confound -- a short course cannot "
        "accumulate zero windows")

    build = (REPO / "code" / "02_descriptive_cohort.R").read_text()
    states = (REPO / "code" / "utils" / "states.R").read_text()
    live = {k for k in spec if not k.startswith("_")}
    unread = {k for k in live if k not in build and k not in states}
    assert not unread, (
        f"predominant_intensity keys declared but read by nothing: {sorted(unread)}")


def test_the_stratification_is_never_called_a_trajectory():
    """Modal time ignores order -- only 44.9% of episodes end in the state they
    began in. Calling these strata trajectories or latent classes is exactly the
    claim the tabled gbmt work could not support, so the words are banned from
    the code that builds them."""
    for rel in ("code/utils/states.R", "code/02_descriptive_cohort.R"):
        src = (REPO / rel).read_text()
        for line in src.splitlines():
            low = line.lower()
            if "predominant" not in low and "modal" not in low:
                continue
            if line.strip().startswith("#"):
                continue          # comments may discuss why it is NOT one
            assert "trajectory" not in low and "latent class" not in low, (
                f"{rel} names the stratification a trajectory: {line.strip()[:70]}")


# ----------------------------------------------------------- unit variation
# Splitting the titration metric by ICU unit and calendar year. The attribution
# rule and the resample count are protocol: two sites crediting units by
# different rules, or drawing different numbers of resamples, produce intervals
# that cannot be set beside each other.

def test_the_unit_variation_rule_is_declared_and_consumed():
    """Two of these keys are read by Phase 0, which emits care_intervals.parquet
    and anchor_year, and the rest by the R script that splits the metric. Every
    declared key must reach the code that applies it -- the same check that
    caught `require_extubated_by_hours` declared and never read."""
    spec = COV.get("unit_variation")
    assert spec is not None, "covariates.json must declare a `unit_variation` block"

    build = (REPO / "code" / "01_build_cohort.py").read_text()
    figure = (REPO / "code" / "06_unit_variation.R").read_text()
    live = {k for k in spec if not k.startswith("_")}
    unread = {k for k in live if k not in build and k not in figure}
    assert not unread, (
        f"unit_variation keys declared but read by nothing: {sorted(unread)}. "
        f"Either consume them or delete them.")

    assert spec["bootstrap_resamples"] > 0
    assert spec["min_events_per_unit"] >= TEMPLATE["reporting"]["small_cell_min_den"], (
        "min_events_per_unit must clear reporting.small_cell_min_den, which is "
        "the disclosure floor, not a suggestion")
    assert "care_setting" in spec["unit_keys"], (
        "care_setting (location_type, else location_category) is the only unit "
        "key another CLIF site can reproduce -- both inputs are required mCIDE "
        "columns, where location_name is optional free text. Dropping it would "
        "make this analysis unpoolable.")
    assert spec["year_from"] in ("anchor_dttm", "block_admission_dttm"), (
        "year_from must name a datetime the cohort frame actually carries")


def test_anchor_year_is_declared_so_the_build_must_produce_it():
    """assert_config_is_honoured() fails the build when a declared time_invariant
    key reaches no column. Declaring anchor_year there is what makes the column a
    guarded contract rather than an incidental side effect of Phase 0."""
    ti = COV["time_invariant"]
    assert "anchor_year" in ti, (
        "anchor_year must be declared in time_invariant, or nothing checks that "
        "Phase 0 actually built it")
    assert ti["anchor_year"]["source"]["column"] == COV["unit_variation"]["year_from"], (
        "the declared source column and unit_variation.year_from have drifted")


def test_the_unit_split_is_never_called_a_quality_measure():
    """No evidence establishes that pairing a bolus with an up-titration is
    better care, so a unit with a lower rate is not worse. The words are banned
    from the code that builds the split -- the same guard, for the same reason,
    as test_the_stratification_is_never_called_a_trajectory."""
    banned = ("quality", "performance", "benchmark")
    src = (REPO / "code" / "06_unit_variation.R").read_text()
    for line in src.splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            continue              # comments may discuss why it is NOT one
        low = line.lower()
        hit = [w for w in banned if w in low]
        assert not hit, (
            f"06_unit_variation.R calls the split a {hit[0]} measure: "
            f"{stripped[:70]}")


def test_the_funnel_limits_claim_no_adjustment_they_do_not_make():
    """Until 2026-10-05 the envelope was binomial limits * sqrt(deff). deff
    measures WITHIN-episode clustering and is grouping-invariant; Spiegelhalter's
    phi measures BETWEEN-unit over-dispersion. Same algebra, different quantity,
    under a citation that implied the other one -- and with nothing
    risk-adjusted, widening absorbed the variation the figure exists to show.
    Both are now measured and applied nowhere. This guards the exact regression:
    re-introducing the factor, or a caption that claims it."""
    src = (REPO / "code" / "06_unit_variation.R").read_text()
    code = [ln for ln in src.splitlines() if not ln.strip().startswith("#")]

    applied = [ln.strip() for ln in code if "sqrt(deff" in ln.replace(" ", "")]
    assert not applied, (
        "the funnel envelope must be unadjusted binomial, but something "
        f"applies sqrt(deff) again: {applied[:2]}")

    low = "\n".join(code).lower()
    for claim in ("widened by", "same clustering correction"):
        assert claim not in low, (
            f"06_unit_variation.R says the limits are {claim!r}; they are "
            "unadjusted binomial. Fix the code or the caption, not both ways.")


def test_over_dispersion_is_measured_declared_and_owned():
    """The funnel's limits assume independence between units, so whether that
    holds has to be reported rather than assumed. A diagnostic that is computed
    but never written cannot be read; an output written but not OWNED survives
    clear_owned_outputs() and goes stale inside the upload set."""
    src = (REPO / "code" / "06_unit_variation.R").read_text()

    owned = src.split("OWNED <- list(")[1].split("))")[0]
    assert "overdispersion.csv" in owned, (
        "overdispersion.csv must be in OWNED, or the first re-run leaves a "
        "stale copy in output/final_no_phi/")
    assert 'file.path(dirs$phase, "overdispersion.csv")' in src, (
        "overdispersion.csv is OWNED but nothing writes it")

    trim = COV["unit_variation"]["overdispersion_trim"]
    assert 0 <= trim < 0.5, (
        "overdispersion_trim is a Winsorisation fraction PER END. At 0.5 every "
        "unit is pulled to the median and phi reads 0 whatever the data says.")


def test_the_caterpillar_caption_refuses_to_be_read_as_a_ranking():
    """Ordering units by rate invites reading the rank as the result. The
    intervals are what say how much of the ordering is resolvable -- at UCMC
    4 of 6 ordered-adjacent pairs overlap on co-administration -- so the caption
    has to state both the count and that overlap is not a pairwise test.

    The count must be COMPUTED by the run that draws the panel. A hardcoded
    number goes stale the first time a site has different units, and would be
    the worst kind of stale: still plausible, and about the figure's own
    interpretation."""
    src = (REPO / "code" / "06_unit_variation.R").read_text()

    assert "n_adj_overlap <- " in src, (
        "the caterpillar must compute its ordered-adjacent overlap count")
    assert "ord$ci_lo[-1] <= ord$ci_hi[-n_units]" in src, (
        "the overlap count must come from the ordered intervals actually drawn")
    assert "n_adj_overlap," in src, (
        "n_adj_overlap is computed but never reaches the caption -- the same "
        "declared-but-unconsumed failure this suite exists to catch")

    low = src.lower()
    assert "not a pairwise test" in low, (
        "the caterpillar caption must say overlap is not a pairwise test: "
        f"{COV['unit_variation']['_INTERVALS_DESCRIBE_THEY_DO_NOT_ADJUST'][:60]}...")
    assert "sampling only" in low, (
        "the caption must say the interval covers sampling only -- it is silent "
        "on whether the charting behind the numerator is accurate, which the ED "
        "and procedural findings showed is the larger doubt")


# ------------------------------------------------------------- indication
# Was the dose change prompted by anything documented? Thresholds and the
# window are protocol: two sites cutting NVPS at different values, or looking
# over different windows, produce numbers that cannot be pooled.

# Published ranges of each instrument. A threshold outside its own scale is
# silently unreachable -- the analysis would run and report a clean zero.
INSTRUMENT_RANGE = {
    "NVPS": (0, 10),      # Non-Verbal Pain Scale
    "CPOT": (0, 8),       # Critical-Care Pain Observation Tool
    "RASS": (-5, 4),      # Richmond Agitation-Sedation Scale
    "PAINAD": (0, 10),    # Pain Assessment in Advanced Dementia
    "BPS": (3, 12),       # Behavioural Pain Scale
}


def test_the_indication_rule_is_declared_and_consumed():
    """Phase 0 decides what to read and export; the R script applies the
    thresholds. Every declared key must reach one of them."""
    spec = COV.get("indication")
    assert spec is not None, "covariates.json must declare an `indication` block"

    build = (REPO / "code" / "01_build_cohort.py").read_text()
    figure = (REPO / "code" / "05_titration.R").read_text()
    live = {k for k in spec if not k.startswith("_")}
    unread = {k for k in live if k not in build and k not in figure}
    assert not unread, (
        f"indication keys declared but read by nothing: {sorted(unread)}. "
        f"Either consume them or delete them.")


def test_indication_thresholds_are_inside_each_instrument_range():
    spec = COV["indication"]
    both = {**spec["pain_instruments"], **spec["sedation_instruments"]}
    assert both, "at least one instrument must be declared"
    for name, thr in both.items():
        assert name in INSTRUMENT_RANGE, (
            f"{name} has no published range in this test's table. Add it, with "
            f"a source, before declaring a threshold on it.")
        lo, hi = INSTRUMENT_RANGE[name]
        assert lo <= thr <= hi, (
            f"{name} threshold {thr} is outside its scale [{lo}, {hi}], so it "
            f"can never fire and the analysis would report a clean zero")


def test_the_primary_indication_window_is_in_the_sensitivity_sweep():
    spec = COV["indication"]
    assert spec["indication_window_hours"] in spec["indication_window_sensitivity_hours"], (
        "the primary window must appear in the sweep, or the curve cannot be "
        "read against the headline number")
    assert max(spec["indication_window_sensitivity_hours"]) <= \
        COV["windows"]["granular"]["extent_hours"], (
        "a window wider than the observation grid cannot be evaluated")


def test_an_optional_instrument_is_never_declared_as_a_covariate():
    """assessment_covariates() raises SystemExit when a declared category
    matches no rows -- correct for a required covariate, fatal for an optional
    site-specific instrument. CPOT is not charted at UCMC, so declaring it in
    time_varying would kill the build here."""
    spec = COV["indication"]
    both = {**spec["pain_instruments"], **spec["sedation_instruments"]}
    declared_as_covariate = {
        v["source"]["category"]
        for k, v in COV["time_varying"].items()
        if not k.startswith("_") and (v.get("source") or {}).get("table") == "assessments"
    }
    optional = set(both) - declared_as_covariate
    assert optional, (
        "every indication instrument is also a required time_varying covariate, "
        "so the tolerant path is never exercised -- check this is intended")
    # RASS and NVPS are legitimately both; the point is that the set is not
    # forced to overlap, and that anything site-specific stays out.
    assert "CPOT" not in declared_as_covariate, (
        "CPOT is declared as a time_varying covariate. assessment_covariates() "
        "raises on a category matching no rows, and CPOT is not charted at "
        "every site -- this would fail the build wherever it is absent.")


def test_the_renamed_unit_figures_are_retired_not_orphaned():
    """clear_owned_outputs() only removes what a script NAMES. The five figures
    renamed on 2026-09-30 (unit_caterpillar_*, unit_funnel_*, year_trend.png)
    would otherwise sit in the shareable tree beside their replacements looking
    current -- the failure the RETIRED mechanism exists to prevent."""
    src = (REPO / "code" / "06_unit_variation.R").read_text()
    block = src[src.index("RETIRED <-"):src.index("n_cleared <-")]
    # Strip comments first. Checking raw text would pass on a RETIRED that had
    # been commented out, since the old names would survive in the prose --
    # a guard that cannot fail is not a guard (lessons.md #17).
    code = "\n".join(l.split("#", 1)[0] for l in block.splitlines())
    for old in ("caterpillar_", "funnel_", "trend.png"):
        assert old in code, (
            f"the RETIRED block no longer names {old!r} in CODE; a renamed "
            f"figure will be left behind in output/final_no_phi/06_unit_variation/")
    assert '"unit_' in code, "RETIRED must name the old `unit_` prefix explicitly"


def test_both_titration_metrics_are_split_by_unit():
    """06 splits two metrics. If one is dropped from METRICS its figures stop
    being owned, and a stale copy survives the next run."""
    src = (REPO / "code" / "06_unit_variation.R").read_text()
    assert 'METRICS <- c("coadmin", "indication")' in src, (
        "the METRICS vector changed; check FIGS, OWNED and RETIRED together")
    for m in ("coadmin", "indication"):
        assert f"{m} = list(" in src, f"ANALYSES has no {m} entry"


if __name__ == "__main__":
    import sys, traceback

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


# --------------------------------------------------------------- dose states
# exposure.dose_states is the newest config block, and a config block nothing
# reads is the failure mode this repo has hit six times: it looks like a policy
# declaration, so the next person trusts it, and editing it is a no-op. These
# tests assert the block REACHES the operations it names.
