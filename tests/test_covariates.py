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
NOTES = (REPO / "docs" / "design_notes.md").read_text()
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


def _parse_section_11_table() -> dict[str, dict[str, str]]:
    """Pull the §11 time-varying table out of design_notes.md as {var: {col: cell}}."""
    import re

    start = NOTES.index("### Time-varying covariates")
    end = NOTES.index("### The three missingness classes")
    rows = {}
    for line in NOTES[start:end].split("\n"):
        if not line.startswith("|") or line.startswith("|---"):
            continue
        cells = [c.replace("*", "").strip() for c in line.strip("|").split("|")]
        if len(cells) != 6 or cells[0] == "Variable":
            continue
        names = re.findall(r"`([^`]+)`", cells[0])
        for n in names:
            rows[n] = {"summary": cells[2], "locf": cells[3], "cap": cells[4], "class": cells[5]}
    return rows


def test_design_notes_section_11_matches_the_config():
    """§11 is a hand-written mirror of a machine-readable file, which is exactly
    the pair that drifts. CRRT-dose-lmtp's own build notes are stamped 0.2.0
    against a 0.13.0 config and are wrong on three counts as a result. This makes
    the mirror checkable instead of aspirational."""
    table = _parse_section_11_table()

    missing = set(TIME_VARYING) - set(table)
    assert not missing, f"declared in covariates.json but absent from §11's table: {sorted(missing)}"

    unknown = {v for v in table if v not in TIME_VARYING and v not in
               ("inf_dose", "bolus_dose", "total_dose")}
    assert not unknown, f"named in §11's table but not declared in covariates.json: {sorted(unknown)}"

    for name, spec in TIME_VARYING.items():
        row = table[name]

        # missingness class
        cls = spec["missing_class"]
        if cls != "not_applicable":
            assert f"`{cls}`" in row["class"], (
                f"{name}: §11 says class {row['class']!r}, config says {cls!r}"
            )

        # LOCF cap
        locf = spec.get("locf") or {}
        cap = locf.get("cap_hours")
        if cap is None:
            assert row["cap"] in ("—", "-", ""), (
                f"{name}: config declares no LOCF cap but §11 shows {row['cap']!r}"
            )
        else:
            assert row["cap"] == f"{cap}h", (
                f"{name}: config cap is {cap}h but §11 shows {row['cap']!r}"
            )

        # LOCF eligibility
        if "eligible" in locf:
            shown_yes = "yes" in row["locf"].lower()
            assert shown_yes == locf["eligible"], (
                f"{name}: config LOCF eligible={locf['eligible']} but §11 shows {row['locf']!r}"
            )


def test_repeat_block_decision_carries_its_required_diagnostics():
    """Keeping every encounter block is defensible only if the dependence it
    admits is measured. The decision and the diagnostics that bound it must travel
    together: a stated choice with no reporting obligation attached is how a known
    limitation becomes an unreported one."""
    eb = COV["encounter_blocks"]
    assert eb.get("blocks_per_patient") == "all", (
        "blocks_per_patient changed; the diagnostics below were written for 'all' "
        "and the limitations text in design_notes.md §11 assumes it"
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


def test_design_notes_section_11_lists_every_time_invariant_variable():
    """The time-varying table is mirror-checked above; the time-invariant one was
    not, and it is the half that just grew from 4 variables to 7. A mirror that
    covers only some of the config is a mirror nobody can rely on."""
    import re

    start = NOTES.index("### Time-invariant covariates")
    end = NOTES.index("### Time-varying covariates")
    block = NOTES[start:end]

    labels = {
        name: spec["label"].lower()
        for name, spec in COV["time_invariant"].items()
        if not name.startswith("_")
    }
    lowered = block.lower()
    missing = [
        n for n, lab in labels.items()
        if n not in block and lab.split(",")[0] not in lowered
    ]
    assert not missing, (
        f"declared in covariates.json time_invariant but absent from §11's table: {missing}"
    )


def test_design_notes_section_11_reports_the_right_number_of_checks():
    """§11 claims a count of the checks in this file. Keep it true."""
    import re

    n = len([k for k in globals() if k.startswith("test_")])
    claimed = re.search(r"`tests/test_covariates\.py`\*\* — (\d+) static checks", NOTES)
    assert claimed, "§11 no longer states how many checks this file holds"
    assert int(claimed.group(1)) == n, (
        f"§11 claims {claimed.group(1)} static checks; this file defines {n}"
    )


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
