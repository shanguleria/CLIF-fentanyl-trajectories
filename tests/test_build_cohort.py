"""Tests for the pure logic in code/01_build_cohort.py, on synthetic frames.

No CLIF data is read. The stages that need real tables are exercised in a
smoke run; what is tested here is the logic that is wrong-or-right independent
of the data: LOCF caps, absence-means-zero, the consumption assertion, unit
conversion, and the Severinghaus transform.
"""
from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "code"))
_spec = importlib.util.spec_from_file_location("build", REPO / "code" / "01_build_cohort.py")
B = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(B)


def _long(n_blocks=2, n_win=None):
    n_win = n_win or B.N_WINDOWS
    rows = []
    for b in range(n_blocks):
        for w in range(n_win):
            rows.append({"encounter_block": f"b{b}", "window_idx": w,
                         "window_start_hr": w * B.WINDOW_H, "alive_admitted": True})
    return pd.DataFrame(rows)


# ----------------------------------------------------- module-level contract
def test_every_module_constant_the_loaders_need_exists():
    """A NameError inside a loader only fires at call time, so importing the
    module is not enough to catch a deleted constant -- this has bitten three
    times during editing, each costing a full run to discover."""
    required = {
        "WINDOW_H": int, "N_WINDOWS": int, "EXTENT_H": int, "STITCH_H": int,
        "EPISODE_GAP_H": (int, float), "MIN_IMV_H": (int, float),
        "INF_HOLD_H": (int, float), "NEE_HOLD_H": (int, float),
        "FIO2_LOOKBACK_H": (int, float), "SPO2_CEILING": (int, float),
        "RA_FIO2": float,
    }
    for name, kind in required.items():
        assert hasattr(B, name), f"module constant {name} is missing"
        assert isinstance(getattr(B, name), kind), f"{name} has the wrong type"

    for name in ("LAB_NEEDED", "VITAL_NEEDED", "ASSESS_NEEDED", "LAB_VARS",
                 "ZERO_VARS", "NOT_VENT_VARS", "SOFA_PRESSORS", "SOFA_INPUT_CAPS",
                 "NEE_COEF", "OWNED", "EXPOSURE", "RETIRED_OUTPUTS"):
        assert hasattr(B, name), f"module constant {name} is missing"
        assert len(getattr(B, name)), f"{name} is empty"


def test_the_loaders_reference_only_names_that_exist():
    """Compile-time check on every global a loader touches."""
    import inspect

    for fn_name in ("load_core", "load_cohort_tables", "_mac_categories",
                    "assert_categories_present"):
        fn = getattr(B, fn_name)
        for name in fn.__code__.co_names:
            if name.isupper() and "_" in name or name.isupper():
                assert hasattr(B, name) or name in dir(__builtins__), (
                    f"{fn_name} references {name}, which does not exist"
                )


def test_every_lab_the_analysis_needs_is_in_the_read_filter():
    """The category filter is pushed down to the parquet read, so a lab absent
    from it is absent from the data, silently."""
    needed = set(B.LAB_VARS.values()) | {"po2_arterial", "creatinine",
                                         "platelet_count"}
    assert needed <= set(B.LAB_NEEDED), (
        f"not read from disk: {sorted(needed - set(B.LAB_NEEDED))}"
    )



# --------------------------------------------------------------- LOCF caps
def test_locf_respects_the_declared_cap_in_windows_not_rows():
    """A 24h cap on a 4h grid is 6 windows, not 24."""
    df = _long(n_blocks=1)
    df["lactate"] = np.nan
    df.loc[0, "lactate"] = 4.2
    out, _, _ = B.apply_missingness(df)
    filled = out["lactate"].notna().sum()
    cap_h = B.COV["time_varying"]["lactate"]["locf"]["cap_hours"]
    assert filled == 1 + cap_h // B.WINDOW_H == 7, f"carried {filled} windows, expected 7"


def test_locf_flag_marks_only_filled_cells():
    df = _long(n_blocks=1)
    df["bun"] = np.nan
    df.loc[0, "bun"] = 30.0
    out, _, _ = B.apply_missingness(df)
    assert not out.loc[0, "bun_locf"], "an observed value must not be flagged as carried"
    assert out.loc[1, "bun_locf"], "a carried value must be flagged"
    assert out["bun_locf"].sum() == 6


def test_locf_does_not_cross_encounter_blocks():
    df = _long(n_blocks=2)
    df["lactate"] = np.nan
    df.loc[df["encounter_block"].eq("b0") & df["window_idx"].eq(B.N_WINDOWS - 1),
           "lactate"] = 9.9
    out, _, _ = B.apply_missingness(df)
    b1 = out[out["encounter_block"] == "b1"]
    assert b1["lactate"].isna().all(), "a value must never carry into the next block"


# ------------------------------------------------------- absence means zero
def test_absence_means_zero_variables_are_zeroed_not_carried():
    df = _long(n_blocks=1)
    df["nee"] = np.nan
    df.loc[0, "nee"] = 0.4
    out, _, _ = B.apply_missingness(df)
    assert out.loc[1, "nee"] == 0.0, "an absent vasopressor record means 0, not the prior dose"
    assert "nee_locf" not in out.columns, "absence-means-zero variables are not LOCF-eligible"


def test_imv_status_absence_means_not_ventilated():
    df = _long(n_blocks=1)
    df["imv_status"] = np.nan
    df.loc[0, "imv_status"] = 1.0
    out, _, _ = B.apply_missingness(df)
    assert out.loc[5, "imv_status"] == 0.0


def test_missingness_is_counted_before_the_fill():
    df = _long(n_blocks=1)
    df["lactate"] = np.nan
    df.loc[0, "lactate"] = 4.2
    _, rep, _ = B.apply_missingness(df)
    row = rep[rep["variable"] == "lactate"].iloc[0]
    assert row["n_missing_pre_locf"] == B.N_WINDOWS - 1, (
        "counting after the fill would hide how much was carried"
    )


def test_missingness_denominator_is_alive_admitted_rows_only():
    df = _long(n_blocks=1)
    df.loc[df["window_idx"] >= 10, "alive_admitted"] = False
    df["lactate"] = np.nan
    _, rep, _ = B.apply_missingness(df)
    assert rep[rep["variable"] == "lactate"].iloc[0]["n_alive_admitted"] == 10


def test_report_separates_pre_locf_from_final_missingness():
    """The point of the report: how much was carried, and what is still absent."""
    df = _long(n_blocks=1)
    df["lactate"] = np.nan
    df.loc[0, "lactate"] = 4.2
    _, rep, _ = B.apply_missingness(df)
    r = rep[rep["variable"] == "lactate"].iloc[0]
    assert r["n_observed"] == 1
    assert r["n_missing_pre_locf"] == B.N_WINDOWS - 1
    assert r["n_filled_by_locf"] == 6
    assert r["n_missing_final"] == B.N_WINDOWS - 7
    assert r["n_observed"] + r["n_filled_by_locf"] + r["n_missing_final"] == r["n_alive_admitted"]


def test_report_shows_zero_by_rule_separately_from_missing():
    """An absent vasopressor record is a zero, not a gap, and must not read as one."""
    df = _long(n_blocks=1)
    df["nee"] = np.nan
    df.loc[0, "nee"] = 0.4
    _, rep, _ = B.apply_missingness(df)
    r = rep[rep["variable"] == "nee"].iloc[0]
    assert r["n_zero_by_rule"] == B.N_WINDOWS - 1
    assert r["n_missing_final"] == 0
    assert r["pct_missing_final"] == 0.0


def test_report_percentages_use_the_alive_admitted_denominator():
    df = _long(n_blocks=1)
    df.loc[df["window_idx"] >= 10, "alive_admitted"] = False
    df["bun"] = np.nan
    _, rep, _ = B.apply_missingness(df)
    r = rep[rep["variable"] == "bun"].iloc[0]
    assert r["n_alive_admitted"] == 10 and r["pct_missing_final"] == 100.0


def test_pattern_table_pools_small_cells():
    min_cell = B.COV["missing_values"]["reporting"]["pattern_min_cell"]
    df = _long(n_blocks=1)
    for v in ("bun", "lactate"):
        df[v] = 1.0
    df.loc[0, "bun"] = np.nan          # a pattern seen once
    _, _, pat = B.apply_missingness(df)
    assert len(pat), "a pattern table must be produced"
    pooled = pat[pat["pattern"].astype(str).str.startswith("<other")]
    assert len(pooled) == 1, f"rare patterns must be pooled below {min_cell}"


# ------------------------------------------------- oxygenation absence reasons
def _oxy(n_pao2, n_spo2, n_usable, oxy=np.nan):
    return pd.DataFrame([{"alive_admitted": True, "oxygenation": oxy,
                          "_n_pao2": n_pao2, "_n_spo2": n_spo2,
                          "_n_spo2_usable": n_usable}])


def _reason_counts(df):
    r = B.oxygenation_absence_reasons(df)
    return {row["variable"].split(": ", 1)[1]: row["n_missing_pre_locf"]
            for _, row in r.iterrows()}


def test_nothing_measured_is_attributed_to_the_measurement():
    c = _reason_counts(_oxy(0, 0, 0))
    assert c["no PaO2 and no SpO2 measured"] == 1
    assert sum(c.values()) == 1, "the three reasons must be mutually exclusive"


def test_spo2_on_the_plateau_is_attributed_to_the_ceiling_not_the_lookback():
    """SpO2 measured but all >= 97: the transform is undefined, not the FiO2 absent."""
    c = _reason_counts(_oxy(0, 5, 0))
    assert c[f"SpO2 present but all >= {B.SPO2_CEILING} (plateau)"] == 1
    assert c["usable measurement but no FiO2 within the lookback"] == 0


def test_a_usable_measurement_with_no_value_is_attributed_to_fio2():
    """A PaO2 exists and no oxygenation came out, so the FiO2 pairing failed."""
    c = _reason_counts(_oxy(3, 0, 0))
    assert c["usable measurement but no FiO2 within the lookback"] == 1
    c2 = _reason_counts(_oxy(0, 5, 5))
    assert c2["usable measurement but no FiO2 within the lookback"] == 1


def test_a_window_with_a_value_is_not_counted_as_absent():
    c = _reason_counts(_oxy(3, 5, 5, oxy=250.0))
    assert sum(c.values()) == 0


def test_plateau_is_distinguishable_from_no_measurement():
    """The flags must come from the UNFILTERED SpO2 series. Deriving both from the
    already-ceiling-filtered series makes n_spo2 == n_spo2_usable, so the plateau
    condition can never fire and every plateau window is misreported as 'nothing
    measured'. That produced a spurious 0.00% plateau on real data."""
    plateau_only = _oxy(0, 5, 0)          # SpO2 present, none below the ceiling
    nothing = _oxy(0, 0, 0)
    c1, c2 = _reason_counts(plateau_only), _reason_counts(nothing)
    key_p = f"SpO2 present but all >= {B.SPO2_CEILING} (plateau)"
    key_n = "no PaO2 and no SpO2 measured"
    assert c1[key_p] == 1 and c1[key_n] == 0
    assert c2[key_n] == 1 and c2[key_p] == 0, (
        "these two must not collapse into one another"
    )


def test_the_three_reasons_partition_the_absent_windows():
    df = pd.concat([_oxy(0, 0, 0), _oxy(0, 5, 0), _oxy(3, 0, 0),
                    _oxy(0, 5, 5), _oxy(2, 2, 2, oxy=300.0)], ignore_index=True)
    c = _reason_counts(df)
    assert sum(c.values()) == 4, "four absent windows, each counted exactly once"


# --------------------------------------------------------- outcome ascertainment
def test_extubation_is_a_gap_longer_than_the_success_window():
    """The rule is 'not followed by reintubation within 72h', so an IMV record
    whose next record is >72h later marks an extubation."""
    succ_h = B.CONFIG["outcomes"]["successful_extubation_hours"]
    hr = pd.Series([0.0, 4.0, 8.0, 200.0])
    nxt = hr.shift(-1)
    is_ext = nxt.isna() | ((nxt - hr) > succ_h)
    assert is_ext.tolist() == [False, False, True, True], (
        "only the record before the long gap, and the last record, are extubations"
    )


def test_a_brief_gap_is_reintubation_not_extubation():
    succ_h = B.CONFIG["outcomes"]["successful_extubation_hours"]
    hr = pd.Series([80.0, 80.0 + succ_h - 1])
    nxt = hr.shift(-1)
    is_ext = nxt.isna() | ((nxt - hr) > succ_h)
    assert is_ext.tolist() == [False, True], (
        "a reintubation inside the window means the first record is not an extubation"
    )


def test_outcome_cannot_be_read_from_the_window_grid():
    """The grid ends at the landmark, so nothing in it is post-landmark. This is
    the bug that made every block censored on the first real run."""
    T = B.CONFIG["cohort"]["landmark_hours"]
    last_window_start = (B.N_WINDOWS - 1) * B.WINDOW_H
    assert last_window_start < T, (
        f"the grid's last window starts at {last_window_start}h and the landmark is "
        f"{T}h, so `window_start_hr >= T` selects nothing -- outcomes must come "
        f"from the raw IMV series"
    )


# ------------------------------------------------------------- bolus units
def test_bolus_bound_nulls_an_implausible_mcg_per_kg_dose():
    """445 UCMC rows are charted mcg/kg with a median of 54.5, which converts to
    ~3,800 mcg for a 70 kg adult. Rather than guess the label is wrong, convert
    per the charted unit and let the bound decide."""
    lo, hi = B.OUTLIERS["med_bolus_mcg"]["fentanyl"]
    assert 54.5 * 70 > hi, "the implausible case must exceed the bound"
    assert 0.5 * 70 < hi, "a genuine low mcg/kg dose must survive it"


def test_bolus_bound_catches_the_observed_mcg_outlier():
    hi = B.OUTLIERS["med_bolus_mcg"]["fentanyl"][1]
    assert 22025.0 > hi, "the observed 22,025 mcg maximum must be nulled"
    assert 50.0 < hi, "the median 50 mcg bolus must survive"


# ------------------------------------------------------- clifpy waterfall casing
def test_waterfall_lowercasing_is_restored_to_mcide_casing():
    """clifpy's waterfall returns device_category as 'imv', not 'IMV'. Comparing
    against the canonical value then matches nothing and empties the cohort
    without raising -- which is exactly what happened on the first UCMC run."""
    df = pd.DataFrame({"device_category": ["imv", "room air", "nasal cannula",
                                           "high flow nc", None],
                       "mode_category": ["pressure control", None, None, None, None]})
    out = B._canonicalise_devices(df)
    assert out["device_category"].tolist()[:4] == [
        "IMV", "Room Air", "Nasal Cannula", "High Flow NC"]
    assert out["mode_category"].iloc[0] == "Pressure Control"
    assert pd.isna(out["device_category"].iloc[4])


def test_an_unknown_device_value_raises_rather_than_becoming_null():
    df = pd.DataFrame({"device_category": ["imv", "jet ventilator"]})
    try:
        B._canonicalise_devices(df)
    except SystemExit as e:
        assert "jet ventilator" in str(e)
    else:
        raise AssertionError("a value outside the schema must raise, not silently null")


def test_empty_cohort_raises_instead_of_crashing_on_nan():
    """Realistic path: every block is ventilated for less than the minimum, so the
    duration filter empties the cohort."""
    blocks = pd.DataFrame([{"encounter_block": "b0", "patient_id": "p0", "age": 60,
                            "block_discharge_dttm": pd.Timestamp("2026-01-05")}])
    anchor = pd.DataFrame([{"encounter_block": "b0",
                            "anchor_dttm": pd.Timestamp("2026-01-01"),
                            "first_imv_episode_hours": B.MIN_IMV_H - 1,
                            "first_episode_records": 2, "n_imv_episodes": 1}])
    try:
        B.build_cohort(blocks, anchor)
    except SystemExit as e:
        assert "empty" in str(e).lower() and "bug, not a" in str(e)
    else:
        raise AssertionError("an empty cohort must raise")


def test_the_anchor_merge_mismatch_asserts_rather_than_excluding():
    """The anchor comes from the same IMV series that selected the blocks, so a
    mismatch is impossible rather than an exclusion, and must not silently drop."""
    blocks = pd.DataFrame([{"encounter_block": "b0", "patient_id": "p0", "age": 60,
                            "block_discharge_dttm": pd.Timestamp("2026-01-05")}])
    empty = pd.DataFrame(columns=["encounter_block", "anchor_dttm",
                                  "first_imv_episode_hours",
                                  "first_episode_records", "n_imv_episodes"])
    try:
        B.build_cohort(blocks, empty)
    except AssertionError as e:
        assert "should be impossible" in str(e)
    except SystemExit:
        raise AssertionError("a merge mismatch must assert, not read as an exclusion")


# ------------------------------------------------------------- IMV episodes
_MAP = pd.DataFrame([{"hospitalization_id": "h0", "encounter_block": "b0"}])
_BASE = pd.Timestamp("2026-01-01")


def _resp(entries):
    """entries: (hour, device_category) -- the waterfalled series."""
    return pd.DataFrame({
        "hospitalization_id": ["h0"] * len(entries),
        "recorded_dttm": [_BASE + pd.Timedelta(hours=h) for h, _ in entries],
        "device_category": [d for _, d in entries]})


def _raw(hours):
    """The raw IMV record timestamps."""
    return pd.DataFrame({
        "hospitalization_id": ["h0"] * len(hours),
        "recorded_dttm": [_BASE + pd.Timedelta(hours=h) for h in hours]})


def _ep(resp_entries, raw_hours):
    return B.imv_episodes(_resp(resp_entries), _raw(raw_hours), _MAP)


def test_a_waterfall_transition_ends_the_episode():
    """100% of waterfalled IMV -> non-IMV transitions go to a real device, so a
    transition is a true extubation and ends the episode precisely."""
    out = _ep([(h, "IMV") for h in (0, 2, 4, 6)] + [(7, "Nasal Cannula")],
              [0, 2, 4, 6])
    assert out["first_imv_episode_hours"].iloc[0] == 6
    assert out["episode_ended_by"].iloc[0] == "observed transition to another device"


def test_a_raw_gap_ends_the_episode_when_no_transition_is_charted():
    """The waterfall carries IMV forward when nothing breaks the segment --
    measured at p90 30h and up to 1,137h past the last raw record, in 19.3% of
    hospitalizations. The raw-gap rule is what catches those."""
    g = B.EPISODE_GAP_H
    out = _ep([(h, "IMV") for h in (0, 2, 4, 4 + g + 5)], [0, 2, 4, 4 + g + 5])
    assert out["first_imv_episode_hours"].iloc[0] == 4
    assert out["n_imv_episodes"].iloc[0] == 2
    assert out["episode_ended_by"].iloc[0] == "no transition charted"


def test_normal_charting_gaps_do_not_fragment_an_episode():
    """Raw IMV gaps at UCMC are p95 5.00h and p98 7.07h, so a threshold at or
    below those would split continuously-ventilated patients."""
    assert B.EPISODE_GAP_H > 7.0, (
        f"episode gap {B.EPISODE_GAP_H}h sits inside the normal charting "
        f"distribution (p98 = 7.07h)"
    )
    out = _ep([(h, "IMV") for h in (0, 4.3, 8.6, 12)], [0, 4.3, 8.6, 12])
    assert out["n_imv_episodes"].iloc[0] == 1


def test_first_to_last_span_is_not_used_as_duration():
    """Span merges separate intubations weeks apart; the raw maximum was 24,063h."""
    out = _ep([(h, "IMV") for h in (0, 1, 500, 501)], [0, 1, 500, 501])
    assert out["first_imv_episode_hours"].iloc[0] == 1
    assert out["n_imv_episodes"].iloc[0] == 2


def test_a_single_imv_record_has_zero_duration():
    out = _ep([(0, "IMV")], [0])
    assert out["first_imv_episode_hours"].iloc[0] == 0
    assert 0 < B.MIN_IMV_H


def test_the_minimum_duration_is_one_analysis_window():
    assert B.MIN_IMV_H == B.WINDOW_H, (
        "the minimum is deliberately window_hours: a block that cannot fill one "
        "analysis window has no trajectory to model"
    )


# --------------------------------------------------------------------- SOFA
def _sofa_row(**kw):
    base = dict(map=85.0, platelet_count=250.0, bilirubin_total=0.5,
                creatinine=0.8, gcs_total=15.0, oxygenation=450.0, imv_status=1.0)
    base.update(kw)
    return B.score_sofa(pd.DataFrame([base])).iloc[0]


def test_healthy_patient_scores_zero():
    r = _sofa_row()
    assert r["sofa_total"] == 0 and r["sofa_n_components"] == 6


def test_absent_pressors_must_be_null_not_zero():
    """The trap: `epi <= 0.1` fires on a charted 0, giving every unpressored
    patient a cardiovascular score of 3."""
    assert _sofa_row()["sofa_cv"] == 0, "no pressor columns at all -> score on MAP"
    on_zero = _sofa_row(epinephrine_mcg_kg_min=0.0)
    assert on_zero["sofa_cv"] == 3, (
        "a charted 0 DOES trigger the <= 0.1 clause -- which is why the builder "
        "must leave absent pressors NULL"
    )


def test_low_dose_norepinephrine_scores_three():
    assert _sofa_row(norepinephrine_mcg_kg_min=0.05)["sofa_cv"] == 3


def test_high_dose_norepinephrine_scores_four():
    assert _sofa_row(norepinephrine_mcg_kg_min=0.5)["sofa_cv"] == 4


def test_respiratory_scores_off_the_ventilator_instead_of_returning_null():
    """clifpy and the epi repo both return NULL for P/F < 200 off IMV/NIPPV/CPAP."""
    off = _sofa_row(oxygenation=150.0, imv_status=0.0)
    assert not pd.isna(off["sofa_resp"]), "severe hypoxaemia off the vent must score"
    assert off["sofa_resp"] == 2
    on = _sofa_row(oxygenation=150.0, imv_status=1.0)
    assert on["sofa_resp"] == 3, "the same P/F on the vent scores higher"


def test_total_is_null_only_when_no_component_is_available():
    empty = B.score_sofa(pd.DataFrame([{}])).iloc[0]
    assert empty["sofa_n_components"] == 0
    assert pd.isna(empty["sofa_total"])


def test_partial_sofa_totals_what_it_has_and_says_how_many():
    r = _sofa_row(platelet_count=np.nan, gcs_total=np.nan)
    assert r["sofa_n_components"] == 4
    assert not pd.isna(r["sofa_total"]), (
        "a partial SOFA is summed, not nulled -- but n_components records the gap"
    )


def test_each_component_hits_its_documented_cutpoints():
    assert _sofa_row(platelet_count=19.0)["sofa_coag"] == 4
    assert _sofa_row(platelet_count=149.0)["sofa_coag"] == 1
    assert _sofa_row(bilirubin_total=12.0)["sofa_liver"] == 4
    assert _sofa_row(bilirubin_total=1.2)["sofa_liver"] == 1
    assert _sofa_row(creatinine=5.0)["sofa_renal"] == 4
    assert _sofa_row(creatinine=1.2)["sofa_renal"] == 1
    assert _sofa_row(gcs_total=5.0)["sofa_cns"] == 4
    assert _sofa_row(gcs_total=14.0)["sofa_cns"] == 1
    assert _sofa_row(oxygenation=399.0)["sofa_resp"] == 1


# ------------------------------------------------- the consumption assertion
def test_consumption_assertion_fires_on_a_declared_but_absent_column():
    df = _long(n_blocks=1)
    try:
        B.assert_config_is_honoured(df)
    except SystemExit as e:
        assert "silent no-op" in str(e)
    else:
        raise AssertionError("a declared covariate with no column must fail the build")


def test_consumption_assertion_passes_when_every_declared_column_exists():
    df = _long(n_blocks=1)
    for k in list(B.COV["time_varying"]) + list(B.COV["time_invariant"]):
        if not k.startswith("_"):
            df[k] = 0.0
    for c in B.EXPOSURE["columns"]:
        df[c] = 0.0
    B.assert_config_is_honoured(df)


# ------------------------------------------------------- units and transforms
def test_infusion_unit_conversion_covers_the_charted_units():
    """The analysis unit is mcg/hr (SG 2026-09-07). Every arm hand-computed,
    because a factor wrong by 60 or 1000 still produces a plausible dose."""
    m = pd.DataFrame({
        "med_dose": [100.0, 2.0, 0.1, 0.05, 2.0],
        "med_dose_unit": ["mcg/hr", "mcg/kg/hr", "mg/hr", "mcg/kg/min", "mcg/min"],
        "weight_kg": [50.0] * 5,
    })
    out = B._to_mcg_hr(m)
    assert out.iloc[0] == 100.0           # already the analysis unit: identity
    assert out.iloc[1] == 100.0           # 2 mcg/kg/hr * 50 kg
    assert out.iloc[2] == 100.0           # 0.1 mg/hr -> 100 mcg/hr, no weight
    assert abs(out.iloc[3] - 150.0) < 1e-9  # 0.05 mcg/kg/min * 60 * 50 kg
    assert abs(out.iloc[4] - 120.0) < 1e-9  # 2 mcg/min * 60


def test_the_dominant_charted_unit_needs_no_weight():
    """99.5% of UCMC fentanyl infusion rows are charted mcg/hr. That arm must
    survive a missing weight -- under the former mcg/kg/hr unit it did not, and
    a null weight silently removed the exposure."""
    m = pd.DataFrame({"med_dose": [75.0], "med_dose_unit": ["mcg/hr"],
                      "weight_kg": [float("nan")]})
    assert B._to_mcg_hr(m).iloc[0] == 75.0


def test_an_unhandled_infusion_unit_raises_rather_than_dropping():
    m = pd.DataFrame({"med_dose": [1.0], "med_dose_unit": ["mcg/kg/day"],
                      "weight_kg": [70.0]})
    try:
        B._to_mcg_hr(m)
    except SystemExit as e:
        assert "mcg/kg/day" in str(e)
    else:
        raise AssertionError("an unhandled dose unit must raise, not silently drop exposure")


def test_severinghaus_matches_the_reference_value():
    """Severinghaus(97) = 90.6 mmHg, per the reference implementation."""
    v = B._severinghaus(pd.Series([97.0]))
    assert abs(v.iloc[0] - 90.6) < 0.2, f"got {v.iloc[0]}"


def test_severinghaus_is_undefined_on_the_plateau():
    assert pd.isna(B._severinghaus(pd.Series([100.0])).iloc[0])


def test_a_record_charted_after_discharge_reaches_no_window():
    """A post-event window is structurally empty, not missing. Charting lag put
    428 IMV records into windows past discharge before this bound existed, which
    made `imv_status == 1` possible where alive_admitted was false -- an invariant
    every downstream consumer would otherwise have to re-derive for itself."""
    anchor = pd.Timestamp("2026-01-01 00:00", tz="UTC")
    cohort = pd.DataFrame([{"encounter_block": 1, "anchor_dttm": anchor,
                            "followup_end_dttm": anchor + pd.Timedelta(hours=10)}])
    rec = pd.DataFrame([
        {"encounter_block": 1, "recorded_dttm": anchor + pd.Timedelta(hours=1)},
        {"encounter_block": 1, "recorded_dttm": anchor + pd.Timedelta(hours=9.9)},
        {"encounter_block": 1, "recorded_dttm": anchor + pd.Timedelta(hours=10)},
        {"encounter_block": 1, "recorded_dttm": anchor + pd.Timedelta(hours=20)},
    ])
    out = B._to_windows(rec, cohort, "recorded_dttm")
    assert len(out) == 2, f"expected the two in-follow-up records, got {len(out)}"
    assert out["recorded_dttm"].max() < cohort["followup_end_dttm"].iloc[0]


def test_the_extubated_gap_rule_zeroes_dose_and_keeps_the_pre_gate_value():
    """design_notes.md §10a(a) was DECLARED in covariates.json and never applied
    until 2026-09-07. The gate must fire, must fire only on alive-admitted
    not-ventilated windows, and must leave total_dose_ungated recoverable --
    that column is what makes the decision reversible without another run."""
    df = pd.DataFrame([
        # ventilated, on a drip -> untouched
        {"encounter_block": 1, "alive_admitted": True,  "imv_status": 1.0,
         "inf_dose": 2.0, "bolus_dose": 0.0, "n_bolus": 0.0, "total_dose": 2.0},
        # extubated but still admitted, bolus given -> zeroed
        {"encounter_block": 1, "alive_admitted": True,  "imv_status": 0.0,
         "inf_dose": 0.0, "bolus_dose": 1.5, "n_bolus": 1.0, "total_dose": 1.5},
        # after discharge -> not alive_admitted, gate must not claim it
        {"encounter_block": 1, "alive_admitted": False, "imv_status": float("nan"),
         "inf_dose": 0.0, "bolus_dose": 0.0, "n_bolus": 0.0, "total_dose": 0.0},
    ])
    out = B.gate_dose_on_ventilation(df.copy())
    assert out.loc[0, "total_dose"] == 2.0, "a ventilated window must be untouched"
    assert out.loc[1, "total_dose"] == 0.0 and out.loc[1, "bolus_dose"] == 0.0
    assert out.loc[1, "total_dose_ungated"] == 1.5, (
        "the pre-gate value must survive, or reversing §10a(a) costs a re-run"
    )
    assert (out["total_dose_ungated"] == df["total_dose"]).all()


def test_the_gate_does_not_reach_windows_with_no_ventilation_record():
    """imv_status is NaN outside the at-risk period and 0 where the waterfall
    knows the patient was off the vent. Only alive_admitted rows may be gated,
    so a NaN outside follow-up cannot silently zero a real dose."""
    df = pd.DataFrame([
        {"encounter_block": 1, "alive_admitted": False, "imv_status": float("nan"),
         "inf_dose": 3.0, "bolus_dose": 0.0, "n_bolus": 0.0, "total_dose": 3.0},
    ])
    out = B.gate_dose_on_ventilation(df.copy())
    assert out.loc[0, "total_dose"] == 3.0


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
