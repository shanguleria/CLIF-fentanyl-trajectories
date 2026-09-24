"""Phase 0 -- build the analytic tables from CLIF.

Outputs (PHI, written to output/intermediate_phi/):
    trajectory_long.parquet     one row per encounter block per window
    time_to_event.parquet       one row per encounter block
    hospital_intervals.parquet  one row per ADT interval

Protocol: config/covariates.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "code"))

from clifpy import (  # noqa: E402
    Adt, CrrtTherapy, Hospitalization, HospitalDiagnosis, Labs,
    MedicationAdminContinuous, MedicationAdminIntermittent, Patient,
    PatientAssessments, RespiratorySupport, Vitals, stitch_encounters,
)
from clifpy.utils.comorbidity import calculate_cci  # noqa: E402
from utils.doses import convert as convert_doses  # noqa: E402
from utils.fio2 import normalize_fio2  # noqa: E402
from utils.outliers import (  # noqa: E402
    apply_long, apply_med_converted, apply_med_raw, fentanyl_charted_max_mcg_hr,
    fentanyl_sanity_ceiling,
    load_config as load_outliers, nee_sanity_ceiling,
)
from utils.waterfall_cache import (  # noqa: E402
    cache_key, describe as describe_cache, load as cache_load,
    store as cache_store,
)
from utils.strobe import render_png, render_text  # noqa: E402
from utils.paths import (  # noqa: E402
    clear_owned_outputs, phase_dir, provenance, site_dirs, write_manifest,
)

CONFIG = json.loads((REPO / "config" / "config.json").read_text())
COV = json.loads((REPO / "config" / "covariates.json").read_text())
OUTLIERS = load_outliers()

GRID = COV["windows"]["granular"]
WINDOW_H = GRID["width_hours"]
EXTENT_H = GRID["extent_hours"]
N_WINDOWS = GRID["n_windows"]
STITCH_H = COV["encounter_blocks"]["stitch_time_interval_hours"]
MIN_IMV_H = CONFIG["cohort"]["min_imv_hours"]

EXPOSURE = COV["exposure"]
INF_HOLD_H = EXPOSURE["infusion"]["hold_hours"]

NEE = COV["time_varying"]["nee"]
NEE_COEF = NEE["coefficients"]
NEE_HOLD_H = NEE["hold_hours"]

OXY = COV["time_varying"]["oxygenation"]
FIO2_LOOKBACK_H = OXY["fio2_lookback_hours"]
SPO2_CEILING = OXY["spo2_ceiling"]
RA_FIO2 = OXY["fio2_scale"]["fraction_band"][0]
EPISODE_GAP_H = CONFIG["cohort"]["imv_episode_gap_hours"]

LAB_VARS = {
    k: v["source"]["category"]
    for k, v in COV["time_varying"].items()
    if not k.startswith("_") and (v.get("source") or {}).get("table") == "labs"
}
# Same shape as LAB_VARS, and for the same reason: the set of assessment
# covariates is whatever the config declares, so adding a score is a config edit
# rather than an edit in two places that can disagree.
ASSESS_VARS = {
    k: v["source"]["category"]
    for k, v in COV["time_varying"].items()
    if not k.startswith("_") and (v.get("source") or {}).get("table") == "assessments"
}
ZERO_VARS = COV["missing_values"]["absence_means_zero"]["members"]
NOT_VENT_VARS = COV["missing_values"]["absence_means_not_ventilated"]["members"]
SOFA_INPUT_CAPS = COV["missing_values"]["sofa_inputs"]["variables"]

# Everything this script owns. Cleared before it runs so a crash cannot leave a stale file
# The shareable tree is subdivided by the script that produced each file; the
# phaseN_ prefixes stay, so the folder says which script and the prefix says
# which phase. manifest.json is the exception and lives at the ROOT of
# out_final: it is not a Phase 0 result but the pipeline's staleness marker,
# which every later phase consults through require_manifest().
PHASE_DIR = "01_cohort"

OWNED = {
    "out_phi": ["trajectory_long.parquet", "trajectory_long.csv",
                "time_to_event.parquet", "time_to_event.csv",
                "hospital_intervals.parquet",
                "exemplar_series.parquet", "exemplar_meta.json",
                "exemplar_id.txt"],
    "out_final": ["manifest.json"],
    "phase": ["strobe.csv", "strobe.txt", "strobe.png",
              "provenance.json", "exemplar_selection.csv"],
    "diagnostics": ["missingness.csv", "missingness_patterns.csv",
                    "diagnostics.csv"],
}

# Paths this script used to write and no longer does. Cleared so a relocation
# cannot leave a stale twin the owned list no longer names.
# Two rounds of retirement. The prefixes are written out with an f-string rather
# than inline so a future find-and-replace on output names cannot silently strip
# them -- which is exactly what happened here on 2026-09-23.
_P0 = "phase0_"
RETIRED_OUTPUTS = [
    # flat locations retired 2026-09-08 when out_final was subdivided by script
    f"output/final_no_phi/{_P0}missingness.csv",
    f"output/final_no_phi/{_P0}missingness_patterns.csv",
    f"output/final_no_phi/{_P0}diagnostics.csv",
    f"output/final_no_phi/{_P0}strobe.csv",
    f"output/final_no_phi/{_P0}strobe.txt",
    f"output/final_no_phi/{_P0}strobe.png",
    f"output/final_no_phi/{_P0}provenance.json",
    f"output/final_no_phi/diagnostics/{_P0}missingness.csv",
    f"output/final_no_phi/diagnostics/{_P0}missingness_patterns.csv",
    f"output/final_no_phi/diagnostics/{_P0}diagnostics.csv",
    # the phase0_ prefix itself retired 2026-09-23 -- the phases are gone
    f"output/final_no_phi/{_P0}manifest.json",
    f"output/final_no_phi/01_cohort/{_P0}strobe.csv",
    f"output/final_no_phi/01_cohort/{_P0}strobe.txt",
    f"output/final_no_phi/01_cohort/{_P0}strobe.png",
    f"output/final_no_phi/01_cohort/{_P0}provenance.json",
    f"output/final_no_phi/01_cohort/diagnostics/{_P0}missingness.csv",
    f"output/final_no_phi/01_cohort/diagnostics/{_P0}missingness_patterns.csv",
    f"output/final_no_phi/01_cohort/diagnostics/{_P0}diagnostics.csv",
]

STROBE: list[tuple[str, int]] = []
FLOW: list[dict] = []


def note(label: str, n: int) -> None:
    """A diagnostic count. Does not enter the cohort flow."""
    STROBE.append((label, n))
    print(f"  {label:.<52} {n:,}")


def flow(label: str, n_after: int, reason: str = "") -> None:
    """One row of the STROBE cohort flow. Exclusions are derived, not asserted."""
    n_before = FLOW[-1]["n_after"] if FLOW else None
    n_excl = None if n_before is None else n_before - n_after
    FLOW.append({"step": label, "n_before": n_before, "n_excluded": n_excl,
                 "reason": reason, "n_after": n_after})
    if n_excl:
        print(f"  {label:.<52} {n_after:>9,}   (-{n_excl:,}: {reason})")
    else:
        print(f"  {label:.<52} {n_after:>9,}")


# --------------------------------------------------------------------- loading
# Only the categories the analysis actually uses are read. Pushed down to the
# parquet read, this is the difference between 9 lab categories and ~50.
LAB_NEEDED = sorted(set(LAB_VARS.values()) | {"po2_arterial", "creatinine",
                                              "platelet_count"})
VITAL_NEEDED = ["spo2", "map", "weight_kg", "height_cm"]
# gcs_total is unioned in explicitly: it is a SOFA INPUT carried by
# _sofa_inputs() with its own cap, not a declared time_varying covariate, so
# it is absent from ASSESS_VARS and would be dropped from the read filter.
ASSESS_NEEDED = sorted(set(ASSESS_VARS.values()) | {"gcs_total"})


def _kw(**extra) -> dict:
    return dict(
        data_directory=CONFIG["data_directory"],
        filetype=CONFIG["filetype"],
        timezone=CONFIG["timezone"],
        output_directory=str(REPO / "logs"),
        **extra,
    )


def _mac_categories() -> list[str]:
    meds = CONFIG["medications"]
    return sorted(set(meds["opioid_infusion_categories"]) | set(meds["other_sedative_categories"])
                  | set(NEE_COEF) | set(SOFA_PRESSORS))


def load_core() -> dict:
    """Tables small enough to load whole, plus the IMV screen."""
    t = {
        "patient": Patient.from_file(**_kw()).df,
        "hospitalization": Hospitalization.from_file(**_kw()).df,
        "adt": Adt.from_file(**_kw()).df,
    }
    for name, df in t.items():
        print(f"  {name:.<52} {len(df):,} rows")

    screen = RespiratorySupport.from_file(
        **_kw(columns=["hospitalization_id", "device_category", "recorded_dttm"])).df
    imv = screen[screen["device_category"] == "IMV"]
    raw_anchor = (imv.groupby("hospitalization_id")["recorded_dttm"].min()
                     .rename("imv_dttm").reset_index())
    print(f"  {'resp_support IMV screen':.<52} {len(raw_anchor):,} hospitalizations")
    # The IMV series is kept whole: outcomes are ascertained after the trajectory
    # window, so they cannot use the span-trimmed waterfall.
    return t, raw_anchor, imv[["hospitalization_id", "recorded_dttm"]].copy()


def load_cohort_tables(hosp_ids: list[str]) -> dict:
    """Everything else, filtered to the cohort's hospitalizations and categories.

    Pushed down to the read: an unfiltered site-wide load makes the respiratory
    waterfall alone take hours.
    """
    hid = {"hospitalization_id": hosp_ids}
    t = {
        "hospital_diagnosis": HospitalDiagnosis.from_file(**_kw(filters=hid)).df,
        "crrt": CrrtTherapy.from_file(**_kw(filters=hid)).df,
        "labs": Labs.from_file(
            **_kw(filters={**hid, "lab_category": LAB_NEEDED})).df,
        "vitals": Vitals.from_file(
            **_kw(filters={**hid, "vital_category": VITAL_NEEDED})).df,
        "assessments": PatientAssessments.from_file(
            **_kw(filters={**hid, "assessment_category": ASSESS_NEEDED})).df,
        "mac": MedicationAdminContinuous.from_file(
            **_kw(filters={**hid, "med_category": _mac_categories()})).df,
        "mai": MedicationAdminIntermittent.from_file(
            **_kw(filters={**hid, "med_category": CONFIG["medications"]["opioid_bolus_categories"]})).df,
    }
    key = cache_key(CONFIG, [_canonicalise_devices])
    print(f"  {describe_cache(key)}")
    cached, missing = cache_load(key, hosp_ids)
    n_hit = 0 if cached is None else int(cached["hospitalization_id"].nunique())
    print(f"  waterfall: {n_hit:,} hospitalizations from cache, "
          f"{len(missing):,} to compute")

    frames = [] if cached is None else [cached]
    if missing:
        rs = RespiratorySupport.from_file(
            **_kw(filters={"hospitalization_id": missing}))
        print(f"  {'respiratory_support to waterfall':.<52} {len(rs.df):,} rows")
        if len(rs.df):
            fresh = _canonicalise_devices(
                rs.waterfall(verbose=False, return_dataframe=True))
            frames.append(fresh)
        else:
            # These hospitalizations sit in an IMV block because a sibling was
            # ventilated; they have no respiratory_support rows of their own.
            fresh = pd.DataFrame(columns=["hospitalization_id"])
            print(f"  {len(missing):,} hospitalization(s) have no "
                  f"respiratory_support rows; recorded as covered")
        total = cache_store(key, fresh, attempted=missing)
        print(f"  waterfall cached: {total:,} rows total")
    t["resp"] = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]

    for name, df in t.items():
        if name != "resp_raw":
            print(f"  {name:.<52} {len(df):,} rows")
    return t


def _canonicalise_devices(df: pd.DataFrame) -> pd.DataFrame:
    """Restore mCIDE casing for the category columns the waterfall lowercases.

    clifpy's waterfall returns device_category as 'imv'/'room air' rather than
    the schema's 'IMV'/'Room Air'. Comparing against the canonical value then
    matches nothing and empties the cohort without raising.
    """
    import yaml

    sch = yaml.safe_load(
        (Path(__import__("clifpy").__file__).parent / "schemas"
         / "respiratory_support_schema.yaml").read_text())
    out = df.copy()
    for col in ("device_category", "mode_category"):
        if col not in out.columns:
            continue
        allowed = next((c.get("permissible_values", []) for c in sch["columns"]
                        if c["name"] == col), [])
        lut = {v.lower(): v for v in allowed}
        lowered = out[col].astype("string").str.lower().str.strip()
        mapped = lowered.map(lut)
        unmapped = lowered.notna() & mapped.isna()
        if unmapped.any():
            raise SystemExit(
                f"{col}: values the schema does not list: "
                f"{sorted(lowered[unmapped].unique())[:10]}"
            )
        out[col] = mapped
    return out


def prepare_long_tables(t: dict, mapping: pd.DataFrame) -> dict:
    """Attach encounter_block and apply outlier bounds once per table.

    Every consumer previously did this itself: vitals was merged and bounded five
    times over 10.5M rows, labs three times over 3.8M. apply_long is idempotent,
    so the repeats were pure cost.
    """
    spec = {
        "labs": ("labs", "lab_category", "lab_value_numeric"),
        "vitals": ("vitals", "vital_category", "vital_value"),
        "assessments": ("patient_assessments", "assessment_category",
                        "numerical_value"),
    }
    for name, (table, cat, val) in spec.items():
        d = t[name].merge(mapping, on="hospitalization_id", how="inner")
        d, rep = apply_long(d, table, cat, val, config=OUTLIERS)
        print(rep)
        t[name] = d
    for name in ("resp", "crrt"):
        t[name] = t[name].merge(mapping, on="hospitalization_id", how="inner")
    return t


def assert_categories_present(t: dict) -> None:
    """Fail loudly if a configured med_category is absent; print what is there."""
    meds = CONFIG["medications"]
    for table, key in (("mac", "opioid_infusion_categories"),
                       ("mai", "opioid_bolus_categories")):
        present = set(t[table]["med_category"].dropna().unique())
        wanted = set(meds[key])
        if not wanted <= present:
            raise SystemExit(
                f"{key}: {sorted(wanted - present)} absent from {table}. "
                f"Present: {sorted(present)[:40]}"
            )
    others = set(meds["other_sedative_categories"]) - set(t["mac"]["med_category"].dropna())
    if others:
        print(f"  NOTE other_sedative_categories absent from mac: {sorted(others)}")


# ------------------------------------------------------------ blocks + cohort
def build_blocks(t: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    _, _, mapping = stitch_encounters(t["hospitalization"], t["adt"], time_interval=STITCH_H)
    flow("Hospitalizations in the CLIF extract", len(mapping))
    flow("Encounter blocks after stitching", mapping["encounter_block"].nunique(),
         "merged into an existing block within "
         f"{STITCH_H}h")

    hb = t["hospitalization"].merge(mapping, on="hospitalization_id", how="left",
                                    validate="one_to_one")
    assert hb["encounter_block"].notna().all(), "hospitalization with no block"

    hb = hb.sort_values(["encounter_block", "admission_dttm"], kind="stable")
    blocks = hb.groupby("encounter_block").agg(
        patient_id=("patient_id", "first"),
        n_hospitalizations=("hospitalization_id", "size"),
        block_admission_dttm=("admission_dttm", "min"),
        block_discharge_dttm=("discharge_dttm", "max"),
        age=("age_at_admission", "first"),
    ).reset_index()
    assert blocks["n_hospitalizations"].sum() == len(hb), "rows lost stitching"

    last = (hb.sort_values(["encounter_block", "discharge_dttm"], kind="stable")
              .groupby("encounter_block").last()[["discharge_category"]].reset_index())
    blocks = blocks.merge(last, on="encounter_block", how="left", validate="one_to_one")
    assert blocks["discharge_category"].notna().all(), "block lost its disposition"
    return blocks, mapping


def hospital_intervals(t: dict, mapping: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in ("hospital_id", "hospital_type", "in_dttm", "out_dttm")
            if c in t["adt"].columns]
    if "hospital_id" not in cols:
        raise SystemExit("adt has no hospital_id; stitch_encounters requires it")
    hi = (t["adt"][["hospitalization_id"] + cols]
          .merge(mapping, on="hospitalization_id", how="inner")
          .sort_values(["encounter_block", "in_dttm"], kind="stable")
          .reset_index(drop=True))
    note("distinct hospital_id", hi["hospital_id"].nunique())
    return hi


def hospital_endpoints(hi: pd.DataFrame) -> pd.DataFrame:
    first = (hi.sort_values(["encounter_block", "in_dttm"], kind="stable")
               .groupby("encounter_block")["hospital_id"].first()
               .rename("hospital_id_admission"))
    last = (hi.sort_values(["encounter_block", "out_dttm"], kind="stable")
              .groupby("encounter_block")["hospital_id"].last()
              .rename("hospital_id_discharge"))
    ends = pd.concat([first, last], axis=1).reset_index()
    n_moved = int((ends["hospital_id_admission"] != ends["hospital_id_discharge"]).sum())
    note("blocks transferring between hospitals", n_moved)
    return ends


def imv_episodes(resp: pd.DataFrame, imv_raw: pd.DataFrame,
                 mapping: pd.DataFrame) -> pd.DataFrame:
    """Continuous IMV episodes per encounter block.

    An episode ends at whichever comes first:
      - a waterfalled IMV -> non-IMV transition. Measured on real data, 100% of
        these go to a real device (NIPPV, nasal cannula, trach collar, ...) and
        none to a null device, so a transition is a true extubation.
      - a gap between RAW IMV records longer than imv_episode_gap_hours. The
        waterfall alone is not sufficient: where no subsequent device is ever
        charted, nothing breaks the segment and it carries IMV forward -- measured
        at p90 30h and up to 1,137h past the last raw record, in 19.3% of
        hospitalizations. Those extubations are invisible to it.
    """
    wf = (resp[["hospitalization_id", "recorded_dttm", "device_category"]]
          .merge(mapping, on="hospitalization_id", how="inner")
          .sort_values(["encounter_block", "recorded_dttm"], kind="stable"))
    wf["is_imv"] = wf["device_category"].eq("IMV")

    imv = wf[wf["is_imv"]].copy()
    # end of the current IMV run in the waterfalled series
    nxt_non_imv = (wf[~wf["is_imv"]][["encounter_block", "recorded_dttm"]]
                   .rename(columns={"recorded_dttm": "transition_dttm"})
                   .sort_values("transition_dttm", kind="stable"))
    imv = pd.merge_asof(imv.sort_values("recorded_dttm"), nxt_non_imv,
                        left_on="recorded_dttm", right_on="transition_dttm",
                        by="encounter_block", direction="forward")

    raw = (imv_raw.merge(mapping, on="hospitalization_id", how="inner")
                  .sort_values(["encounter_block", "recorded_dttm"], kind="stable"))
    raw_gap = (raw.groupby("encounter_block")["recorded_dttm"].diff()
                  .dt.total_seconds() / 3600)
    raw["_break_before"] = raw_gap.isna() | (raw_gap > EPISODE_GAP_H)
    raw["_episode"] = raw.groupby("encounter_block")["_break_before"].cumsum().astype(int)

    first = raw[raw["_episode"] == 1]
    out = first.groupby("encounter_block").agg(
        anchor_dttm=("recorded_dttm", "min"),
        raw_episode_end=("recorded_dttm", "max"),
        first_episode_records=("recorded_dttm", "size"),
    ).reset_index()
    out["n_imv_episodes"] = (raw.groupby("encounter_block")["_episode"].max()
                                .reindex(out["encounter_block"]).to_numpy())

    # the waterfall transition, if one occurs before the raw gap closes the episode
    trans = (imv.groupby("encounter_block")["transition_dttm"].min()
                .rename("wf_transition").reset_index())
    out = out.merge(trans, on="encounter_block", how="left")
    # The episode ENDS at the last IMV record; the transition is what tells us it
    # ended rather than continued. A transition charted within one gap-threshold of
    # that last record means the extubation was observed.
    out["first_episode_end"] = out["raw_episode_end"]
    out["first_imv_episode_hours"] = (
        (out["first_episode_end"] - out["anchor_dttm"]).dt.total_seconds() / 3600)
    lag_h = ((out["wf_transition"] - out["raw_episode_end"]).dt.total_seconds() / 3600)
    out["episode_ended_by"] = np.where(
        out["wf_transition"].notna() & (lag_h >= 0) & (lag_h <= EPISODE_GAP_H),
        "observed transition to another device", "no transition charted")

    note("blocks with more than one IMV episode", int((out["n_imv_episodes"] > 1).sum()))
    for k, v in out["episode_ended_by"].value_counts().items():
        note(f"  first episode ended by {k}", int(v))
    q = out["first_imv_episode_hours"].quantile([.25, .5, .75])
    print(f"  first IMV episode hours: median {q[.5]:.1f} "
          f"(IQR {q[.25]:.1f}-{q[.75]:.1f}), max {out['first_imv_episode_hours'].max():.1f}")
    return out[["encounter_block", "anchor_dttm", "first_imv_episode_hours",
                "first_episode_records", "n_imv_episodes", "episode_ended_by"]]


def build_cohort(blocks: pd.DataFrame, anchor: pd.DataFrame) -> pd.DataFrame:
    c = blocks.merge(anchor, on="encounter_block", how="inner", validate="one_to_one")
    # Not a flow step: the anchor comes from the same IMV series that selected these
    # blocks, so this can only ever be equal. It is a consistency check.
    assert len(c) == len(blocks), (
        f"anchor merge lost {len(blocks) - len(c)} blocks; the IMV screen and the "
        f"episode builder disagree, which should be impossible"
    )
    if c.empty:
        raise SystemExit(
            "no block has an intubation anchor. An empty cohort is a bug, not a "
            "finding -- check that device_category still carries its mCIDE casing."
        )
    c = c[c["first_imv_episode_hours"] >= MIN_IMV_H]
    flow(f"Ventilated at least {MIN_IMV_H}h (one analysis window)", len(c),
         f"first continuous IMV episode shorter than {MIN_IMV_H}h")

    c = c[c["age"] >= CONFIG["cohort"]["min_age"]]
    flow(f"Adult blocks (age >= {CONFIG['cohort']['min_age']})", len(c),
         f"age < {CONFIG['cohort']['min_age']} or age missing")

    c["followup_end_dttm"] = c[["block_discharge_dttm"]].min(axis=1)
    c = c[c["followup_end_dttm"] > c["anchor_dttm"]]
    flow("Blocks with follow-up after the anchor", len(c),
         "discharged or died at or before the anchor")

    c = c.sort_values("encounter_block", kind="stable").reset_index(drop=True)
    c["id_num"] = np.arange(1, len(c) + 1)
    if c.empty:
        raise SystemExit(
            "cohort is empty. Check the STROBE counts above for the step that "
            "dropped everything -- an empty cohort is a bug, not a finding."
        )
    per_patient = c.groupby("patient_id").size()
    note("patients contributing more than one block", int((per_patient > 1).sum()))
    note("blocks per patient (max)", int(per_patient.max()))
    return c


def window_grid(cohort: pd.DataFrame) -> pd.DataFrame:
    idx = np.arange(N_WINDOWS)
    g = cohort.loc[:, ["encounter_block", "patient_id", "id_num", "anchor_dttm",
                       "followup_end_dttm"]]
    g = g.loc[g.index.repeat(N_WINDOWS)].reset_index(drop=True)
    g["window_idx"] = np.tile(idx, len(cohort))
    g["window_start_hr"] = g["window_idx"] * WINDOW_H
    g["win_start"] = g["anchor_dttm"] + pd.to_timedelta(g["window_start_hr"], unit="h")
    g["win_end"] = g["win_start"] + pd.Timedelta(hours=WINDOW_H)
    g["alive_admitted"] = g["win_start"] < g["followup_end_dttm"]
    note("patient-window rows", len(g))
    note("alive-admitted rows", int(g["alive_admitted"].sum()))
    return g


# -------------------------------------------------------------------- exposure
def _hourly_scaffold(cohort: pd.DataFrame) -> pd.DataFrame:
    hours = np.arange(0, EXTENT_H)
    s = cohort.loc[:, ["encounter_block", "anchor_dttm", "followup_end_dttm"]]
    s = s.loc[s.index.repeat(len(hours))].reset_index(drop=True)
    s["hr"] = np.tile(hours, len(cohort))
    s["cell_dttm"] = s["anchor_dttm"] + pd.to_timedelta(s["hr"], unit="h")
    return s[s["cell_dttm"] < s["followup_end_dttm"]].copy()


def infusion_grid(t: dict, mapping: pd.DataFrame, cohort: pd.DataFrame,
                  cats: list[str], rate_fn) -> pd.DataFrame:
    """Hourly grid of infusion rate for `cats`. Config: exposure.infusion.

    rate_fn maps the charted (dose, unit, weight) to the target rate; fentanyl
    uses _to_mcg_hr, the sedatives go through the dose_units table.
    """
    m = t["mac"][t["mac"]["med_category"].isin(cats)].merge(
        mapping, on="hospitalization_id", how="inner")
    m = m[m["encounter_block"].isin(cohort["encounter_block"])]

    m, rep = apply_med_raw(m, "med_category", "med_dose", "med_dose_unit", config=OUTLIERS)
    print(rep)

    m = m.merge(cohort[["encounter_block", "anchor_dttm", "weight_kg"]],
                on="encounter_block", how="left")
    m["rate"] = rate_fn(m)
    # A charted stop is a rate of zero, not a missing value -- and so is anything
    # the nurse recorded as NOT ADMINISTERED. Every `stop` at UCMC is already
    # `not_administered`, but the two are separate CLIF columns and a site may
    # not pair them, so both are checked. Measured 2026-09-24: 802 rows carried a
    # positive rate while charted not_administered and were being counted as drug
    # given at full rate.
    act = m.get("mar_action_category", pd.Series(index=m.index, dtype=object))
    grp = m.get("mar_action_group", pd.Series(index=m.index, dtype=object))
    not_given = ((act.astype("string").str.lower() == "stop")
                 | (grp.astype("string").str.lower() == "not_administered"))
    n_zeroed = int((not_given & (m["rate"] > 0)).sum())
    m.loc[not_given, "rate"] = 0.0
    if n_zeroed:
        print(f"  infusion: {n_zeroed:,} record(s) charted stop/not_administered "
              f"carried a positive rate and are zeroed")

    m["hr"] = ((m["admin_dttm"] - m["anchor_dttm"]).dt.total_seconds() // 3600).astype("Int64")
    # A record inside the window whose converted rate is NaN disappears here.
    # It is not necessarily wrong -- a dose charted with no value cannot be used
    # -- but a silent drop of exposure is exactly what this pipeline reports
    # everywhere else, so it is counted out loud.
    in_window = (m["hr"] >= 0) & (m["hr"] < EXTENT_H)
    n_nan = int((in_window & m["rate"].isna()).sum())
    if n_nan:
        print(f"  infusion: {n_nan:,} of {int(in_window.sum()):,} in-window record(s) "
              f"carry no convertible rate and are dropped")
    m = m[in_window & m["rate"].notna()]

    # THE LAST charted rate in the hour, which means sorting on TIME. Until
    # 2026-09-24 the sort key was ["encounter_block", "hr", "rate"] -- no time
    # component at all -- so `.last()` returned the HIGHEST rate in the hour.
    # The two rules agree on a single-record hour and on an uptitration, and
    # diverge wherever the rate falls inside one. A charted stop is exactly that
    # case, so a stop was invisible whenever anything else was charted in the
    # same hour. Measured before the fix: 9,312 cells (21.3% of multi-record
    # cells) disagreed, the error was ALWAYS upward, median 50 mcg/hr, and 35.0%
    # of episodes had at least one affected window. covariates.json ->
    # exposure.infusion.algorithm had declared the correct rule all along.
    last = (m.sort_values(["encounter_block", "hr", "admin_dttm", "rate"], kind="stable")
              .groupby(["encounter_block", "hr"], as_index=False)["rate"].last())

    grid = _hourly_scaffold(cohort).merge(last, on=["encounter_block", "hr"], how="left")
    grid = grid.sort_values(["encounter_block", "hr"], kind="stable")
    grid["rate"] = grid.groupby("encounter_block")["rate"].ffill(limit=INF_HOLD_H)

    # Never carry past the block's last charted record.
    span = m.groupby("encounter_block")["hr"].max().rename("last_hr")
    grid = grid.merge(span, on="encounter_block", how="left")
    grid.loc[grid["hr"] > grid["last_hr"].fillna(-1), "rate"] = np.nan
    grid["rate"] = grid["rate"].fillna(0.0)          # absence_means_zero
    return grid[["encounter_block", "hr", "rate"]]


def _to_mcg_hr(m: pd.DataFrame) -> pd.Series:
    """Convert charted fentanyl infusion dose to mcg/hr using the charted unit.

    mcg/hr is how fentanyl is ordered, and it is what 99.5% of UCMC rows are
    already charted in, so the dominant path is an identity and weight never
    touches it. Weight is needed only for the weight-based minority.
    """
    unit = m["med_dose_unit"].astype("string").str.lower().str.strip()
    dose = pd.to_numeric(m["med_dose"], errors="coerce")
    w = pd.to_numeric(m["weight_kg"], errors="coerce")
    out = pd.Series(np.nan, index=m.index, dtype="float64")
    out[unit == "mcg/hr"] = dose
    out[unit == "mcg/kg/hr"] = dose * w
    out[unit == "mg/hr"] = dose * 1000.0
    out[unit == "mcg/kg/min"] = dose * 60.0 * w
    out[unit == "mcg/min"] = dose * 60.0
    unresolved = dose.notna() & out.isna()
    if unresolved.any():
        bad = sorted(unit[unresolved].dropna().unique())
        raise SystemExit(
            f"fentanyl infusion: {int(unresolved.sum()):,} rows in unhandled units {bad}. "
            f"Add them to _to_mcg_hr rather than dropping -- a dropped unit removes "
            f"exposure for whichever patients were charted that way."
        )
    return out


def gate_dose_on_ventilation(long: pd.DataFrame) -> pd.DataFrame:
    """Force dose to 0 in alive-admitted windows the patient was not ventilated in.

    Applied after status_covariates because it needs
    imv_status. total_dose_ungated preserves the pre-gate value, so reversing the
    decision is a column swap rather than another run.
    """
    off = long["alive_admitted"] & (long["imv_status"].fillna(0) == 0)
    long["total_dose_ungated"] = long["total_dose"]
    hit = off & (long["total_dose"] > 0)
    note("windows zeroed by the extubated-gap rule", int(hit.sum()))
    if hit.any():
        print(f"    across {long.loc[hit, 'encounter_block'].nunique():,} blocks; "
              f"median {long.loc[hit, 'total_dose'].median():.2f} mcg/hr, "
              f"{int((long.loc[hit, 'inf_dose'] == 0).sum()):,} bolus-only")
    for c in ("inf_dose", "inf_mcg", "bolus_dose", "bolus_mcg", "n_bolus",
              "total_dose", "window_mcg"):
        long.loc[off, c] = 0.0
    return long


def _to_target_rate(m: pd.DataFrame) -> pd.Series:
    """Convert a sedative infusion to its own target unit. Config: dose_units."""
    rate, rep = convert_doses(m, "med_category", "med_dose", "med_dose_unit")
    print(f"  {rep}")
    return rate


def sedative_exposure(t: dict, mapping: pd.DataFrame, cohort: pd.DataFrame,
                      windows: pd.DataFrame) -> pd.DataFrame:
    """propofol_dose and midazolam_dose per window. Config: exposure.sedatives."""
    spec = COV["exposure"]["sedatives"]
    out = windows[["encounter_block", "window_idx"]].copy()
    for cat in CONFIG["medications"]["other_sedative_categories"]:
        col = f"{cat}_dose"
        if col not in spec["columns"]:
            raise SystemExit(
                f"{cat} is in config.json other_sedative_categories but "
                f"{col} is not declared in covariates.json exposure.sedatives.columns"
            )
        print(f"  {cat} -> {col} ({spec['units'][col]})")
        grid = infusion_grid(t, mapping, cohort, [cat], _to_target_rate)
        grid["window_idx"] = grid["hr"] // WINDOW_H
        w = (grid.groupby(["encounter_block", "window_idx"], as_index=False)["rate"]
                 .mean().rename(columns={"rate": col}))
        out = out.merge(w, on=["encounter_block", "window_idx"], how="left")
        out[col] = out[col].fillna(0.0)          # absence_means_zero
        vent = out[col] > 0
        note(f"windows with any {cat}", int(vent.sum()))
        if vent.any():
            print(f"    median {out.loc[vent, col].median():.2f}  "
                  f"max {out[col].max():.2f} {spec['units'][col]}")
    return out


def _bolus_events(t: dict, mapping: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    """Individual bolus administrations: exact admin_dttm, converted mcg, bounded.

    Split out of bolus_doses() on 2026-09-24 so the exemplar export can reach the
    per-administration grain that the 4h aggregation throws away. ONE copy of the
    unit conversion and the outlier bound live here -- a second copy is how two
    places quietly start disagreeing about what a mg is.
    """
    cats = CONFIG["medications"]["opioid_bolus_categories"]
    b = t["mai"][t["mai"]["med_category"].isin(cats)].merge(
        mapping, on="hospitalization_id", how="inner")
    b = b.merge(cohort[["encounter_block", "anchor_dttm", "weight_kg"]],
                on="encounter_block", how="inner")

    unit = (b["med_dose_unit"].astype("string").str.lower().str.strip()
            .replace({"mcg of opiate": "mcg", "": pd.NA}))
    dose = pd.to_numeric(b["med_dose"], errors="coerce")
    w = pd.to_numeric(b["weight_kg"], errors="coerce")
    mcg = pd.Series(np.nan, index=b.index, dtype="float64")
    mcg[unit == "mcg"] = dose
    mcg[unit == "mg"] = dose * 1000.0
    mcg[unit == "mcg/kg"] = dose * w
    mcg[unit == "mg/kg"] = dose * 1000.0 * w

    # Unconvertible rows are reported, not raised: at UCMC they are 7.6% of rows
    # and carry no dose at all, so raising would block the run over empty records.
    unresolved = dose.notna() & mcg.isna()
    if unresolved.any():
        counts = unit[unresolved].fillna("<no unit>").value_counts()
        print(f"  bolus: {int(unresolved.sum()):,} of {len(b):,} rows carry a dose in an "
              f"unhandled unit and are NULLED")
        for u, n in counts.items():
            print(f"    {u}: {n:,}")
    b["mcg"] = mcg

    lo, hi = OUTLIERS["med_bolus_mcg"]["fentanyl"]
    over = b["mcg"].notna() & ((b["mcg"] < lo) | (b["mcg"] > hi))
    if over.any():
        by_unit = unit[over].fillna("<no unit>").value_counts()
        print(f"  bolus: {int(over.sum()):,} converted doses outside [{lo:g}, {hi:g}] mcg, "
              f"nulled")
        for u, n in by_unit.items():
            print(f"    {u}: {n:,}")
    b.loc[over, "mcg"] = np.nan

    b["window_idx"] = ((b["admin_dttm"] - b["anchor_dttm"]).dt.total_seconds()
                       // (WINDOW_H * 3600)).astype("Int64")
    b = b[(b["window_idx"] >= 0) & (b["window_idx"] < N_WINDOWS) & b["mcg"].notna()]
    return b


def bolus_doses(b: pd.DataFrame) -> pd.DataFrame:
    """Window bolus total / weight / window_hours. Summed, never carried forward.

    Takes the per-administration frame from _bolus_events() rather than building
    it, so the unit conversion and its console report happen once per run even
    though two consumers now need the result.
    """
    agg = b.groupby(["encounter_block", "window_idx"], as_index=False).agg(
        bolus_mcg=("mcg", "sum"), n_bolus=("mcg", "size"))
    # mcg delivered in the window, spread over its hours -> mcg/hr, the same scale
    # as the infusion arm, so the two are additive. bolus_mcg travels alongside
    # for the amount scale: a bolus is charted as an amount already, so unlike
    # the infusion arm it needs no reconstruction and is exact in a short window.
    agg["bolus_dose"] = agg["bolus_mcg"] / WINDOW_H
    return agg[["encounter_block", "window_idx", "bolus_dose", "bolus_mcg",
                "n_bolus"]]


def window_exposure(grid: pd.DataFrame, bolus: pd.DataFrame,
                    windows: pd.DataFrame) -> pd.DataFrame:
    g = grid.copy()
    g["window_idx"] = g["hr"] // WINDOW_H
    # mean AND sum in one pass. The mean is the time-weighted RATE (mcg/hr); the
    # sum over hourly cells, each an hour wide, is the AMOUNT delivered (mcg).
    #
    # These are NOT interchangeable by a factor of WINDOW_H. _hourly_scaffold
    # drops cells past followup_end_dttm, so a window straddling discharge or
    # death holds fewer than WINDOW_H cells -- 632 ventilated windows at UCMC
    # (151 with 1 hour, 215 with 2, 266 with 3; measured 2026-09-09). There the
    # mean is still the right rate but mean * WINDOW_H overstates the amount.
    # Taking the sum directly is exact in both cases; reconstructing it from the
    # mean is not, which is why both are computed here rather than derived later.
    # Note extubation does NOT shorten a window: the scaffold ends at discharge
    # or death, and an extubated patient keeps all WINDOW_H cells.
    inf = g.groupby(["encounter_block", "window_idx"], as_index=False)["rate"].agg(
        inf_dose="mean", inf_mcg="sum")

    out = windows.merge(inf, on=["encounter_block", "window_idx"], how="left")
    out = out.merge(bolus, on=["encounter_block", "window_idx"], how="left")
    for c in ("inf_dose", "inf_mcg", "bolus_dose", "bolus_mcg", "n_bolus"):
        out[c] = out[c].fillna(0.0)
    out["total_dose"] = out["inf_dose"] + out["bolus_dose"]        # mcg/hr
    out["window_mcg"] = out["inf_mcg"] + out["bolus_mcg"]          # mcg

    # The derived ceiling is exact for the infusion arm only: it is
    # max(mcg/hr) / min(weight). A window's bolus SUM has no principled ceiling,
    # since several bounded boluses can stack, so that arm is reported not asserted.
    ceiling = fentanyl_sanity_ceiling(OUTLIERS)
    over_inf = out.loc[out["alive_admitted"], "inf_dose"] > ceiling
    if over_inf.any():
        raise SystemExit(
            f"{int(over_inf.sum()):,} windows have an INFUSION rate above the derived "
            f"ceiling of {ceiling:g} mcg/hr, which is proof the bounds did not run."
        )
    # Between the charted mcg/hr bound and the derived ceiling sits a real region:
    # a weight-based arm at its own bound, times a large weight. Bound-consistent,
    # clinically extreme, and reported rather than silently accepted.
    charted_max = fentanyl_charted_max_mcg_hr(OUTLIERS)
    high = out.loc[out["alive_admitted"], "inf_dose"] > charted_max
    if high.any():
        note(f"windows whose INFUSION rate exceeds the charted mcg/hr bound "
             f"({charted_max:g})", int(high.sum()))
        print(f"    max {out.loc[out['alive_admitted'], 'inf_dose'].max():.0f} mcg/hr; "
              f"these come from the weight-based arm and are bound-consistent, "
              f"not proof of a fault")
    over_total = out.loc[out["alive_admitted"], "total_dose"] > ceiling
    if over_total.any():
        top = out.loc[out["alive_admitted"] & (out["total_dose"] > ceiling), "total_dose"]
        note(f"windows whose total_dose exceeds {ceiling:g} mcg/hr (bolus stacking)",
             int(over_total.sum()))
        print(f"    max {top.max():.1f} mcg/hr; these are extreme but not "
              f"proof of a bounds failure")
    note("alive-admitted windows with any fentanyl", int((out.loc[out['alive_admitted'], 'total_dose'] > 0).sum()))
    note("alive-admitted windows with a bolus", int((out.loc[out['alive_admitted'], 'n_bolus'] > 0).sum()))
    return out


# ------------------------------------------------------------------ covariates
def _window_of(dttm: pd.Series, anchor: pd.Series) -> pd.Series:
    return ((dttm - anchor).dt.total_seconds() // (WINDOW_H * 3600)).astype("Int64")


def _to_windows(df: pd.DataFrame, cohort: pd.DataFrame, time_col: str) -> pd.DataFrame:
    """Map records to windows. Bounded by the window grid AND by follow-up.

    A record charted at or after block_discharge_dttm belongs to no window: a
    post-event window is structurally empty, not missing (covariates.json
    windows._anchor_note). Without the follow-up bound, charting lag put 428
    IMV records and 44 CRRT records into windows the patient had already been
    discharged from, so `imv_status == 1` could be true where alive_admitted
    was false. Same bound the hourly scaffold already applies.
    """
    d = df.merge(cohort[["encounter_block", "anchor_dttm", "followup_end_dttm"]],
                 on="encounter_block", how="inner")
    d["window_idx"] = _window_of(d[time_col], d["anchor_dttm"])
    d = d[d[time_col] < d["followup_end_dttm"]]

    return d[(d["window_idx"] >= 0) & (d["window_idx"] < N_WINDOWS)]


def lab_covariates(t: dict, cohort: pd.DataFrame) -> pd.DataFrame:
    labs = _to_windows(t["labs"], cohort, "lab_result_dttm")

    out = None
    for var, cat in LAB_VARS.items():
        how = COV["time_varying"][var]["summary"]
        sub = labs[labs["lab_category"] == cat]
        agg = (sub.groupby(["encounter_block", "window_idx"], as_index=False)
                  ["lab_value_numeric"].agg(how).rename(columns={"lab_value_numeric": var}))
        out = agg if out is None else out.merge(
            agg, on=["encounter_block", "window_idx"], how="outer")
    return out if out is not None else pd.DataFrame()


def assessment_covariates(t: dict, cohort: pd.DataFrame) -> pd.DataFrame:
    """Per-window sedation and pain scores. Config: time_varying.rass / .nvps.

    Deliberately the same shape as lab_covariates(): variables, source categories
    and summary rules all come from covariates.json, so adding a score is a config
    edit. gcs_total is NOT here -- it is a SOFA input with its own cap, carried by
    _sofa_inputs().
    """
    if not ASSESS_VARS:
        return pd.DataFrame()

    asm = _to_windows(t["assessments"], cohort, "recorded_dttm")
    vocabulary = {k for k in COV["summary_rules"] if not k.startswith("_")}

    out = None
    for var, cat in ASSESS_VARS.items():
        spec = COV["time_varying"][var]
        how = spec["summary"]
        # summary_rules says any rule outside its vocabulary must raise rather
        # than default. That matters most HERE: pandas .agg would happily accept
        # "median", so an undeclared rule would work in code while violating the
        # protocol, and nothing downstream would ever say so.
        if how not in vocabulary:
            raise SystemExit(
                f"time_varying.{var}.summary = {how!r} is not in "
                f"covariates.json summary_rules ({sorted(vocabulary)}). "
                f"Add the rule to the vocabulary first."
            )
        val = spec["source"]["value"]
        sub = asm[asm["assessment_category"] == cat]
        if sub.empty:
            raise SystemExit(
                f"assessment_category {cat!r} (declared as time_varying.{var}) "
                f"matched no rows in the cohort. Check the category spelling and "
                f"its CASE -- the filter matches literally."
            )
        agg = (sub.groupby(["encounter_block", "window_idx"], as_index=False)
                  [val].agg(how).rename(columns={val: var}))
        print(f"    {var:.<24} {cat} ({how}) -> {len(agg):,} block-windows")
        out = agg if out is None else out.merge(
            agg, on=["encounter_block", "window_idx"], how="outer")
    return out if out is not None else pd.DataFrame()


def nee_covariate(t: dict, mapping: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    """Max of the summed vasopressor step function. Config: time_varying.nee."""
    m = t["mac"][t["mac"]["med_category"].isin(NEE_COEF)].merge(
        mapping, on="hospitalization_id", how="inner")
    m = m[m["encounter_block"].isin(cohort["encounter_block"])]
    if m.empty:
        return pd.DataFrame(columns=["encounter_block", "window_idx", "nee"])

    m, rep = apply_med_raw(m, "med_category", "med_dose", "med_dose_unit", config=OUTLIERS)
    print(rep)
    m = _attach_current_weight(m, t)
    conv = m.copy()
    conv["dose_std"], drep = convert_doses(conv, "med_category", "med_dose",
                                           "med_dose_unit")
    print(drep)
    conv, rep = apply_med_converted(conv, "med_category", "dose_std", config=OUTLIERS)
    print(rep)

    stopped = conv.get("mar_action_category", pd.Series(index=conv.index, dtype=object))
    conv.loc[stopped.astype("string").str.lower() == "stop", "dose_std"] = 0.0

    conv = conv.merge(cohort[["encounter_block", "anchor_dttm"]], on="encounter_block")
    conv["hr"] = ((conv["admin_dttm"] - conv["anchor_dttm"]).dt.total_seconds() // 3600
                  ).astype("Int64")
    conv = conv[(conv["hr"] >= 0) & (conv["hr"] < EXTENT_H) & conv["dose_std"].notna()]

    last = (conv.sort_values(["encounter_block", "med_category", "hr", "dose_std"],
                             kind="stable")
                .groupby(["encounter_block", "med_category", "hr"], as_index=False)
                ["dose_std"].last())

    scaffold = _hourly_scaffold(cohort)[["encounter_block", "hr"]]
    drugs = pd.DataFrame({"med_category": list(NEE_COEF)})
    full = scaffold.merge(drugs, how="cross").merge(
        last, on=["encounter_block", "med_category", "hr"], how="left")
    full = full.sort_values(["encounter_block", "med_category", "hr"], kind="stable")
    full["dose_std"] = (full.groupby(["encounter_block", "med_category"])["dose_std"]
                            .ffill(limit=NEE_HOLD_H))
    full["contrib"] = full["dose_std"].fillna(0.0) * full["med_category"].map(NEE_COEF)

    hourly = full.groupby(["encounter_block", "hr"], as_index=False)["contrib"].sum()
    hourly["window_idx"] = hourly["hr"] // WINDOW_H
    nee = (hourly.groupby(["encounter_block", "window_idx"], as_index=False)["contrib"]
                 .max().rename(columns={"contrib": "nee"}))

    ceiling = nee_sanity_ceiling(NEE_COEF, OUTLIERS)
    if (nee["nee"] > ceiling).any():
        raise SystemExit(f"nee exceeds its arithmetic ceiling of {ceiling:g}")
    return nee


def oxygenation_covariate(t: dict, cohort: pd.DataFrame) -> pd.DataFrame:
    """One column on the P/F scale, Severinghaus fallback. Config: time_varying.oxygenation."""
    rs = t["resp"][t["resp"]["encounter_block"].isin(cohort["encounter_block"])].copy()
    rs["fio2_set"], rep = normalize_fio2(rs["fio2_set"])
    print(rep)
    fio2 = (rs.loc[rs["fio2_set"].notna(),
                   ["encounter_block", "recorded_dttm", "fio2_set"]]
              .sort_values("recorded_dttm", kind="stable"))
    any_row = (rs[["encounter_block", "recorded_dttm"]]
                 .assign(_covered=1).sort_values("recorded_dttm", kind="stable"))

    def pair(readings: pd.DataFrame, time_col: str) -> pd.DataFrame:
        r = readings.sort_values(time_col, kind="stable")
        m = pd.merge_asof(r, fio2, left_on=time_col, right_on="recorded_dttm",
                          by="encounter_block", direction="backward",
                          tolerance=pd.Timedelta(hours=FIO2_LOOKBACK_H))
        m = pd.merge_asof(m.sort_values(time_col), any_row, left_on=time_col,
                          right_on="recorded_dttm", by="encounter_block",
                          direction="backward",
                          tolerance=pd.Timedelta(hours=FIO2_LOOKBACK_H),
                          suffixes=("", "_any"))
        assumed = m["fio2_set"].isna() & m["_covered"].isna()
        m.loc[assumed, "fio2_set"] = RA_FIO2
        m["assumed"] = assumed
        return m[m["fio2_set"].notna() & (m["fio2_set"] > 0)]

    labs = t["labs"]
    pao2 = labs.loc[labs["lab_category"] == "po2_arterial",
                    ["encounter_block", "lab_result_dttm", "lab_value_numeric"]]
    pf = pair(pao2, "lab_result_dttm")
    pf["ratio"] = pf["lab_value_numeric"] / pf["fio2_set"]

    vit = t["vitals"]
    # Keep the unfiltered series: the ceiling filter is what separates a plateau
    # window from one with no measurement at all, so the flags must see both.
    spo2_all = vit.loc[vit["vital_category"] == "spo2",
                       ["encounter_block", "recorded_dttm", "vital_value"]]
    spo2 = spo2_all[spo2_all["vital_value"] < SPO2_CEILING]
    sf = pair(spo2.rename(columns={"recorded_dttm": "obs_dttm"}), "obs_dttm")
    sf["ratio"] = _severinghaus(sf["vital_value"]) / sf["fio2_set"]

    pf = _to_windows(pf.rename(columns={"lab_result_dttm": "obs_dttm"}), cohort, "obs_dttm")
    sf = _to_windows(sf, cohort, "obs_dttm")

    def summarise(d: pd.DataFrame, name: str) -> pd.DataFrame:
        g = d.groupby(["encounter_block", "window_idx"])
        return pd.DataFrame({
            name: g["ratio"].min(),
            f"{name}_all_assumed": g["assumed"].all(),
        }).reset_index()

    pf_w, sf_w = summarise(pf, "pf_ratio"), summarise(sf, "sf_ratio")
    out = pf_w.merge(sf_w, on=["encounter_block", "window_idx"], how="outer")

    # Raw availability BEFORE pairing, so a missing oxygenation can be attributed
    # to the measurement, the SpO2 ceiling, or the FiO2 lookback.
    def flag(d: pd.DataFrame, time_col: str, name: str) -> pd.DataFrame:
        w = _to_windows(d, cohort, time_col)
        return (w.groupby(["encounter_block", "window_idx"], as_index=False)
                 .size().rename(columns={"size": name})[
                     ["encounter_block", "window_idx", name]])

    # Plateau readings paired to FiO2: Severinghaus is undefined at or above the
    # ceiling, but the FiO2 still bounds what the P/F could be.
    plateau = spo2_all[spo2_all["vital_value"] >= SPO2_CEILING]
    if len(plateau):
        pl = pair(plateau.rename(columns={"recorded_dttm": "obs_dttm"}), "obs_dttm")
        pl = _to_windows(pl, cohort, "obs_dttm")
        pl_w = (pl.groupby(["encounter_block", "window_idx"], as_index=False)
                  ["fio2_set"].max().rename(columns={"fio2_set": "_plateau_fio2"}))
    else:
        pl_w = pd.DataFrame(columns=["encounter_block", "window_idx", "_plateau_fio2"])

    avail = flag(pao2, "lab_result_dttm", "_n_pao2")
    for d, col in ((spo2_all, "_n_spo2"), (spo2, "_n_spo2_usable")):
        avail = avail.merge(flag(d, "recorded_dttm", col),
                            on=["encounter_block", "window_idx"], how="outer")
    out = out.merge(avail, on=["encounter_block", "window_idx"], how="outer")
    out = out.merge(pl_w, on=["encounter_block", "window_idx"], how="left")

    use_pf = out["pf_ratio"].notna()
    out["oxygenation"] = out["pf_ratio"].where(use_pf, out["sf_ratio"])
    out["oxygenation_source"] = np.select(
        [use_pf & ~out["pf_ratio_all_assumed"].fillna(False),
         use_pf & out["pf_ratio_all_assumed"].fillna(False),
         out["sf_ratio"].notna() & ~out["sf_ratio_all_assumed"].fillna(False),
         out["sf_ratio"].notna() & out["sf_ratio_all_assumed"].fillna(False)],
        ["pf", "pf_room_air", "sf", "sf_room_air"], default="none")
    return out[["encounter_block", "window_idx", "oxygenation", "oxygenation_source",
                "pf_ratio", "sf_ratio", "_n_pao2", "_n_spo2", "_n_spo2_usable",
                "_plateau_fio2"]]


def _severinghaus(spo2: pd.Series) -> pd.Series:
    s = pd.to_numeric(spo2, errors="coerce") / 100.0
    s = s.where((s > 0) & (s < 1))
    a = 11700.0 / ((1.0 / s) - 1.0)
    b = np.sqrt(50.0 ** 3 + a ** 2)
    return np.cbrt(b + a) - np.cbrt(b - a)


SOFA_PRESSORS = ["norepinephrine", "epinephrine", "dopamine", "dobutamine"]


def _sofa_inputs(t: dict, cohort: pd.DataFrame) -> pd.DataFrame:
    """Per-window SOFA components not already carried as covariates."""
    labs = _to_windows(t["labs"], cohort, "lab_result_dttm")
    plt_ = (labs[labs["lab_category"] == "platelet_count"]
            .groupby(["encounter_block", "window_idx"], as_index=False)["lab_value_numeric"]
            .min().rename(columns={"lab_value_numeric": "platelet_count"}))
    creat = (labs[labs["lab_category"] == "creatinine"]
             .groupby(["encounter_block", "window_idx"], as_index=False)["lab_value_numeric"]
             .max().rename(columns={"lab_value_numeric": "creatinine"}))

    vit = _to_windows(t["vitals"], cohort, "recorded_dttm")
    mp = (vit[vit["vital_category"] == "map"]
          .groupby(["encounter_block", "window_idx"], as_index=False)["vital_value"]
          .min().rename(columns={"vital_value": "map"}))

    asm = _to_windows(t["assessments"], cohort, "recorded_dttm")
    gcs = (asm[asm["assessment_category"] == "gcs_total"]
           .groupby(["encounter_block", "window_idx"], as_index=False)["numerical_value"]
           .min().rename(columns={"numerical_value": "gcs_total"}))

    out = plt_
    for part in (creat, mp, gcs):
        out = out.merge(part, on=["encounter_block", "window_idx"], how="outer")
    return out


def _attach_current_weight(med: pd.DataFrame, t: dict) -> pd.DataFrame:
    """Attach the most recent charted weight at each admin time.

    clifpy demands a weight whenever the PREFERRED unit is weight-based, and that
    branch is first in its CASE, so a missing weight masks every other cause of a
    conversion failure. Attaching the column here also makes clifpy skip its own
    vitals lookup. NEE follows CURRENT weight, unlike the dose denominator, which
    is fixed at the anchor -- see covariates.json weight._DO_NOT_UNIFY.
    """
    w = (t["vitals"].loc[t["vitals"]["vital_category"] == "weight_kg",
                 ["encounter_block", "recorded_dttm", "vital_value"]]
            .dropna().sort_values("recorded_dttm", kind="stable"))

    out = med.sort_values("admin_dttm", kind="stable")
    out = pd.merge_asof(out, w, left_on="admin_dttm", right_on="recorded_dttm",
                        by="encounter_block", direction="backward",
                        suffixes=("", "_w"))
    out = out.rename(columns={"vital_value": "weight_kg"})

    first = w.groupby("encounter_block")["vital_value"].first().rename("_first_w")
    out = out.merge(first, on="encounter_block", how="left")
    filled = out["weight_kg"].isna() & out["_first_w"].notna()
    out.loc[filled, "weight_kg"] = out.loc[filled, "_first_w"]
    if int(filled.sum()):
        print(f"    weight backfilled from the block's first charted value: "
              f"{int(filled.sum()):,} rows administered before any weight")
    still = int(out["weight_kg"].isna().sum())
    if still:
        print(f"    {still:,} rows have no weight anywhere in the block")
    return out.drop(columns=[c for c in ("_first_w", "recorded_dttm") if c in out.columns])


def _sofa_pressors(t: dict, mapping: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    """Max mcg/kg/min per window for the four SOFA vasopressors.

    Absent drugs stay NULL. Coding them 0 makes the `<= 0.1` clause in the
    cardiovascular score fire, giving every unpressored patient a score of 3.
    """
    m = t["mac"][t["mac"]["med_category"].isin(SOFA_PRESSORS)].merge(
        mapping, on="hospitalization_id", how="inner")
    m = m[m["encounter_block"].isin(cohort["encounter_block"])]
    cols = ["encounter_block", "window_idx"] + [f"{d}_mcg_kg_min" for d in SOFA_PRESSORS]
    if m.empty:
        return pd.DataFrame(columns=cols)

    m, _ = apply_med_raw(m, "med_category", "med_dose", "med_dose_unit", config=OUTLIERS)
    m = _attach_current_weight(m, t)
    conv = m.copy()
    conv["dose_std"], _ = convert_doses(conv, "med_category", "med_dose", "med_dose_unit")
    conv = conv[conv["dose_std"] > 0]
    conv = _to_windows(conv, cohort, "admin_dttm")
    if conv.empty:
        return pd.DataFrame(columns=cols)

    wide = (conv.groupby(["encounter_block", "window_idx", "med_category"], as_index=False)
                ["dose_std"].max()
                .pivot(index=["encounter_block", "window_idx"],
                       columns="med_category", values="dose_std").reset_index())
    wide.columns.name = None
    for d in SOFA_PRESSORS:
        wide[f"{d}_mcg_kg_min"] = wide[d] if d in wide.columns else np.nan
    return wide[cols]


SOFA_PARTS = ["sofa_cv", "sofa_coag", "sofa_liver", "sofa_resp", "sofa_cns", "sofa_renal"]


def _severinghaus_at(spo2: float) -> float:
    s = spo2 / 100.0
    a = 11700.0 / ((1.0 / s) - 1.0)
    b = np.sqrt(50.0 ** 3 + a ** 2)
    return float(np.cbrt(b + a) - np.cbrt(b - a))


_SEVERINGHAUS_AT_CEILING = _severinghaus_at(SPO2_CEILING - 0.01)


def oxygenation_absence_reasons(long: pd.DataFrame) -> pd.DataFrame:
    """Why an alive-admitted window carries no oxygenation, before any carry-forward.

    Three mutually exclusive causes: nothing was measured; SpO2 was measured but
    sat on the plateau where the Severinghaus transform is undefined; or a usable
    measurement existed but no FiO2 could be paired to it within the lookback.
    """
    a = long[long["alive_admitted"]].copy()
    n_alive_admitted = len(a)
    pre_na = a["oxygenation"].isna()
    if "oxygenation_locf" in a.columns:
        pre_na = pre_na | a["oxygenation_locf"].fillna(False)

    n_pao2 = a["_n_pao2"].fillna(0)
    n_spo2 = a["_n_spo2"].fillna(0)
    n_usable = a["_n_spo2_usable"].fillna(0)

    nothing = pre_na & (n_pao2 == 0) & (n_spo2 == 0)
    plateau = pre_na & (n_pao2 == 0) & (n_spo2 > 0) & (n_usable == 0)
    no_fio2 = pre_na & ~nothing & ~plateau

    rows = []
    for label, mask, why in (
        ("no PaO2 and no SpO2 measured", nothing,
         "nothing to pair; not recoverable from this data"),
        (f"SpO2 present but all >= {SPO2_CEILING} (plateau)", plateau,
         "Severinghaus is undefined on the plateau; P/F is right-censored, not high"),
        ("usable measurement but no FiO2 within the lookback", no_fio2,
         f"the {FIO2_LOOKBACK_H}h fio2 pairing window is the binding constraint here"),
    ):
        k = int(mask.sum())
        rows.append({"variable": f"oxygenation absent: {label}",
                     "kind": "absence_reason", "class": why, "locf_cap_hours": "",
                     "n_alive_admitted": n_alive_admitted, "n_observed": 0, "n_zero_by_rule": 0,
                     "n_missing_pre_locf": k,
                     "pct_missing_pre_locf": round(100.0 * k / n_alive_admitted, 2),
                     "n_filled_by_locf": 0, "pct_filled_by_locf": 0.0,
                     "n_missing_final": k,
                     "pct_missing_final": round(100.0 * k / n_alive_admitted, 2)})
    if "_plateau_fio2" in a.columns:
        f = a.loc[plateau, "_plateau_fio2"].dropna()
        n_pl = int(plateau.sum())
        print(f"    plateau windows: {n_pl:,}; FiO2 pairable for {len(f):,} "
              f"({100*len(f)/max(n_pl,1):.1f}%)")
        if len(f):
            bound = _SEVERINGHAUS_AT_CEILING / f
            print(f"      FiO2 among them: median {f.median():.2f} "
                  f"(IQR {f.quantile(.25):.2f}-{f.quantile(.75):.2f})")
            print(f"      implied P/F lower bound: median {bound.median():.0f}, "
                  f"and {100*(bound < 300).mean():.1f}% are below 300")
            for lo, hi in ((0, 200), (200, 300), (300, 400), (400, 1e9)):
                k = int(((bound >= lo) & (bound < hi)).sum())
                lbl = f"{lo}-{hi}" if hi < 1e9 else f"{lo}+"
                print(f"        bound {lbl:>9}: {k:>7,}  ({100*k/len(f):5.1f}%)")

    total = int(pre_na.sum())
    print(f"    oxygenation absent before LOCF: {total:,} of {n_alive_admitted:,} "
          f"({100*total/n_alive_admitted:.1f}%)")
    for r in rows:
        print(f"      {r['variable'][20:]:<52} {r['n_missing_pre_locf']:>8,}"
              f" {r['pct_missing_pre_locf']:>6.2f}%")
    return pd.DataFrame(rows)


def _sofa_report_row(long: pd.DataFrame) -> pd.DataFrame:
    """Per-component coverage and the sofa_total row of the missingness report."""
    a = long[long["alive_admitted"]]
    n = len(a)
    print(f"    {'component':<14}{'scored':>10}{'pct':>8}")
    for c in SOFA_PARTS:
        k = int(a[c].notna().sum())
        print(f"    {c:<14}{k:>10,}{100*k/n:>7.1f}%")
    dist = a["sofa_n_components"].value_counts().sort_index()
    shown = "  ".join(f"{int(k)}:{100*v/n:.1f}%" for k, v in dist.items())
    print(f"    components per window -> {shown}")
    print(f"    sofa_total median {a['sofa_total'].median():.0f} "
          f"(IQR {a['sofa_total'].quantile(.25):.0f}-{a['sofa_total'].quantile(.75):.0f})")
    miss = int(a["sofa_total"].isna().sum())
    return pd.DataFrame([{
        "variable": "sofa_total", "kind": "derived",
        "class": "scored after all inputs are filled", "locf_cap_hours": "",
        "n_alive_admitted": n, "n_observed": n - miss, "n_zero_by_rule": 0,
        "n_missing_pre_locf": miss,
        "pct_missing_pre_locf": round(100.0 * miss / n, 2) if n else float("nan"),
        "n_filled_by_locf": 0, "pct_filled_by_locf": 0.0,
        "n_missing_final": miss,
        "pct_missing_final": round(100.0 * miss / n, 2) if n else float("nan"),
    }])


def score_sofa(df: pd.DataFrame) -> pd.DataFrame:
    """Six SOFA components and their total. Vincent 1996."""
    g = lambda c: df[c] if c in df.columns else pd.Series(np.nan, index=df.index)
    dopa, epi = g("dopamine_mcg_kg_min"), g("epinephrine_mcg_kg_min")
    norepi, dobu = g("norepinephrine_mcg_kg_min"), g("dobutamine_mcg_kg_min")
    mp, plt_, bili = g("map"), g("platelet_count"), g("bilirubin_total")
    creat, gcs, pf = g("creatinine"), g("gcs_total"), g("oxygenation")
    on_vent = g("imv_status").fillna(0).astype(bool)

    out = pd.DataFrame(index=df.index)
    out["sofa_cv"] = np.select(
        [(dopa > 15) | (epi > 0.1) | (norepi > 0.1),
         (dopa > 5) | (epi.notna() & (epi <= 0.1)) | (norepi.notna() & (norepi <= 0.1)),
         (dopa.notna() & (dopa <= 5)) | (dobu > 0),
         mp < 70, mp >= 70],
        [4, 3, 2, 1, 0], default=np.nan)
    out["sofa_coag"] = np.select(
        [plt_ < 20, plt_ < 50, plt_ < 100, plt_ < 150, plt_ >= 150],
        [4, 3, 2, 1, 0], default=np.nan)
    out["sofa_liver"] = np.select(
        [bili >= 12, bili >= 6, bili >= 2, bili >= 1.2, bili < 1.2],
        [4, 3, 2, 1, 0], default=np.nan)
    out["sofa_resp"] = np.select(
        [(pf < 100) & on_vent, (pf < 100) & ~on_vent,
         (pf < 200) & on_vent, (pf < 200) & ~on_vent,
         pf < 300, pf < 400, pf >= 400],
        [4, 3, 3, 2, 2, 1, 0], default=np.nan)
    out["sofa_cns"] = np.select(
        [gcs < 6, gcs <= 9, gcs <= 12, gcs <= 14, gcs == 15],
        [4, 3, 2, 1, 0], default=np.nan)
    out["sofa_renal"] = np.select(
        [creat >= 5, creat >= 3.5, creat >= 2, creat >= 1.2, creat < 1.2],
        [4, 3, 2, 1, 0], default=np.nan)

    parts = ["sofa_cv", "sofa_coag", "sofa_liver", "sofa_resp", "sofa_cns", "sofa_renal"]
    out["sofa_n_components"] = out[parts].notna().sum(axis=1)
    out["sofa_total"] = out[parts].sum(axis=1, min_count=1)
    return out


def status_covariates(t: dict, cohort: pd.DataFrame) -> pd.DataFrame:
    rs = _to_windows(t["resp"], cohort, "recorded_dttm")
    imv = (rs.assign(v=rs["device_category"].eq("IMV"))
             .groupby(["encounter_block", "window_idx"], as_index=False)["v"].max()
             .rename(columns={"v": "imv_status"}))

    crrt = _to_windows(t["crrt"], cohort, "recorded_dttm")
    cr = (crrt.assign(v=1).groupby(["encounter_block", "window_idx"], as_index=False)["v"]
              .max().rename(columns={"v": "crrt_status"}))
    return imv.merge(cr, on=["encounter_block", "window_idx"], how="outer")


def attach_weight(t: dict, cohort: pd.DataFrame) -> pd.DataFrame:
    """Dose denominator, fixed at the anchor. Backward match preferred; report the lag."""
    w = t["vitals"].loc[t["vitals"]["vital_category"] == "weight_kg",
                ["encounter_block", "recorded_dttm", "vital_value"]].dropna()
    w = w.sort_values("recorded_dttm", kind="stable")
    anchors = cohort[["encounter_block", "anchor_dttm"]].sort_values("anchor_dttm",
                                                                     kind="stable")

    back = pd.merge_asof(anchors, w, left_on="anchor_dttm", right_on="recorded_dttm",
                         by="encounter_block", direction="backward")
    fwd = pd.merge_asof(anchors, w, left_on="anchor_dttm", right_on="recorded_dttm",
                        by="encounter_block", direction="forward")
    out = back.rename(columns={"vital_value": "weight_kg", "recorded_dttm": "weight_dttm"})
    need = out["weight_kg"].isna()
    out.loc[need, "weight_kg"] = fwd.loc[need, "vital_value"].to_numpy()
    out.loc[need, "weight_dttm"] = fwd.loc[need, "recorded_dttm"].to_numpy()
    out["weight_direction"] = np.where(need, "forward", "backward")
    out["weight_lag_hours"] = (
        (out["anchor_dttm"] - out["weight_dttm"]).dt.total_seconds() / 3600).abs()

    note("blocks with no weight at all", int(out["weight_kg"].isna().sum()))
    note("blocks whose weight came from after the anchor", int(need.sum()))
    lag = out["weight_lag_hours"].dropna()
    if len(lag):
        print(f"  weight lag hours: median {lag.median():.1f}, "
              f"p95 {lag.quantile(0.95):.1f}, max {lag.max():.1f}")
    return out[["encounter_block", "weight_kg", "weight_dttm", "weight_direction",
                "weight_lag_hours"]]


def bmi_admission(t: dict, cohort: pd.DataFrame) -> pd.DataFrame:
    vit = t["vitals"].merge(cohort[["encounter_block", "block_admission_dttm"]],
                    on="encounter_block", how="inner")
    vit["lag_h"] = ((vit["recorded_dttm"] - vit["block_admission_dttm"])
                    .dt.total_seconds() / 3600)

    def first_of(cat: str) -> pd.DataFrame:
        d = vit[(vit["vital_category"] == cat) & vit["vital_value"].notna()].copy()
        early = d[(d["lag_h"] >= 0) & (d["lag_h"] <= 24)]
        rest = d[~d.index.isin(early.index)]
        pick = pd.concat([early.assign(_p=0), rest.assign(_p=1)])
        pick = pick.sort_values(["encounter_block", "_p", "lag_h"], kind="stable")
        return (pick.groupby("encounter_block", as_index=False)
                    .first()[["encounter_block", "vital_value", "lag_h"]]
                    .rename(columns={"vital_value": cat, "lag_h": f"{cat}_lag_h"}))

    wt, ht = first_of("weight_kg"), first_of("height_cm")
    b = wt.merge(ht, on="encounter_block", how="outer")
    b["bmi_admission"] = b["weight_kg"] / (b["height_cm"] / 100.0) ** 2
    b["bmi_lag_hours"] = b[["weight_kg_lag_h", "height_cm_lag_h"]].max(axis=1)
    note("blocks with a BMI", int(b["bmi_admission"].notna().sum()))
    return b[["encounter_block", "bmi_admission", "bmi_lag_hours"]]


def time_invariant(t: dict, mapping: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    pat = t["patient"][["patient_id", "sex_category", "race_category"]].rename(
        columns={"sex_category": "sex", "race_category": "race"})
    ti = cohort[["encounter_block", "patient_id", "age"]].merge(
        pat, on="patient_id", how="left")

    # CCI is per hospitalization; a block may hold several, so take the max.
    cci = calculate_cci(t["hospital_diagnosis"], hierarchy=True)
    cci = cci.to_pandas() if hasattr(cci, "to_pandas") else cci
    score_col = next((c for c in cci.columns if "cci" in c.lower()
                      and c != "hospitalization_id"), None)
    if score_col is None:
        raise SystemExit(f"no CCI score column in {list(cci.columns)}")
    cci = cci.merge(mapping, on="hospitalization_id", how="inner")
    cci = (cci.groupby("encounter_block", as_index=False)[score_col].max()
              .rename(columns={score_col: "cci"}))

    n_dx = t["hospital_diagnosis"].merge(mapping, on="hospitalization_id", how="inner")
    with_dx = set(n_dx["encounter_block"].unique())
    ti = ti.merge(cci, on="encounter_block", how="left")
    note("blocks with no diagnosis rows (CCI unknown, not 0)",
         int((~ti["encounter_block"].isin(with_dx)).sum()))
    ti.loc[~ti["encounter_block"].isin(with_dx), "cci"] = np.nan
    return ti


# -------------------------------------------------------------------- assembly
def apply_missingness(long: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Zero the absence-means-something variables, count, then LOCF. Order matters.

    Returns (long, per_variable, per_pattern). Counts are over alive-admitted rows only.
    """
    df = long.copy()
    alive_admitted = df["alive_admitted"]
    n_alive_admitted = int(alive_admitted.sum())

    tv = [k for k in COV["time_varying"] if not k.startswith("_") and k in df.columns]
    ti = [k for k in COV["time_invariant"] if not k.startswith("_") and k in df.columns]
    exposure = [c for c in EXPOSURE["columns"] if c in df.columns]

    observed = {v: int(df.loc[alive_admitted, v].notna().sum()) for v in tv}

    zeroed = {}
    for v in ZERO_VARS + NOT_VENT_VARS:
        if v in df.columns:
            blank = alive_admitted & df[v].isna()
            zeroed[v] = int(blank.sum())
            df.loc[blank, v] = 0.0

    pre = {v: int(df.loc[alive_admitted, v].isna().sum()) for v in tv}

    locf_caps = {
        k: (v["locf"] or {}).get("cap_hours")
        for k, v in COV["time_varying"].items()
        if not k.startswith("_") and (v.get("locf") or {}).get("eligible") and k in df.columns
    }
    df = df.sort_values(["encounter_block", "window_idx"], kind="stable")
    for v, cap in SOFA_INPUT_CAPS.items():
        if v in df.columns:
            df[v] = df.groupby("encounter_block")[v].ffill(
                limit=max(int(cap // WINDOW_H), 1))

    filled = {}
    for v, cap in locf_caps.items():
        limit = max(int(cap // WINDOW_H), 1)
        before = df[v].isna()
        df[v] = df.groupby("encounter_block")[v].ffill(limit=limit)
        df[f"{v}_locf"] = before & df[v].notna()
        filled[v] = int((df[f"{v}_locf"] & alive_admitted).sum())

    post = {v: int(df.loc[alive_admitted, v].isna().sum()) for v in tv}

    def pct(n: int) -> float:
        return round(100.0 * n / n_alive_admitted, 2) if n_alive_admitted else float("nan")

    rows = []
    for v in tv:
        cap = locf_caps.get(v)
        rows.append({
            "variable": v,
            "kind": "time_varying",
            "class": COV["time_varying"][v]["missing_class"],
            "locf_cap_hours": cap if cap is not None else "",
            "n_alive_admitted": n_alive_admitted,
            "n_observed": observed[v],
            "n_zero_by_rule": zeroed.get(v, 0),
            "n_missing_pre_locf": pre[v],
            "pct_missing_pre_locf": pct(pre[v]),
            "n_filled_by_locf": filled.get(v, 0),
            "pct_filled_by_locf": pct(filled.get(v, 0)),
            "n_missing_final": post[v],
            "pct_missing_final": pct(post[v]),
        })
    for v in exposure + ti:
        n_miss = int(df.loc[alive_admitted, v].isna().sum()) if v in exposure \
            else int(df[v].isna().sum())
        denom = n_alive_admitted if v in exposure else len(df)
        rows.append({
            "variable": v,
            "kind": "exposure" if v in exposure else "time_invariant",
            "class": "", "locf_cap_hours": "",
            "n_alive_admitted": denom,
            "n_observed": denom - n_miss,
            "n_zero_by_rule": 0,
            "n_missing_pre_locf": n_miss,
            "pct_missing_pre_locf": round(100.0 * n_miss / denom, 2) if denom else float("nan"),
            "n_filled_by_locf": 0, "pct_filled_by_locf": 0.0,
            "n_missing_final": n_miss,
            "pct_missing_final": round(100.0 * n_miss / denom, 2) if denom else float("nan"),
        })
    per_variable = pd.DataFrame(rows)

    if "oxygenation_source" in df.columns:
        src = (df.loc[alive_admitted, "oxygenation_source"].fillna("none")
                 .value_counts().rename_axis("oxygenation_source")
                 .reset_index(name="n"))
        src["pct"] = (100.0 * src["n"] / n_alive_admitted).round(2)
        print("\n  oxygenation provenance (alive-admitted windows)")
        print(src.to_string(index=False))

    for v in ZERO_VARS + NOT_VENT_VARS:
        if v in df.columns:
            df[v] = pd.to_numeric(df[v], errors="coerce").astype("float64")

    per_pattern = _missingness_patterns(df[alive_admitted], tv)
    return df, per_variable, per_pattern


def _missingness_patterns(df: pd.DataFrame, variables: list[str]) -> pd.DataFrame:
    """Which variables go missing together. Cells below the threshold are pooled."""
    min_cell = COV["missing_values"]["reporting"]["pattern_min_cell"]
    if not variables or df.empty:
        return pd.DataFrame()
    key = df[variables].isna()
    pattern = key.apply(lambda r: "".join("1" if x else "0" for x in r), axis=1)
    counts = pattern.value_counts().reset_index()
    counts.columns = ["pattern", "n"]
    small = counts["n"] < min_cell
    pooled = int(counts.loc[small, "n"].sum())
    n_pooled = int(small.sum())
    out = counts[~small].copy()
    if n_pooled:
        out = pd.concat([out, pd.DataFrame([{
            "pattern": f"<other: {n_pooled} patterns each below {min_cell}>",
            "n": pooled}])], ignore_index=True)
    out["pct"] = (100.0 * out["n"] / len(df)).round(2)
    out.insert(0, "variables", "|".join(variables))
    return out


def assert_config_is_honoured(long: pd.DataFrame) -> None:
    """Every declared covariate must reach a column, or the build fails."""
    declared = [k for k in COV["time_varying"] if not k.startswith("_")]
    declared += [k for k in COV["time_invariant"] if not k.startswith("_")]
    declared += EXPOSURE["columns"]
    missing = [c for c in declared if c not in long.columns]
    if missing:
        raise SystemExit(
            f"declared in config/covariates.json but absent from trajectory_long: {missing}. "
            f"Adding a key to that config must never be a silent no-op."
        )


def build_time_to_event(cohort: pd.DataFrame, long: pd.DataFrame,
                        imv_records: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    """One row per block, origin = landmark T.

    Successful extubation = extubation not followed by reintubation within
    successful_extubation_hours. Ventilation status comes from the RAW IMV series,
    not the window grid: the grid stops at T, so outcomes ascertained from it are
    all censored by construction.
    """
    T = CONFIG["cohort"]["landmark_hours"]
    succ_h = CONFIG["outcomes"]["successful_extubation_hours"]
    death_cats = set(CONFIG["outcomes"]["mortality_categories"])

    at_T = long[(long["window_start_hr"] == T - WINDOW_H) & long["alive_admitted"]]
    eligible = set(at_T.loc[at_T["imv_status"] == 1, "encounter_block"])
    flow(f"LANDMARK COHORT -- time_to_event (ventilated at T={T}h)", len(eligible),
         f"extubated, died or discharged before T={T}h")

    tte = cohort[cohort["encounter_block"].isin(eligible)].copy()
    tte["landmark_eligible"] = True
    tte["hours_to_discharge"] = (
        (tte["block_discharge_dttm"] - tte["anchor_dttm"]).dt.total_seconds() / 3600)
    tte["died"] = tte["discharge_category"].isin(death_cats)

    imv = (imv_records.merge(mapping, on="hospitalization_id", how="inner")
                      .merge(tte[["encounter_block", "anchor_dttm"]],
                             on="encounter_block", how="inner"))
    imv["hr"] = (imv["recorded_dttm"] - imv["anchor_dttm"]).dt.total_seconds() / 3600
    imv = imv[imv["hr"] >= 0].sort_values(["encounter_block", "hr"], kind="stable")

    # An extubation is an IMV record whose next IMV record is more than succ_h
    # later, or which has none. The first such event at or after T is the outcome.
    imv["next_hr"] = imv.groupby("encounter_block")["hr"].shift(-1)
    gap = imv["next_hr"] - imv["hr"]
    imv["is_extubation"] = imv["next_hr"].isna() | (gap > succ_h)
    ext = imv[imv["is_extubation"] & (imv["hr"] >= T)]
    first_ext = ext.groupby("encounter_block")["hr"].min().rename("ext_hr")
    tte = tte.merge(first_ext, on="encounter_block", how="left")

    note("blocks with an extubation at or after T", int(tte["ext_hr"].notna().sum()))

    # Death inside the succ_h window after extubation is the competing event,
    # not a success, per the VFD convention.
    off = tte["hours_to_discharge"] - tte["ext_hr"]
    died_in_window = tte["died"] & (off < succ_h)
    success = tte["ext_hr"].notna() & ~died_in_window

    tte["extubation_event"] = np.select([success, tte["died"]], [1, 2], default=0)
    tte["extubation_time"] = np.where(
        tte["extubation_event"] == 1, tte["ext_hr"] - T,
        (tte["hours_to_discharge"] - T).clip(lower=0))

    tte["mortality_event"] = np.where(tte["died"], 1, 2)
    tte["mortality_time"] = (tte["hours_to_discharge"] - T).clip(lower=0)

    unresolved = set(CONFIG["outcomes"]["unresolved_discharge_categories"])
    n_unres = int(tte["discharge_category"].isin(unresolved).sum())
    if n_unres:
        note("blocks with an unresolved discharge_category", n_unres)
    for code, label in ((1, "successful extubation"), (2, "death"), (0, "censored")):
        note(f"  extubation outcome = {label}",
             int((tte["extubation_event"] == code).sum()))
    note("  30-day analysis: died", int((tte["mortality_event"] == 1).sum()))

    # Tracheostomy is event code 3 and its rule is still undecided.
    tte["tracheostomy_pending"] = True
    return tte


def exemplar_export(cohort: pd.DataFrame, long: pd.DataFrame, grid: pd.DataFrame,
                    bolus_events: pd.DataFrame, t: dict, dirs: dict) -> None:
    """Select ONE episode by the declared rule and export its sub-hourly series.

    F1 draws a single patient in detail. Baker et al. Sci Rep 2020;10:10718
    Figure 1, which F1 is modelled on, is captioned only "Representative ICU
    admission" and states NO selection rule -- the one thing F1 has to improve
    on. So the criteria are declared in covariates.json, applied here, and the
    episode is DRAWN AT RANDOM from those that qualify. Nothing about the choice
    depends on anyone looking at a patient, which is also what lets the figure be
    designed without patient data ever being displayed.

    If the drawn episode reads badly, change the RULE and re-run. Picking a
    different one from the eligible set by eye reintroduces exactly the bias the
    rule exists to remove.
    """
    spec = COV.get("exemplar")
    if spec is None:
        raise SystemExit("config/covariates.json declares no `exemplar` block; "
                         "F1 has no selection rule to apply.")
    if spec["draw"] != "uniform_random_from_eligible":
        raise SystemExit(f"unsupported exemplar draw: {spec['draw']!r}")

    print("\nExemplar (F1) selection")
    anchor_n = long["encounter_block"].nunique()

    # ---- ventilation shape, from the SAME window grid the states use ----------
    # extubated = alive, admitted, and not on IMV, after having been on it. Read
    # off imv_status so this figure and code/utils/states.R cannot disagree about
    # what extubation is.
    w = long.sort_values(["encounter_block", "window_idx"])
    w = w.assign(_vent=(w["imv_status"].notna() & (w["imv_status"] == 1)))
    g = w.groupby("encounter_block", sort=False)
    ever_vent = g["_vent"].any()
    # Block-level death, from discharge_category -- it already covers Expired
    # AND Hospice, so both forms of comfort care are caught.
    died_block = g["died"].any()
    # first window that is alive-admitted and NOT ventilated, after ventilation
    def _extub_window(d):
        v = d["_vent"].to_numpy()
        if not v.any():
            return pd.NA
        after = np.flatnonzero(~v & d["alive_admitted"].to_numpy())
        after = after[after > np.flatnonzero(v)[0]]
        if not len(after):
            return pd.NA
        first = after[0]
        # no reintubation inside the window: the trace ends at extubation, and a
        # patient who goes back on the vent would make that ending a lie.
        if v[first:].any():
            return pd.NA
        return int(d["window_idx"].to_numpy()[first])
    extub_w = g.apply(_extub_window, include_groups=False)
    extub_w.name = "extub_window"

    # ---- fentanyl shape, from the hourly grid --------------------------------
    gr = grid.sort_values(["encounter_block", "hr"])
    def _fent(d):
        r = d["rate"].to_numpy()
        nz = np.flatnonzero(r > 0)
        if not len(nz):
            return pd.Series({"cont_hours": 0, "gap_hours": 0})
        inner = r[nz[0]:nz[-1] + 1] == 0          # zeros BETWEEN infusions only
        best = 0
        run = 0
        for z in inner:
            run = run + 1 if z else 0
            best = max(best, run)
        return pd.Series({"cont_hours": int((r > 0).sum()), "gap_hours": int(best)})
    fent = gr.groupby("encounter_block", sort=False).apply(_fent, include_groups=False)

    n_bolus = bolus_events.groupby("encounter_block").size().rename("n_bolus")

    # ---- assessment density, from the RAW timestamped records ----------------
    asm = t["assessments"].merge(cohort[["encounter_block", "anchor_dttm"]],
                                 on="encounter_block", how="inner")
    asm["t_hr"] = ((asm["recorded_dttm"] - asm["anchor_dttm"]).dt.total_seconds()
                   / 3600.0)
    asm = asm[(asm["t_hr"] >= 0) & (asm["t_hr"] <= EXTENT_H)]
    dens = (asm[asm["assessment_category"].isin(["RASS", "NVPS"])]
            .groupby(["encounter_block", "assessment_category"]).size()
            .unstack(fill_value=0))
    for c in ("RASS", "NVPS"):
        if c not in dens.columns:
            dens[c] = 0

    f = (pd.DataFrame(index=pd.Index(sorted(long["encounter_block"].unique()),
                                     name="encounter_block"))
         .join(extub_w).join(fent).join(n_bolus).join(dens[["RASS", "NVPS"]])
         .join(died_block.rename("died")))
    f = f.fillna({"cont_hours": 0, "gap_hours": 0, "n_bolus": 0,
                  "RASS": 0, "NVPS": 0, "died": False})

    # ---- the funnel. Printed so a threshold that empties the pool is visible --
    # immediately, and loosening it is a protocol change rather than a surprise.
    # The declared deadline is APPLIED, not assumed. It happens to equal the
    # window extent today, so the check is inert on this config -- but a site
    # that sets it to 48 must get 48, and a threshold that is read by nothing is
    # the failure mode tests/test_covariates.py exists to catch.
    extub_hr_col = f["extub_window"].astype("Float64") * WINDOW_H
    tests = [
        (f"extubated by {spec['require_extubated_by_hours']}h, no reintubation",
         f["extub_window"].notna()
         & (extub_hr_col <= spec["require_extubated_by_hours"]).fillna(False)),
        ("survived the hospitalisation (excludes comfort care)",
         ~f["died"].astype(bool)) if spec.get("require_survived_hospitalization")
        else ("survival not required", pd.Series(True, index=f.index)),
        (f"continuous infusion >= {spec['min_continuous_hours']}h",
         f["cont_hours"] >= spec["min_continuous_hours"]),
        (f">= {spec['min_boluses']} boluses", f["n_bolus"] >= spec["min_boluses"]),
        (f">= {spec['min_rass_observations']} RASS",
         f["RASS"] >= spec["min_rass_observations"]),
        (f">= {spec['min_nvps_observations']} NVPS",
         f["NVPS"] >= spec["min_nvps_observations"]),
    ]
    keep = pd.Series(True, index=f.index)
    counts = []
    for label, ok in tests:
        alone = int(ok.sum())
        keep = keep & ok
        print(f"    {label:<48s} {alone:>6,} alone, {int(keep.sum()):>6,} cumulative")
        counts.append({"criterion": label, "n_passing_alone": alone,
                       "n_passing_cumulative": int(keep.sum())})
    n_eligible = int(keep.sum())
    print(f"    {'ELIGIBLE':<48s} {n_eligible:>6,} of {anchor_n:,} episodes")
    counts.append({"criterion": "ELIGIBLE", "n_passing_alone": n_eligible,
                   "n_passing_cumulative": n_eligible})

    # MEASURED, not filtered on. A >= 6h interruption was briefly a criterion and
    # was removed (SG, 2026-09-24): it cut the pool by 90% and would have made
    # the exemplar unrepresentative of how fentanyl is actually delivered here.
    # The number is reported because it is a finding in its own right.
    for n_b in (2, 6, 10, 20):
        n_hit = int((f["n_bolus"] >= n_b).sum())
        print(f"    (not a criterion) >= {n_b:>2d} boluses: {n_hit:,} of "
              f"{anchor_n:,} episodes, {100 * n_hit / anchor_n:.1f}%")
        counts.append({"criterion": f"NOT A CRITERION: >= {n_b} boluses",
                       "n_passing_alone": n_hit, "n_passing_cumulative": pd.NA})
    for g_h in (6, 12):
        n_gap = int((f["gap_hours"] >= g_h).sum())
        print(f"    (not a criterion) ever off fentanyl >= {g_h}h: "
              f"{n_gap:,} of {anchor_n:,} episodes, {100 * n_gap / anchor_n:.1f}%")
        counts.append({"criterion": f"NOT A CRITERION: off-fentanyl gap >= {g_h}h",
                       "n_passing_alone": n_gap, "n_passing_cumulative": pd.NA})
    pd.DataFrame(counts).to_csv(dirs["phase"] / "exemplar_selection.csv", index=False)

    if not n_eligible:
        print("  NO episode meets every criterion. F1 cannot be drawn; loosen the "
              "thresholds in covariates.json `exemplar` and re-run.")
        return

    rng = np.random.default_rng(spec["draw_seed"])
    # NOT str(): encounter_block is numeric at this site, and coercing the
    # drawn key to text made every downstream .loc and == miss.
    chosen = f.index[keep][rng.integers(n_eligible)]
    extub_window = int(f.loc[chosen, "extub_window"])
    extub_hr = float(extub_window * WINDOW_H)
    # first_imv_episode_hours is computed from the raw records, so it locates
    # extubation more precisely than the 4h grid can. Use it when the two agree;
    # the assertion is what stops a silent disagreement becoming a wrong figure.
    fine = long.loc[long["encounter_block"] == chosen, "first_imv_episode_hours"]
    fine = float(fine.iloc[0]) if len(fine) and pd.notna(fine.iloc[0]) else np.nan
    if np.isfinite(fine) and extub_hr <= fine <= extub_hr + WINDOW_H:
        extub_hr = fine
    print(f"  drawn at random (seed {spec['draw_seed']}) from {n_eligible:,}; "
          f"extubated at {extub_hr:.1f}h")

    # ---- the series, de-identified at construction ---------------------------
    inf = grid[grid["encounter_block"] == chosen][["hr", "rate"]]
    inf = pd.DataFrame({"t_hr": inf["hr"].astype(float), "series": "infusion",
                        "value": inf["rate"].astype(float)})

    bo = bolus_events[bolus_events["encounter_block"] == chosen]
    bo = pd.DataFrame({
        "t_hr": ((bo["admin_dttm"] - bo["anchor_dttm"]).dt.total_seconds() / 3600.0),
        "series": "bolus", "value": bo["mcg"].astype(float)})

    a = asm[asm["encounter_block"] == chosen]
    ord_ = pd.DataFrame({
        "t_hr": a["t_hr"].astype(float),
        "series": a["assessment_category"].str.lower(),
        "value": pd.to_numeric(a["numerical_value"], errors="coerce")})
    ord_ = ord_[ord_["series"].isin(["rass", "nvps"]) & ord_["value"].notna()]

    series = pd.concat([inf, bo, ord_], ignore_index=True)
    series = series[(series["t_hr"] >= 0) & (series["t_hr"] <= EXTENT_H)]
    series = series.sort_values(["series", "t_hr"]).reset_index(drop=True)

    # Asserted, not assumed. Same contract as the 100-episode raster in
    # 04_delivery_states.R: relative hours only, no dates, no identifiers.
    ident = {"encounter_block", "patient_id", "hospitalization_id", "anchor_dttm"}
    assert not (ident & set(series.columns)), "exemplar series carries an identifier"
    assert not any(pd.api.types.is_datetime64_any_dtype(series[c])
                   for c in series.columns), "exemplar series carries a datetime"
    assert set(series.columns) == {"t_hr", "series", "value"}
    assert series["t_hr"].min() >= 0

    out = dirs["out_phi"]
    series.to_parquet(out / "exemplar_series.parquet", index=False)
    with open(out / "exemplar_meta.json", "w") as fh:
        json.dump({"extent_h": EXTENT_H, "window_h": WINDOW_H,
                   "extubation_hr": extub_hr, "n_eligible": n_eligible,
                   "criteria": {k: v for k, v in spec.items()
                                if not k.startswith("_")
                                and k not in ("draw", "draw_seed")}},

                  fh, indent=2)
    # The chosen id stays on the PHI side so the site can reproduce the figure;
    # it never reaches the shareable tree, where only the rule and the N appear.
    (out / "exemplar_id.txt").write_text(f"{chosen}\n")
    print(f"  exemplar series: {len(series):,} rows across "
          f"{series['series'].nunique()} series")


def main() -> None:
    dirs = site_dirs(REPO)
    dirs["phase"] = phase_dir(dirs, PHASE_DIR)
    dirs["diagnostics"] = dirs["phase"] / "diagnostics"
    dirs["diagnostics"].mkdir(parents=True, exist_ok=True)
    n_cleared = clear_owned_outputs(dirs, OWNED, retired=RETIRED_OUTPUTS)
    if n_cleared:
        print(f"cleared {n_cleared} output(s) from a previous run")
    prov = provenance(CONFIG)
    print(f"site {prov['site_name']}  clif {prov['clif_version']}  code {prov['code_version']}")
    print(f"grid {WINDOW_H}h x {N_WINDOWS} to {EXTENT_H}h   stitch {STITCH_H}h   "
          f"infusion hold {INF_HOLD_H}h\n")

    print("Loading core tables")
    t, raw_anchor, imv_records = load_core()

    print("\nEncounter blocks")
    blocks, mapping = build_blocks(t)
    hi = hospital_intervals(t, mapping)
    ends = hospital_endpoints(hi)

    # Restrict before the waterfall: it is the expensive step and only the
    # cohort's rows can affect the result.
    imv_ids = set(raw_anchor["hospitalization_id"])
    imv_blocks = set(mapping.loc[mapping["hospitalization_id"].isin(imv_ids), "encounter_block"])
    blocks = blocks[blocks["encounter_block"].isin(imv_blocks)]
    flow("Blocks with any IMV record", len(blocks), "no invasive ventilation recorded")
    mapping = mapping[mapping["encounter_block"].isin(imv_blocks)]
    hosp_ids = sorted(mapping["hospitalization_id"].astype(str).unique())
    note("hospitalizations to load", len(hosp_ids))

    print("\nLoading cohort tables")
    t.update(load_cohort_tables(hosp_ids))
    assert_categories_present(t)
    t = prepare_long_tables(t, mapping)

    print("\nCohort")
    anchor = imv_episodes(t["resp"], imv_records, mapping)
    cohort = build_cohort(blocks, anchor)

    print("\nWeight and BMI")
    cohort = cohort.merge(attach_weight(t, cohort), on="encounter_block", how="left")
    cohort = cohort[cohort["weight_kg"].notna()]
    flow("ANALYTIC COHORT -- trajectory_long", len(cohort),
         "no weight charted anywhere in the block")
    bmi = bmi_admission(t, cohort)

    print("\nWindows")
    windows = window_grid(cohort)

    print("\nExposure")
    grid = infusion_grid(t, mapping, cohort,
                         CONFIG["medications"]["opioid_infusion_categories"],
                         _to_mcg_hr)
    bolus_events = _bolus_events(t, mapping, cohort)
    bolus = bolus_doses(bolus_events)
    long = window_exposure(grid, bolus, windows)

    print("\n  Sedatives (descriptive companions, infusions only)")
    long = long.merge(sedative_exposure(t, mapping, cohort, windows),
                      on=["encounter_block", "window_idx"], how="left")

    print("\nCovariates")
    for part in (lab_covariates(t, cohort),
                 assessment_covariates(t, cohort),
                 nee_covariate(t, mapping, cohort),
                 oxygenation_covariate(t, cohort),
                 status_covariates(t, cohort)):
        if len(part):
            long = long.merge(part, on=["encounter_block", "window_idx"], how="left")

    long = gate_dose_on_ventilation(long)

    hi = hi[hi["encounter_block"].isin(cohort["encounter_block"])]
    ti = time_invariant(t, mapping, cohort).merge(bmi, on="encounter_block", how="left")
    ti = ti.merge(ends, on="encounter_block", how="left")
    long = long.merge(ti.drop(columns=["patient_id"]), on="encounter_block", how="left")
    long = long.merge(
        cohort[["encounter_block", "weight_kg", "first_imv_episode_hours",
                "n_imv_episodes", "discharge_category"]],
        on="encounter_block", how="left")

    # died is derived here rather than left to each consumer, so the mortality
    # rule lives in one place. time_to_event applies the identical rule to the
    # landmark subset; this carries it for the WHOLE analytic cohort, which the
    # Phase 1 state description needs in order to tell died from discharged.
    long["died"] = long["discharge_category"].isin(
        set(CONFIG["outcomes"]["mortality_categories"]))

    print("\n  SOFA")
    sofa_in = _sofa_inputs(t, cohort)
    press = _sofa_pressors(t, mapping, cohort)
    if len(press):
        sofa_in = sofa_in.merge(press, on=["encounter_block", "window_idx"], how="outer")
    long = long.merge(sofa_in, on=["encounter_block", "window_idx"], how="left")

    print("\nMissingness and LOCF")
    long, per_variable, per_pattern = apply_missingness(long)

    # SOFA is scored only once every input has been carried forward. Scoring it
    # earlier used raw bilirubin and oxygenation while the other four components
    # were filled, which is what held 6-component coverage at 9.3%.
    print("\n  SOFA")
    long = pd.concat([long.reset_index(drop=True),
                      score_sofa(long).reset_index(drop=True)], axis=1)
    per_variable = pd.concat(
        [per_variable, _sofa_report_row(long),
         oxygenation_absence_reasons(long)], ignore_index=True)
    # Emit the plateau as a COVARIATE before the internals are dropped (SG,
    # 2026-09-09). A window where SpO2 was measured, every reading sat at or above
    # the ceiling, and no PaO2 was drawn carries real information: the patient is
    # oxygenating well and the P/F is right-censored HIGH, not unknown. Filling
    # such a window with a median P/F imputes moderate hypoxaemia for exactly the
    # patients doing best. This flag lets a model use the fact directly.
    #
    # Note the condition does NOT include "oxygenation is missing": it is a
    # statement about what was measured in the window, so it stays true and
    # interpretable even where LOCF later carried a value in.
    long["spo2_plateau"] = (
        (long["_n_pao2"].fillna(0) == 0)
        & (long["_n_spo2"].fillna(0) > 0)
        & (long["_n_spo2_usable"].fillna(0) == 0)
    ).astype(int)
    note(f"windows on the SpO2 plateau (>= {SPO2_CEILING}, no PaO2)",
         int(long.loc[long["alive_admitted"], "spo2_plateau"].sum()))
    pl = long.loc[long["alive_admitted"] & long["spo2_plateau"].astype(bool)]
    if len(pl):
        print(f"    of those, oxygenation still missing after LOCF: "
              f"{int(pl['oxygenation'].isna().sum()):,} "
              f"({100 * pl['oxygenation'].isna().mean():.1f}%)")

    long = long.drop(columns=[c for c in long.columns if c.startswith("_n_")])
    cols = ["variable", "kind", "locf_cap_hours", "n_observed", "n_zero_by_rule",
            "pct_missing_pre_locf", "pct_filled_by_locf", "pct_missing_final"]
    print("\n  missingness, alive-admitted rows only")
    print(per_variable[cols].to_string(index=False))

    assert_config_is_honoured(long)

    print("\nTime to event")
    tte = build_time_to_event(cohort, long, imv_records, mapping)

    out = dirs["out_phi"]
    # Parquet is what later phases read; CSV is the same table for human review.
    long.to_parquet(out / "trajectory_long.parquet", index=False)
    long.to_csv(out / "trajectory_long.csv", index=False)
    tte.to_parquet(out / "time_to_event.parquet", index=False)
    tte.to_csv(out / "time_to_event.csv", index=False)
    hi.to_parquet(out / "hospital_intervals.parquet", index=False)

    # Deliberately AFTER the analytic tables are on disk. F1 is a figure input,
    # not a core table: a fault in the exemplar selection must not cost the whole
    # Phase 0 rebuild, which is exactly what it did on 2026-09-24.
    exemplar_export(cohort, long, grid, bolus_events, t, dirs)
    diag = dirs["diagnostics"]
    per_variable.to_csv(diag / "missingness.csv", index=False)
    if len(per_pattern):
        per_pattern.to_csv(diag / "missingness_patterns.csv", index=False)
    pd.DataFrame(STROBE, columns=["step", "n"]).to_csv(
        diag / "diagnostics.csv", index=False)

    pd.DataFrame(FLOW).to_csv(dirs["phase"] / "strobe.csv", index=False)
    (dirs["phase"] / "strobe.txt").write_text(render_text(FLOW) + "\n")
    render_png(FLOW, dirs["phase"] / "strobe.png",
               title=f"Phase 0 cohort flow -- {CONFIG['site_name']}")
    print()
    print(render_text(FLOW))
    (dirs["phase"] / "provenance.json").write_text(json.dumps(prov, indent=2))

    # Written last: its presence is what marks these outputs complete and current.
    write_manifest(dirs, CONFIG, REPO, {
        "trajectory_long": len(long),
        "time_to_event": len(tte),
        "hospital_intervals": len(hi),
    })

    print(f"\nwritten to {out}")
    print(f"  trajectory_long.parquet  {len(long):,} rows x {long.shape[1]} cols")
    print(f"  time_to_event.parquet    {len(tte):,} rows")
    print(f"  hospital_intervals.parquet {len(hi):,} rows")
    for name in ("trajectory_long.csv", "time_to_event.csv"):
        mb = (out / name).stat().st_size / 1e6
        print(f"  {name:<26} {mb:,.1f} MB  (review copy)")
    print("  manifest.json written -- outputs are marked complete")


if __name__ == "__main__":
    main()
