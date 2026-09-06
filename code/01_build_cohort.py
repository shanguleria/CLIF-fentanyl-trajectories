"""Phase 0 -- build the analytic tables from CLIF.

Outputs (both PHI, written to data/intermediate_phi/):
    trajectory_long.parquet     one row per encounter block per window
    time_to_event.parquet       one row per encounter block
    hospital_intervals.parquet  one row per ADT interval

Protocol: config/covariates.json. Column spec: docs/design_notes.md §10, §11.
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
from utils.fio2 import normalize_fio2  # noqa: E402
from utils.outliers import (  # noqa: E402
    apply_long, apply_med_converted, apply_med_raw, fentanyl_sanity_ceiling,
    load_config as load_outliers, nee_sanity_ceiling,
)
from utils.paths import provenance, site_dirs  # noqa: E402

CONFIG = json.loads((REPO / "config" / "config.json").read_text())
COV = json.loads((REPO / "config" / "covariates.json").read_text())
OUTLIERS = load_outliers()

GRID = COV["windows"]["granular"]
WINDOW_H = GRID["width_hours"]
EXTENT_H = GRID["extent_hours"]
N_WINDOWS = GRID["n_windows"]
STITCH_H = COV["encounter_blocks"]["stitch_time_interval_hours"]

EXPOSURE = COV["exposure"]
INF_HOLD_H = EXPOSURE["infusion"]["hold_hours"]
GRID_MIN = EXPOSURE["grid_resolution_minutes"]

NEE = COV["time_varying"]["nee"]
NEE_COEF = NEE["coefficients"]
NEE_HOLD_H = NEE["hold_hours"]
NEE_PREFERRED = {**{d: "mcg/kg/min" for d in NEE_COEF}, "vasopressin": "u/min"}

OXY = COV["time_varying"]["oxygenation"]
FIO2_LOOKBACK_H = OXY["fio2_lookback_hours"]
SPO2_CEILING = OXY["spo2_ceiling"]
RA_FIO2 = OXY["fio2_scale"]["fraction_band"][0]
SPAN_LEAD_H = OXY["waterfall_span"]["lead_hours"]
SPAN_TRAIL_H = OXY["waterfall_span"]["trail_hours"]

LAB_VARS = {
    k: v["source"]["category"]
    for k, v in COV["time_varying"].items()
    if not k.startswith("_") and (v.get("source") or {}).get("table") == "labs"
}
ZERO_VARS = COV["missing_values"]["absence_means_zero"]["members"]
NOT_VENT_VARS = COV["missing_values"]["absence_means_not_ventilated"]["members"]

STROBE: list[tuple[str, int]] = []


def note(label: str, n: int) -> None:
    STROBE.append((label, n))
    print(f"  {label:.<52} {n:,}")


# --------------------------------------------------------------------- loading
LAB_NEEDED = sorted(set(LAB_VARS.values()) | {"po2_arterial", "creatinine", "platelet_count"})
VITAL_NEEDED = ["spo2", "map", "weight_kg", "height_cm"]
ASSESS_NEEDED = ["gcs_total"]


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
    return t, raw_anchor


def load_cohort_tables(hosp_ids: list[str], span: pd.DataFrame | None = None) -> dict:
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
    rs = RespiratorySupport.from_file(**_kw(filters=hid))
    t["resp_raw"] = rs.df.copy()
    n_raw = len(rs.df)
    if span is not None:
        d = rs.df.merge(span, on="hospitalization_id", how="inner")
        keep = (d["recorded_dttm"] >= d["span_lo"]) & (d["recorded_dttm"] <= d["span_hi"])
        rs.df = d[keep].drop(columns=["span_lo", "span_hi"])
        print(f"  {'respiratory_support':.<52} {n_raw:,} rows -> "
              f"{len(rs.df):,} in span ({100*len(rs.df)/max(n_raw,1):.1f}%)")
    else:
        print(f"  {'respiratory_support':.<52} {n_raw:,} rows  (untrimmed)")
    t["resp"] = _canonicalise_devices(rs.waterfall(verbose=False, return_dataframe=True))
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
    note("hospitalizations", len(mapping))
    note("encounter blocks after stitching", mapping["encounter_block"].nunique())

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


def find_anchor(t: dict, mapping: pd.DataFrame) -> pd.DataFrame:
    """First IMV record per encounter block."""
    rs = t["resp"].merge(mapping, on="hospitalization_id", how="inner")
    imv = rs[rs["device_category"] == "IMV"]
    anchor = (imv.groupby("encounter_block")["recorded_dttm"].min()
                 .rename("anchor_dttm").reset_index())
    note("blocks with any IMV record", len(anchor))
    return anchor


def build_cohort(blocks: pd.DataFrame, anchor: pd.DataFrame) -> pd.DataFrame:
    c = blocks.merge(anchor, on="encounter_block", how="inner", validate="one_to_one")
    note("blocks with an intubation anchor", len(c))
    if c.empty:
        raise SystemExit(
            "no block has an intubation anchor. An empty cohort is a bug, not a "
            "finding -- check that device_category still carries its mCIDE casing."
        )
    c = c[c["age"] >= CONFIG["cohort"]["min_age"]]
    note(f"adult blocks (age >= {CONFIG['cohort']['min_age']})", len(c))

    c["followup_end_dttm"] = c[["block_discharge_dttm"]].min(axis=1)
    c = c[c["followup_end_dttm"] > c["anchor_dttm"]]
    note("blocks with follow-up after the anchor", len(c))

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
    g["at_risk"] = g["win_start"] < g["followup_end_dttm"]
    note("patient-window rows", len(g))
    note("at-risk rows", int(g["at_risk"].sum()))
    return g


# -------------------------------------------------------------------- exposure
def _hourly_scaffold(cohort: pd.DataFrame) -> pd.DataFrame:
    hours = np.arange(0, EXTENT_H)
    s = cohort.loc[:, ["encounter_block", "anchor_dttm", "followup_end_dttm"]]
    s = s.loc[s.index.repeat(len(hours))].reset_index(drop=True)
    s["hr"] = np.tile(hours, len(cohort))
    s["cell_dttm"] = s["anchor_dttm"] + pd.to_timedelta(s["hr"], unit="h")
    return s[s["cell_dttm"] < s["followup_end_dttm"]].copy()


def infusion_grid(t: dict, mapping: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    """Hourly grid of fentanyl infusion rate, mcg/kg/hr. Config: exposure.infusion."""
    cats = CONFIG["medications"]["opioid_infusion_categories"]
    m = t["mac"][t["mac"]["med_category"].isin(cats)].merge(
        mapping, on="hospitalization_id", how="inner")
    m = m[m["encounter_block"].isin(cohort["encounter_block"])]

    m, rep = apply_med_raw(m, "med_category", "med_dose", "med_dose_unit", config=OUTLIERS)
    print(rep)

    m = m.merge(cohort[["encounter_block", "anchor_dttm", "weight_kg"]],
                on="encounter_block", how="left")
    m["rate"] = _to_mcg_kg_hr(m)
    # A charted stop is a rate of zero, not a missing value.
    stopped = m.get("mar_action_category", pd.Series(index=m.index, dtype=object))
    m.loc[stopped.astype("string").str.lower() == "stop", "rate"] = 0.0

    m["hr"] = ((m["admin_dttm"] - m["anchor_dttm"]).dt.total_seconds() // 3600).astype("Int64")
    m = m[(m["hr"] >= 0) & (m["hr"] < EXTENT_H) & m["rate"].notna()]
    last = (m.sort_values(["encounter_block", "hr", "rate"], kind="stable")
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


def _to_mcg_kg_hr(m: pd.DataFrame) -> pd.Series:
    """Convert charted infusion dose to mcg/kg/hr using the charted unit."""
    unit = m["med_dose_unit"].astype("string").str.lower().str.strip()
    dose = pd.to_numeric(m["med_dose"], errors="coerce")
    w = pd.to_numeric(m["weight_kg"], errors="coerce")
    out = pd.Series(np.nan, index=m.index, dtype="float64")
    out[unit == "mcg/kg/hr"] = dose
    out[unit == "mcg/hr"] = dose / w
    out[unit == "mg/hr"] = dose * 1000.0 / w
    out[unit == "mcg/kg/min"] = dose * 60.0
    out[unit == "mcg/min"] = dose * 60.0 / w
    unresolved = dose.notna() & out.isna()
    if unresolved.any():
        bad = sorted(unit[unresolved].dropna().unique())
        raise SystemExit(
            f"fentanyl infusion: {int(unresolved.sum()):,} rows in unhandled units {bad}. "
            f"Add them to _to_mcg_kg_hr rather than dropping -- a dropped unit removes "
            f"exposure for whichever patients were charted that way."
        )
    return out


def bolus_doses(t: dict, mapping: pd.DataFrame, grid: pd.DataFrame,
                cohort: pd.DataFrame) -> pd.DataFrame:
    """Window bolus total / weight / window_hours. Summed, never carried forward."""
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
    agg = b.groupby(["encounter_block", "window_idx"], as_index=False).agg(
        bolus_mcg=("mcg", "sum"), n_bolus=("mcg", "size"))
    agg = agg.merge(cohort[["encounter_block", "weight_kg"]], on="encounter_block", how="left")
    agg["bolus_dose"] = agg["bolus_mcg"] / agg["weight_kg"] / WINDOW_H
    return agg[["encounter_block", "window_idx", "bolus_dose", "n_bolus"]]


def window_exposure(grid: pd.DataFrame, bolus: pd.DataFrame,
                    windows: pd.DataFrame) -> pd.DataFrame:
    g = grid.copy()
    g["window_idx"] = g["hr"] // WINDOW_H
    inf = g.groupby(["encounter_block", "window_idx"], as_index=False)["rate"].mean()
    inf = inf.rename(columns={"rate": "inf_dose"})

    out = windows.merge(inf, on=["encounter_block", "window_idx"], how="left")
    out = out.merge(bolus, on=["encounter_block", "window_idx"], how="left")
    for c in ("inf_dose", "bolus_dose", "n_bolus"):
        out[c] = out[c].fillna(0.0)
    out["total_dose"] = out["inf_dose"] + out["bolus_dose"]

    ceiling = fentanyl_sanity_ceiling(OUTLIERS)
    over = out.loc[out["at_risk"], "total_dose"] > ceiling
    if over.any():
        raise SystemExit(
            f"{int(over.sum()):,} windows exceed the derived ceiling of {ceiling:g} "
            f"mcg/kg/hr, which is proof the bounds did not run, not a clinical finding."
        )
    note("at-risk windows with any fentanyl", int((out.loc[out['at_risk'], 'total_dose'] > 0).sum()))
    note("at-risk windows with a bolus", int((out.loc[out['at_risk'], 'n_bolus'] > 0).sum()))
    return out


# ------------------------------------------------------------------ covariates
def _window_of(dttm: pd.Series, anchor: pd.Series) -> pd.Series:
    return ((dttm - anchor).dt.total_seconds() // (WINDOW_H * 3600)).astype("Int64")


def _to_windows(df: pd.DataFrame, cohort: pd.DataFrame, time_col: str) -> pd.DataFrame:
    d = df.merge(cohort[["encounter_block", "anchor_dttm"]], on="encounter_block", how="inner")
    d["window_idx"] = _window_of(d[time_col], d["anchor_dttm"])
    return d[(d["window_idx"] >= 0) & (d["window_idx"] < N_WINDOWS)]


def lab_covariates(t: dict, mapping: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    labs = t["labs"].merge(mapping, on="hospitalization_id", how="inner")
    labs, rep = apply_long(labs, "labs", "lab_category", "lab_value_numeric", config=OUTLIERS)
    print(rep)
    labs = _to_windows(labs, cohort, "lab_result_dttm")

    out = None
    for var, cat in LAB_VARS.items():
        how = COV["time_varying"][var]["summary"]
        sub = labs[labs["lab_category"] == cat]
        agg = (sub.groupby(["encounter_block", "window_idx"], as_index=False)
                  ["lab_value_numeric"].agg(how).rename(columns={"lab_value_numeric": var}))
        out = agg if out is None else out.merge(
            agg, on=["encounter_block", "window_idx"], how="outer")
    return out if out is not None else pd.DataFrame()


def nee_covariate(t: dict, mapping: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    """Max of the summed vasopressor step function. Config: time_varying.nee."""
    from clifpy.utils.unit_converter import convert_dose_units_by_med_category

    m = t["mac"][t["mac"]["med_category"].isin(NEE_COEF)].merge(
        mapping, on="hospitalization_id", how="inner")
    m = m[m["encounter_block"].isin(cohort["encounter_block"])]
    if m.empty:
        return pd.DataFrame(columns=["encounter_block", "window_idx", "nee"])

    m, rep = apply_med_raw(m, "med_category", "med_dose", "med_dose_unit", config=OUTLIERS)
    print(rep)
    conv, _ = convert_dose_units_by_med_category(
        m, vitals_df=t["vitals"], preferred_units=NEE_PREFERRED)

    # clifpy leaves the raw value in place when it cannot convert; check the unit string.
    want = conv["med_category"].map(NEE_PREFERRED).astype("string").str.lower()
    got = conv["med_dose_unit_converted"].astype("string").str.lower()
    bad = conv["med_dose"].notna() & (got.isna() | (got != want))
    if bad.any():
        counts = conv.loc[bad].groupby(["med_category", "med_dose_unit_converted"]).size()
        raise SystemExit(f"unit conversion failed and was not nulled:\n{counts}")

    conv = conv.rename(columns={"med_dose_converted": "dose_std"})
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


def oxygenation_covariate(t: dict, mapping: pd.DataFrame,
                          cohort: pd.DataFrame) -> pd.DataFrame:
    """One column on the P/F scale, Severinghaus fallback. Config: time_varying.oxygenation."""
    rs = t["resp"].merge(mapping, on="hospitalization_id", how="inner")
    rs = rs[rs["encounter_block"].isin(cohort["encounter_block"])].copy()
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

    labs = t["labs"].merge(mapping, on="hospitalization_id", how="inner")
    labs, _ = apply_long(labs, "labs", "lab_category", "lab_value_numeric", config=OUTLIERS)
    pao2 = labs.loc[labs["lab_category"] == "po2_arterial",
                    ["encounter_block", "lab_result_dttm", "lab_value_numeric"]]
    pf = pair(pao2, "lab_result_dttm")
    pf["ratio"] = pf["lab_value_numeric"] / pf["fio2_set"]

    vit = t["vitals"].merge(mapping, on="hospitalization_id", how="inner")
    vit, _ = apply_long(vit, "vitals", "vital_category", "vital_value", config=OUTLIERS)
    spo2 = vit.loc[(vit["vital_category"] == "spo2")
                   & (vit["vital_value"] < SPO2_CEILING),
                   ["encounter_block", "recorded_dttm", "vital_value"]]
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

    use_pf = out["pf_ratio"].notna()
    out["oxygenation"] = out["pf_ratio"].where(use_pf, out["sf_ratio"])
    out["oxygenation_source"] = np.select(
        [use_pf & ~out["pf_ratio_all_assumed"].fillna(False),
         use_pf & out["pf_ratio_all_assumed"].fillna(False),
         out["sf_ratio"].notna() & ~out["sf_ratio_all_assumed"].fillna(False),
         out["sf_ratio"].notna() & out["sf_ratio_all_assumed"].fillna(False)],
        ["pf", "pf_room_air", "sf", "sf_room_air"], default="none")
    return out[["encounter_block", "window_idx", "oxygenation", "oxygenation_source",
                "pf_ratio", "sf_ratio"]]


def _severinghaus(spo2: pd.Series) -> pd.Series:
    s = pd.to_numeric(spo2, errors="coerce") / 100.0
    s = s.where((s > 0) & (s < 1))
    a = 11700.0 / ((1.0 / s) - 1.0)
    b = np.sqrt(50.0 ** 3 + a ** 2)
    return np.cbrt(b + a) - np.cbrt(b - a)


SOFA_PRESSORS = ["norepinephrine", "epinephrine", "dopamine", "dobutamine"]


def _sofa_inputs(t: dict, mapping: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    """Per-window SOFA components not already carried as covariates."""
    labs = t["labs"].merge(mapping, on="hospitalization_id", how="inner")
    labs, _ = apply_long(labs, "labs", "lab_category", "lab_value_numeric", config=OUTLIERS)
    labs = _to_windows(labs, cohort, "lab_result_dttm")
    plt_ = (labs[labs["lab_category"] == "platelet_count"]
            .groupby(["encounter_block", "window_idx"], as_index=False)["lab_value_numeric"]
            .min().rename(columns={"lab_value_numeric": "platelet_count"}))
    creat = (labs[labs["lab_category"] == "creatinine"]
             .groupby(["encounter_block", "window_idx"], as_index=False)["lab_value_numeric"]
             .max().rename(columns={"lab_value_numeric": "creatinine"}))

    vit = t["vitals"].merge(mapping, on="hospitalization_id", how="inner")
    vit, _ = apply_long(vit, "vitals", "vital_category", "vital_value", config=OUTLIERS)
    vit = _to_windows(vit, cohort, "recorded_dttm")
    mp = (vit[vit["vital_category"] == "map"]
          .groupby(["encounter_block", "window_idx"], as_index=False)["vital_value"]
          .min().rename(columns={"vital_value": "map"}))

    asm = t["assessments"].merge(mapping, on="hospitalization_id", how="inner")
    asm, _ = apply_long(asm, "patient_assessments", "assessment_category",
                        "numerical_value", config=OUTLIERS)
    asm = _to_windows(asm, cohort, "recorded_dttm")
    gcs = (asm[asm["assessment_category"] == "gcs_total"]
           .groupby(["encounter_block", "window_idx"], as_index=False)["numerical_value"]
           .min().rename(columns={"numerical_value": "gcs_total"}))

    out = plt_
    for part in (creat, mp, gcs):
        out = out.merge(part, on=["encounter_block", "window_idx"], how="outer")
    return out


def _sofa_pressors(t: dict, mapping: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    """Max mcg/kg/min per window for the four SOFA vasopressors.

    Absent drugs stay NULL. Coding them 0 makes the `<= 0.1` clause in the
    cardiovascular score fire, giving every unpressored patient a score of 3.
    """
    from clifpy.utils.unit_converter import convert_dose_units_by_med_category

    m = t["mac"][t["mac"]["med_category"].isin(SOFA_PRESSORS)].merge(
        mapping, on="hospitalization_id", how="inner")
    m = m[m["encounter_block"].isin(cohort["encounter_block"])]
    cols = ["encounter_block", "window_idx"] + [f"{d}_mcg_kg_min" for d in SOFA_PRESSORS]
    if m.empty:
        return pd.DataFrame(columns=cols)

    m, _ = apply_med_raw(m, "med_category", "med_dose", "med_dose_unit", config=OUTLIERS)
    conv, _ = convert_dose_units_by_med_category(
        m, vitals_df=t["vitals"],
        preferred_units={d: "mcg/kg/min" for d in SOFA_PRESSORS})
    got = conv["med_dose_unit_converted"].astype("string").str.lower()
    conv = conv[got.eq("mcg/kg/min")].rename(columns={"med_dose_converted": "dose_std"})
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


def score_sofa(df: pd.DataFrame) -> pd.DataFrame:
    """Six SOFA components and their total. Vincent 1996; see design_notes.md §11."""
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


def status_covariates(t: dict, mapping: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    rs = _to_windows(t["resp"].merge(mapping, on="hospitalization_id", how="inner"),
                     cohort, "recorded_dttm")
    imv = (rs.assign(v=rs["device_category"].eq("IMV"))
             .groupby(["encounter_block", "window_idx"], as_index=False)["v"].max()
             .rename(columns={"v": "imv_status"}))

    crrt = t["crrt"].merge(mapping, on="hospitalization_id", how="inner")
    crrt = _to_windows(crrt, cohort, "recorded_dttm")
    cr = (crrt.assign(v=1).groupby(["encounter_block", "window_idx"], as_index=False)["v"]
              .max().rename(columns={"v": "crrt_status"}))
    return imv.merge(cr, on=["encounter_block", "window_idx"], how="outer")


def attach_weight(t: dict, mapping: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    """Dose denominator, fixed at the anchor. Backward match preferred; report the lag."""
    vit = t["vitals"].merge(mapping, on="hospitalization_id", how="inner")
    vit, _ = apply_long(vit, "vitals", "vital_category", "vital_value", config=OUTLIERS)
    w = vit.loc[vit["vital_category"] == "weight_kg",
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


def bmi_admission(t: dict, mapping: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    vit = t["vitals"].merge(mapping, on="hospitalization_id", how="inner")
    vit, _ = apply_long(vit, "vitals", "vital_category", "vital_value", config=OUTLIERS)
    vit = vit.merge(cohort[["encounter_block", "block_admission_dttm"]],
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

    Returns (long, per_variable, per_pattern). Counts are over at-risk rows only.
    """
    df = long.copy()
    at_risk = df["at_risk"]
    n_at_risk = int(at_risk.sum())

    tv = [k for k in COV["time_varying"] if not k.startswith("_") and k in df.columns]
    ti = [k for k in COV["time_invariant"] if not k.startswith("_") and k in df.columns]
    exposure = [c for c in EXPOSURE["columns"] if c in df.columns]

    observed = {v: int(df.loc[at_risk, v].notna().sum()) for v in tv}

    zeroed = {}
    for v in ZERO_VARS + NOT_VENT_VARS:
        if v in df.columns:
            blank = at_risk & df[v].isna()
            zeroed[v] = int(blank.sum())
            df.loc[blank, v] = 0.0

    pre = {v: int(df.loc[at_risk, v].isna().sum()) for v in tv}

    locf_caps = {
        k: (v["locf"] or {}).get("cap_hours")
        for k, v in COV["time_varying"].items()
        if not k.startswith("_") and (v.get("locf") or {}).get("eligible") and k in df.columns
    }
    df = df.sort_values(["encounter_block", "window_idx"], kind="stable")
    filled = {}
    for v, cap in locf_caps.items():
        limit = max(int(cap // WINDOW_H), 1)
        before = df[v].isna()
        df[v] = df.groupby("encounter_block")[v].ffill(limit=limit)
        df[f"{v}_locf"] = before & df[v].notna()
        filled[v] = int((df[f"{v}_locf"] & at_risk).sum())

    post = {v: int(df.loc[at_risk, v].isna().sum()) for v in tv}

    def pct(n: int) -> float:
        return round(100.0 * n / n_at_risk, 2) if n_at_risk else float("nan")

    rows = []
    for v in tv:
        cap = locf_caps.get(v)
        rows.append({
            "variable": v,
            "kind": "time_varying",
            "class": COV["time_varying"][v]["missing_class"],
            "locf_cap_hours": cap if cap is not None else "",
            "n_at_risk": n_at_risk,
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
        n_miss = int(df.loc[at_risk, v].isna().sum()) if v in exposure \
            else int(df[v].isna().sum())
        denom = n_at_risk if v in exposure else len(df)
        rows.append({
            "variable": v,
            "kind": "exposure" if v in exposure else "time_invariant",
            "class": "", "locf_cap_hours": "",
            "n_at_risk": denom,
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
        src = (df.loc[at_risk, "oxygenation_source"].fillna("none")
                 .value_counts().rename_axis("oxygenation_source")
                 .reset_index(name="n"))
        src["pct"] = (100.0 * src["n"] / n_at_risk).round(2)
        print("\n  oxygenation provenance (at-risk windows)")
        print(src.to_string(index=False))

    per_pattern = _missingness_patterns(df[at_risk], tv)
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


def build_time_to_event(cohort: pd.DataFrame, long: pd.DataFrame) -> pd.DataFrame:
    """One row per block, origin = landmark T. Two event codings; see §10, §10a."""
    T = CONFIG["cohort"]["landmark_hours"]
    succ_h = CONFIG["outcomes"]["successful_extubation_hours"]
    death_cats = set(CONFIG["outcomes"]["mortality_categories"])

    at_T = long[(long["window_start_hr"] == T - WINDOW_H) & long["at_risk"]]
    eligible = set(at_T.loc[at_T["imv_status"] == 1, "encounter_block"])
    note(f"landmark-eligible blocks (ventilated at T={T}h)", len(eligible))

    tte = cohort[cohort["encounter_block"].isin(eligible)].copy()
    tte["landmark_eligible"] = True

    post = long[(long["window_start_hr"] >= T) & long["at_risk"]]
    vent = post[post["imv_status"] == 1]
    last_vent = vent.groupby("encounter_block")["window_start_hr"].max()
    tte["last_vent_hr"] = tte["encounter_block"].map(last_vent)

    tte["died"] = tte["discharge_category"].isin(death_cats)
    tte["hours_to_discharge"] = (
        (tte["block_discharge_dttm"] - tte["anchor_dttm"]).dt.total_seconds() / 3600)

    ext_hr = tte["last_vent_hr"] + WINDOW_H
    off_long_enough = (tte["hours_to_discharge"] - ext_hr) >= succ_h
    tte["extubation_event"] = np.select(
        [tte["died"] & ~off_long_enough.fillna(False), off_long_enough.fillna(False)],
        [2, 1], default=0)
    tte["extubation_time"] = np.where(
        tte["extubation_event"] == 1, ext_hr - T,
        tte["hours_to_discharge"] - T)

    tte["mortality_event"] = np.where(tte["died"], 1, 2)
    tte["mortality_time"] = (tte["hours_to_discharge"] - T).clip(lower=0)

    unresolved = set(CONFIG["outcomes"]["unresolved_discharge_categories"])
    n_unres = int(tte["discharge_category"].isin(unresolved).sum())
    if n_unres:
        note("blocks with an unresolved discharge_category", n_unres)
    for code, label in ((1, "successful extubation"), (2, "death")):
        note(f"  extubation outcome = {label}", int((tte["extubation_event"] == code).sum()))
    return tte


def main() -> None:
    dirs = site_dirs(REPO)
    prov = provenance(CONFIG)
    print(f"site {prov['site_name']}  clif {prov['clif_version']}  code {prov['code_version']}")
    print(f"grid {WINDOW_H}h x {N_WINDOWS} to {EXTENT_H}h   stitch {STITCH_H}h   "
          f"infusion hold {INF_HOLD_H}h\n")

    print("Loading core tables")
    t, raw_anchor = load_core()

    print("\nEncounter blocks")
    blocks, mapping = build_blocks(t)
    hi = hospital_intervals(t, mapping)
    ends = hospital_endpoints(hi)

    # Restrict before the waterfall: it is the expensive step and only the
    # cohort's rows can affect the result.
    imv_ids = set(raw_anchor["hospitalization_id"])
    imv_blocks = set(mapping.loc[mapping["hospitalization_id"].isin(imv_ids), "encounter_block"])
    blocks = blocks[blocks["encounter_block"].isin(imv_blocks)]
    note("blocks containing any IMV record", len(blocks))
    mapping = mapping[mapping["encounter_block"].isin(imv_blocks)]
    hosp_ids = sorted(mapping["hospitalization_id"].astype(str).unique())
    note("hospitalizations to load", len(hosp_ids))

    # The waterfall is the whole cost of a run; trim it to the analysis span.
    # Equivalence measured in validation/waterfall_span_equivalence.py.
    blk_anchor = (raw_anchor.merge(mapping, on="hospitalization_id", how="inner")
                            .groupby("encounter_block", as_index=False)["imv_dttm"].min())
    span = mapping.merge(blk_anchor, on="encounter_block", how="inner")
    span["span_lo"] = span["imv_dttm"] - pd.Timedelta(hours=SPAN_LEAD_H)
    span["span_hi"] = span["imv_dttm"] + pd.Timedelta(hours=EXTENT_H + SPAN_TRAIL_H)
    span = span[["hospitalization_id", "span_lo", "span_hi"]]

    print("\nLoading cohort tables")
    t.update(load_cohort_tables(hosp_ids, span=span))
    assert_categories_present(t)

    print("\nCohort")
    anchor = find_anchor(t, mapping)
    cohort = build_cohort(blocks, anchor)

    print("\nWeight and BMI")
    cohort = cohort.merge(attach_weight(t, mapping, cohort), on="encounter_block", how="left")
    cohort = cohort[cohort["weight_kg"].notna()]
    note("blocks with a usable weight", len(cohort))
    bmi = bmi_admission(t, mapping, cohort)

    print("\nWindows")
    windows = window_grid(cohort)

    print("\nExposure")
    grid = infusion_grid(t, mapping, cohort)
    bolus = bolus_doses(t, mapping, grid, cohort)
    long = window_exposure(grid, bolus, windows)

    print("\nCovariates")
    for part in (lab_covariates(t, mapping, cohort),
                 nee_covariate(t, mapping, cohort),
                 oxygenation_covariate(t, mapping, cohort),
                 status_covariates(t, mapping, cohort)):
        if len(part):
            long = long.merge(part, on=["encounter_block", "window_idx"], how="left")

    hi = hi[hi["encounter_block"].isin(cohort["encounter_block"])]
    ti = time_invariant(t, mapping, cohort).merge(bmi, on="encounter_block", how="left")
    ti = ti.merge(ends, on="encounter_block", how="left")
    long = long.merge(ti.drop(columns=["patient_id"]), on="encounter_block", how="left")
    long = long.merge(cohort[["encounter_block", "weight_kg"]], on="encounter_block",
                      how="left")

    print("\n  SOFA")
    sofa_in = _sofa_inputs(t, mapping, cohort)
    press = _sofa_pressors(t, mapping, cohort)
    if len(press):
        sofa_in = sofa_in.merge(press, on=["encounter_block", "window_idx"], how="outer")
    long = long.merge(sofa_in, on=["encounter_block", "window_idx"], how="left")
    long = pd.concat([long.reset_index(drop=True),
                      score_sofa(long).reset_index(drop=True)], axis=1)
    comp = long.loc[long["at_risk"], "sofa_n_components"]
    note("at-risk windows with a complete 6-component SOFA", int((comp == 6).sum()))
    note("at-risk windows with no SOFA component at all", int((comp == 0).sum()))

    print("\nMissingness and LOCF")
    long, per_variable, per_pattern = apply_missingness(long)
    cols = ["variable", "kind", "locf_cap_hours", "n_observed", "n_zero_by_rule",
            "pct_missing_pre_locf", "pct_filled_by_locf", "pct_missing_final"]
    print("\n  missingness, at-risk rows only")
    print(per_variable[cols].to_string(index=False))

    assert_config_is_honoured(long)

    print("\nTime to event")
    tte = build_time_to_event(cohort, long)

    out = dirs["data_phi"]
    long.to_parquet(out / "trajectory_long.parquet", index=False)
    tte.to_parquet(out / "time_to_event.parquet", index=False)
    hi.to_parquet(out / "hospital_intervals.parquet", index=False)
    per_variable.to_csv(dirs["out_final"] / "phase0_missingness.csv", index=False)
    if len(per_pattern):
        per_pattern.to_csv(dirs["out_final"] / "phase0_missingness_patterns.csv", index=False)
    pd.DataFrame(STROBE, columns=["step", "n"]).to_csv(
        dirs["out_final"] / "phase0_strobe.csv", index=False)
    (dirs["out_final"] / "phase0_provenance.json").write_text(json.dumps(prov, indent=2))

    print(f"\nwritten to {out}")
    print(f"  trajectory_long.parquet  {len(long):,} rows x {long.shape[1]} cols")
    print(f"  time_to_event.parquet    {len(tte):,} rows")
    print(f"  hospital_intervals.parquet {len(hi):,} rows")


if __name__ == "__main__":
    main()
