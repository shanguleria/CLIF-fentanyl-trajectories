# CLIF-fentanyl-trajectories

Group-based trajectory modeling of fentanyl dosing in mechanically ventilated
ICU patients, using the Common Longitudinal ICU data Format (CLIF).

## Objective

Identify and characterize distinct trajectories of fentanyl exposure over the
first days of invasive mechanical ventilation (IMV), and test whether trajectory
class predicts liberation from ventilation and survival.

The pipeline is both **descriptive** (cohort-level dosing curves with
balanced-panel overlays) and **model-based** (group-based trajectory models via
`gbmt`, latent class mixed models via `lcmm`, followed by competing-risks
outcome analysis). Design rationale, methodological decisions, and the evidence
behind them are in [`docs/design_notes.md`](docs/design_notes.md).

## Required CLIF tables and fields

Built against the **CLIF 2.1.0** specification. Field names and permissible
values below were verified against `clifpy`'s bundled schemas.

### Shared by all phases

| Table | Fields | Notes |
|---|---|---|
| `patient` | `patient_id`, `sex_category`, `race_category`, `ethnicity_category`, `birth_date` | `death_dttm` exists but is **unstable across CLIF sites** and is deliberately not used; mortality comes from `hospitalization.discharge_category` |
| `hospitalization` | `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission`, `admission_type_category`, `discharge_category` | mortality from `discharge_category` ∈ {`Expired`, `Hospice`} |
| `adt` | `hospitalization_id`, `in_dttm`, `out_dttm`, `location_category` | ICU location, length of stay |

### Exposure — fentanyl dosing

| Table | Fields | Notes |
|---|---|---|
| `medication_admin_continuous` | `hospitalization_id`, `admin_dttm`, `med_category`, `med_dose`, `med_dose_unit`, `med_route_category`, `mar_action_category` | infusion rates; LOCF between charted changes |
| `medication_admin_intermittent` | `hospitalization_id`, `admin_dttm`, `med_category`, `med_dose`, `med_dose_unit`, `mar_action_category` | boluses; **summed within window, never carried forward** |
| `vitals` | `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value` | weight, for mcg/kg/hr normalization |

Propofol and midazolam are extracted into their own columns from the same
medication tables, for later use.

### Ventilation status and outcomes

| Table | Fields | Categorical values used |
|---|---|---|
| `respiratory_support` | `hospitalization_id`, `recorded_dttm`, `device_category`, `mode_category`, `tracheostomy` | `device_category` = `IMV` defines ventilated windows; `Trach Collar` and the `tracheostomy` (0/1) flag identify tracheostomy |
| `patient_procedures` | `hospitalization_id`, `procedure_code`, `procedure_code_format`, `procedure_billed_dttm` | `ICD10PCS` / `CPT` cross-check on tracheostomy timing |

### Covariates

Every covariate's window-aggregation rule and missingness class is declared in
**[`config/covariates.json`](config/covariates.json)** and mirrored in design
notes §11. That file is the source of truth; `tests/test_covariates.py` enforces
its internal consistency.

| Table | Fields | Used for |
|---|---|---|
| `labs` | `hospitalization_id`, `lab_result_dttm`, `lab_category`, `lab_value_numeric` | `bun`, `bicarbonate`, `pco2_arterial`, `lactate`, `inr`, `bilirubin_total`, `po2_arterial`, `creatinine`, `platelet_count` |
| `vitals` | `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value` | `spo2`, `map` (SOFA); `weight_kg`, `height_cm` (BMI, dose denominator) |
| `medication_admin_continuous` | `hospitalization_id`, `admin_dttm`, `med_category`, `med_dose`, `med_dose_unit`, `mar_action_category` | NEE: norepinephrine, epinephrine, phenylephrine, dopamine, vasopressin, angiotensin |
| `respiratory_support` | **all columns** | FiO₂ waterfall + `device_category` for `imv_status`. Load every column — the waterfall needs `device_name`, `mode_category`, `lpm_set`, `peep_set` to build its blocks |
| `crrt_therapy` | `hospitalization_id`, `recorded_dttm`, `crrt_mode_category` | `crrt_status`. **Point-in-time table, no start/stop** — interval reconstruction, not a lookup |
| `hospital_diagnosis` | `hospitalization_id`, `diagnosis_code`, `diagnosis_code_format`, `poa_present` | Charlson Comorbidity Index |
| `patient_assessments` | `hospitalization_id`, `recorded_dttm`, `assessment_category`, `numerical_value` | `gcs_total` (SOFA CNS); RASS is ordinal, so **not** averaged |
| `patient` | `patient_id`, `sex_category`, `race_category` | sex; **race — Table 1 reporting only, not a model covariate** |
| `adt` | `hospitalization_id`, `hospital_id`, `hospital_type`, `in_dttm`, `out_dttm`, `location_category` | `hospital_id_admission` / `_discharge`; **`hospital_id` is required by `stitch_encounters`** |

**Encounter blocks are the unit of analysis.** Hospitalizations are stitched with
clifpy `stitch_encounters(..., time_interval=6)`: a `hospitalization_id` is one
encounter, not one clinical course, and a patient intubated, transferred, and
still intubated has one ventilation episode. Without stitching that trajectory is
truncated at the transfer. See design notes §11.

**Two clifpy cautions**, both verified against 0.3.8 in this repo's `.venv` and
documented in `covariates.json`:

- `compute_sofa_polars` must be called with `fill_na_scores_with_zero=False`, and
  its cardiovascular component ignores vasopressin, phenylephrine, angiotensin
  and milrinone. SOFA is the convenience summary here; the explicit markers are
  the severity measure.
- `convert_dose_units_by_med_category` does **not** null a dose it cannot
  convert — it returns the raw value and reports the failure in the unit string.
  Guard the unit string and raise. See design notes §11.

## Cohort identification

- **Population**: adults (≥18) receiving invasive mechanical ventilation.
- **Time anchor**: intubation (first `device_category == "IMV"`).
- **Unit of analysis**: one row per hospitalization per time window.
- **Descriptive phase**: all intubated patients, contributing while ventilated.
- **Modeling phases**: landmark cohort — alive and still ventilated at `T`
  (default 72h). Patients extubated and reintubated before `T` remain in the
  cohort; dose in extubated windows is `0`. Rationale in design notes §8/§10a.
- **Outputs**: parquet (patient-level, PHI) and CSV (aggregate, shareable).

## Configuration

Copy the template and edit it for your site. `config/config.json` is
**gitignored** because it contains a local path.

```bash
cp config/config_template.json config/config.json
python3 code/check_config.py      # confirms it is usable before you run anything
```

One config file, one site, one data source — the standard CLIF layout. To run
against a different source (a MIMIC-to-CLIF conversion, say), edit `site_name`,
`data_directory` and `dataset_version` in this file; there is no second code path.

**`config/covariates.json` is a different kind of file and is committed.** It
holds the covariate protocol — every variable's source, window-aggregation rule,
LOCF cap and missingness class — and carries its own `definition_version`, which
is stamped onto shareable outputs. The test for what belongs there: *if two sites
set this differently, is the pooled result still meaningful?* If no, it is
protocol. Nothing estimand-defining may live only in the gitignored
`config.json`.

```bash
.venv/bin/python tests/test_covariates.py   # 16 checks on the covariate protocol
.venv/bin/python tests/test_fio2.py         # 10 checks on FiO2 unit handling
.venv/bin/python tests/test_outliers.py     # 14 checks on the outlier bounds
.venv/bin/python tests/test_build_cohort.py # 17 checks on Phase 0 logic
```

### Key settings

The first block is the consortium-standard CLIF header. `data_directory`,
`filetype` and `timezone` are the names **clifpy itself requires**, so this file
is handed to the library with no translation layer.

| Key | Meaning |
|---|---|
| `site_name` | site identifier stamped on every shareable output |
| `clif_version` | the CLIF **spec** version your tables implement |
| `dataset_version` | conversion/ETL release, if the tables came from one; empty for a site's own data |
| `data_directory` | absolute path to the directory holding `clif_<table>.parquet` |
| `filetype` | `parquet` or `csv`; passed straight to clifpy |
| `timezone` | your site's local zone. Every window is a wall-clock interval, so this shifts results if wrong |
| `tables_in_use` | CLIF tables loaded by Phase 0, and the set `check_config.py` verifies |
| `medications.*` | `med_category` values in **your** vocabulary — probe your data first |
| `cohort.window_hours` | aggregation window (default 4) |
| `cohort.landmark_hours` | landmark `T` (default 72) |
| `cohort.balanced_panel_hours` | thresholds for the Phase 1 overlays |
| `outcomes.mortality_source` | `discharge_category` (`death_dttm` is unstable) |
| `outcomes.mortality_categories` | values counted as death; sensitivity set drops `Hospice` |
| `outcomes.unresolved_discharge_categories` | censored, **not** counted alive |
| `model.scaling` | `gbmt` normalization — **0**, see design notes §3 |
| `model.nstart` | EM random restarts — ≥50, see design notes §4/E3 |

**`clif_version` is always the CLIF *specification* version.** A conversion
release such as the MIMIC CLIF conversion 1.1.0 belongs in `dataset_version` —
there is no CLIF spec version 1.1.0.

**`cohort`, `outcomes` and `model` define the estimand.** They are the same at
every site by design: if two sites set them differently the pooled result is not
meaningful. Do not edit them to make a run work. `site_name`, `data_directory`,
`filetype`, `timezone` and `medications` are the parts that legitimately differ.

### Where outputs land

| Path | Contents | Shareable? |
|---|---|---|
| `data/intermediate_phi/` | analytic tables, one row per patient-window | **No — PHI** |
| `output/intermediate_phi/` | model objects, class assignments | **No — PHI** |
| `output/final_no_phi/` | aggregate tables and figures — the coordinating-center upload set, PHI-checked at assembly | **Yes** |

`data/`, `output/intermediate_phi/`, and `logs/` are gitignored. Both PHI
directories get a `README.md` warning label written at runtime by
`site_dirs()`. Shareable outputs carry a provenance block from `provenance()`:
`site_name`, `clif_version`, `dataset_version`, `code_version` (git describe),
`generated`.

## Prerequisites

**Python** 3.11+ for data preparation (Phase 0), **R** 4.4+ for modeling
(Phases 1–6). `gbmt` and `lcmm` are R-only, which is why the pipeline is split.

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

Rscript -e 'install.packages("renv"); renv::restore()'
```

## Running the pipeline

**macOS / Linux**

```bash
./run_pipeline.sh
```

**Windows (PowerShell)**

```powershell
.\run_pipeline.ps1
```

Both runners run the preflight (`code/check_config.py`) first, then execute the
steps below in order, stopping on first error. The preflight confirms the config
is complete and that `data_directory` actually holds the tables named in
`tables_in_use` — it checks by filename and opens nothing, so a bad path fails in
a second rather than an hour into Phase 0. Run it on its own at any time:

```bash
python3 code/check_config.py
```

## Pipeline steps

| Step | Language | Script | Description |
|---|---|---|---|
| 0 | Python | `code/01_build_cohort.py` | Load CLIF tables via clifpy, build the windowed trajectory table and the time-to-event table |
| 1 | R | `code/02_descriptive_trajectory.R` | Cohort dose curves, balanced-panel overlays, retention table, zero fraction |
| 2 | R | `code/03_landmark_cohort.R` | Apply landmark `T`, report retention and failed-extubation counts |
| 3 | R | `code/04_gbmt_classes.R` | `gbmt` on combined dose, `ng` sweep, class enumeration |
| 4 | R | `code/05_lcmm_classes.R` | `lcmm::hlme` on the same data; compare partitions by ARI |
| 5 | R | `code/06_two_indicator.R` | Two-indicator infusion-vs-bolus strategy model |
| 6 | R | `code/07_outcomes.R` | Competing-risks outcome models from the landmark |

## Project structure

```
CLIF-fentanyl-trajectories/
├── README.md
├── requirements.txt              # Python, pinned
├── renv.lock                     # R, pinned
├── run_pipeline.sh / .ps1
├── config/
│   ├── covariates.json           # committed -- covariate PROTOCOL (definition_version)
│   ├── outlier_config.json       # committed -- the only place bounds are written
│   ├── config_template.json      # committed
│   └── config.json               # gitignored, site-local
├── code/
│   ├── check_config.py           # preflight: is config.json usable?
│   ├── 01_build_cohort.py        # Phase 0  (Python / clifpy)
│   ├── 02_descriptive_trajectory.R
│   ├── 03_landmark_cohort.R
│   ├── 04_gbmt_classes.R
│   ├── 05_lcmm_classes.R
│   ├── 06_two_indicator.R
│   ├── 07_outcomes.R
│   └── utils/
│       ├── fio2.py               # FiO2 scale detection + normalisation
│       ├── outliers.py           # applies config/outlier_config.json
│       ├── paths.R               # output dirs + provenance
│       ├── paths.py              #   (the two must agree)
│       └── dependencies.R        # package list for renv's scanner
├── tests/
│   ├── test_covariates.py        # integrity checks on config/covariates.json
│   ├── test_fio2.py              # FiO2 must be a fraction; enforced, not assumed
│   ├── test_outliers.py          # bounds are applied, and gaps are reported
│   └── test_build_cohort.py      # Phase 0 logic on synthetic frames
├── validation/                   # methodological evidence, synthetic data
│   ├── scaling_experiments.R
│   └── composition_bias_demo.R
├── docs/
│   └── design_notes.md           # protocol, decisions, evidence
├── data/intermediate_phi/        # gitignored -- PHI
├── output/
│   ├── intermediate_phi/         # gitignored -- PHI
│   └── final_no_phi/             # the shareable set
└── logs/                         # gitignored
```

## Definitions and provenance

- **Dose**: within-window total fentanyl (infusion + bolus), expressed
  **mcg/kg/hr**. Infusion uses LOCF then a time-weighted mean; boluses are summed
  and never carried forward.
- **Successful extubation**: extubation not followed by reintubation within 72h.
  Competes with **death** and **tracheostomy**. A failed extubation is not an
  event. Late extubations that cannot be confirmed are censored.
- **Mortality**: `hospitalization.discharge_category` ∈ {`Expired`, `Hospice`},
  sensitivity analysis with `Expired` alone. `patient.death_dttm` is **not** used
  — it is unstable across CLIF sites. Without linked vital status this is
  therefore *in-hospital mortality censored at day 30*, not 30-day all-cause
  mortality, and should be named that way. `Still Admitted` / `Missing` /
  `Other` are censored, never counted as alive-and-discharged.
- **Scaling**: `gbmt` normalizes **within unit**; any `scaling ≥ 1` erases
  absolute dose level. This pipeline uses `scaling = 0`. Evidence:
  `validation/scaling_experiments.R`.

Shareable outputs carry a provenance block from `provenance()`: `site_name`,
`clif_version`, `dataset_version`, `code_version` (`git describe --always
--dirty`, so a dirty tree is visible), `generated`.

## Contributing

R scripts follow the house style in
`~/Desktop/Research/CLIF/R_setup.R` — a Purpose/Author/Created/Inputs/Outputs
banner, the `requireNamespace` + `library(p, character.only = TRUE)` package
loop, `here()` for all paths (never `setwd()`), and a `sessionInfo()` provenance
write at the end of each script.

**Any generated code must be read by a human before it is run or committed.** See
[`docs/code_review_checklist.md`](docs/code_review_checklist.md).

## Data safety

Patient-level data never leaves `data/intermediate_phi/` and
`output/intermediate_phi/`, both gitignored and both carrying a warning label
written at runtime. Only `output/final_no_phi/` is intended for transfer; its
contents are aggregate by construction and PHI-checked at assembly.

## Onboarding a new site

1. Clone the repo; `cp config/config_template.json config/config.json`.
2. Set `site_name`, `clif_version`, `data_directory`, `filetype` and `timezone`.
3. Install pinned dependencies (see Prerequisites).
4. Run `python3 code/check_config.py` until it prints `READY.`
5. Probe your `med_category` values and update `config.medications`.
6. Run the pipeline.
7. Review `output/final_no_phi/` and send only that.
