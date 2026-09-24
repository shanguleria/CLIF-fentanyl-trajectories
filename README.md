# CLIF-fentanyl-trajectories

Group-based trajectory modeling of fentanyl dosing in mechanically ventilated
ICU patients, using the Common Longitudinal ICU data Format (CLIF).

## Objective

Describe **how fentanyl is delivered** to adults over the first days of invasive
mechanical ventilation (IMV) — continuous infusion, intermittent bolus, both, or
neither — and how that delivery pattern changes over the ventilation course.

The pipeline is **descriptive**: cohort-level dosing curves with balanced-panel
overlays, delivery-state prevalence and transitions, and a landmarked view of the
first 72 hours. The question, the measurement frame, the standing principles and
the options considered and rejected are in
[`docs/principles.md`](docs/principles.md); literature in
[`docs/references.md`](docs/references.md). **Definitions live in `config/`**, not
in either document.

> **Scope change, 2026-09-23.** This project began as a trajectory-modelling
> study: find latent classes of fentanyl exposure with `gbmt`/`lcmm`, then use
> class membership as an exposure. That approach was **tested and retired** —
> dose level proved continuous with no subpopulation structure (`gbmt` classes
> predicted an episode's mean dose at R² = 0.82; `lcmm` converged only at two
> classes; the two methods agreed at an adjusted Rand index of 0.037). Steps 3–6
> below remain in the repo and still run, but they are **tabled**, not part of the
> current analysis, and their outputs on disk predate their own scripts. Do not
> cite a number from them without re-running. See `.claude/claude-todo.md`.

## Required CLIF tables and fields

Built against the **CLIF 2.1.0** specification. Field names and permissible
values below were verified against `clifpy`'s bundled schemas.

### Shared by all phases

| Table | Fields | Notes |
|---|---|---|
| `patient` | `patient_id`, `sex_category`, `race_category`, `ethnicity_category`, `birth_date` | `death_dttm` exists but is **unstable across CLIF sites** and is deliberately not used; mortality comes from `hospitalization.discharge_category` |
| `hospitalization` | `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission`, `admission_type_category`, `discharge_category` | mortality from `discharge_category` ∈ {`Expired`, `Hospice`} |
| `adt` | `hospitalization_id`, `hospital_id`, `hospital_type`, `in_dttm`, `out_dttm` | encounter stitching and `hospital_id_admission` / `_discharge`. **`location_category` is not read** — the pipeline has no ICU flag |

### Exposure — fentanyl dosing

| Table | Fields | Notes |
|---|---|---|
| `medication_admin_continuous` | `hospitalization_id`, `admin_dttm`, `med_category`, `med_dose`, `med_dose_unit`, `med_route_category`, `mar_action_category` | infusion rates; LOCF between charted changes |
| `medication_admin_intermittent` | `hospitalization_id`, `admin_dttm`, `med_category`, `med_dose`, `med_dose_unit`, `mar_action_category` | boluses; **summed within window, never carried forward** |
| `vitals` | `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value` | weight, for `nee` and BMI normalization |

Propofol, midazolam and dexmedetomidine are extracted into their own columns
from `medication_admin_continuous` as **descriptive companions** to the fentanyl
exposure — infusions only, each in the unit it is ordered in (mcg/kg/min, mg/hr
and mcg/kg/hr respectively). They are not `gbmt` indicators: mixed units across
drugs would force `scaling >= 2`, which erases the absolute dose level the design
rests on. Add a drug by adding its category to
`medications.other_sedative_categories` **and** its column, unit, conversion and
outlier bound — a test fails if any of the four is missing.

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
| `patient_assessments` | `hospitalization_id`, `recorded_dttm`, `assessment_category`, `numerical_value` | `gcs_total` (SOFA CNS); `RASS` → `rass` (**min**, deepest sedation in the window); `NVPS` → `nvps` (**max**, worst pain). Both are ordinal, so neither is averaged — see design notes §11. **Category names are matched literally and are case-sensitive**: UCMC charts `RASS` and `NVPS` capitalised but `gcs_total` lowercase. A site charting CPOT rather than NVPS needs its own entry in `covariates.json` |
| `patient` | `patient_id`, `sex_category`, `race_category` | sex; **race — Table 1 reporting only, not a model covariate** |
| `adt` | `hospitalization_id`, `hospital_id`, `hospital_type`, `in_dttm`, `out_dttm` | `hospital_id_admission` / `_discharge` and encounter stitching; **`hospital_id` is required by `stitch_encounters`**. `location_category` is **not** read — there is no ICU flag in the pipeline, and "in the ICU" is proxied by `imv_status` wherever an at-risk denominator needs it |

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

The suite is the preflight in `run_pipeline.sh` / `.ps1`, which runs these in
order and stops on the first failure. `pytest tests/` runs them all at once.

```bash
.venv/bin/python code/check_config.py        # the site config is usable at all
.venv/bin/python tests/test_covariates.py    # the covariate protocol is internally consistent
.venv/bin/python tests/test_fio2.py          # FiO2 unit handling (percent vs fraction)
.venv/bin/python tests/test_outliers.py      # outlier bounds, both layers
.venv/bin/python tests/test_waterfall_cache.py  # the cache key covers every input
.venv/bin/python tests/test_doses.py         # dose unit conversion
.venv/bin/python tests/test_pooling.py       # federated-pooling exports are exactly poolable
.venv/bin/python tests/test_paths.py         # output locations and staleness guards
.venv/bin/python tests/test_build_cohort.py  # Phase 0 logic
```

Counts are deliberately not stated here. The numbers that used to be — 16, 14,
30, 11, 12 — were wrong by 2026-09-23 for five of the six files listed (the true
figures were 24, 15, 57, 14 and 9), and two suites were missing from the list
entirely. A hand-maintained count with nothing enforcing it is a claim that goes
stale silently; `pytest tests/` reports the real number.

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
| `cohort.extended_*` | **declared, not consumed.** The 12h/7d descriptive view is deferred; see `covariates.json` `windows.extended._STATUS` |
| `outcomes.mortality_source` | `discharge_category` (`death_dttm` is unstable) |
| `outcomes.mortality_categories` | values counted as death; sensitivity set drops `Hospice` |
| `outcomes.unresolved_discharge_categories` | censored, **not** counted alive |
| `model.scaling` | `gbmt` normalization — **0**, see design notes §3 |
| `model.nstart` | EM random restarts — ≥50, see design notes §7 |

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
| `output/intermediate_phi/` | analytic tables (one row per block-window), model objects, class assignments | **No — PHI** |
| `output/final_no_phi/` | aggregate tables and figures — the coordinating-center upload set, PHI-checked at assembly | **Yes** |

**Nothing under `output/` is tracked** — every site generates its own, PHI
artifacts never leave the site, and the PHI-free set reaches the coordinating
centre by upload rather than by git. `logs/` is gitignored too. The PHI directory
gets a `README.md` warning label written at runtime by `site_dirs()`.

`output/final_no_phi/` is subdivided by the **script** that produced each file.
The `phaseN_` file prefixes were removed on 2026-09-23 — the folder already says
which script produced a file, and the phases themselves are gone. `manifest.json`
sits at the root because it is the pipeline's staleness marker rather than one
script's result. `output/intermediate_phi/` stays flat; it is the machine handoff
between scripts.

What each script writes:

| Script | Files |
|---|---|
| `01_build_cohort.py` | `01_cohort/` — `strobe.{csv,txt,png}`, `provenance.json`, `diagnostics/{missingness,missingness_patterns,diagnostics}.csv`, `exemplar_selection.csv`, plus `manifest.json` at the root of `final_no_phi/` |
| `02_descriptive_cohort.R` | `02_descriptive/` — `baseline_characteristics.csv`, `retention.csv`, `dose_summary.csv`, `balanced_panels.csv`, `imv_encounter_duration.csv`, `pooling_{continuous,categorical}.csv`, `provenance.json`, `captions.md`; figures `fentanyl_{curves,balanced_panels}.png`, `sedative_curves.png` |
| `03_delivery_states.R` | `03_states/` — `state_{prevalence,prevalence_at_risk,transitions}.csv`, `dose_state_{prevalence,transitions}.csv`, `state_{prevalence,prevalence_at_risk,alluvial,raster}.png`, `dose_state_{prevalence,alluvial}.png`, `provenance.json`, `captions.md` |
| `05_exemplar.R` | `05_exemplar/` — `exemplar.png`, `provenance.json`, `captions.md`; PHI handoff `output/intermediate_phi/exemplar_{series.parquet,meta.json,id.txt}` (the chosen episode id stays PHI-side and never reaches the shareable tree) |
| `04_landmark_cohort.R` | `04_landmark/` — `landmark_flow.{csv,txt}`, `T_sensitivity.csv`, `failed_extubation.csv`, `dose_curve.{csv,png}`, `dependence.csv`, `pooling_{continuous,categorical}.csv`, `provenance.json`, `captions.md`; PHI handoff `output/intermediate_phi/landmark_cohort.parquet` |

The two `pooling_*.csv` files exist for **federated pooling**: they carry
`n`, `mean`, `sd`, `sum` and `sum_sq` per variable per stratum (and per window),
so a coordinating centre can compute an exact pooled mean and SD without a
median, which cannot be pooled. Cells below `reporting.small_cell_min_den` are
suppressed.

`manifest.json` is written **last** and is what marks the Phase 0 outputs
complete and current. Every R script calls `require_manifest()` first, which
refuses to read tables produced by a different code version or a drifted config.
Note that it digests whole config files, so editing even a comment in
`config/covariates.json` invalidates the build.

Shareable outputs carry a provenance block from `provenance()`: `site_name`,
`clif_version`, `dataset_version`, `code_version` (git describe), `generated`.

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
| 1 | R | `code/02_descriptive_cohort.R` | Cohort description: Table 1, retention (the at-risk denominator), dose summaries, balanced panels, pooling exports, dose-curve figures |
| 2 | R | `code/03_delivery_states.R` | The seven delivery states: prevalence per window (on both the whole-cohort and the still-ventilated denominator), transitions, the alluvial, and a 100-episode per-patient raster. Every stacked area carries a numbers-at-risk table for window 0 and the windows **ending** at 24 / 48 / 72h. Figures are journal-style — no title or subtitle on the panel; captions are written to `captions.md`. Run twice — delivery **route** and dose **intensity band** |
| 3 | R | `code/04_landmark_cohort.R` | Apply landmark `T`, report retention and failed-extubation counts |
| 4 | R | `code/05_exemplar.R` | **F1** — one ventilation course at sub-hourly resolution: infusion rate, boluses, RASS and NVPS. The episode is **drawn at random from those meeting the pre-specified rule** in `covariates.json` → `exemplar`, never chosen by inspection |

These five steps are the current analysis and are exactly what the runners
execute. The **tabled** trajectory-modelling scripts moved to `code/tabled/` on
2026-09-23 (`04_gbmt_classes.R`, `05_lcmm_classes.R`, `06_transition_model.R`,
`07_outcomes.R`) — see **Objective**. They still run, but by hand:
`Rscript code/tabled/04_gbmt_classes.R`. They were also removed from both runners,
because `07_outcomes.R` is an unconditional `stop()` and while it sat in the step
list `./run_pipeline.sh` exited non-zero on every run — the documented
one-command build could never succeed.

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
│   ├── 01_build_cohort.py        # build the analytic tables (Python / clifpy)
│   ├── 02_descriptive_cohort.R   # cohort description + dose curves
│   ├── 03_delivery_states.R      # delivery states: prevalence, transitions, figures
│   ├── 04_landmark_cohort.R      # landmark T (available, not the current analysis)
│   ├── 05_exemplar.R             # F1: one ventilation course in detail
│   ├── tabled/                   # trajectory modelling, parked -- not in the runners
│   │   ├── 04_gbmt_classes.R
│   │   ├── 05_lcmm_classes.R
│   │   ├── 06_transition_model.R
│   │   └── 07_outcomes.R
│   └── utils/
│       ├── doses.py              # dose unit conversion (not clifpy's)
│       ├── fio2.py               # FiO2 scale detection + normalisation
│       ├── outliers.py           # applies config/outlier_config.json
│       ├── paths.R               # output dirs + provenance
│       ├── paths.py              #   (the two must agree)
│       ├── states.R              # the delivery-state definitions, once
│       ├── figures.R             # shared theme, palettes, panel stacking, captions
│       ├── exemplar.R            # F1 shape, gap handling, the synthetic patient
│       ├── pooling.R             # federated-pooling exports
│       └── dependencies.R        # package list for renv's scanner
├── tests/
│   ├── test_covariates.py        # integrity checks on config/covariates.json
│   ├── test_fio2.py              # FiO2 must be a fraction; enforced, not assumed
│   ├── test_outliers.py          # bounds are applied, and gaps are reported
│   ├── test_build_cohort.py      # Phase 0 logic on synthetic frames
│   ├── test_paths.py             # the PHI boundary, and paths.R == paths.py
│   ├── test_waterfall_cache.py   # the cache key covers every input
│   └── test_doses.py             # every charted dose unit converts correctly
├── docs/
│   ├── principles.md            # the question, the frame, principles, rejected options
│   └── references.md            # literature, with DOIs
├── output/
│   ├── intermediate_phi/         # ALL patient-level artifacts
│   │                             #   .parquet = pipeline input, .csv = review copy
│   └── final_no_phi/             # PHI-free: STROBE flow (.csv/.txt/.png), provenance
│       └── diagnostics/          #   missingness, patterns, counts
└── logs/                         # gitignored
                                  # NOTE: nothing under output/ is tracked;
                                  # every site generates its own
```

## Definitions and provenance

- **Dose**: within-window total fentanyl (infusion + bolus), expressed
  **mcg/hr** — the unit fentanyl is ordered in, and the unit 99.5% of UCMC rows
  are already charted in, so weight does not enter the dominant path. Infusion
  uses LOCF then a time-weighted mean; boluses are summed over the window and
  divided by its hours, never carried forward. (Was mcg/kg/hr until 2026-09-07.)
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
  absolute dose level, and `scaling = 2` additionally divides by a within-patient
  SD that is exactly zero for a near-all-zero indicator, without raising. This
  pipeline uses `scaling = 0`; the mechanism is in design notes §3.

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
