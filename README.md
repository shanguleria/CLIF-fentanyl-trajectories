# CLIF-fentanyl-trajectories

How fentanyl is delivered to adults on invasive mechanical ventilation —
continuous infusion, intermittent bolus, both, or neither — and how that
changes over the ventilation course. Computed from a site's own CLIF 2.1
tables, inside that site.

Everything runs on one machine in your institution. Nothing reaches the network.
The only thing that leaves is a folder of PHI-free aggregate CSVs and figures
that you look at first and send deliberately.

---

## Run it at your site

Eight steps. Expect step 6 to take 20–40 minutes and to stop at least once on
something real about your data — see [Troubleshooting](#troubleshooting).

### 1. Install the tools

- **Python 3.11+** — builds the cohort (step 6, Phase 0)
- **R 4.4+** — every analysis script after that
- **git** — clone rather than download a ZIP, or the run cannot record which
  code version produced it

### 2. Get the code and its dependencies

Clone **at the release tag**, not the default branch:

```bash
git clone --branch v0.7.1 https://github.com/shanguleria/CLIF-fentanyl-trajectories.git
cd CLIF-fentanyl-trajectories

python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

Rscript -e 'install.packages("renv"); renv::restore()'
```

Use `pip install -r requirements.txt`, not `pip install clifpy`. The versions are
pinned exactly because clifpy's minor releases have changed CLIF datetime
timezone handling, and that silently changes results.

### 3. Check your CLIF tables

Eleven tables, all CLIF 2.1, as `clif_<table>.parquet` (or `.csv`) in one
directory:

| table | used for |
|---|---|
| `patient`, `hospitalization`, `adt` | cohort, encounter stitching, care location |
| `respiratory_support` | the ventilation anchor and episode boundaries |
| `medication_admin_continuous` | fentanyl infusions, sedatives, vasopressors |
| `medication_admin_intermittent` | fentanyl boluses |
| `vitals`, `labs` | severity covariates, weight for dosing |
| `patient_assessments` | RASS, and your pain instrument — see [What differs between sites](#what-differs-between-sites) |
| `crrt_therapy`, `hospital_diagnosis` | renal replacement, comorbidity index |

Step 5 verifies they are all present before anything is read. Full field-level
requirements are in [Required CLIF tables and fields](#required-clif-tables-and-fields).

### 4. Write your site config

```bash
cp config/config_template.json config/config.json
```

Edit five things, and nothing else:

| key | what to set |
|---|---|
| `site_name` | your site identifier; stamped on every output |
| `clif_version` | the CLIF **spec** version your tables implement (`2.1.0`) |
| `data_directory` | absolute path to the directory holding the tables |
| `filetype` | `parquet` or `csv` |
| `timezone` | your site's local zone — every analysis window is a wall-clock interval |

Leave `tables_in_use` alone. Set `dose_unit_aliases` only if step 6 reports a
`med_dose_unit` string it cannot convert.

**`config/protocol.json` is a different kind of file. Do not edit it.** It holds
everything that defines what is measured, it is tracked in git, and it is
compared byte-for-byte across sites — a change there makes your results
unpoolable with everyone else's. If a run fails and the fix looks like
changing a number in `protocol.json`, that is a conversation with the study
team, not a local edit.

### 5. Check the config before you run anything

```bash
python3 code/check_config.py
```

Runs in a second, opens no table. It confirms the config is complete, that
`data_directory` really holds the tables named in `tables_in_use`, and that every
category declared in `protocol.json` is a real CLIF mCIDE value. Repeat until it
prints `READY.`

### 6. Run the pipeline

```bash
./run_pipeline.sh          # macOS / Linux
.\run_pipeline.ps1         # Windows PowerShell
```

Both runners do the same thing: preflight, then the six analysis steps in order,
stopping on the first error, then a verification pass over the shareable outputs.
If that last pass fails the outputs exist but **must not be sent**.

### 7. Look at the results, then send them

```
output/final_no_phi/     <- send this, and only this
```

Read `output/final_no_phi/01_cohort/strobe.txt` first — it is the cohort flow,
and the step that excluded the most patients is usually where a data problem
shows up. Then skim `01_cohort/diagnostics/missingness.csv`.

Send `output/final_no_phi/` entire. **Do not send** `output/intermediate_phi/`,
which is patient-level and must never leave your institution.

### 8. Update to a newer version

```bash
git fetch --tags && git checkout v0.8.0
pip install -r requirements.txt
Rscript -e 'renv::restore()'
./run_pipeline.sh
```

A protocol change bumps `definition_version` in `protocol.json`, and every site
has to re-run before results can be combined. The study team will say when.

---

## Keep patient data inside your site

Two output directories, and the split is the whole safety model:

| directory | contents | leaves the site? |
|---|---|---|
| `output/intermediate_phi/` | one row per patient per window | **never** |
| `output/final_no_phi/` | aggregate counts, rates, figures | yes, deliberately |

`.gitignore` excludes all of `output/` with a bare rule, so nothing from either
can be committed — verified by test, not by inspection
(`tests/test_paths.py`). The pipeline writes a warning README into the PHI
directory at runtime so it survives `git clean`.

No step in the pipeline contacts the network.

---

## What differs between sites

The first thing that differs, and the one most likely to be wrong, is **which
pain instrument you chart**.

| `assessment_category` | scale | threshold | covariate |
|---|---|---|---|
| `NVPS` | 0–10 | ≥ 4 | `nvps` |
| `cpot_total` | 0–8 | ≥ 3 | `cpot_total` |
| `BPS` | 3–12 | ≥ 6 | `bps` |

**You should not need to edit anything.** All three are declared; each is cut at
its own threshold; a site that charts none is refused rather than quietly
producing an empty pain panel. Phase 0 prints a zero-record line for every
instrument you do not chart, and that line is expected.

Watch the case: the CLIF mCIDE spells `NVPS` and `BPS` upper-case and the CPOT
family lower-case (`cpot_total`). Step 5 validates every declared category
against clifpy's packaged schema, so a misspelling fails in a second rather than
silently matching nothing for the whole study.

Raw scores are never pooled across instruments — the scales are not comparable.
Dichotomising each at its own threshold is what makes them combinable.

---

## Troubleshooting

### It stops with an error

Expect this once or twice. Each guard below exists because the alternative was a
plausible wrong number. Fix at the named line; do not drop the rows.

| message | meaning |
|---|---|
| `rows in unhandled units` | a fentanyl dose unit `_to_mcg_hr` does not convert, and the run **stops**. Add it to `_to_mcg_hr`; if it truly cannot be converted (e.g. a bare `mcg` amount), list it in `config.json` `accepted_unhandled_fentanyl_units` only after confirming the share of rows with your data owner. Accepted rows are counted on the console and excluded |
| `med_category ... matched no rows` | fentanyl is missing from a medication table; the message prints what *is* there |
| `device_category value not in clifpy's schema` | a casing or vocabulary mismatch in `respiratory_support` |
| `adt is missing hospital_id` | encounter stitching cannot run |
| `not one pain instrument matched a row` | see [What differs between sites](#what-differs-between-sites) |
| `an empty cohort is a bug, not a finding` | an earlier filter removed everything; read `strobe.txt` |
| `dose above the derived ceiling` | proof the outlier bounds did not run |

### It finishes, but a number looks wrong

These degrade quietly, so check them by hand:

1. **The Phase 0 bolus report.** A dose in an unmapped `med_dose_unit` is
   nulled, and because absence of a bolus means zero, a nulled row reads
   downstream as *no bolus given*. This is the most likely cause of a wrong
   co-administration rate. Add your spellings to `dose_unit_aliases`.
2. **`imv_episode_gap_hours = 8`** is tuned to one site's charting cadence and
   nothing checks it locally. It sets the ventilation anchor, so it sets every
   window. Print your own IMV inter-record gap quantiles before trusting it.
3. **`care_setting` can collapse.** It is `location_type` falling back to
   `location_category`. If your ADT populates `location_type` on non-ICU rows,
   the ED merges into the ICU with no error. Check the level counts in
   `care_intervals.parquet` before reading the unit figures.
4. **Small sites suppress more.** A caterpillar with two or three units is not a
   variation estimate. Report the suppression count.

---

## Try it without real data

The single-patient figure renders from a synthetic trace, so the figure code can
be checked on a machine with no CLIF tables at all:

```bash
EXEMPLAR_SOURCE=synthetic Rscript code/03_exemplar.R
```

Writes to `output/dev_synthetic/` — never to the shareable tree — and labels the
panel `GENERATED DATA, NOT A PATIENT`.

---

## Tests

```bash
.venv/bin/python -m pytest tests/ -q
```

195 tests. They check config integrity, unit conversion, the PHI boundary, the
pooling arithmetic and the heterogeneity statistics — on synthetic frames and
config files, never on patient data, so they run anywhere. Both runners execute
them as the preflight.

---

# Reference

Everything below is detail. You do not need it to run the pipeline.

---

## Objective

Describe **how fentanyl is delivered** to adults over the first days of invasive
mechanical ventilation (IMV) — continuous infusion, intermittent bolus, both, or
neither — and how that delivery pattern changes over the ventilation course.

The pipeline is **descriptive**: cohort-level dosing curves with balanced-panel
overlays, delivery-state prevalence and transitions, and a landmarked view of the
first 72 hours.

**Definitions live in `config/`.** `protocol.json` is the protocol — the
cohort rules, the exposure, the covariates, the medication vocabulary, the
outcome definitions and the reporting parameters — and `outlier_config.json` the
bounds; both are tracked, both
carry their rationale in `_why` keys, and both are the source of truth for any
question about what a number means. The project's guiding principles, the
literature and the design notes are kept outside the repository.

> **Scope change, 2026-09-23.** This project began as a trajectory-modelling
> study: find latent classes of fentanyl exposure with `gbmt`/`lcmm`, then use
> class membership as an exposure. That approach was **tested and retired** —
> dose level proved continuous with no subpopulation structure (`gbmt` classes
> predicted an episode's mean dose at R² = 0.82; `lcmm` converged only at two
> classes; the two methods agreed at an adjusted Rand index of 0.037). Steps 3–6
> below remain in the repo and still run, but they are **tabled**, not part of the
> current analysis, and their outputs on disk predate their own scripts. Do not
> cite a number from them without re-running. See `.claude/claude-todo.md`.

## Pipeline steps

| Step | Language | Script | Description |
|---|---|---|---|
| 0 | Python | `code/01_build_cohort.py` | Load CLIF tables via clifpy, build the windowed trajectory table and the time-to-event table |
| 1 | R | `code/02_descriptive_cohort.R` | Cohort description: Table 1 — **stratified on predominant fentanyl intensity**, an ordered declared rule so the p-value is a trend test — retention (the at-risk denominator), dose summaries, balanced panels, pooling exports, dose-curve figures |
| 2 | R | `code/03_exemplar.R` | **F1** — one ventilation course at sub-hourly resolution: infusion rate, boluses, RASS and the site's pain score. The episode is **drawn at random from those meeting the pre-specified rule** in `protocol.json` → `exemplar`, never chosen by inspection |
| 3 | R | `code/04_delivery_states.R` | The seven delivery states: prevalence per window (on both the whole-cohort and the still-ventilated denominator), transitions, the alluvial, and a 100-episode per-patient raster. Every stacked area carries a numbers-at-risk table for window 0 and the windows **ending** at 24 / 48 / 72h. Figures are journal-style — no title or subtitle on the panel; captions are written to `captions.md`. Run twice — delivery **route** and dose **intensity band** |
| 4 | R | `code/05_titration.R` | **Titration practice, both halves.** *Response*: how often a rate increase is accompanied by a bolus within ±30 min, on raw charted timestamps rather than the hourly grid, with rate decreases as a negative control. *Indication*: how often an increase or a bolus had a qualifying pain score or agitation within ±1 h, swept to 4 h — a level at each instrument's own threshold, never a change |
| 5 | R | `code/06_unit_variation.R` | **Variation in both step-4 metrics** — co-administration and documented indication, each split by the ICU unit the patient was in when the event was charted and by calendar year. Caterpillar and funnel views per unit key. Caterpillar intervals are a cluster bootstrap resampling episodes rather than events; funnel limits are unadjusted binomial, with over-dispersion measured into `overdispersion.csv` and applied nowhere. Describes variation; models no association |

These six steps are the current analysis and are exactly what the runners
execute. The **tabled** scripts moved to `code/tabled/` on 2026-09-23
(`gbmt_classes.R`, `lcmm_classes.R`, `transition_model.R`, `outcomes.R`), joined
by `landmark_cohort.R` on 2026-09-30 — see **Objective**. The landmark script's
only consumer of `landmark_cohort.parquet` was the already-tabled gbmt/lcmm
pair, so a live step existed largely to feed tabled work; `06_unit_variation.R`
took the slot it vacated and nothing else renumbered. `time_to_event.parquet` is
unaffected — it is a Phase 0 product that steps 1 and 3 still read, so the STROBE
landmark row survives. They still run, but by hand:
`Rscript code/tabled/gbmt_classes.R`. Their names and output folders carry **no
step number**: the number line belongs to the live pipeline, and while these
scripts kept their old numbers two different scripts claimed 04 and two claimed
05. They were also removed from both runners, because `outcomes.R` is an
unconditional `stop()` and while it sat in the step
list `./run_pipeline.sh` exited non-zero on every run — the documented
one-command build could never succeed.

## Required CLIF tables and fields

Built against the **CLIF 2.1.0** specification. Field names and permissible
values below were verified against `clifpy`'s bundled schemas.

### Shared by all phases

| Table | Fields | Notes |
|---|---|---|
| `patient` | `patient_id`, `sex_category`, `race_category`, `ethnicity_category`, `birth_date` | `death_dttm` exists but is **unstable across CLIF sites** and is deliberately not used; mortality comes from `hospitalization.discharge_category` |
| `hospitalization` | `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission`, `admission_type_category`, `discharge_category` | mortality from `discharge_category` ∈ {`Expired`, `Hospice`} |
| `adt` | `hospitalization_id`, `hospital_id`, `hospital_type`, `in_dttm`, `out_dttm`, `location_category`, `location_name`, `location_type` | encounter stitching, `hospital_id_admission` / `_discharge`, and the ICU stays `06_unit_variation.R` attributes titration events to. `location_category` selects ICU rows; `location_type` is the mCIDE unit taxonomy and is the only unit key another site can pool on; `location_name` is optional site-local free text |

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
| `patient_procedures` | — | **Not read, and not required.** Named in the *tabled* tracheostomy rule as a supporting identification route. It was listed in `tables_in_use` until 2026-10-26 with no consumer anywhere in `code/`, which failed the preflight at a site that lacks the table for no analytic reason. Add it back when that rule is implemented |

### Covariates

Every covariate's window-aggregation rule and missingness class is declared in
**[`config/protocol.json`](config/protocol.json)** and mirrored in design
notes §11. That file is the source of truth; `tests/test_protocol.py` enforces
its internal consistency.

| Table | Fields | Used for |
|---|---|---|
| `labs` | `hospitalization_id`, `lab_result_dttm`, `lab_category`, `lab_value_numeric` | `bun`, `bicarbonate`, `pco2_arterial`, `lactate`, `inr`, `bilirubin_total`, `po2_arterial`, `creatinine`, `platelet_count` |
| `vitals` | `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value` | `spo2`, `map` (SOFA); `weight_kg`, `height_cm` (BMI, dose denominator) |
| `medication_admin_continuous` | `hospitalization_id`, `admin_dttm`, `med_category`, `med_dose`, `med_dose_unit`, `mar_action_category` | NEE: norepinephrine, epinephrine, phenylephrine, dopamine, vasopressin, angiotensin |
| `respiratory_support` | **all columns** | FiO₂ waterfall + `device_category` for `imv_status`. Load every column — the waterfall needs `device_name`, `mode_category`, `lpm_set`, `peep_set` to build its blocks |
| `crrt_therapy` | `hospitalization_id`, `recorded_dttm`, `crrt_mode_category` | `crrt_status`. **Point-in-time table, no start/stop** — interval reconstruction, not a lookup |
| `hospital_diagnosis` | `hospitalization_id`, `diagnosis_code`, `diagnosis_code_format`, `poa_present` | Charlson Comorbidity Index |
| `patient_assessments` | `hospitalization_id`, `recorded_dttm`, `assessment_category`, `numerical_value` | `gcs_total` (SOFA CNS); `RASS` → `rass` (**min**, deepest sedation in the window); and **whichever pain instrument the site charts** → `nvps` / `cpot_total` / `bps` (**max**, worst pain). All are ordinal, so none is averaged. **Also the indication instruments** in `protocol.json` → `indication`: `NVPS` (≥4), `cpot_total` (≥3), `BPS` (≥6) for pain and `RASS` (≥+1) for agitation, used to ask whether a dose change had a documented trigger — see [Pain instruments across sites](#pain-instruments-across-sites) |
| `patient` | `patient_id`, `sex_category`, `race_category` | sex; **race — Table 1 reporting only, not a model covariate** |
| `adt` | `hospitalization_id`, `hospital_id`, `hospital_type`, `in_dttm`, `out_dttm`, `location_category`, `location_name`, `location_type` | `hospital_id_admission` / `_discharge` and encounter stitching; **`hospital_id` is required by `stitch_encounters`**. The `location_*` trio was added 2026-09-30 and feeds `icu_intervals.parquet`, which credits each titration event to the unit the patient was in when it was charted. **It is a STRATIFIER, never a cohort filter**: "in the ICU" is still proxied by `imv_status` wherever an at-risk denominator needs it, and no episode is excluded for its location |

**Encounter blocks are the unit of analysis.** Hospitalizations are stitched with
clifpy `stitch_encounters(..., time_interval=6)`: a `hospitalization_id` is one
encounter, not one clinical course, and a patient intubated, transferred, and
still intubated has one ventilation episode. Without stitching that trajectory is
truncated at the transfer. See design notes §11.

**Two clifpy cautions**, both verified against 0.3.8 in this repo's `.venv` and
documented in `protocol.json`:

- `compute_sofa_polars` must be called with `fill_na_scores_with_zero=False`, and
  its cardiovascular component ignores vasopressin, phenylephrine, angiotensin
  and milrinone. SOFA is the convenience summary here; the explicit markers are
  the severity measure.
- `convert_dose_units_by_med_category` does **not** null a dose it cannot
  convert — it returns the raw value and reports the failure in the unit string.
  Guard the unit string and raise. See design notes §11.

### Pain instruments across sites

**UCMC charts NVPS; most other CLIF sites chart CPOT.** The pipeline does not
assume either. Three behavioural pain instruments are declared, each cut at its
own conventional threshold, and a site populates whichever it uses:

| `assessment_category` | scale | threshold | covariate |
|---|---|---|---|
| `NVPS` | 0–10 | ≥ 4 | `nvps` |
| `cpot_total` | 0–8 | ≥ 3 | `cpot_total` |
| `BPS` | 3–12 | ≥ 6 (the convention is "> 5") | `bps` |

Three things make this work, and all three were broken before 2026-10-26:

1. **The names are CLIF mCIDE values, matched literally.** The mCIDE does not
   capitalise uniformly: `NVPS`, `RASS`, `PAINAD` and `BPS` are upper case while
   the CPOT family is lower case (`cpot_total`, plus four `cpot_*` sub-scores).
   `protocol.json` declared `"CPOT"` — not a permissible value — for four weeks,
   which at a CPOT site would have matched zero rows, printed a reassuring *"not
   charted here"*, passed the ≥1-instrument guard **because RASS was present**,
   and reported the pain half of the indication metric as a silent zero.
   `check_config.py` now validates every declared category against clifpy's
   packaged schema in the preflight, so this cannot survive a run.
2. **The covariates are `site_optional`.** `assessment_covariates()` raises on a
   declared category matching no rows — right for RASS or lactate, which means a
   broken extract, and fatal for one of three interchangeable instruments. The
   pain covariates carry `site_optional: true`: an absent one prints a loud
   zero-record line and becomes an **all-NA column**, which must still exist
   because `assert_config_is_honoured()` requires every declared key to reach a
   column. **At least one must be present or the build refuses** — a site with no
   analgesia signal at all must not produce a silent all-NA pain panel.
3. **The figure takes its scale from config.** `03_exemplar.R` reads the
   instrument's range from `outlier_config.json` → `analysis_unit_bounds`, so a
   CPOT trace is drawn on a 0–8 axis. It was hardcoded to NVPS's 0–10, which
   would have labelled ticks at 9 and 10 that no CPOT value can reach — a wrong
   figure with no error. The axis title comes from
   `indication.instrument_display_names` ("CPOT", not `cpot_total`).

`PAINAD` is charted at UCMC on 45.2% of episodes and is **deliberately excluded**:
it scores pain in dementia, and mixing it into the pain class mixes two
constructs. Raw scores are never pooled or averaged across instruments — the
scales are not comparable. Dichotomising each at its own threshold is exactly
what makes them combinable: after the cut, every instrument contributes the same
yes/no.

**Adding a fourth instrument** is a config change in three places and no code:
`indication.pain_instruments` (threshold), `indication.instrument_display_names`
(panel label), a `time_varying` entry with `site_optional: true` plus membership
in `missing_values.time_varying.members`, and an `analysis_unit_bounds` entry.
`tests/test_protocol.py` fails if you miss any of them.

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
`data_directory` in this file; there is no second code path.

**`config/protocol.json` is a different kind of file and is committed.** It
holds the whole protocol: every variable's source, window-aggregation rule, LOCF
cap and missingness class, **and** the cohort rules, outcome definitions,
medication vocabulary, model and reporting parameters. It carries its own
`definition_version`, which is stamped onto shareable outputs. The test for what
belongs there: *if two sites set this differently, is the pooled result still
meaningful?* If no, it is protocol.

Those last five blocks — `cohort`, `outcomes`, `medications`, `model`,
`reporting` — lived in `config.json` until 2026-10-26. Because that file is
gitignored and is deliberately excluded from `PROTOCOL_FILES`, the parameters
that *define the cohort* were invisible to the only mechanism that proves two
sites ran the same protocol: a site could change `imv_episode_gap_hours`, build a
different cohort, and `check_provenance.py` would still report agreement. The
file said so itself and then held the estimand. Now `config.json` is site
identity only, and **nothing estimand-defining can live there** because there is
nothing estimand-defining left in it.

The suite is the preflight in `run_pipeline.sh` / `.ps1`, which runs these in
order and stops on the first failure. `pytest tests/` runs them all at once.

```bash
.venv/bin/python code/check_config.py        # the site config is usable at all
.venv/bin/python tests/test_protocol.py    # the covariate protocol is internally consistent
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

**`config/config.json` — site identity. Edit all of this.**

| Key | Meaning |
|---|---|
| `site_name` | site identifier stamped on every shareable output |
| `clif_version` | the CLIF **spec** version your tables implement |
| `data_directory` | absolute path to the directory holding `clif_<table>.parquet` |
| `filetype` | `parquet` or `csv`; passed straight to clifpy |
| `timezone` | your site's local zone. Every window is a wall-clock interval, so this shifts results if wrong |
| `tables_in_use` | CLIF tables loaded by Phase 0, and the set `check_config.py` verifies by filename |
| `accepted_unhandled_fentanyl_units` | fentanyl dose units you have chosen to leave unconverted (default `[]`, so any unconvertible unit stops the run). Counted and excluded, not guessed. UCMC lists `mcg` (90 rows site-wide, pending confirmation). Digested into `site_config_digest`, so the centre can see it |
| `dose_unit_aliases` | **your** free-text `med_dose_unit` strings → the canonical bolus units. `med_dose_unit` has no mCIDE, so this is the one estimand-affecting thing that is irreducibly site-specific. UCMC maps `"mcg of opiate"` → `mcg`. An unmapped unit carrying a dose is reported and nulled, and a nulled bolus reads downstream as *no bolus given* — **so read the Phase 0 bolus report before trusting the co-administration numbers at a new site** |

**`config/protocol.json` — the estimand. Do not edit any of this to make a run
work.** If two sites set these differently the pooled result is not meaningful,
and `check_provenance.py` exits non-zero on divergence.

| Key | Meaning |
|---|---|
| `medications.*` | `med_category` values. These **are** mCIDE-controlled (75 permissible values for continuous, 165 for intermittent, and the two lists differ), so they are the same at every site — what varies is whether your site populates them. Phase 0 asserts they are present and prints what it found; `check_config.py` validates each against the schema for the table that role is read from |
| `cohort.window_hours` | aggregation window (default 4) |
| `cohort.landmark_hours` | landmark `T` (default 72) |
| `cohort.imv_episode_gap_hours` | gap that separates two IMV episodes (default 8) |
| `cohort.min_imv_hours`, `cohort.min_age` | duration and age floors |
| `cohort.balanced_panel_hours` | thresholds for the Phase 1 overlays |
| `cohort.extended_*` | **declared, not consumed.** The 12h/7d descriptive view is deferred; see `windows.extended._STATUS` |
| `outcomes.mortality_categories` | `discharge_category` values counted as death (`death_dttm` is unstable); sensitivity set drops `Hospice` |
| `outcomes.unresolved_discharge_categories` | censored, **not** counted alive |
| `outcomes.successful_extubation_hours` | gap after which an extubation counts as successful |
| `reporting.small_cell_min_den` | the disclosure floor; `unit_variation.min_events_per_unit` must clear it |
| `indication.*` | pain and sedation instruments, thresholds and the pairing window — see [Pain instruments across sites](#pain-instruments-across-sites) |
| `model.scaling` | `gbmt` normalization — **0**, see design notes §3 |
| `model.seed`, `model.ng_range`, `model.polynomial_degree` | trajectory-model settings; used only by the **tabled** scripts in `code/tabled/` |

**`clif_version` is always the CLIF *specification* version** — `2.1.0` here.
A conversion release is **not** a CLIF version: the MIMIC-to-CLIF build 1.1.0
implements the CLIF 2.x spec, and there is no CLIF spec version 1.1.0. If you are
running against a conversion, put its release in `site_name`
(`MIMIC-IV-Ext-CLIF-1.1.0`), which is free text and is stamped on every output.
A separate `dataset_version` key was removed on 2026-10-26: it applied only to
MIMIC, was empty at every real site, and an always-empty key in a file every site
must edit is noise.

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
| `01_build_cohort.py` | `01_cohort/` — `strobe.{csv,txt,png}`, `provenance.json`, `diagnostics/{missingness,missingness_patterns,diagnostics}.csv`, `exemplar_selection.csv`, plus `manifest.json`; PHI handoffs include `titration_{rate,bolus}_events.parquet` (raw charted timestamps — the titration analysis never uses the hourly grid) and `care_intervals.parquet` (attributed care-location stays — ICU and ED — on the relative clock, so the attribution rule can change without a Phase 0 rebuild) at the root of `final_no_phi/` |
| `02_descriptive_cohort.R` | `02_descriptive/` — `baseline_characteristics.csv`, `retention.csv`, `dose_summary.csv`, `balanced_panels.csv`, `imv_encounter_duration.csv`, `pooling_{continuous,categorical}.csv`, `provenance.json`, `captions.md`; figures `fentanyl_{curves,balanced_panels}.png`, `sedative_curves.png` |
| `03_exemplar.R` | `03_exemplar/` — `exemplar.png`, `provenance.json`, `captions.md`; PHI handoff `output/intermediate_phi/exemplar_{series.parquet,meta.json,id.txt}` (the chosen episode id stays PHI-side and never reaches the shareable tree) |
| `04_delivery_states.R` | `04_states/` — `state_{prevalence,prevalence_at_risk,transitions}.csv`, `dose_state_{prevalence,transitions}.csv`, `state_{prevalence,prevalence_at_risk,alluvial,raster}.png`, `dose_state_{prevalence,alluvial}.png`, `provenance.json`, `captions.md` |
| `05_titration.R` | `05_titration/` — `coadministration.csv`, `window_sensitivity.csv`, `charting_precision.csv`, `rate_change_magnitude.csv`, `charting_agreement.csv`, `adherence_distribution.csv`, `provenance.json`, `captions.md`; figures `coadministration.png`, `window_sensitivity.png`; the indication analysis `indication.csv`, `indication_window_sensitivity.csv`, `indication.png`; PHI handoffs `output/intermediate_phi/titration_adherence.parquet` (per-encounter adherence; only its distribution ships), `titration_events_classified.parquet` and `indication_events.parquet` (the classified events and per-event indication flags `06_unit_variation.R` splits, so each definition has one home) |
| `06_unit_variation.R` | `06_unit_variation/` — `unit_adherence.csv`, `year_adherence.csv`, `attribution_funnel.csv`, `overdispersion.csv`, `provenance.json`, `captions.md`; ten figures, named by metric: `coadmin_*` (was a bolus given with the increase?) and `indication_*` (was there a documented pain score or agitation?), each as `_caterpillar_<key>.png`, `_funnel_<key>.png` and `_year.png` — one pair of caterpillar/funnel per `unit_variation.unit_keys` (`care_setting`, `location_name`). Writes no PHI handoff |
| `tabled/landmark_cohort.R` | **Tabled 2026-09-30, not run by the runners.** `landmark/` — `landmark_flow.{csv,txt}`, `T_sensitivity.csv`, `failed_extubation.csv`, `dose_curve.{csv,png}`, `dependence.csv`, `pooling_{continuous,categorical}.csv`, `provenance.json`, `captions.md`; PHI handoff `output/intermediate_phi/landmark_cohort.parquet` |

The two `pooling_*.csv` files exist for **federated pooling**: they carry
`n`, `mean`, `sd`, `sum` and `sum_sq` per variable per stratum (and per window),
so a coordinating centre can compute an exact pooled mean and SD without a
median, which cannot be pooled. Cells below `reporting.small_cell_min_den` are
suppressed. Only `02_descriptive_cohort.R` emits them now — the landmark script
that was the second producer was tabled on 2026-09-30, and
`tests/test_pooling.py` was narrowed to match rather than left to skip
vacuously.

`manifest.json` is written **last** and is what marks the Phase 0 outputs
complete and current. Every R script calls `require_manifest()` first, which
refuses to read tables produced by a different code version or a drifted config.
Note that it digests whole config files, so editing even a comment in
`config/protocol.json` invalidates the build.

Shareable outputs carry a provenance block from `provenance()` — see
[Definitions and provenance](#definitions-and-provenance) for the fields and
how a coordinating centre confirms that several sites ran the same analysis.

## Project structure

```
CLIF-fentanyl-trajectories/
├── README.md
├── requirements.txt              # Python, pinned
├── renv.lock                     # R, pinned
├── run_pipeline.sh / .ps1
├── config/
│   ├── protocol.json             # committed -- THE PROTOCOL (definition_version):
│   │                             #   cohort, outcomes, medications, model, reporting,
│   │                             #   covariates, exposure, indication, unit_variation
│   ├── outlier_config.json       # committed -- the only place bounds are written
│   ├── config_template.json      # committed
│   └── config.json               # gitignored -- SITE IDENTITY ONLY: paths, timezone,
│                                 #   tables_in_use, dose_unit_aliases
├── code/
│   ├── check_config.py           # preflight: is config.json usable?
│   ├── check_provenance.py       # coordinating centre: did all sites run the same code?
│   ├── 01_build_cohort.py        # build the analytic tables (Python / clifpy)
│   ├── 02_descriptive_cohort.R   # cohort description + dose curves
│   ├── 03_exemplar.R             # F1: one ventilation course in detail
│   ├── 04_delivery_states.R      # delivery states: prevalence, transitions, figures
│   ├── 05_titration.R            # bolus co-administration at uptitration
│   ├── 06_unit_variation.R       # variation by ICU unit and calendar year
│   ├── tabled/                   # parked -- not in the runners
│   │   ├── gbmt_classes.R        # -> output/parked/gbmt/
│   │   ├── lcmm_classes.R        # -> parked/lcmm/
│   │   ├── transition_model.R    # -> parked/transitions/
│   │   ├── landmark_cohort.R     # -> parked/landmark/   (tabled 2026-09-30)
│   │   └── outcomes.R            #    unnumbered: the number line is the live pipeline's
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
│   ├── test_protocol.py        # integrity checks on config/protocol.json
│   ├── test_fio2.py              # FiO2 must be a fraction; enforced, not assumed
│   ├── test_outliers.py          # bounds are applied, and gaps are reported
│   ├── test_build_cohort.py      # Phase 0 logic on synthetic frames
│   ├── test_paths.py             # the PHI boundary, and paths.R == paths.py
│   ├── test_waterfall_cache.py   # the cache key covers every input
│   ├── test_pooling.py           # pooling exports are exactly poolable
│   └── test_doses.py             # every charted dose unit converts correctly
├── multisite/                    # gitignored -- OTHER SITES' bundles as received,
│   ├── ucmc/                     #   one dir per site, each an unmodified copy of
│   ├── rush/                     #   that site's output/final_no_phi/. Read by
│   └── northwestern/             #   check_provenance.py and pool_sites.R
├── output/
│   ├── pooled/                   # pool_sites.R output: one CSV per input file
│   │                             #   with a POOLED row + heterogeneity.csv
│   ├── intermediate_phi/         # ALL patient-level artifacts
│   │                             #   .parquet = pipeline input, .csv = review copy
│   ├── final_no_phi/             # THE UPLOAD SET. PHI-free; one dir per script,
│   │   └── diagnostics/          #   each with a provenance.json
│   ├── parked/                   # results of code/tabled/ -- NOT uploaded (2026-10-26)
│   └── dev_synthetic/            # EXEMPLAR_SOURCE=synthetic renders; never real data
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

### The provenance block

Every directory under `output/final_no_phi/` carries a `provenance.json` from
`provenance()` (`code/utils/paths.py`, mirrored in `paths.R`):

| Field | Meaning |
|---|---|
| `site_name` | from `config.json` |
| `clif_version` | the CLIF **spec** version |
| `definition_versions` | one entry per protocol file — they are bumped independently |
| `code_version` | `git describe --always --dirty`; **advisory only** |
| `code_digest` | SHA-256 over the code the runner executes — the authoritative "same code?" answer |
| `protocol_digests` | `protocol.json` + `outlier_config.json`; **must be identical across sites** |
| `site_config_digest` | `config.json`; **expected to differ** — it holds `site_name`, `data_directory` and `dose_unit_aliases`. Recorded so the centre can *see* a site's unit normalisation, never compared for equality |
| `file_digests` | one digest per covered file, so a divergence can be named |
| `generated` | local timestamp |

`code_version` is advisory because `git describe` answers the question only
where a checkout exists, carries a tag and is clean — it returns `"unknown"` for
a ZIP download and a bare `-dirty` for an uncommitted edit. `code_digest` is
content-based, so it works in all three cases and compares by equality. Both
languages must produce the same digest; `tests/test_paths.py` asserts it.

### Confirming that several sites ran the same analysis

At the coordinating centre, over the site bundles:

```bash
python3 code/check_provenance.py multisite/site_a multisite/site_b multisite/site_c
```

It reports one row per site and then separates the two kinds of divergence:

- **Protocol** (`protocol_digests`, `definition_versions`, `clif_version`) —
  exits non-zero. A pooled estimate across different protocol versions is not
  an estimate of one quantity.
- **Code** (`code_digest`) — reported, never fatal, and named down to the file.
  A site may legitimately have had to extend `_to_mcg_hr` with a locally charted
  dose unit to run at all; this pipeline's own error message instructs it to.

## Data safety

Patient-level data never leaves `data/intermediate_phi/` and
`output/intermediate_phi/`, both gitignored and both carrying a warning label
written at runtime. Only `output/final_no_phi/` is intended for transfer; its
contents are aggregate by construction and PHI-checked at assembly.

## Contributing

R scripts follow a house style, which any script in `code/` demonstrates — read
`code/05_titration.R` as the reference. It is: a Purpose/Author/Created/Inputs/
Outputs banner; the `requireNamespace` + `library(p, character.only = TRUE)`
package loop; `here()` for every path and never `setwd()`; and a `sessionInfo()`
provenance write at the end.

(This previously pointed at `R_setup.R` under a home directory, which resolves
only on the original author's machine — a documented path that cannot work at
another site is worse than no pointer, since the reader cannot tell whether they
are missing a file or the reference is stale.)

**Any generated code must be read by a human before it is run or committed.**
