# Fentanyl trajectory modeling — design notes

**Status:** working draft for review · **Date:** 2026-09-04 · **Author:** Shan Guleria
**Scope:** how to apply group-based trajectory modeling to ICU fentanyl dosing
data. Decisions, the evidence behind them, and open questions.

All numeric claims carry a source pointer. Everything cited comes from synthetic
or public data — no patient data was used to produce any figure in this file.

---

## 1. Two models, two questions

| | **Model A — exposure** | **Model B — management strategy** |
|---|---|---|
| Question | What are the trajectories of total opioid exposure? | Is the patient managed by escalating a drip, or by intermittent boluses? |
| Indicators | **1**: total fentanyl mcg/hr per 4h window | **2**: infusion rate + bolus-equivalent rate |
| `scaling` | `0` | `0` (required — see §4) |
| Status | **Primary.** Do this first. | **Secondary.** Do after A is settled. |

Model A is the cleaner and more defensible primary analysis. Model B answers a
genuinely different and interesting question, but it is the more fragile model
and should not be the headline result.

### Model A indicator

Combining infusion and push into one quantity resolves a units mismatch:
infusion rate is an *intensity* (mcg/hr), cumulative bolus is a *quantity*
(mcg/kg). Converting the window's bolus total to an equivalent rate makes them
commensurable and additive:

```
total_fent_mcg_kg_hr =
    ( infusion_mcg_delivered_in_window + sum_of_push_mcg_in_window )
    / weight_kg / 4
```

This is what the patient actually received, it stays in clinically interpretable
units, and it keeps the model at a single indicator.

### Model B indicators

Keep the two streams separate, both expressed as mcg/hr over the window:

- `inf_rate`  — time-weighted mean infusion rate
- `push_rate` — window bolus total / weight / 4

---

## 2. A correction to record

It is **not** GBTM that fails on zero-heavy dosing data. Model A — a single
combined-exposure indicator — is a perfectly appropriate GBTM application, and
zeros are harmless there once `scaling = 0` removes all division.

What fails is specifically a **near-all-zero second indicator under within-unit
normalisation** (Model B run at the package default `scaling = 2`) — see §3 for
the mechanism.

---

## 3. What `scaling` actually does

`gbmt` normalises each indicator **within each unit** — every patient against
their *own* mean and SD across the observation period, not against the cohort.
(source: `references/gbmt_R.pdf`, `gbmt` Details section)

| `scaling` | Formula | Erases | Requires > 0 |
|---|---|---|---|
| **0** none | `y` | nothing | no |
| **1** centering | `y - mean_i` | the patient's own level | no |
| **2** standardization *(package default)* | `(y - mean_i) / sd_i` | level **and** variability | no |
| **3** ratio to mean | `y / mean_i` | level (centers at 1) | **yes** |
| **4** log ratio to mean | `log(y / mean_i)` | level (centers at 0) | **yes** |

**The consequence that matters clinically:** with any `scaling >= 1`, absolute
dose cannot define a group. A patient flat at 200 mcg/hr and one flat at 50
become the *same trajectory* after within-patient normalisation. If
"persistently high-dose" is a phenotype we want to find, `scaling >= 1` destroys
it by construction.

`scaling = 3` and `4` are unusable for infusion data regardless — they require
strictly positive values, and the data has structural zeros whenever the drip is
off.

**The failure mode to know about, if anyone is ever tempted to change the
setting.** A patient who never received a bolus has a within-patient SD of
exactly zero for that indicator, and `scaling = 2` divides by it. `gbmt` raises
no error and issues no warning — it returns a confidently wrong partition. That
is a property of the formula, not a finding about any particular dataset, and it
is why `scaling = 0` is mandatory rather than preferred.

---

## 4. Evidence — removed 2026-09-07 (SG)

This section held four synthetic-data experiments (E1–E4) demonstrating what each
`scaling` value erases, and the scripts that produced them
(validation/scaling_experiments.R, validation/composition_bias_demo.R) have
been deleted along with it. Two further scripts it cited were never in this repo
at all.

The decisions those experiments supported are unchanged and are stated as
mechanisms rather than as measurements: `scaling = 0` and the divide-by-zero-SD
failure in §3, `nstart >= 50` and the polynomial-degree cap in §7, the
balanced-panel rationale in §10 Phase 1. The section number is left in place so
every later cross-reference (§5–§12) keeps its number.

---

## 5. Data preparation spec

| Step | Rule | Why |
|---|---|---|
| Time anchor | intubation (preferred) or ICU admission | trajectories need a common origin |
| Window | **4h** (granular, to 72h) or 12h (extended, to 7d) | hourly over 7 days = 168 points, far more resolution than a degree-2/3 polynomial can use |
| Truncation | fixed window (e.g. 72h or 7d), stated explicitly | polynomial basis needs common support |
| **Infusion rate** | **LOCF**, then **time-weighted mean** within the window | a rate persists until changed; a plain mean of records is wrong whenever rate changes are unevenly spaced |
| **Bolus doses** | **SUM within window. Never LOCF.** | a bolus is an event, not a state — carrying it forward replicates one dose across every later window and massively overcounts |
| Empty windows | true `0` for both streams | drip off is a real zero, not missing |
| Weight | no longer the fentanyl denominator (§5); still fixes `nee` and BMI | a normaliser, not a covariate. Two different weight rules are required — see §11 *Weight* |

### The hourly grid — RESOLVED 2026-09-05 (SG)

The infusion rate is materialised as an **hourly waterfall grid**: every 1h
timepoint carries the rate in effect, so a 4h window total is a sum over its four
cells. Same structure as `nee` (§11), which keeps the two exposures consistent.

| Step | Rule |
|---|---|
| Bin | to hours; take the **last** charted rate per (block, hour), stable-sorted on (time, value) |
| Stops | `mar_action_category == "stop"` is a rate of **0**, not a missing value |
| Fill | forward-fill at most `hold_hours` |
| Remainder | **0** — `absence_means_zero`, per the table above |
| Bounds | never past the block's first/last charted medication record, or the `alive_admitted` span |
| Extubated windows | `imv_status == 0` forces dose to 0, applied **after** the grid (§10a(a)) |

**This resolves the open "are stops charted explicitly?" question by making the
answer not matter for correctness** — but only because of the hold cap. Without
it, LOCF on an infusion rate runs to the end of follow-up: a drip that stopped at
hour 10 with no stop row charted would still show a rate at hour 72. Because this
is the **exposure**, that is not noise around the result, it *is* the result.

> ⚠ `hold_hours = 4` is **inherited from `nee`**, where the CRRT coordinating
> site's median inter-record interval was 57 min (p90 68 min). **The fentanyl
> charting interval at UCMC is the single parameter most worth measuring before
> the first real run.** If fentanyl is charted less densely than vasopressors, 4h
> under-fills and dose is understated; if stops are reliably charted, the cap
> rarely binds at all.

**The grid is a discretised time-weighted mean, not an exact one.** A rate change
at 10:30 is attributed to whichever cell the last-rate-per-hour rule picks, so the
error is bounded by the within-hour timing of rate changes — at most one hour of
one rate per change. At 4h windows that is small and uniform across patients, but
it is an approximation and belongs in the Methods. `grid_resolution_minutes` in
`config/covariates.json` lowers it to 15 or 30 for a sensitivity analysis; nothing
else changes.

### Worked example of the LOCF trap

If a patient receives one 50 mcg bolus in window 3 and the pipeline applies LOCF
to the bolus column, windows 4, 5, 6 … all report 50 mcg. A 50 mcg exposure
becomes several hundred. Both streams come out of the same pipeline, so the rule
must be applied per-column, not globally.

---

## 6. Where the covariates go

`gbmt` has **no covariate argument** — a time-varying variable is either an
indicator that defines the groups, or it is not in the model.

| Role | Variables | Where |
|---|---|---|
| Trajectory indicator | fentanyl dose (± the companion sedatives — see §9) | `x.names` |
| Normaliser | none for fentanyl as of 2026-09-07; weight still normalises `nee` | data prep — mcg/hr (§11) |
| Membership predictor | age, sex, CCI, BMI at admission; plus **window-0** SOFA, NEE, oxygenation | stage 2: `multinom(group ~ ...)` |
| Distal outcome | successful extubation, 30-day mortality, VFD-28 | stage 2 — competing risks, see Phase 6 |
| Group descriptor | bolus counts, % windows with a bolus, % windows on CRRT | descriptive table only — characterises groups without letting them drive the grouping |
| **Reported, not modelled** | race, hospital at admission / at discharge | Table 1 baseline characteristics only (§11) |

The exact variable list is in [`config/covariates.json`](../config/covariates.json)
and §11; this table only says what each is *for*. Two notes on what changed:

- **RASS is not collected.** It appeared here as a candidate indicator, but it is
  in no other section, in no phase, and in the config. Adding an ordinal
  assessment as a second indicator runs into the same mixed-units problem as
  the companion sedatives (§3, §9) — if it is wanted, it goes through that decision,
  not in by default.
- **Severity enters as a membership predictor at window 0, not as "admission
  SOFA".** SOFA, NEE and oxygenation are time-varying here (§11), so the baseline
  value is a window-0 slice of a trajectory column rather than a separate
  admission-time variable. There is no separate `admission_type` covariate; it
  was listed here and in Phase 0 but is collected by neither. **`race` IS now
  collected** (SG, 2026-09-05) — for Table 1 reporting only, not as a membership
  predictor.

The stage-2 approach ("classify–analyze") treats estimated group membership as
observed and therefore understates standard errors. Acceptable when APPA > 0.9;
below that, a proper 3-step / BCH correction is needed, which means `lcmm` or
Mplus rather than `gbmt`.

---

## 7. Class enumeration protocol

Fixed in advance, applied as a conjunction — never minimum BIC alone.

| Criterion | Threshold |
|---|---|
| BIC elbow (not global minimum) | ΔBIC flattens |
| APPA, every group | ≥ 0.70 |
| OCC, every group | > 5.0 |
| Smallest group | ≥ 5% of patients |
| Entropy | > 0.80 |
| Interpretability | distinct, clinically meaningful shapes |
| Stability | same solution across random restarts |

BIC magnitudes are **not comparable across `scaling` settings**, because each
setting changes the values the likelihood is computed on. Comparing them is a
category error, not a close call.

Nor across different indicator sets, different rows, or different time windows.

---

## 8. Landmark design and cohort definition

### The problem

Extubation and death both drive fentanyl dosing to zero. If patients who are
extubated early are carried forward as `dose = 0`, the exposure becomes a
deterministic function of the outcome: a "rapid taper to zero" trajectory group
*is* the early-extubation group. Reporting that this group has superior
liberation outcomes would restate how the variable was constructed, not
discover anything.

The mirror-image problem is immortal time bias. A patient with a complete 72h
de-escalating trajectory necessarily stayed intubated for 72 hours. Comparing
them against patients without such a trajectory partly compares "stayed
ventilated" against "did not." This is the structure Anderson, Cain & Gelber
(1983) identified when comparing chemotherapy responders to non-responders:
responding takes time, so anyone who died early was automatically a
non-responder.

### The design

```
hour 0 ──── intubation (time anchor)
   │
   │   [0, T] ── trajectory measured here; class assigned from THIS window only
   │
hour T ──── LANDMARK
   │         Cohort = alive AND still mechanically ventilated at T
   │         Follow-up clock RESETS to zero here
   │
   └──────── outcomes from T onward (competing risks)
```

Four rules, all load-bearing:

1. The risk set is defined **at T**; anyone with the outcome before T is excluded.
2. Class assignment uses **only** [0, T] data — no information from after T.
3. Survival time is measured **from T**, not from intubation. Forgetting this
   step reintroduces exactly the bias the design removes.
4. Every analysed patient survived to T by construction, so class membership
   cannot act as a proxy for having survived.

### The estimand

> The association between fentanyl trajectory class and subsequent outcome,
> among **encounter blocks** (ventilation episodes, not patients),
> **conditional on being alive and mechanically ventilated at T.**

The conditioning is part of the estimand, not a limitation for the discussion
section. It belongs in the Methods, the abstract, and arguably the title.

**The unit is the encounter block** *(SG, 2026-09-05)*. A patient admitted twice
for different reasons may have entirely different fentanyl trajectories, and each
is a distinct clinical episode worth describing. Say "ventilation episodes", not
"patients", in every count, table and figure legend — including the CONSORT/STROBE
flow. See §11.

### Why not simply zero-fill (evidence)

| Alternative | Why it fails | Source |
|---|---|---|
| Carry `dose = 0` forward after extubation | Outcome becomes embedded in the exposure; also often factually wrong, since many patients receive fentanyl post-extubation for pain | §2, §8 above |
| Keep short trajectories, unbalanced panel | `gbmt` **silently caps** the polynomial degree at (shortest unit's windows − 1), by warning only, so one briefly ventilated patient constrains the trajectory shape available to the entire cohort | `gbmt` source, degree check |


### Justification language (Methods-ready draft)

> Patients were required to be alive and receiving mechanical ventilation at the
> landmark time of T hours after intubation. This requirement follows from the
> definition of the exposure rather than from analytic convenience: a sedation
> *trajectory* is a pattern observed over time, and a patient ventilated for
> fewer than T hours has no opportunity to exhibit one. Assigning such patients
> a trajectory — by carrying dose forward as zero after extubation — would make
> the exposure a deterministic function of the outcome, since extubation itself
> drives sedative dosing to zero. Follow-up for all outcomes therefore begins at
> T, and estimates apply to patients still receiving mechanical ventilation at
> that time: the population in whom sedation strategy remains a live clinical
> decision.

### Reporting requirements

- Flow diagram: N intubated → N alive and ventilated at T → N analysed.
- **A table comparing included versus excluded patients.** Not optional — it is
  the disclosure that makes the rest credible, and the first thing a careful
  reviewer looks for.
- Sensitivity analyses at multiple T (24 / 48 / 72h), reported together.
- Explicit estimand language wherever results are stated.
- Do **not** generalise to all ventilated patients.

### Choosing T

Run this before any modelling — it is a two-line tabulation and it makes the
decision for you:

| T | N alive & ventilated | % of intubated cohort | Windows available (4h) |
|---|---|---|---|
| 12h | ? | ? | 3 |
| 24h | ? | ? | 6 |
| 48h | ? | ? | 12 |
| 72h | ? | ? | 18 |

A degree-2 polynomial needs only 3 points, so 24h (6 windows) is already ample.
The binding constraint is cohort retention, not trajectory resolution — which
argues for the shortest T that still captures clinically meaningful variation.

### The alternative: joint modelling

`lcmm::Jointlcmm` models the trajectory and the time-to-event **simultaneously**,
with latent classes explaining both. Verified against the package: it supports
competing risks with cause-specific baseline hazards, `hazardtype = "Specific"`
for class-specific hazards, and `cause(X)` syntax for cause-specific covariate
effects. Extubation-versus-death is therefore natively handled.

Because informative dropout is modelled rather than assumed away, **no landmark
is required** — patients contribute whatever follow-up they have.

**The catch that changes the paper.** The classes are constructed *using* the
outcome. Class-specific hazards are the model output, not a downstream test, so
"trajectory class predicts extubation" cannot be reported as an independent
finding. It is a different inferential structure — descriptive of joint patterns
rather than exposure → outcome.

**Practical costs:** many more parameters and frequent convergence failure
(requires `gridsearch()` for starting values); class enumeration must balance
longitudinal against survival fit, muddying the clean §7 protocol; and it is
markedly harder to explain to clinical co-authors and reviewers.

**Decision:** landmark design for the primary analysis; joint model as a
sensitivity analysis demonstrating the conclusions survive relaxing the landmark
assumption, or as a separate paper.

---

## 9. Open questions

- [x] **RESOLVED 2026-09-05 (SG).** Confirm infusion stops are charted as
      explicit `0` rows (§5). Superseded by the hourly grid + `hold_hours` cap,
      which makes the answer non-load-bearing for correctness. The **fentanyl
      charting interval** is now the thing to measure instead, since it sets the
      cap.
- [ ] Decide the time anchor: intubation vs ICU admission.
- [ ] Decide the truncation window, and quantify how many patients it excludes —
      informative censoring (death, early extubation) is the largest validity
      threat in this design.
- [ ] Check the zero fraction of the Model A combined indicator once built. Heavy
      zero-inflation may argue for `crimCV`'s zero-inflated Poisson instead.
- [ ] Script the `scaling` × BIC comparison so §7's numbers are reproducible.
- [ ] Decide whether the companion sedatives join as additional indicators — note
      that mixing units would force `scaling >= 2`, which §3 says destroys the
      level information we care about. Fentanyl-equivalents may be the escape.
- [ ] Fix the **landmark time** for outcome ascertainment (§10 Phase 4). Measuring
      trajectories over 0–72h and then counting extubations from hour 0 is
      immortal time bias — a patient cannot be in the cohort unless they stayed
      intubated long enough to have a trajectory.
- [ ] Confirm whether successful extubation is itself a competing-risks problem
      (extubation vs death vs tracheostomy) or can be treated as binary.
- [x] **RESOLVED 2026-09-05.** Write the per-variable within-window aggregation
      rules (§11). Now declared in `config/covariates.json` and mirrored in §11;
      enforced by `tests/test_covariates.py`.
- [x] **RESOLVED 2026-09-05 (SG). How many encounter blocks may one patient
      contribute?** **All of them** — the unit of analysis is the ventilation
      episode, not the patient, because two admissions for different reasons are
      genuinely different clinical courses. Residual within-patient correlation
      remains and cannot be absorbed by `gbmt` or `crr`; §11 makes the two
      bounding counts required outputs. Original framing: `encounter_block` is now the analysis key for Phases 0-6, and
      stitching does not merge a January admission with a June one. A patient
      intubated twice enters as two units that are not independent. Neither
      `gbmt::gbmt` nor `cmprsk::crr` takes a clustering argument — verified by
      reading their formals — so the correlation cannot be absorbed at the model
      stage the way `CRRT-dose-lmtp` absorbs it in `lmtp`. It must be removed at
      the cohort stage or accepted and bounded. See §11.
- [ ] **Specify the adjustment sets.** §6 names the *roles* — age, sex, CCI, BMI
      and window-0 severity as membership predictors — but no section states which
      covariates enter `multinom(group ~ ...)` or `glm(outcome ~ group + ...)`.
      Surfaced by the 2026-09-05 consistency audit, which found the covariate list
      defined and the model formulas not. The two are different decisions and only
      the first has been made.
- [ ] **Resolve §10(a)**: choose grid + landmark (3h/[0,48]/T=48 or 12h/[0,72]/T=72).
- [ ] **Resolve §10(c)**: confirm "cumulative dose" means within-window total,
      not a running sum since intubation.
- [ ] Measure the zero fraction of `total_dose` per window before Phase 2. If
      high, `crimCV`'s zero-inflated Poisson may fit better than `gbmt`.
- [ ] Decide whether the companion sedatives get their own dose columns in Table 1
      now (cheap to add, expensive to backfill), even if unused until later.

---

## 10. Analysis roadmap

Revised 2026-09-05 per SG. Seven gated phases; each produces something
reviewable before the next begins.

### Phase 0 — data structures

**Table 1 — trajectory (long).** One row per **encounter block** per window.

Column set and per-variable rules are declared in
[`config/covariates.json`](../config/covariates.json); §11 is the mirror. This
table is the block-level shape only.

| Block | Columns |
|---|---|
| Keys | `patient_id` (chr), `encounter_block`, `id_num` (int), `window_idx`, `window_start_hr` |
| Fentanyl dose — **three columns, kept separate** | `inf_dose`, `bolus_dose`, `total_dose` (within-window, mcg/hr; rules in §5) |
| Other sedatives | `propofol_dose`, `midazolam_dose`, `dexmedetomidine_dose` — infusions only, each in its own clinical unit (§11) |
| Normaliser | `weight_kg` (fixed at the anchor), `weight_lag_hours` |
| Time-invariant | `age`, `sex`, `race`, `cci`, `bmi_admission`, `bmi_lag_hours` |
| Site / hospital | `hospital_id_admission`, `hospital_id_discharge` (§11) |
| Time-varying severity | `sofa_total`, `nee`, `oxygenation`, `oxygenation_source` |
| Time-varying labs | `bun`, `bicarbonate`, `pco2_arterial`, `lactate`, `inr`, `bilirubin_total` |
| **Status flags** | `imv_status` — required to identify extubated windows and count failed extubations. `crrt_status` |
| Provenance | `<var>_locf` per LOCF-eligible variable; `alive_admitted` per window |

`alive_admitted` is not optional bookkeeping: it is the denominator for every
missingness count, and a post-event window is structurally empty rather than
missing (§11).

**⚠ "At risk" means two different things in this document, so the column does
not use the phrase.** `alive_admitted` (renamed from `at_risk` on 2026-09-07) is
TRUE while the patient is alive, admitted, and inside the extent — it says
nothing about ventilation. Phase 1's number-at-risk row means **still
ventilated**, which is `imv_status == 1` and a very different denominator: at
window 17 (68–72h) `alive_admitted` counts 12,884 episodes and still-ventilated
counts 6,778. Reading the column as the Phase 1 denominator would nearly double
it at 72h and move the dose curve by roughly the +50% that §10 works out below.
A record charted at or after discharge reaches no window at all, so
`imv_status == 1` cannot be true where `alive_admitted` is false.

**Table 1b — hospital intervals.** One row per ADT interval: `encounter_block`,
`hospitalization_id`, `hospital_id`, `hospital_type`, `in_dttm`, `out_dttm`.
Written because `hospital_id` is time-varying and the two columns above are only
its endpoints (§11).

Keeping all three dose columns defers the Model A / Model B choice to analysis
time rather than baking it into the pipeline.

**Table 2 — time-to-event.** One row per **encounter block**. The two planned outcome
analyses need **different event codings**, so a single `event`/`time` pair will
not serve both:

| Analysis | Time origin | Event coding |
|---|---|---|
| Extubation competing with death | landmark T | 0 = censored, 1 = extubation, 2 = death |
| 30-day mortality competing with discharge alive | landmark T | 0 = censored, 1 = death, 2 = discharge alive |

Carry both pairs (or a tidy long form), plus the landmark eligibility flag.

### Phase 1 — descriptive cohort dose trajectory
**No landmark, no exposure window, no outcome model.** All intubated encounter
blocks, contributing for as long as the patient remains ventilated.

| View | Window | Extent | Points | Status |
|---|---|---|---|---|
| Granular | **4h** | 72h | 18 | built |
| Extended | 12h | 7 days | 14 | **deferred, not built** — SG 2026-09-07; see `covariates.json` `windows.extended._STATUS` |

The extended view is descriptive only and gates nothing: T is chosen from the
retention table below, and §8's "Choosing T" tabulation tops out at 72h while
arguing for the *shortest* viable T. Building it means re-running the hourly
waterfall over a longer horizon, not re-aggregating the finished table.

Median/IQR (lead with these — dose is right-skewed), mean/SD alongside, and a
**number-at-risk row** (at risk = still intubated, i.e. `imv_status == 1` — NOT
the `alive_admitted` column; see §10 Phase 0 above). Stopping fentanyl is not
exclusionary: `dose = 0` is a real observation for a ventilated patient.

**Overlay balanced-panel curves** at ≥24h, ≥48h, ≥72h. The all-ventilated curve
answers a different question at every timepoint because its denominator keeps
changing; a balanced panel freezes the denominator so movement reflects real
within-patient change. Crossing = composition; parallel = real change.

**Also plot three curves, not one:** (1) % receiving any fentanyl, (2) median
across all still ventilated, (3) median among those receiving any. Weaning-to-zero among
the still-ventilated and dropout of low-dose patients push the overall median in
opposite directions; a single curve hides both.

**Denominator must be stated in Methods.** A still-ventilated denominator and a
full-cohort-with-zeros denominator answer different questions, and the latter
mostly measures extubation rate rather than dosing. Use the still-ventilated
denominator, and say so.

**The age filter, resolved.** `min_age = 18` excluded exactly zero blocks on
every cohort derivation, which looked like a bug. It is not: `age_at_admission`
is populated for **all 166,814** hospitalizations at UCMC (zero nulls), the
minimum is exactly 18 and 1,143 patients are aged 18. **The extract is adult-only
by construction.** The filter is correct and inert here and will bind at a site
whose extract includes children — keep it, and do not read a zero exclusion as
evidence the field is unpopulated.

#### Federated pooling exports

*(SG, 2026-09-07.)* Phase 1 emits two files whose only purpose is to be pooled
across sites, because **a median cannot be pooled and a display string cannot be
pooled at all**:

| File | Grain |
|---|---|
| `phase1_pooling_continuous.csv` | one row per variable × stratum, and per variable × window × denominator |
| `phase1_pooling_categorical.csv` | one row per variable × level × stratum |

Each continuous row carries `n`, `mean`, `sd`, **`sum`, `sum_sq`**, `min`, `max`
and the median/IQR alongside. Carrying the two sums rather than only mean and SD
is what makes the pooled figures **exact** rather than an approximation that
assumes equal variances:

```
mean_pooled = Σ(sum) / Σ(n)
var_pooled  = (Σ(sum_sq) − Σ(sum)² / Σ(n)) / (Σ(n) − 1)
```

Verified on real numbers rather than asserted: pooling the landmark-eligible and
not-eligible strata as if they were two sites reproduces the overall row to
within the 6-decimal rounding, for every baseline variable
(`tests/test_pooling.py`).

**Means are supplied, not preferred.** Dose here is heavily right-skewed and
about half of ventilated windows are exactly zero, so a mean misrepresents a
site's typical patient in the other direction from the median. Pool the means
because they can be pooled; report the medians because they describe. Both are
in the file.

**Cells below `reporting.small_cell_min_den` (10) are suppressed** — every
statistic blanked and `n_suppressed_small_cell` set — because a mean over n = 1
is that patient's value. No cell at UCMC currently trips it.

**Race is collapsed to Black / White / Other for Table 1 display only** (SG). The
map is in `covariates.json` `time_invariant.race.reporting_collapse`, not in the
R. `Unknown` stays its own row: it means the question was asked and not answered,
not a small race group, and at 8.7% folding it into "Other" would misstate what
that category contains. The **full seven-category distribution is still exported
at full granularity** in the categorical pooling file, so collapsing loses
nothing.

#### What Phase 1 measured at UCMC (2026-09-07)

Run: `code/02_descriptive_trajectory.R`. Every figure below is reproducible from
`output/final_no_phi/phase1_*.csv`.

**Retention, and the landmark choice.** 14,897 ventilation episodes at hour 0;
still ventilated 96.0% at 12h, 78.4% at 24h, 57.1% at 48h, **45.2% (6,728) at
72h** (source: `phase1_choosing_T.csv`). The 72h row equals the landmark cohort
exactly, by construction — both are "ventilated in window 17".

**The three curves move in opposite directions, which is the finding.** The
median across all ventilated falls 19 → 0 mcg/hr and is pinned at zero from hour
24, while the median among receivers *rises* 50 → 75 mcg/hr (IQR 25–100 → 25–150)
and the proportion receiving any falls 63.5% → 35.8% (source:
`phase1_dose_summary.csv`). Read as one curve this looks like steady weaning to
nothing. What is actually happening is that fentanyl exposure **narrows to fewer
episodes rather than falling within them** — those still on it at 72h are on
slightly more than at intubation. This is precisely the artefact the three-curve
rule exists to catch.

**The mean is the summary that survives, and it shows a genuine de-escalation.**
Mean dose across all ventilated falls 47.1 → 32.8 mcg/hr over 72h; within the
frozen ≥72h panel it falls 49.1 → 32.8 (source: `phase1_dose_summary.csv`,
`phase1_balanced_panels.csv`). The median cannot show this because over half of
ventilated windows are exactly zero, so from hour 24 it sits on the floor and
every balanced panel collapses onto the same line — the comparison stops
discriminating precisely where the cohort starts shrinking fastest. **Phase 1's
figures therefore report the balanced panels on the mean and the proportion, not
the median.** This bears directly on Phase 3: a model fitted to a quantity that
is zero in the majority of windows is fitting the zero process as much as the
dose process.

**The decline is real, not compositional.** The ≥24 / ≥48 / ≥72h balanced panels
run parallel to the all-ventilated curve and slightly above it — a mean **+2.66
mcg/hr** on the dose scale (range 0.00 to +5.58) and **+0.72 percentage points**
on the proportion scale (range −0.9 to +1.8) (source: `phase1_balanced_panels.csv`
against `phase1_dose_summary.csv`). Composition contributes a little; the
narrowing of exposure is overwhelmingly within-patient de-escalation.

**Zero fraction — the `gbmt`-versus-`crimCV` number.** **51.7%** of ventilated
windows in the landmark cohort carry `total_dose == 0` (50.8% across the whole
cohort) (source: `phase1_zero_fraction.csv`). The non-zero part is unimodal and
right-skewed with no second mode, deciles 12 → 200 mcg/hr with a median of 75
(source: `phase1_dose_distribution.csv`). That is a zero-inflated continuous
distribution. §9 raised `crimCV`'s zero-inflated Poisson as the alternative;
**SG confirmed 2026-09-07 that Phase 3 remains `gbmt` and Phase 4 remains
`lcmm`**, so the roadmap is unchanged and `crimCV` is not pursued. Record the
zero fraction in Methods regardless — a model fitted to a quantity that is zero
in half of windows is fitting the zero process as much as the dose process, and a
reviewer will ask.

**Why so many zeros — investigated 2026-09-07, they are real.** A median of zero
from hour 24 looks implausible for patients still ventilated at 72h, so the
grid was audited against the raw records rather than assumed correct:

- **Charting cadence is not the cause.** Fentanyl `medication_admin_continuous`
  records are charted hourly at UCMC — inter-record gap median 1.00h, p95 2.58h,
  and only **2.6%** of gaps exceed the 4h `hold_hours` cap. The LOCF cap is not
  manufacturing zeros. This closes the standing "VERIFY THE FENTANYL CHARTING
  INTERVAL" item in `covariates.json`, which had been the single most
  load-bearing unmeasured parameter.
- **The grid is faithful to the records.** Of blocks ventilated in window 6 with
  a raw fentanyl record inside h24–28, only **2.5%** are zero in the grid, and
  every one of those records is a `stop` action carrying a dose of 0 or null.
  Where an infusion was running, the grid shows it.
- **13.9% of the ≥72h panel never receive fentanyl at all** across the whole 72h,
  and at hour 24 about half of the panel is genuinely off it.
- **It is mostly NOT substitution onto another drip.** Dexmedetomidine was added
  as a collected sedative on 2026-09-07 partly to test this, and the measurement
  contradicts the obvious guess. Of the 93,929 ventilated windows with zero
  fentanyl, only **36.1%** carry propofol, midazolam or dexmedetomidine, and that
  share *falls* with time: 65.5% at hour 0, 32.0% at 24h, 25.9% at 48h, **23.3%
  at 68h**. By the end of the window three-quarters of zero-fentanyl episodes are
  on none of the four collected drugs. Propofol accounts for 26.0% of
  zero-fentanyl windows and dexmedetomidine for 13.2%.
- **The pattern is de-escalation, not swapping.** Between hour 0 and hour 68 the
  share receiving fentanyl falls 63.5% → 35.8% and propofol 73.8% → 29.1%, while
  dexmedetomidine rises only 13.0% → 16.5% and midazolam 1.5% → 2.8%. The rise in
  dexmedetomidine (+3.5pp) comes nowhere near offsetting the fall in the other two
  (−72pp combined). Continuous sedation is being withdrawn, with a modest shift
  toward dexmedetomidine among those who stay on something.
- **Residual uncollected agents are small but real.** UCMC also charts
  hydromorphone (27,218 rows), ketamine (23,835), remifentanil (10,050) and
  morphine (8,677) infusions, none collected — together under 11% of fentanyl's
  659,338 rows, so they cannot account for the bulk of the zero-fentanyl windows.
  State as a limitation nonetheless: the study describes *fentanyl* trajectories,
  not total analgosedation intensity.

**Companion sedatives.** Three are collected, infusions only, each in the unit it
is ordered in. Share of ventilated windows carrying any, hour 0 → 68h (source:
`phase1_dose_summary.csv`, `phase1_zero_fraction.csv`):

| Drug | Unit | h0 | h68 | Windows with any |
|---|---|---:|---:|---:|
| propofol | mcg/kg/min | 73.8% | 29.1% | 47.2% |
| dexmedetomidine | mcg/kg/hr | 13.0% | 16.5% | 16.9% |
| midazolam | mg/hr | 1.5% | 2.8% | 1.9% |

Median dose among receivers is flat-to-rising for all three: propofol 27.5 → 30.0
mcg/kg/min, dexmedetomidine 0.50 → 0.80 mcg/kg/hr, midazolam 1.5 → 3.75 mg/hr —
the same narrowing-not-lowering shape as fentanyl. **Midazolam is not a candidate
second indicator here** on prevalence grounds alone (98.1% of windows zero),
independently of the units problem.

**Ventilation episodes and the repeat-patient dependence.** 14,897 episodes from
13,627 patients; **932 patients (6.8%) contribute more than one**, max 16 blocks
(source: `phase1_imv_episodes.csv` and the run log). First-episode duration is
median 27.7h (IQR 13.8–66.1, p95 216.9, max 1,369h). `n_imv_episodes` per block
is median 1 (p95 6, **max 146**) — the tail is long-stay patients repeatedly
coming on and off the ventilator across a single admission, and it is worth
confirming that the 8h `imv_episode_gap_hours` rule is not fragmenting one
course into many. These are the first two of the three dependence diagnostics
§11 requires; the third (classes containing two episodes from the same patient)
cannot be produced until Phase 3 assigns classes.

### Phase 2 — dose trajectory in the landmarked window
Trajectories over [0, T] among patients alive and ventilated at T, with T chosen
from the Phase 1 retention table. Anticipating T = 72h → 18 windows at 4h.

Failed extubations count toward the landmark (§10a). Report the % of the cohort
that would be lost if they were excluded, then make the judgement call with that
number in hand — excluding them is also defensible.

**This is the same cohort as the Phase 1 ≥T balanced panel** (modulo failed
extubations), so the curve should closely reproduce it. Phase 2 is therefore not new analysis so much as promoting that
panel to the analytic cohort. It serves as the cohort-selection tool for
Phases 3–6. Rationale in §8.

#### What Phase 2 produced at UCMC (2026-09-07)

Run: `code/03_landmark_cohort.R`. **T = 72h is the primary analysis** (SG,
2026-09-07), with 24 and 48h swept as the sensitivity §8 requires.

**Landmark flow.** 14,897 analytic episodes → **6,728 alive and ventilated at
T = 72h (45.2%)** → 6,728 analysed. Nothing is dropped after the landmark:
failed extubations are retained per §10a (source: `phase2_landmark_flow.csv`).

**Failed extubations inside [0, T], the §10 disclosure.** Two definitions bound
the same quantity, because a single non-ventilated window may be a charting gap
rather than an extubation (source: `phase2_failed_extubation.csv`):

| Definition | n | % of landmark cohort | Cohort if excluded |
|---|---:|---:|---:|
| any non-ventilated window in [0,T) | 445 | 6.6% | 6,283 |
| a run of ≥ 2 windows (≥ 8h, the episode gap rule) | 343 | 5.1% | 6,385 |

Excluding them would cost about 5%. §10 says make the call with that number in
hand; they are retained, because cohort membership is "ventilated at T", which
needs no look-ahead.

**T sensitivity** (source: `phase2_T_sensitivity.csv`):

| T | n eligible | % of intubated | Windows | Failed extub. | Mean dose | Median dose | % windows zero |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 24h | 11,683 | 78.4% | 6 | 0.2% | 50.6 | **18.8** | 41.1% |
| 48h | 8,504 | 57.1% | 12 | 2.0% | 47.8 | 10.0 | 47.6% |
| **72h** | **6,728** | **45.2%** | **18** | **5.1%** | **45.5** | **0.0** | **51.7%** |

Two things move together as T lengthens and both are worth stating in Methods:
the cohort halves, and the outcome becomes more zero-inflated — **at T = 24h the
median dose is 18.8 mcg/hr, at T = 72h it is exactly zero.** Baseline composition
barely shifts (age 60/60/59, 40.8/40.4/40.0% female, SOFA 8 throughout), so the
retained cohort is not a different kind of patient, only a smaller and
longer-ventilated one.

**Consistency with Phase 1.** §10 says Phase 2 is the Phase 1 ≥T balanced panel
promoted to the analytic set rather than new analysis. The script asserts it:
maximum absolute difference in mean dose and in n across all 18 windows is
**0.0000 and 0** respectively, and a discrepancy raises rather than warns.

**The panel is balanced, deliberately.** Every one of the 6,728 episodes
contributes all 18 windows. Re-filtering on `imv_status == 1` inside the landmark
window would unbalance it — and `gbmt` silently caps the polynomial degree at
(shortest unit's windows − 1), so an unbalanced panel would quietly constrain the
trajectory shape. A transiently extubated window carries dose 0 by §10a(a), which
is the settled treatment and not missingness. `n_ventilated` and
`mean_ventilated_only` ride alongside in the CSV so the difference stays visible.

**Repeat-episode dependence in the landmark cohort** (§11's first two required
diagnostics; the third needs Phase 3): 6,728 episodes from **6,357 patients**;
**266 patients (4.2%) contribute more than one**, accounting for 637 episodes
(9.5%), maximum 9. Milder than the full analytic cohort (6.8%), because
contributing twice requires two separate 72h ventilation courses.

**Handoff.** `output/intermediate_phi/landmark_cohort.parquet` — 121,104 rows
(6,728 × 18), with `id_num` a dense integer rank of `encounter_block` because
`gbmt` and `lcmm` want a numeric unit. Both id columns ride along so the
blocks-per-patient choice stays reversible.

### Phase 3 — `gbmt`, single indicator
`total_dose`, `scaling = 0`, sweep `ng = 1:6`, select by the full §7 conjunction.
Expect difficulty if the zero fraction is high — measure it in Phase 1 first.

### Phase 4 — `lcmm`, same data
`hlme()` off the identical table. Expect **fewer** classes: random effects absorb
heterogeneity `gbmt` can only handle by adding groups. Compare partitions by ARI.

### Phase 5 — two indicators: infusion vs bolus strategy

**Specification**

| Setting | Value | Why |
|---|---|---|
| `x.names` | `c("inf_dose", "bolus_dose")` | both in mcg/hr |
| `scaling` | **`0` — mandatory** | any `scaling >= 1` normalises within unit, erasing absolute dose level; `scaling = 2` additionally divides by a within-patient SD that is exactly zero for a near-all-zero indicator, with no error raised (§3) |
| `nstart` | **≥ 50** | `gbmt` initialises EM from a Ward hierarchical clustering, which is scale-sensitive, so random restarts are needed rather than one deterministic start |
| `d` | 2 (3 if windows allow) | balanced panel, so the degree cap below does not bind |

**What the model now estimates.** Σ_j becomes 2×2 per class, and its off-diagonal
is the within-class infusion–bolus covariance — arguably *the* parameter of
interest for a strategy question. Classes are shapes in two dimensions, e.g.
"steady drip, few boluses" / "low drip, bolus-driven" / "escalating both."

**Checks specific to this phase**

1. Fraction of patients with an all-zero `bolus_dose` across every window. At
   `scaling = 0` these are numerically harmless but contribute a degenerate
   dimension; if it is most of the cohort, the second indicator is buying little.
2. Confirm both indicators actually separate the classes. A class distinguished
   only along the near-zero bolus axis is an artifact, not a phenotype.
3. **ARI against the Phase 3/4 single-indicator classes.** High ARI means the
   second indicator added nothing and Model A already captured the structure.

**Alternative formulation worth testing.** Instead of two dose columns, use
`total_dose` plus **`bolus_fraction`** = bolus_dose / total_dose. This encodes
strategy directly, is bounded [0, 1], and separates "how much" from "delivered
how." Caveat: undefined when `total_dose = 0` — code as `NA`, which `gbmt`
accepts as missing rather than as a zero.

### Phase 6 — class membership as predictor

Classes from Phase 5 (fall back to Phase 4) → competing-risks models, follow-up
starting at T. Classification-uncertainty caveat from §6 applies: classify–analyze
treats estimated membership as observed and understates SEs; acceptable at
APPA > 0.9, otherwise a 3-step/BCH correction is needed. `Jointlcmm` in reserve
as sensitivity analysis only (§8).

#### Outcome set

| | Event of interest | Competing events | Method |
|---|---|---|---|
| **Primary** | successful extubation (≥72h off) | death, **tracheostomy** | CIF (Aalen–Johansen); Gray's test |
| **Secondary** | in-hospital death by day 30 | discharge alive | CIF |
| **Supportive** | ventilator-free days at 28d | — | median (IQR), rank-based |

Event coding for the primary outcome: `0 = censored, 1 = successful extubation,
2 = death, 3 = tracheostomy`.

#### The extubated-day-3, died-day-20 case

The two metrics **disagree on this patient, by design**:

| Metric | Value | Why |
|---|---|---|
| Time-to-extubation, death competing | **event = successful extubation at day 3** | competing risks apply only to events that *prevent* the event of interest. Death on day 20 is after the fact and does not retroactively undo the day-3 extubation |
| Ventilator-free days (28d) | **0** | Schoenfeld convention assigns 0 VFDs to death before day 28 regardless of ventilator status |
| 30-day mortality | **event** | died in hospital |

This is not an inconsistency to resolve — the metrics answer different
questions. Time-to-extubation asks *how quickly are patients liberated*; VFD asks
*how much alive-and-ventilator-free time did they accrue*. An exposure that
liberates fast but kills patients looks good on the first and bad on the second,
and reporting both is the honest presentation. **State this exact example in
Methods.**

#### Reporting best practice

- Use the **cumulative incidence function** (Aalen–Johansen), never 1 − Kaplan–
  Meier, which overestimates incidence in the presence of competing risks.
- Report **both** hazard types (Austin, Lee & Fine 2016, *Circulation*):
  cause-specific hazard (Cox) for aetiology; subdistribution hazard (Fine–Gray)
  for absolute risk / prediction. They answer different questions and
  disagreement between them is informative.
- R: `tidycmprsk::cuminc()` / `crr()`, `cmprsk::cuminc()`, or
  `survival::finegray()` + `coxph()`.
- VFDs are known to be non-normal and to conflate mortality with duration
  (Yehya, Harhay et al. 2019, *AJRCCM*). Report median (IQR) and use rank-based
  tests — never a t-test on means. Consider a win-ratio / hierarchical composite
  (death first, then duration) as a modern alternative.

#### ⚠ Event-time subtlety for "successful extubation"

The ≥72h rule means success cannot be declared at the moment of extubation.
Convention: **event time = time of extubation**, with success determined
retrospectively. Two consequences to specify:

- A patient extubated close to the end of follow-up (e.g. day 28 with data
  ending day 30) has **indeterminate** success. **DECIDED (2026-09-05, SG):
  censor** at last observation. Report how many patients this affects.
- This is a look-ahead in *event classification*, not in cohort selection, so it
  does not compromise the landmark (§10a). Standard for VFD-style definitions,
  but worth a Methods sentence.

#### Tracheostomy as a competing event

Decided 2026-09-05 per SG. A patient who receives a tracheostomy can no longer
experience "successful extubation" as defined, so tracheostomy **competes** with
it alongside death.

**Identification in CLIF — verified against clifpy schemas:**

| Route | Field | Notes |
|---|---|---|
| **Primary** | `respiratory_support.tracheostomy` (0/1) | dedicated flag; trach date = first row with `tracheostomy == 1` |
| Supporting | `respiratory_support.device_category == "Trach Collar"` | trach present and spontaneously breathing |
| Cross-check | `patient_procedures.procedure_code` (`ICD10PCS` / `CPT`) with `procedure_billed_dttm` | billing-derived; useful to validate timing |

Note `device_category` also carries `IMV`, which is what defines the
`imv_status` column in Table 1 and identifies extubated windows.

**Ordering rule.** Events are assessed in time order, but **a failed extubation
is not an event** — it is part of the ongoing ventilation course. So a patient
extubated day 3, reintubated day 4, tracheostomised day 8 has their first
qualifying event at day 8: tracheostomy. Only successful extubation, death, or
tracheostomy terminate follow-up.

**⚠ Tracheostomy is a treatment decision, not a biological event.** Timing is
strongly practice-variable — some centres tracheostomise at day 7, others past
day 14. In a multi-site CLIF study, a between-class difference in the
tracheostomy CIF may reflect local practice rather than the patient's course.
Report tracheostomy timing by site, and name this in limitations.

**The alternative framing a reviewer may raise.** If the question were
*ventilator liberation* (off MV ≥72h by any airway) rather than *successful
extubation* (ETT removed, no reintubation), tracheostomy would not compete — a
tracheostomised patient later weaned off the ventilator would count as
liberated. These are different clinical questions; the chosen outcome is
extubation, and that choice should be stated rather than left implicit.

#### Mortality definition in CLIF

`discharge_category` (hospitalization table). Verified permissible values from
`clifpy/schemas/hospitalization_schema.yaml`:

> Home · Skilled Nursing Facility (SNF) · **Expired** · Acute Inpatient Rehab
> Facility · **Hospice** · Long Term Care Hospital (LTACH) · Acute Care Hospital ·
> Group Home · Chemical Dependency · Against Medical Advice (AMA) · Assisted
> Living · **Still Admitted** · **Missing** · **Other** · Psychiatric Hospital ·
> Shelter · Jail

**Proposed:** death = `Expired` OR `Hospice`. Defensible and common in ICU
cohorts — hospice discharge is functionally a death outcome — but it is a choice
that must be stated and sensitivity-analysed with `Expired` alone.

**Three values need explicit rules:**

- **`Still Admitted`** — outcome unresolved at extraction. Censor at last
  observation; do **not** treat as alive-and-discharged.
- **`Missing`** / **`Other`** — decide and report how many patients are affected.

**Naming precision.** `discharge_category` records disposition *at hospital
discharge*, so it cannot see deaths after discharge. Without linked vital-status
data this is **in-hospital mortality censored at day 30**, not 30-day
all-cause mortality. Framing it as *"death competing with discharge alive"* is
the honest construction and already handles this correctly — discharge alive
removes the patient from observation. Name it accordingly in the manuscript.

---

### Resolved

- ~~Grid/landmark mismatch~~ — Phase 1 descriptive at 4h/72h and 12h/7d;
  Phase 2 exposure [0, T=72h] at 4h (18 windows).
- ~~7-day grid as exposure~~ — descriptive only.
- ~~"Cumulative dose" ambiguity~~ — **within-window total**, expressed mcg/hr.
- ~~Window width~~ — **4h**.
- ~~Other sedatives~~ — propofol, midazolam and dexmedetomidine columns in Table 1.

---

## 10a. Reintubation and the definition of successful extubation

Decided 2026-09-05 per SG, with specifications noted.

### Outcome definition

**Successful extubation = extubation not followed by reintubation within 72
hours**, mirroring ventilator-free-day methodology, which gives no credit for a
failed extubation.

### Cohort inclusion

A patient extubated before hour 72 but reintubated before hour 72 **remains in
the landmark cohort**.

**This is sound, and simpler than it looks.** The rule reduces to *"intubated at
hour T"* regardless of any intervening extubation — which is exactly the standard
landmark risk-set definition. Critically, it requires **no look-ahead**: cohort
membership is determined entirely by the patient's state at hour 72, using no
information from after T. A rule that instead asked "was the hour-30 extubation
successful?" would need data through hour 102 to classify a patient at hour 72,
smuggling post-landmark information into cohort selection. The chosen rule
avoids that.

### ⚠ Two clocks, same number

The 72 in each definition measures a **different interval**:

| Use | Clock starts at | 72h means |
|---|---|---|
| Cohort landmark T | **intubation** | still ventilated 72h after intubation |
| Successful extubation | **extubation** | not reintubated within 72h of coming off |

Consistent thresholds are fine, but they are logically independent choices and
must be described separately in Methods or a reader will conflate them.

### Still to specify

**(a) Dose in the extubated gap — DECIDED: `0`. CLOSED 2026-09-07.**
*(SG: "Boluses of sedation while not ventilated should not count for this
study.")* Implemented in `gate_dose_on_ventilation()` and no longer open. The
estimand is sedation and analgesia **during mechanical ventilation**, and every
affected window is bolus-only with no infusion running — post-extubation PRN
analgesia, a different clinical process. Measured inside the landmark cohort at
T = 72h: **317 windows (0.26%) across 144 episodes (2.1%)**, all bolus-only,
moving any cohort mean by 0.22%. No sensitivity analysis on the ungated column is
planned; `total_dose_ungated` is carried so the decision could be revisited
without another Phase 0 run, not because a planned analysis needs it. Original
statement follows. (2026-09-05, SG.) A patient
extubated at 30h and reintubated at 40h gets `dose = 0` for the intervening
windows: that is what actually happened. The `imv_status` column in Table 1
identifies these windows. **Track and report the number and % of landmark-cohort
patients with a failed extubation inside [0, T]** — it bounds the entanglement in
(b) and supports the Phase 2 judgement call on whether to exclude them.

**(b) Mild entanglement to disclose.** A trajectory class characterised by
"dips to zero mid-window then returns" would, by construction, be largely the
failed-extubation group. This is not the fatal circularity of §8 — failed
extubation is not the outcome, successful extubation is — but the two are
related, and it should be named in the limitations rather than discovered by a
reviewer. Quantifying (a) also bounds this.

**(c) Death and tracheostomy need rules in the outcome definition.** A patient
extubated at 30h who dies at 50h without reintubation must **not** count as a
successful extubation — assign to the competing death event, per VFD convention.
Decide separately whether tracheostomy-and-liberation is a success, a competing
event, or a censoring event.

---

### Still unresolved

**(a) See §10a "Still to specify"** — gap coding, entanglement disclosure, and
death/tracheostomy rules in the outcome definition.

---

## 11. Within-window aggregation rules

**Status: RESOLVED, 2026-09-05 (SG).** This section is the human-readable mirror
of [`config/covariates.json`](../config/covariates.json) (`definition_version`
0.1.0). **That file is the source of truth.** If a value here disagrees with it,
the file is right and this section is stale.

Every time-varying variable needs an explicit rule. There is no safe default, and
the right rule differs by what the variable *is*. The
state-versus-event-versus-absence distinction is the one that matters: apply a
state rule to an event column and you overcount; apply an event rule to a state
column and you undercount; treat an absence as ignorance and you fabricate
missingness where there was none.

Grid: **4h windows, half-open `[start, end)`, 18 windows to 72h**, anchored at
intubation. The 12h/168h extended grid (Phase 1 only) uses identical rules.

### Time-invariant covariates

| Variable | Source | Rule |
|---|---|---|
| `age` | `hospitalization.age_at_admission` | inclusion criterion; missing → drop |
| `sex` | `patient.sex_category` | missing retained as its own level |
| `race` | `patient.race_category` | **Table 1 only — not a model covariate.** See below |
| `cci` | `hospital_diagnosis` (ICD10CM/ICD9CM + `poa_present`) | clifpy `calculate_cci(hierarchy=True)` |
| `bmi_admission` | `vitals` `weight_kg`, `height_cm` | nearest to admission within `[0, +24h]`, else first available; **report `bmi_lag_hours`** |
| `hospital_id_admission` | `adt.hospital_id`, first interval by `in_dttm` | hospital at the start of the encounter block |
| `hospital_id_discharge` | `adt.hospital_id`, last interval by `out_dttm` | hospital at the end; the pair identifies within-block transfer |

**`race` is collected for reporting, not for modelling** *(SG, 2026-09-05)*.
Baseline characteristics need it whether or not any model uses it. Entering it as
a covariate is a separate decision that has **not** been made — the role is stated
in the config so nobody adds it to a formula assuming it was already agreed. Note
that `Unknown` is a permissible CLIF category, not a missing value: report it as
its own row, and count a null `race_category` separately.

**`ethnicity_category` is deliberately not collected** *(SG, 2026-09-05)*. It is
the companion field — race and Hispanic ethnicity are separate questions under the
OMB categories CLIF follows, and Table 1 often reports both — so its absence is a
decision, not an oversight. Do not add it back without asking.

**The two `hospital_id` columns are endpoints, not the whole variable.** They come
from `adt` intervals rather than from the first/last hospitalization, because
`hospital_id` can change *within* one hospitalization and a hospitalization-grain
rule would miss a mid-stay transfer. They differ exactly when a patient moves
between hospitals inside a block; report the count where they disagree (zero at a
single-hospital site, and the pair costs nothing there).

> ⚠ `CRRT-dose-lmtp` collapsed `hospital_id` to the admitting hospital and
> recorded what that cost: it *"silently discarded the very variation that makes
> it worth carrying, and at a single-hospital site the loss would have been
> invisible"* (`code/01_build_cohort.py:116-135`). Phase 0 therefore **also**
> writes `hospital_intervals.parquet` at ADT-interval grain. Cheap now, expensive
> to backfill — the same argument that put `propofol_dose` in Table 1 before
> anything used it. The two named columns are what Table 1 reports; the intervals
> are what a per-window hospital assignment would need.

### Encounter blocks — the unit of analysis

A `hospitalization_id` is one encounter, not one clinical course. Hospitalizations
are stitched with clifpy `stitch_encounters(hospitalization, adt,
time_interval=6)`, matching that repo (`lmtp_design.json:538`) and clifpy's own
default.

**This matters directly here:** a patient intubated, transferred, and still
intubated has **one** ventilation episode. Without stitching their trajectory is
truncated at the transfer and the remainder is either dropped or treated as a
second patient.

Block-level fields: `block_admission_dttm` = min(`admission_dttm`);
`block_discharge_dttm` = max(`discharge_dttm`); `discharge_category` from the
**last** hospitalization by `discharge_dttm`. Assert no block loses its
disposition — the outcome definition depends on it.

**`stitch_encounters` requires `adt.hospital_id`** (with `in_dttm`, `out_dttm`,
`location_category`). That is not optional, and is a second reason `hospital_id`
must be collected — a site that cannot supply it cannot run this step.

**`encounter_block` is the analysis key throughout, Phases 0–6** *(SG,
2026-09-05)*. `id_num` is a dense integer rank of `encounter_block`, carried
because `gbmt` and `lcmm` want a numeric unit. `patient_id` is carried alongside
it — not as a key, but because the problem below cannot be measured without it.

#### ⚠ The one real problem: a patient can have more than one qualifying block

Stitching merges hospitalizations within 6h. It does **not** merge a January
admission with a June one. A patient intubated on both occasions produces **two
encounter blocks**, and under `encounter_block` as the analysis key they enter the
model as two independent units.

They are not independent. Two ventilation courses in one patient share
comorbidity, physiology, and often the same unit and clinicians. The consequences
run in a known direction:

- **Standard errors on class membership are understated** — the effective sample
  size is smaller than the row count.
- **A spurious class can appear.** GBTM will happily fit a "patient signature"
  group whose coherence is repeat admissions rather than a distinct dosing
  strategy — and that class would look clinically interesting.
- **Phase 6 inherits it.** Competing-risks estimates assume independent
  observations too.

**We cannot fix this the way `CRRT-dose-lmtp` does.** That study also uses
`encounter_block` as the row unit, but hands `lmtp` a separate clustering `id` =
`patient_id` (`03_lmtp_fit.R:1086`), which absorbs the correlation. **No method in
this pipeline accepts one.** Verified by execution, 2026-09-05:

| Function | Arguments | Cluster argument? |
|---|---|---|
| `gbmt::gbmt` | `x.names, unit, time, ng, d, data, scaling, pruning, delete.empty, nstart, tol, maxit, quiet` | **none** |
| `cmprsk::crr` | `ftime, fstatus, cov1, cov2, tf, cengroup, failcode, cencode, subset, na.action, gtol, maxiter, init, variance` | **none** (`cengroup` stratifies the censoring distribution, it does not cluster) |

So the correlation cannot be absorbed at the model stage.

#### First, a correction: `gbmt` runs either way

"No clustering argument" is not "will not work". **Verified by execution,
2026-09-05:** on a synthetic cohort of **168 episodes from 120 patients, 48 of
whom contributed two**, `gbmt` converged in **5.9s** at `ng = 3` and returned
three clean groups (sizes 54 / 44 / 70). Repeat episodes never prevent the model
from fitting.

**The cost is inferential, not functional**, and it lands in one place that
matters more than the others:

- **Class enumeration (§7).** Correlated units inflate the log-likelihood faster
  than BIC's `npar × log(ss)` penalty grows, so the criterion tilts toward
  **more** classes. That is a directional bias on the single most consequential
  decision in the study.
- **Stage-2 standard errors** are understated; CIs too narrow.
- **A "patient signature" class** can appear — coherent because of repeat
  admissions rather than a dosing strategy.

The magnitude is unmeasured. `validation/repeat_encounter_cost.R` is written to
measure it (BIC-selected `ng` first-episodes-only vs all-episodes, over a sweep of
repeat fractions, with ARI between the solutions) but **has not been run to
completion** — the full sweep is slow and was deprioritised.

#### Current default: keep all blocks — and the choice stays open

*(SG, 2026-09-05.)* **A patient admitted twice for different reasons may have
totally different fentanyl trajectories, and each is clinically interesting.** The
study therefore describes **ventilation episodes**, not patients, and a patient
contributing two episodes contributes two genuinely distinct observations rather
than one observation counted twice.

> **This is deliberately reversible.** *"We will have `patient_id` and
> `encounter_block` as columns and can decide later which to use and whether to
> exclude patients with multiple encounters."* Phase 0 carries **both** on every
> row, so switching to first-episode-only is a one-line filter on the finished
> table — sort by block start, group by `patient_id`, keep the first — not a
> pipeline change. Nothing downstream may hardcode the choice; read
> `encounter_blocks.blocks_per_patient` from the config.

This resolves the estimand question rather than the statistical one, and the
difference is worth being precise about:

- **What it settles.** The unit of analysis is now stated and defensible. The two
  episodes are not a duplicate; they are different clinical courses, and a design
  that discarded the second would discard real information.
- **What remains.** Residual within-patient correlation is still present, and no
  method here can absorb it. Its size is an empirical question, not a
  philosophical one, and it depends entirely on **how many patients repeat and how
  alike their episodes are.**

**Therefore these are required outputs, not optional diagnostics:**

1. The number and % of patients contributing **more than one** encounter block,
   and the distribution of blocks per patient. This is the single number that
   bounds everything above — if it is 2%, nothing here matters; if it is 20%, the
   limitations paragraph has to be specific.
2. For any retained class solution, the number of classes containing **two
   episodes from the same patient**, against what independence would predict. A
   "patient signature" class is the concrete failure mode, and this is what
   detects it rather than assuming it away.
3. A limitations sentence naming the dependence and the fact that `gbmt` and
   `crr` cannot cluster on it.

Both counts come free from Phase 0 and cost nothing to carry. They are what makes
the decision above defensible to a reviewer rather than merely stated.

#### What stitching fixes, and what it does not

A benefit worth stating: `discharge_category` includes **`Acute Care Hospital`**,
which is a transfer *out*. Without stitching, a patient transferred to a partner
hospital and dying there reads as "discharged to Acute Care Hospital" — alive —
and the mortality outcome is wrong. Taking `discharge_category` from the **last**
hospitalization in the block fixes that.

It only fixes it **when both hospitalizations are in this CLIF dataset.** A
transfer to an outside system is still lost, and still reads as a discharge alive.
That is a limitation to state, not something stitching solves.

Component column names for CCI are **read from clifpy, never transcribed** — a
hand-written list drifts, and the first draft of exactly such a list in
`CRRT-dose-lmtp` got five of seventeen names wrong.

### Time-varying covariates

| Variable | Type | Summary | LOCF | Cap | Missingness class |
|---|---|---|---|---|---|
| `inf_dose`, `bolus_dose`, `total_dose` | exposure | **see §5** — not restated here | — | — | drip off is a true `0` |
| `sofa_total` | derived score | scored from filled components | **no** | — | `time_varying` |
| `nee` | state (rate) | **max of the summed step function** | **no** | — | `absence_means_zero` |
| `oxygenation` | derived ratio | **min** (worst) | yes | 8h | `time_varying` |
| `oxygenation_source` | provenance | stamped, not summarised | never | — | own level, never imputed |
| `imv_status` | state | **any in window** | **no** | — | `absence_means_not_ventilated` |
| `crrt_status` | state | **any in window** | **no** | — | `absence_means_zero` |
| `bun` | measurement | **max** | yes | 24h | `time_varying` |
| `bicarbonate` | measurement | **min** | yes | 24h | `time_varying` |
| `pco2_arterial` | measurement | **max** | yes | 24h | `time_varying` |
| `lactate` | measurement | **max** | yes | 24h | `time_varying` |
| `inr` | measurement | **max** | yes | 24h | `time_varying` |
| `bilirubin_total` | measurement | **max** | yes | **72h** | `time_varying` |

`bicarbonate` is the one lab summarised by `min`: unlike the other five, its
abnormal direction is down.

`imv_status` uses **any-in-window** rather than `CRRT-dose-lmtp`'s
state-at-window-end-boundary. At 4h resolution a boundary state discards most of
the window. This is a deliberate divergence, not an inheritance.

### The three missingness classes

**1. `time_varying`** — no observation means NA, then LOCF from the most recent
earlier window, **subject to a per-variable cap**. Members: `sofa_total`,
`oxygenation`, and the six labs.

**2. `absence_means_zero`** — set to `0` for every `alive_admitted` window with no record;
**not LOCF-eligible**; the count set to zero is reported. Members: `nee`,
`crrt_status`.

> For a continuously infused medication, the absence of a record does not mean the
> value was not measured — it means the drug was not running. Carrying the previous
> window's vasopressor dose forward into a window with no infusion recorded would
> *invent* pressor exposure for a patient who had been weaned off it, most often in
> exactly the recovering patients. Labs are the opposite case: an unmeasured
> lactate is unknown, not zero.

`CRRT-dose-lmtp`'s first run reported `nee` 17% and `inotrope` 84% "missing" —
which was **absence being mislabelled as ignorance**. *(source:
`CRRT-dose-lmtp/config/lmtp_design.json:468`)*

**3. `absence_means_not_ventilated`** — same treatment, for `imv_status`.
Invasive ventilation is charted continuously, so a window with no device record
was almost certainly not on a ventilator rather than unobserved.

### Why the LOCF caps exist, and why they are per-variable

`CRRT-dose-lmtp` has **no** time cap on its node-level LOCF because its grid is
three 24h nodes — the structure caps the carry at ~48h. **Ours does not.** With
18 windows an uncapped carry reaches 72h, and the sibling repo shipped exactly
that bug: `CLIF-epidemiology-of-CRRT` forward-filled labs with no `limit=`, and
its baseline lactates were **~20% stale carry-ins, median lag 27h, maximum ~1331h
(~55 days)** — and the staleness was **differential by treatment arm (28.7%
low-dose vs 13.7% high-dose)**, which is a bias, not noise. *(source:
`crrt-manuscript-tools/.claude/lessons.md:581`; fix at
`CLIF-epidemiology-of-CRRT/code/02_construct_crrt_tableone.py:487-500`)*

**All six labs carry a 24h cap** (6 windows). *(SG, 2026-09-05.)* This supersedes
the per-variable clinical half-lives first drafted here (lactate 8h, pCO₂ 8h,
bicarbonate 12h).

The reasoning is about the decision being represented, not the kinetics: **a
clinician acts on the last value available to them**, not on the value the
half-life would justify, and an ICU patient gets minimum daily labs. A 24h carry
is therefore both what the sampling supports and a faithful representation of the
information the bedside actually had.

| Cap | Variables | Reasoning |
|---|---|---|
| **24h** (6 windows) | `bun`, `bicarbonate`, `pco2_arterial`, `lactate`, `inr` | Minimum daily labs in the ICU; the clinician acts on the last available value. |
| **72h** (18 windows) | `bilirubin_total` | *(SG, 2026-09-06.)* The one lab that departs from the uniform rule. It is the slowest-moving of the six — bilirubin changes over days — and is not drawn daily in most ICU patients. **Measured:** at 24h it was present in only **14.1%** of `alive_admitted` windows, which capped 6-component SOFA at roughly that figure regardless of every other component. The kinetics support the longer carry; 24h was set by charting cadence, the wrong constraint for this analyte. |
| **8h** (2 windows) | `oxygenation` | **Not a lab.** SpO₂ is charted at least hourly, so an 8h gap in oxygenation is a data fault rather than a draw-cadence artefact. The daily-labs argument does not extend to it. |

**The cap still binds.** It bounds any carry at 6 windows and rules out the
uncapped failure mode entirely. But 24h is loosest exactly where the biology is
fastest — lactate's plasma half-life is ~20 min, and pCO₂ tracks minute
ventilation, which changes within minutes of a vent adjustment (frequent by
construction in a cohort anchored at intubation). **The `<var>_locf` flags are
what make that checkable**: before leaning on `lactate` or `pco2_arterial` in an
analysis, look at what fraction of their values were carried rather than
measured.

Per-variable caps are retained in the schema, so a single variable can be
tightened later without restructuring anything.

**Whole-gap, not carry-with-expiry.** A gap *longer* than the cap is NA for its
entire length, not carried for `cap_hours` and then dropped. The rule classifies
the gap, and a gap judged too long to bridge was too long throughout it. *(ported
from `CRRT-dose-lmtp/code/02_build_lmtp_df.py:559-621`)*

**Never extrapolate** past the patient's last observation of that variable, or
past the end of their `alive_admitted` period. Extrapolating past the last charted value is
silent and biases in the same direction a real effect would.

**First window** with no earlier observation stays NA. Imputation belongs with the
model, not with dataset construction.

### NEE

Six drugs, summed as a step function, **maximum over the window**:

| Drug | Factor | Preferred unit |
|---|---|---|
| norepinephrine | 1.0 | mcg/kg/min |
| epinephrine | 1.0 | mcg/kg/min |
| phenylephrine | 0.1 | mcg/kg/min |
| dopamine | 0.01 | mcg/kg/min |
| **vasopressin** | **2.5** | **u/min** (not weight-based) |
| angiotensin | 10.0 | mcg/kg/min |

*(source: `CRRT-dose-lmtp/config/lmtp_design.json:193-198`, verified by direct
read 2026-09-05)*

A row-wise maximum is wrong: a patient on moderate doses of three pressors is
sicker than one on a slightly higher dose of a single agent. The step function
also holds a rate forward between records (`hold_hours: 4`), so a drug charted at
01:00 still counts at 01:30 when a second is charted. `mar_action_category ==
"stop"` is a rate of **0**, not a missing value.

> ⚠ **`hold_hours: 4` is inherited, not verified here.** It rests on a median
> inter-record interval of 57 min (p90 68 min) *at the CRRT coordinating site*.
> Confirm the vasopressor charting interval at UCMC before relying on it.

**Source — Goradia et al.** *(supplied by SG, 2026-09-05.)*

> Goradia S, Abu Sardaneh A, Narayan SW, Penm J, Patanwala AE. **Vasopressor dose
> equivalence: A scoping review and suggested formula.** *J Crit Care*
> 2021;61:233–240. [doi:10.1016/j.jcrc.2020.11.002](https://doi.org/10.1016/j.jcrc.2020.11.002).
> PMID **33220576**.

Its suggested formula matches this table exactly — phenylephrine ÷ 10, dopamine ÷
100, vasopressin × 2.5 per u/min, angiotensin II × 10. Verified against Crossref
and PubMed, 2026-09-05.

Worth recording: `CRRT-dose-lmtp` carries these same six factors with **no
citation anywhere in that repo**, its only support being a dimensional check. Cite
Goradia, not the sibling repo.

**The unit guard is mandatory and must raise.** clifpy does not null a dose it
cannot convert — it leaves the *raw* value in `med_dose_converted` and reports the
failure in `med_dose_unit_converted`. A check-for-NA guard therefore sees nothing
wrong, and 20 ng/kg/min of angiotensin becomes "20 mcg/kg/min", which a
coefficient of 10 turns into 200. Measured in `CRRT-dose-lmtp` before its guard
existed: **`nee` reached 8,001 mcg/kg/min-equivalent**. Raise rather than drop —
dropping removes a drug from NEE for the patients who received it, biasing a
confounder in a known direction in the sickest patients. *(source:
`CRRT-dose-lmtp/code/02_build_lmtp_df.py:1135-1179`)*

Bounds must be applied to **raw** doses per (drug, charted unit) *before*
conversion; they are written in the converted unit and are meaningless against a
raw one.

### Oxygenation — one column, not two

`oxygenation` is **one covariate on the P/F scale**: measured P/F where an
arterial gas exists, Severinghaus-derived P/F otherwise. The fallback is **per
window, not per reading**, so a window's value comes from exactly one source, and
`oxygenation_source` names it (`pf`, `pf_room_air`, `sf`, `sf_room_air`, `none`).

This is the correction of a real defect. `CRRT-dose-lmtp` originally carried
`pf_ratio` and `sf_ratio` as separate covariates intending S/F to act as a
fallback, and **nothing implemented the fallback** — both entered as independent
covariates, each median-imputed when absent, while **85% of blocks missing P/F had
S/F observed**. The model was filling a constant into a column whose information
sat in the next one. *(source: `lmtp_design.json:419`)*

**Severinghaus, not Rice.** Rice et al. (Chest 2007) gives a linear
`S/F = 64 + 0.84 × P/F`. Brown et al. (Chest 2016, **PMID 26836924**) showed a
Severinghaus-based nonlinear imputation beats linear and log-linear on both error
and mortality association, largest at low P/F, and prospectively validated it
(**PMID 28538439**). The SpO₂–PaO₂ relationship is sigmoidal, so no linear map
fits across the range. clifpy implements exactly this in its SOFA respiratory
component, so following it also keeps us consistent with the CLIF reference
implementation. Rice was measured against the CRRT cohort before being rejected:
fitted slopes of **0.44 / 0.39 / 0.37 against Rice's 0.84**, degrading with P/F.

Gate: SpO₂ **strictly below 97**. Above it the dissociation curve is flat and an
imputed PaO₂ reports the FiO₂ rather than the patient. Severinghaus(97) = 90.6
mmHg; SpO₂ = 100 returns NaN.

**FiO₂ must be a fraction, and that is enforced rather than assumed** *(SG,
2026-09-05)*. `code/utils/fio2.py`, tested in `tests/test_fio2.py` (10 checks,
wired into both runners). Getting it wrong is silent: clifpy's SOFA gates on
`fio2_set BETWEEN 0.21 AND 1` (clifpy `utils/sofa.py:273`), so at a percent-scale site every
FiO₂ — and every P/F with it — is nulled with no error raised.

The rule that makes this safe:

> **Scale is decided at the column level. Bounds are applied at the value level.**

A column is divided by 100 only when ≥95% of its non-null values sit in `[21, 100]`.
An individual out-of-range value in an otherwise-fractional column is a **data
entry error, not a unit**, and is nulled rather than rescaled — otherwise
`fio2 = 88880` becomes a plausible-looking number instead of the discard it should
be. A column that is neither clearly fraction nor clearly percent **raises**:
mixed units inside one column is a site data problem for a person, not something a
heuristic should resolve. Measured behaviour on the three cases:

| Column | Detected | Action |
|---|---|---|
| 99.5% in `[0.21, 1]`, one value of 88880 | fraction | left alone; **1 value nulled**, and the report names it |
| 100% in `[21, 100]` | percent | **whole column ÷ 100**; 0 nulled |
| 50/50 fraction and percent | — | **raises**, naming both shares |

**Plateau windows are left NA** *(SG, 2026-09-06 — option A)*. Measured before
deciding: **78,187 `alive_admitted` windows (31.1%)**, FiO₂ pairable for **92.8%**, implied
floor `Severinghaus(96.99)/FiO₂` with **median 226, 81.4% below 300, 21.4% below
200**. Three consequences belong in the limitations:

1. These are **not a healthy subgroup** — their true P/F is right-censored, not high.
2. Whatever fills them downstream learns from windows where oxygenation *was*
   measured, and an ABG is drawn when someone is worried — so the filling
   distribution is **sicker** than these patients.
3. A window at FiO₂ 0.30 (floor 302) and one at FiO₂ 1.00 (floor 90) are filled
   **identically**, discarding the FiO₂ information entirely.

Rejected alternatives: **(B)** carry the floor as the value — it systematically
understates and would fabricate a P/F in 31% of windows; **(C)** leave NA but
carry FiO₂ as its own covariate — keeps the information without inventing a
ratio, and remains the natural next step if the plateau proves to matter. The
diagnostic that produced these numbers runs every time, so revisiting is cheap.

**FiO₂ lookback = 4h.** Pair every PaO₂ (and every qualifying SpO₂) with the most
recent non-null `fio2_set` **at or before** it, searching back at most 4h —
`merge_asof(direction="backward", tolerance=4h)` *is* that rule. Never pair to a
future FiO₂.

The eight-step pipeline is written out in one place in `covariates.json`
(`oxygenation.pipeline`) rather than scattered, because each step was agreed
separately and that is what makes the chain impossible to see. Three parts of it
carry warnings worth repeating here:

- **clifpy's waterfall has a known defect we must repair, not accept.** Its
  `fill_block` treats trach collar as a segment breaker and does not confine the
  damage to trach-collar rows — **it wipes `fio2_set` across the whole
  encounter**. At the CRRT coordinating site, 189/2,141 encounters (8.8%) contain
  a trach collar and carry **50.4% of all post-waterfall FiO₂ nulls**; 100% of the
  10,748 room-air rows still missing FiO₂ sit in a trach-collar encounter, even
  though clifpy sets room air to 0.21 — *it sets the value and the fill then
  destroys it*. Repair by re-merging the raw table. **Report upstream.**
- **HFNC, CPAP, NIPPV, Other and Trach Collar are left missing on purpose**, and
  the count is reported. This is the substantive decision, not an omission: on
  those devices flow and FiO₂ are independent settings, so a formula would
  manufacture values that look measured.
- **S/F cannot rescue a missing FiO₂.** FiO₂ is the shared denominator of both
  ratios, so S/F is a fallback for a missing PaO₂ and *never* for a missing FiO₂ —
  when FiO₂ is absent both die together. The `room_air_when_unmonitored` rule
  addresses that case and is deliberately narrow: it fires only where **no**
  `respiratory_support` row exists within the lookback, never where a row exists
  carrying a null `fio2_set` (there the device mix is dominated by nulls and nasal
  cannula, and nasal cannula at unknown flow is not 0.21).

### SOFA — computed here, from our own per-window components

`sofa_total` is **not** taken from clifpy. The six components are scored from the
same per-window aggregates the rest of the table uses, so `sofa_resp` is derived
from the same `oxygenation` column the analysis reports. Cutpoints are Vincent
1996; implementation in `code/01_build_cohort.py:score_sofa`, 9 tests.

Inputs: `map` (min), `platelet_count` (min), `bilirubin_total` (max),
`creatinine` (max), `gcs_total` (min), `oxygenation` (min), `imv_status`, and the
four SOFA vasopressors as max mcg/kg/min.

**A correction to an earlier draft of this section.** It listed clifpy's
cardiovascular component as defective for counting only dopamine, dobutamine,
epinephrine and norepinephrine — omitting vasopressin, phenylephrine and
angiotensin. **That was wrong.** The original SOFA cardiovascular score is
*defined* on those four agents; a patient on vasopressin alone scoring by MAP is a
limitation of the score itself, not of the implementation. clifpy and
`CLIF-epidemiology-of-CRRT` both render it faithfully. Recorded because the
earlier claim would have justified a "fix" that silently redefined the score.

What is genuinely worth changing:

| | clifpy / epi repo | Here |
|---|---|---|
| **Respiratory below P/F 200 off IMV/NIPPV/CPAP** | returns **NULL** | scores **2** (<200) / **3** (<100); on the vent, 3 / 4 |
| **Missing components** | `fill_na_scores_with_zero` defaults to scoring them **0**, i.e. normal | summed with `min_count=1`, and **`sofa_n_components` records how many were available** |
| **FiO₂ scale** | gated `BETWEEN 0.21 AND 1`, nulling everything at a percent-scale site | normalised first (§11 above), so the gate is moot |
| **Renal** | creatinine only | same — CLIF core has no `intake_output`, so urine output is **absent, not proxied**. A stated limitation. |

**The null-versus-zero trap, which is the reason this needed care.** The standard
cardiovascular chain scores 3 on `epinephrine ≤ 0.1`, meaning *receiving* a low
dose. If an absent drug is coded `0` rather than left NULL, `0 ≤ 0.1` is true and
**every unpressored patient scores 3.** Demonstrated: a patient on no pressors
with MAP 85 scores **0** with NULLs and **3** with zeros.

This collides directly with this study's own `absence_means_zero` convention for
`nee` (§11 above). The two are deliberately different: **`nee` zeroes an absent
infusion; the SOFA pressor inputs must not.** `_sofa_pressors()` keeps them NULL,
and a test asserts both halves.

### Weight: two rules, deliberately different

Weight is a **normaliser, not a covariate**, and it needs two rules that must not
be unified:

- **Dose denominator — fixed at the anchor (intubation).** If the denominator
  moves, a weight-normalised dose can change because the *weight* changed rather than
  *dosing*, and the trajectory shape becomes partly an artefact of fluid balance.
  Match backward first (a weight recorded after intubation already reflects
  resuscitation), forward only for the residue, and **report the lag** —
  `CRRT-dose-lmtp` declined to cap this with a `tolerance=`, on the grounds that
  choosing a cutoff is a study decision rather than a coding one, and reported
  median 14h / p95 131h / max 585h instead.
- **NEE denominator — current weight** (clifpy `find_most_recent_weight`). NEE is
  an intensity the clinician is titrating now.

Do not pass the cohort's fixed per-patient weight into clifpy's converter, which
would silently make NEE use the fixed weight too.

### Reporting missingness — the order of operations matters

1. Complete the (patient × window) grid, so "no row" and "row with NA" become the
   same thing and every count uses the same denominator.
2. Apply `absence_means_zero` and `absence_means_not_ventilated`.
3. **Count missingness — `alive_admitted` windows only.**
4. **Then** LOCF.
5. Emit a parallel `<var>_locf` boolean per variable.

Counting *after* the fill makes the extent of filling invisible. The `alive_admitted`
restriction matters too: a post-event window has no covariates because follow-up
had ended, which is **structure, not data quality** — mixing the two would make
late windows look far worse than they are.

Report **per variable** *and* **per pattern**: a per-variable table cannot show
which variables go missing *together*, and that is what determines whether an
imputation model is well posed. Suppress pattern cells below 11, rolling the
remainder into one "other" row that states how many it absorbed.

Phase 0 emits `output/final_no_phi/diagnostics/phase0_missingness.csv` with **both sides of
the fill**, since they answer different questions — how much was carried, and
what is still absent in the analysis data:

| Column | Meaning |
|---|---|
| `n_observed` | actually measured in the window |
| `n_zero_by_rule` | set to 0 by `absence_means_zero` / `absence_means_not_ventilated` — **not** missing |
| `n_missing_pre_locf`, `pct_missing_pre_locf` | before any carry-forward |
| `n_filled_by_locf`, `pct_filled_by_locf` | how much the cap actually carried |
| `n_missing_final`, `pct_missing_final` | what remains for the model to handle |

`n_observed + n_filled_by_locf + n_missing_final == n_alive_admitted` for every
LOCF-eligible variable, which is asserted in `tests/test_build_cohort.py`.
Separating `n_zero_by_rule` is what stops an absent vasopressor record reading as
a data gap — the mislabelling that made the reference repo report `nee` as 17%
"missing" at its first node.

`oxygenation_source` is reported alongside as its own breakdown, so the share of
oxygenation resting on the room-air assumption is visible rather than folded in.

**Oxygenation absence is decomposed by cause** *(SG, 2026-09-06)*, as three
mutually exclusive rows in the same CSV. "Missing" is not one thing here, and the
three have different remedies:

| Reason | What it means | Remedy |
|---|---|---|
| no PaO₂ and no SpO₂ measured | nothing to pair | none — not recoverable from this data |
| SpO₂ present but all ≥ 97 | the Severinghaus transform is undefined on the plateau | the true P/F is **right-censored, not high** — see the limitation above |
| usable measurement but no FiO₂ within the lookback | a PaO₂ or a sub-ceiling SpO₂ existed and could not be paired | the **4h `fio2_lookback_hours`** is the binding constraint, and widening it is a real option |

Only the third is a lookback problem. Reporting them together is what stops a
lookback fix being applied to a plateau problem, or vice versa.

**`pco2_venous` was considered as a fallback for `pco2_arterial` and rejected**
*(SG, 2026-09-06)*. It would have recovered a great deal — 181,540 venous results
against 635,165 arterial, with 1,495 hospitalizations (8.1%) holding venous and no
arterial — and Chong et al. (PMID 33780397), already cited here for the pH offset,
gives arterial pCO₂ = venous pCO₂ − 5 mmHg. It is rejected because **CLIF does not
distinguish central from peripheral venous sampling**, and peripheral venous pCO₂
diverges from arterial most severely in shock and post-arrest patients — precisely
this study's population. A single fixed offset would be least valid exactly where
it matters most. The missingness is accepted and reported instead.

**No variable is dropped for excess missingness.** A variable both frequently
missing and poorly predicted by the others is as likely an extract or mapping
problem as genuine clinical non-measurement — that is a finding, not a reason to
discard it.

### Defining a continuous IMV episode

CLIF has **no intubation event** — only `device_category` transitions — so the
anchor is the start of the first continuous IMV episode, a *proxy* for
intubation. A patient transferred in already ventilated, or admitted with a
tracheostomy, has hour 0 = first observation rather than start of ventilation.
That is a limitation to state, not one we can fix.

An episode ends at whichever comes first:

| signal | why it is needed | measured |
|---|---|---|
| **waterfalled IMV → non-IMV transition** | the waterfall exists to fill charting gaps, so a transition it shows is a real device change | over 300 hospitalizations, **all 433** transitions go to a real device (NIPPV, nasal cannula, trach collar, …) and **none to a null device** |
| **raw IMV gap > `imv_episode_gap_hours`** | where no subsequent device is ever charted, nothing breaks the segment and the waterfall carries IMV forward | it over-extends past the last raw IMV record by **p90 30h, max 1,137h**, in **19.3%** of hospitalizations — those extubations are invisible to it |

Neither is sufficient alone. The waterfall is precise where a transition exists;
the gap rule catches the one-in-five where none is ever charted.

**`imv_episode_gap_hours = 8`** *(SG, 2026-09-06)*, set from the data rather than
convention. Raw IMV inter-record gaps are **p50 1.00h, p75 3.33h, p90 4.37h, p95
5.00h, p98 7.07h** — so only **1.5%** exceed 8h. Far enough past routine charting
not to fragment a continuous episode, close enough to catch a genuine unobserved
extubation, and it accommodates a long off-unit absence such as CT followed by a
lengthy operation.

A 4h threshold — tempting, since it matches the analysis window — would have been
**wrong**: 14.6% of *normal charting* gaps exceed 4h.

The episode **ends at the last IMV record**; the transition is what tells us it
ended rather than continued. `episode_ended_by` records which signal fired.

**`min_imv_hours = 4`** is a separate parameter, deliberately tied to
`window_hours`: a block that cannot fill one analysis window has no trajectory to
model.

### The waterfall runs on the whole encounter block, and is cached

Span trimming to `[anchor − 24h, anchor + 96h]` was **retired 2026-09-06 (SG)**.
It was sized for the 72h trajectory window, but the **30-day** mortality and
extubation outcomes need ventilation status out to day 30 — a 96h window could
never have served them. Removing it also removes a circularity: the trim was
derived *from* the anchor, so the anchor could not be derived from the trimmed
table.

The equivalence measurements that justified trimming remain in the config as the
record of why it was safe for the trajectory window; they no longer describe what
runs.

Full-stay waterfalling is affordable because it is **cached**
(`code/utils/waterfall_cache.py`). The cache is content-addressed and stored **per
hospitalization**, which is sound because the waterfall is verified
**per-encounter independent** — a 24-id run restricted to 12 gives identical rows
and zero differing cells against a 12-id run.

The cohort therefore **never enters the key**:

| cohort change | result |
|---|---|
| narrows (new exclusion) | **full hit** — nothing recomputed |
| widens (relaxed criterion) | computes **only the new patients** |

The key pins the source file's size and mtime, the clifpy version, and a **hash of
our own transform source** — so editing `_canonicalise_devices` invalidates it
automatically, rather than relying on a version constant someone must remember to
bump. A stale cache cannot be used silently: any changed input yields a different
key and a miss.

### Re-runs must not read stale outputs

Phase 0 is re-run often, and several runs during development died part-way. A
crash between two writes leaves a **mismatched pair** from two different code
versions, which is worse than either file being absent — and a stale parquet on
disk looks perfectly valid to whatever reads it next.

Three mechanisms, in order:

1. **`clear_owned_outputs()` runs first.** The script declares the files it owns
   and deletes them before doing any work, so a crash leaves nothing.
2. **`phase0_manifest.json` is written last**, after every output has succeeded.
   It carries the code version, the SHA-256 of `config.json`, `covariates.json`
   and `outlier_config.json`, and the row counts. Its presence *is* the
   completion signal.
3. **`require_manifest()` gates every downstream phase.** An absent manifest, or
   a config digest that no longer matches, stops the phase with a message telling
   the reader to re-run rather than analyse stale tables. Implemented in both
   `paths.py` and `paths.R`, which must agree.

**Both tables are also written as `.csv` beside the `.parquet`** for human review
*(SG, 2026-09-06)*. The parquet is what later phases read; the CSV exists so the
tables can be opened and checked. Both are **inside `output/intermediate_phi/`** —
a per-patient-window CSV is PHI regardless of format, and a test asserts it is
covered by the ignore rule.

### What enforces all of this

`config/covariates.json` is only a policy statement unless something checks it.
Two mechanisms:

- **`tests/test_fio2.py`** — 10 checks on the unit rule above.
- **`tests/test_outliers.py`** — 14 checks on `config/outlier_config.json`, the
  only place a bound is written. Bounds live in **two layers** because a bound is
  only comparable to a value already in its unit: `med_dose_raw` per (drug,
  charted unit) applied *before* conversion, `med_dose_converted` per drug applied
  *after*. Angiotensin is bounded there rather than inherited — clifpy has no
  angiotensin entry, which made its bounds a partial net. Two derived sanity
  ceilings are computed rather than written down, so they cannot drift when a
  bound changes: fentanyl **11,000 mcg/hr** and NEE **17.45 mcg/kg/min-equiv** —
  the latter reproducing exactly the figure `CRRT-dose-lmtp` reports, which is an
  independent check that both ported tables match theirs.
- **`tests/test_covariates.py`** — 20 static checks, all verified to fire by
  breaking them: summary rules are in the dispatch vocabulary; every variable has
  exactly one missingness class; class membership lists agree with the
  per-variable declarations; LOCF-eligible variables have caps and ineligible ones
  do not; no cap is shorter than one window (a guaranteed no-op reads as policy
  and is not); absence-means-zero variables are never LOCF-eligible; NEE
  coefficients match the source category list and carry a citation with a DOI;
  window arithmetic is self-consistent; **the grid declared here matches the one
  in the site config**, which is the classic "same threshold in two places, now
  drifted" failure; **and the table in this section matches the config** —
  variable for variable, cap for cap, class for class.
- **That last check is why this section can be trusted.** A hand-written mirror of
  a machine-readable file is exactly the pair that drifts:
  `CRRT-dose-lmtp/docs/lmtp_df_build_notes.md` is stamped `0.2.0` against a
  `0.13.0` config and is wrong on three counts as a result — the outlier source,
  the pH handling, and the very separate-S/F-column design this section rejects.
  It even states its own precedence rule ("if a number below disagrees with that
  file, the file is right and this document is stale"), which is an admission that
  nothing enforces it. Here something does.
- **A consumption assertion in `01_build_cohort.py`** *(still to write)* — every
  variable the config declares must reach a column in `trajectory_long.parquet`,
  or the build fails loudly. Adding a key to this config must never be a silent
  no-op. In `CRRT-dose-lmtp` this exact check caught `pf_source`: declared from the
  start, never built, and therefore absent from the frame while the config's own
  justification assumed it was present.

---

## 12. Source map

| Claim | Source |
|---|---|
| `scaling` formulas, within-unit normalisation | `references/gbmt_R.pdf`, `gbmt` Details |
| BIC formula `-2*logLik + npar*log(ss)`, lower is better | `gbmt:::icCalc` |
| Landmark / immortal time bias | Anderson, Cain & Gelber (1983) *JCO*; Dafni (2011) *Circ Cardiovasc Qual Outcomes* |
| `Jointlcmm` competing-risks support | `args(lcmm::Jointlcmm)`; `Jointlcmm.Rd` |
| `gbmt` has no covariate argument | `gbmt()` signature |
| APPA ≥ 0.7, OCC > 5 thresholds | Nagin (2005), *Group-based modeling of development* |
| Class enumeration criteria performance | Nylund, Asparouhov & Muthén (2007) |
| NEE coefficients (norepi 1.0, epi 1.0, phenylephrine 0.1, dopamine 0.01, vasopressin 2.5 per u/min, angiotensin 10.0) | Goradia et al., *J Crit Care* 2021;61:233–240, doi:10.1016/j.jcrc.2020.11.002, PMID 33220576. Same factors appear uncited in `CRRT-dose-lmtp/config/lmtp_design.json:193-198` |
| Uniform 24h lab LOCF cap | SG, 2026-09-05: clinicians act on the last available value, and ICU patients get minimum daily labs |
| NEE `hold_hours` = 4; vasopressor median inter-record interval 57 min, p90 68 min | `CRRT-dose-lmtp/docs/lmtp_df_build_notes.md:217-219` (CRRT coordinating site; unverified at UCMC) |
| NEE reached 8,001 mcg/kg/min-equivalent before the unit guard | `CRRT-dose-lmtp/code/02_build_lmtp_df.py:1138-1152` |
| Severinghaus over Rice; slopes 0.44/0.39/0.37 vs Rice's 0.84 | `CRRT-dose-lmtp/config/lmtp_design.json:420`; Brown PMID 26836924, validated PMID 28538439 |
| 85% of blocks missing P/F had S/F observed | `CRRT-dose-lmtp/config/lmtp_design.json:419` |
| FiO₂ lookback 4h; SpO₂ ceiling 97, strict | `lmtp_design.json:134`, `:155`; code `02_build_lmtp_df.py:933`, `:961` |
| clifpy waterfall wipes `fio2_set` across whole trach-collar encounters; 189/2,141 (8.8%) carry 50.4% of nulls | `CRRT-dose-lmtp/config/lmtp_design.json:375` |
| Uncapped lab LOCF: ~20% stale, median lag 27h, max ~1331h, differential 28.7% vs 13.7% | `crrt-manuscript-tools/.claude/lessons.md:581` |
| `nee` 17% / `inotrope` 84% "missing" = absence mislabelled as ignorance | `CRRT-dose-lmtp/config/lmtp_design.json:468` |
| clifpy SOFA defects (vasopressor set, non-IMV P/F NULL, FiO₂ 0.21–1 gate, creatinine-only renal) | clifpy 0.3.8 `utils/sofa.py:16-19`, `:172-176`, `:273`; verified in `.venv` 2026-09-05 |
| Weight lag when matched to an anchor: median 14h, p95 131h, max 585h | `CRRT-dose-lmtp/code/01_build_cohort.py:588-593` |
