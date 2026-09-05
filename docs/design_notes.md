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
| Indicators | **1**: total fentanyl mcg/kg/hr per 4h window | **2**: infusion rate + bolus-equivalent rate |
| `scaling` | `0` | `0` (required — see §4) |
| Status | **Primary.** Do this first. | **Secondary.** Do after A is settled. |

Model A is the cleaner and more defensible primary analysis. Model B answers a
genuinely different and interesting question, but it is the more fragile model
and should not be the headline result.

### Model A indicator

Combining infusion and push into one quantity resolves a units mismatch:
infusion rate is an *intensity* (mcg/kg/hr), cumulative bolus is a *quantity*
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

Keep the two streams separate, both expressed as mcg/kg/hr over the window:

- `inf_rate`  — time-weighted mean infusion rate
- `push_rate` — window bolus total / weight / 4

---

## 2. A correction to record

It is **not** GBTM that fails on zero-heavy dosing data. Model A — a single
combined-exposure indicator — is a perfectly appropriate GBTM application, and
zeros are harmless there once `scaling = 0` removes all division.

What fails is specifically a **near-all-zero second indicator under within-unit
normalisation** (Model B run at the package default `scaling = 2`). See §4/E2.

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
dose cannot define a group. A patient flat at 200 mcg/kg/hr and one flat at 50
become the *same trajectory* after within-patient normalisation. If
"persistently high-dose" is a phenotype we want to find, `scaling >= 1` destroys
it by construction.

`scaling = 3` and `4` are unusable for infusion data regardless — they require
strictly positive values, and the data has structural zeros whenever the drip is
off.

---

## 4. Evidence

All three experiments are reproducible via `code/03_scaling_experiments.R`
(output: `output/intermediate/scaling_experiments.csv`). ARI = adjusted Rand
index against the known simulated truth; 1.0 = perfect recovery, 0 = chance.

### E1 — what each scaling value erases

12 synthetic patients: 4 flat-high, 4 flat-low, 4 rising.

| Setting | ARI vs truth |
|---|---:|
| `scaling = 0` | **1.000** |
| `scaling = 1` | 0.542 |
| `scaling = 2` | 0.505 |

(source: `code/03_scaling_experiments.R` §E1; `scaling_experiments.csv`
rows `E1_what_scaling_erases`)

Only `scaling = 0` recovers the truth. Under `scaling = 2` the flat-high and
flat-low patients merge, exactly as the algebra predicts.

### E2 — a near-all-zero second indicator

12 patients, second indicator (`push`) nonzero for only 2 of them.

| Setting | ARI vs truth |
|---|---:|
| `scaling = 0` | **1.000** |
| `scaling = 2` | **−0.053** |

(source: `code/03_scaling_experiments.R` §E2; `scaling_experiments.csv`
rows `E2_zero_inflated_indicator`)

A negative ARI is *worse than chance*. The mechanism: a patient who never
received a bolus has a within-patient SD of exactly 0 for that column, and
`scaling = 2` divides by it. **The run produced no error and no warning** — it
returned a confidently wrong partition. This is the single most dangerous
finding in this document.

### E3 — is `scaling = 0` sensitive to the units of each indicator?

Theory says no: a Gaussian mixture with a freely estimated full covariance per
group absorbs any linear rescaling of a variable. But `gbmt` initialises EM from
a Ward hierarchical clustering, which **is** scale-sensitive, so the two runs can
converge to different local optima.

Tested by multiplying `push` by 10 and re-fitting, over 6 seeds, on deliberately
poorly separated data:

| Seed | ARI(base, push × 10) |
|---|---:|
| 1, 2, 4, 5 | 1.000 |
| 3 | 0.485 |
| 6 | 0.606 |

Mean 0.849; the partition changed in **2 of 6 runs**.
(source: `code/03_scaling_experiments.R` §E3; `scaling_experiments.csv`
rows `E3_unit_sensitivity_scaling0`)

**Implication:** `scaling = 0` is scale-invariant in principle but not reliably
in practice on realistic (overlapping) data. Mitigations for Model B:

1. Set `nstart = 50` or higher so EM restarts randomly rather than relying on
   the Ward solution.
2. Put both indicators on similar numeric ranges before fitting, even though
   theory says it should not matter.
3. Report assignment stability across restarts and across a units change, using
   the ARI helper in `code/02_indicator_sensitivity.R`.

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
| Weight | denominator of the dose (mcg/kg/hr) | weight is a normaliser, not a covariate |

**Verify before relying on LOCF:** check what fraction of infusion stops are
charted as an explicit `0` row versus an end timestamp with nothing after. The
LOCF approach is only valid if stops are explicit. *(open — see §9)*

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
| Trajectory indicator | fentanyl dose (± propofol, midazolam, RASS) | `x.names` |
| Normaliser | weight | data prep — mcg/kg/hr |
| Membership predictor | age, admission SOFA, race, admission type | stage 2: `multinom(group ~ ...)` |
| Distal outcome | vent-free days, delirium, mortality | stage 2: `glm(outcome ~ group + ...)` |
| Group descriptor | bolus counts, % windows with a bolus | descriptive table only — characterises groups without letting them drive the grouping |

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

BIC magnitudes are **not comparable across `scaling` settings** — the same
`ng = 3` model on identical data ranged from **+18,491** (`scaling = 0`) to
**−3,950** (`scaling = 3`). *(source: session run over `agrisus2`, `scaling` 0–4,
`d = 2`, `ng = 3`; not yet scripted — see §9)*

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
> **conditional on being alive and mechanically ventilated at T.**

The conditioning is part of the estimand, not a limitation for the discussion
section. It belongs in the Methods, the abstract, and arguably the title.

### Why not simply zero-fill (evidence)

| Alternative | Why it fails | Source |
|---|---|---|
| Carry `dose = 0` forward after extubation | Outcome becomes embedded in the exposure; also often factually wrong, since many patients receive fentanyl post-extubation for pain | §2, §8 above |
| Keep short trajectories, unbalanced panel | `gbmt` **silently caps** the polynomial degree at (shortest unit's windows − 1). Requesting `d = 4` with one 4-window patient yields effective `d = 3`, by warning only | `code/03_scaling_experiments.R` §E4 |

E4 result: requested `d = 2` → effective `d = 2`; requested `d = 4` → effective
`d = 3`, warning `'d' was set to the maximum feasible value: 3`. One briefly
ventilated patient therefore constrains the trajectory shape available to the
entire cohort. (source: `output/intermediate/scaling_experiments.csv`, rows
`E4_unbalanced_panel_degree_cap`)

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

- [ ] Confirm infusion stops are charted as explicit `0` rows in our data (§5).
- [ ] Decide the time anchor: intubation vs ICU admission.
- [ ] Decide the truncation window, and quantify how many patients it excludes —
      informative censoring (death, early extubation) is the largest validity
      threat in this design.
- [ ] Check the zero fraction of the Model A combined indicator once built. Heavy
      zero-inflation may argue for `crimCV`'s zero-inflated Poisson instead.
- [ ] Script the `scaling` × BIC comparison so §7's numbers are reproducible.
- [ ] Decide whether propofol / midazolam join as additional indicators — note
      that mixing units would force `scaling >= 2`, which §3 says destroys the
      level information we care about. Fentanyl-equivalents may be the escape.
- [ ] Fix the **landmark time** for outcome ascertainment (§10 Phase 4). Measuring
      trajectories over 0–72h and then counting extubations from hour 0 is
      immortal time bias — a patient cannot be in the cohort unless they stayed
      intubated long enough to have a trajectory.
- [ ] Confirm whether successful extubation is itself a competing-risks problem
      (extubation vs death vs tracheostomy) or can be treated as binary.
- [ ] Write the per-variable within-window aggregation rules (§11). Each
      time-varying covariate needs its own rule; there is no sensible default.
- [ ] **Resolve §10(a)**: choose grid + landmark (3h/[0,48]/T=48 or 12h/[0,72]/T=72).
- [ ] **Resolve §10(c)**: confirm "cumulative dose" means within-window total,
      not a running sum since intubation.
- [ ] Measure the zero fraction of `total_dose` per window before Phase 2. If
      high, `crimCV`'s zero-inflated Poisson may fit better than `gbmt`.
- [ ] Decide whether propofol / midazolam get their own dose columns in Table 1
      now (cheap to add, expensive to backfill), even if unused until later.

---

## 10. Analysis roadmap

Revised 2026-09-05 per SG. Seven gated phases; each produces something
reviewable before the next begins.

### Phase 0 — data structures

**Table 1 — trajectory (long).** One row per patient per window.

| Block | Contents |
|---|---|
| Keys | `patient_id` (chr), `id_num` (int), `window_idx`, `window_start_hr` |
| Fentanyl dose — **three columns, kept separate** | `inf_dose`, `bolus_dose`, `total_dose` (within-window totals, expressed mcg/kg/hr) |
| Other sedatives | `propofol_dose`, `midazolam_dose` — added now; cheap here, expensive to backfill |
| Time-invariant covariates | age, sex, race, admission type, admission SOFA |
| Time-varying covariates | NEE, P/F, labs — one aggregation rule each (§11) |
| **IMV status** | on/off ventilator in the window — required to identify extubated windows and count failed extubations |

Keeping all three dose columns defers the Model A / Model B choice to analysis
time rather than baking it into the pipeline.

**Table 2 — time-to-event.** One row per patient. The two planned outcome
analyses need **different event codings**, so a single `event`/`time` pair will
not serve both:

| Analysis | Time origin | Event coding |
|---|---|---|
| Extubation competing with death | landmark T | 0 = censored, 1 = extubation, 2 = death |
| 30-day mortality competing with discharge alive | landmark T | 0 = censored, 1 = death, 2 = discharge alive |

Carry both pairs (or a tidy long form), plus the landmark eligibility flag.

### Phase 1 — descriptive cohort dose trajectory
**No landmark, no exposure window, no outcome model.** All intubated patients,
contributing for as long as they remain ventilated.

| View | Window | Extent | Points |
|---|---|---|---|
| Granular | **4h** | 72h | 18 |
| Extended | 12h | 7 days | 14 |

Median/IQR (lead with these — dose is right-skewed), mean/SD alongside, and a
**number-at-risk row** (at-risk = still intubated). Stopping fentanyl is not
exclusionary: `dose = 0` is a real observation for a ventilated patient.

**Overlay balanced-panel curves** at ≥24h, ≥48h, ≥72h. The all-at-risk curve
answers a different question at every timepoint because its denominator keeps
changing; a balanced panel freezes the denominator so movement reflects real
within-patient change. Crossing = composition; parallel = real change. Worked
example and figure: `code/04_composition_bias_demo.R` →
`output/intermediate/composition_bias_demo.png`.

**Also plot three curves, not one:** (1) % receiving any fentanyl, (2) median
across all at-risk, (3) median among those receiving any. Weaning-to-zero among
the still-ventilated and dropout of low-dose patients push the overall median in
opposite directions; a single curve hides both.

**Denominator must be stated in Methods.** Same simulated data, same constant
patient doses: at-risk denominator **+50%**, full-cohort-with-zeros denominator
**−95%**. The latter mostly measures extubation rate, not dosing. Use at-risk.

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
| `x.names` | `c("inf_dose", "bolus_dose")` | both in mcg/kg/hr |
| `scaling` | **`0` — mandatory** | §4/E2: `scaling = 2` gave ARI **−0.053**, worse than chance, with no error |
| `nstart` | **≥ 50** | §4/E3: partition changed in 2 of 6 runs from Ward-init local optima |
| `d` | 2 (3 if windows allow) | balanced panel, so no degree cap (§8/E4) |

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
- ~~"Cumulative dose" ambiguity~~ — **within-window total**, expressed mcg/kg/hr.
- ~~Window width~~ — **4h**.
- ~~Other sedatives~~ — propofol and midazolam columns added to Table 1 now.

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

**(a) Dose in the extubated gap — DECIDED: `0`.** (2026-09-05, SG.) A patient
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

Every time-varying variable needs an explicit rule. There is no safe default,
and the right rule differs by what the variable *is*. Rules apply identically to
whichever window width is chosen (§10a). To be filled in and reviewed before
Phase 1 is built:

| Variable | Type | Proposed rule | Rationale |
|---|---|---|---|
| Fentanyl infusion rate | state (rate) | LOCF → time-weighted mean | persists until changed |
| Fentanyl boluses | event | **sum**, never LOCF | discrete events, not a state |
| NEE (norepinephrine equivalents) | state (rate) | LOCF → time-weighted mean | same logic as fentanyl drip |
| P/F ratio | intermittent measurement | *worst* or *nearest-to-window-end* — **decide** | driven by ABG timing, not a continuous state |
| Labs | intermittent measurement | *last* or *worst* in window — **decide** | irregular sampling; "mean" is rarely meaningful |
| RASS | intermittent assessment | *modal* or *worst* — **decide** | ordinal; a mean of ordinal scores is not interpretable |

The state-versus-event-versus-measurement distinction is the one that matters:
apply a state rule to an event column and you overcount; apply an event rule to
a state column and you undercount.

---

## 12. Source map

| Claim | Source |
|---|---|
| `scaling` formulas, within-unit normalisation | `references/gbmt_R.pdf`, `gbmt` Details |
| BIC formula `-2*logLik + npar*log(ss)`, lower is better | `gbmt:::icCalc` |
| E4 polynomial-degree cap on unbalanced panels | `code/03_scaling_experiments.R` §E4 |
| Landmark / immortal time bias | Anderson, Cain & Gelber (1983) *JCO*; Dafni (2011) *Circ Cardiovasc Qual Outcomes* |
| `Jointlcmm` competing-risks support | `args(lcmm::Jointlcmm)`; `Jointlcmm.Rd` |
| E1 / E2 / E3 ARI values | `code/03_scaling_experiments.R` → `output/intermediate/scaling_experiments.csv` |
| ARI method and implementation | `code/02_indicator_sensitivity.R` |
| Indicator-set sensitivity (ARI 0.16–0.41 on `agrisus2`) | `output/intermediate/indicator_sensitivity_ari.csv` |
| BIC by `ng`, scree plot | `code/01_gbmt_example.R` → `output/intermediate/gbmt_ic_comparison.csv` |
| `gbmt` has no covariate argument | `gbmt()` signature |
| APPA ≥ 0.7, OCC > 5 thresholds | Nagin (2005), *Group-based modeling of development* |
| Class enumeration criteria performance | Nylund, Asparouhov & Muthén (2007) |
