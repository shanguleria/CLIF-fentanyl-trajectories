# Fentanyl administration on invasive mechanical ventilation — project principles

**Status:** current · **Date:** 2026-09-23 · **Author:** Shan Guleria

Definitions of terms live in `config/`, which is read by the code.

---

## 1. The question

Use CLIF's granular **medication, vitals, respiratory support and patient
assessment** data to characterize **in detail how fentanyl is administered to
patients on invasive mechanical ventilation** — how much, by what route, and how
that changes over the ventilation course.

The primary purpose is **descriptive**. 

## 2. The measurement frame

- **t = 0 is the start of invasive mechanical ventilation.** CLIF has no
  intubation event, only `device_category` transitions, so the anchor is the
  start of the first continuous IMV episode and is a **proxy** for intubation. A
  patient transferred in already ventilated has t = 0 at first observation
  instead. This is a limitation to state, not to fix.
- **4-hour windows over the first 72 hours** — 18 windows per episode.
- Within each window we record:
  - the **fentanyl infusion rate in mcg/hr** (LOCF between charted rate changes,
    then the time-weighted mean over the window — a rate persists until changed,
    so a plain mean of records is wrong when changes are unevenly spaced);
  - the **number and dose of fentanyl boluses in mcg** (summed within the window,
    **never** carried forward — a bolus is an event, not a state).
- Both streams are a true **`0`** in a window with no record. Drip off is a real
  zero, not a missing value.
- These support **cumulative** and **time-dependent** exposure, computed at
  analysis time rather than baked in.
- **The unit of analysis is the encounter block**, not the hospitalization: a
  patient transferred while ventilated has one ventilation course, and without
  stitching that course is truncated at the transfer.

Every parameter above has a single home in `config/`
(`window_hours`, `granular_extent_hours`, `landmark_hours`,
`imv_episode_gap_hours`, `min_imv_hours`) with its reason beside it. 

Width-sensitivity check on the state definitions still outstanding, but 4h is clinically plausible.

## 3. What this project is, and is not

**It is** a descriptive characterization of fentanyl delivery, delivered as the
five figures in §4.

**It is not, for now, a trajectory-classification study.** Group-based and
multi-trajectory modelling is intended **later work** and is deemed too
complicated to tackle now. `code/tabled/gbmt_classes.R` and
`code/tabled/lcmm_classes.R` remain in the repo and still run, but they are
tabled and are not in either runner. `code/tabled/landmark_cohort.R` joined
them on 2026-09-30 (§5). They carry no step number, and neither do their output
folders (`gbmt/`, `lcmm/`, `transitions/`, `landmark/`): the number line belongs
to the live pipeline. Their outputs on disk predate their own scripts —
do not cite a number from them without re-running. What was learned before
tabling: dose level is continuously distributed with no subpopulation gaps, so a
*latent* class claim was not supportable, while a *declared* band — which asserts
nothing and only labels — remains legitimate.

**It is not, for now, an outcome study.** RASS and NVPS assessments and a
covariate set for illness severity, comorbidity and age **are collected**, so
associating dosing practice with outcomes, or adjusting for covariates, is open
later. Nothing in the current analysis depends on it.

## 4. The deliverable — five figures

| | Figure | Status |
|---|---|---|
| F1 | One example IMV course: continuous fentanyl, boluses, time off fentanyl, with documented RASS and NVPS throughout. Modelled on Baker et al. Figure 1 — **one panel, three y scales** | **exists** — `exemplar` |
| F2 | States of 100 example patients as a per-patient raster, same states as the alluvial. Modelled on Iyer et al. Figure 1B; reads as panel A to F4 | **exists** — `state_raster` |
| F3 | Prevalence of states over time | **exists** — `state_prevalence` |
| F4 | Transitions between states | **exists** — `state_alluvial` |
| F5 | Point prevalence on an **at-risk denominator**: among patients still ventilated and alive in each window, the fraction receiving ≥1 bolus, a continuous infusion, both, or neither | **exists** — `state_prevalence_at_risk` |

States are **infusion only / bolus only / both / no fentanyl / extubated /
discharged alive / died**, defined by delivery **route** rather than by a
threshold, so there are no cut points to defend. `extubated` is **transient**,
not absorbing — patients are reintubated. A second, parallel definition cuts the
same windows into declared **intensity bands**; both run over one cohort to
triangulate. The definitions live once, in `code/utils/states.R`.

**Every figure in the pipeline is drawn journal-style** — no title or subtitle
on the panel, axis titles and a boxed legend at the foot. Each of scripts 02–05
writes its captions to `captions.md` beside its figures, with a guard that no
figure may ship without one. The single exception is the synthetic development
render of F1, which carries a "NOT A PATIENT" stamp on its face.

**How the F1 exemplar is selected** was the last open question here and is now
closed. The criteria are declared in `config/covariates.json` → `exemplar`,
applied by `01_build_cohort.py`, and the episode is **drawn at random** from
those that qualify. Baker's figure, which F1 is modelled on, is captioned only
"Representative ICU admission" and states no rule at all; a declared filter plus
a seeded draw is the one thing F1 set out to improve on it. If a drawn episode
reads badly the **rule** changes and the pipeline re-runs — picking a different
one from the eligible set by eye would reintroduce exactly the bias the rule
removes.

F5 is the figure that fixes a real dilution in F3: F3's denominator is every
episode at every window with terminal states carried forward, so the fentanyl
percentages shrink partly because patients leave rather than because practice
changes.

## 5. Landmark

A cohort defined by "still ventilated at 72h" cannot be assembled by looking
backwards from an outcome: conditioning on survival to a time point and then
measuring from admission credits the exposed group with time in which they could
not have had the event. That is **immortal time bias**, and the solution is to set
a **landmark** — fix the origin at a chosen T, require survival and ventilation
to T for entry, and measure only forward from there.

Everything before T is exposure history; everything after is
outcome. `code/tabled/landmark_cohort.R` implements it, and T-sensitivity is
reported alongside the primary because the choice of T is a judgement.

**TABLED 2026-09-30 (SG).** The landmark was always *available, not current*, and
its only consumer of `landmark_cohort.parquet` was the already-tabled gbmt/lcmm
pair — so a live pipeline step existed largely to feed tabled work. It moved to
`code/tabled/`, lost its step number, and `06_unit_variation.R` took the slot.
The reasoning above is unchanged and the script still runs by hand; per §3, do
not cite a number from it without re-running.

`time_to_event.parquet` is a **Phase 0** product and is untouched — steps 1 and 3
still read it, so the STROBE landmark row (episodes still ventilated at T) is
still produced. It now describes a cohort no live script analyses.

## 5b. Titration practice

A question adjacent to the five figures, added 2026-09-24: **when the fentanyl
infusion rate is increased, is a bolus given with it?** A rate change alone
approaches the new steady state over roughly four to five half-lives; a bolus at
the moment of uptitration gets there in minutes.

Descriptive for now — `code/05_titration.R`. Associations between adherence and
average rate, cumulative dose or time to extubation are deferred, and the
per-encounter adherence is emitted within a 24h landmark as well as overall so
that work needs no rebuild.

Three things make the number interpretable rather than decorative. **Rate
decreases are a negative control**: a bolus at a downtitration has no
pharmacologic rationale, so its pairing rate is the floor for coincidence.
**The window is varied from 5 to 60 minutes**, because a pairing that only
appears at wide windows is bolus frequency, not practice. And the analysis uses
**raw charted timestamps, never the hourly grid**, whose 1h binning would move a
rate change away from the bolus that accompanied it.

## 5c. How Table 1 is stratified

**Predominant fentanyl intensity** (SG, 2026-09-24), replacing landmark-eligible
vs not. That split was a leftover from the modelling design: "eligible vs not"
is "survived ventilated to 72h vs not", a severity contrast rather than a
delivery contrast, and its p-value invited reading as a finding.

Over an episode's at-risk time — every ventilated window — the patient is in
exactly one declared intensity band at all times, so one band holds more time
than the others. That band labels the episode. A **declared** rule, not a latent
one: the same ground on which the bands themselves survived after gbmt was
tabled. It asserts nothing about subpopulations existing; it only labels.

The same rule on delivery **route** was measured first and rejected — it left
93.5% of episodes in two classes and `continuous + bolus` with n = 124. The
intensity bands split 52 / 21 / 11 / 16%, four usable strata, and being
**ordered** they permit a trend test rather than an omnibus one.

**It is not a trajectory and must not be called one.** Modal time ignores order:
only 41.5% of episodes end in the band they began in
(`validation/modal_state_probe.py`), so "started high and
weaned" and "stayed medium throughout" can share a label. Two artifacts are
reported as Table 1 rows rather than argued away — modal share, which says how
decisive each label is, and at-risk hours, which expose a duration confound
(a short course cannot accumulate zero windows, so `low` skews short).

## 5d. Variation by unit and by calendar year

Applies to **both** titration metrics as of 2026-09-30: co-administration
(§5b) and documented indication (§5e). One machine, two metrics — the
attribution, the bootstrap, the suppression and the three figure forms are
written once and driven from a spec list, so a third metric is a config-shaped
addition rather than a copy. Figures are named by metric (`coadmin_*`,
`indication_*`) because `unit_caterpillar` did not say *a caterpillar of what*,
which was terse with one metric and ambiguous with two.

Added 2026-09-30 — `code/06_unit_variation.R`. §5b's pooled number hides two
axes: **which ICU** the patient was in when the rate change was charted, and
**which year** it was. The extract spans 2018–2024 at UCMC, so the pandemic years
are inside it rather than adjacent to it.

**This reverses a prior decision, narrowly.** On 2026-09-23 ICU location was
considered and declined, on the ground that "adding location introduces a new
exclusion that would have to be defined and defended." That was about location as
a **cohort filter**. Here it is a **stratification label**: it excludes nobody,
the cohort stays at 14,897 episodes, and `imv_status` remains the ICU proxy
wherever an at-risk denominator needs one. The 2026-09-23 decision still stands
for the use it actually addressed.

**Which care locations count.** The ICU is not the only place a ventilated
patient is. Measured 2026-09-30, of 64,037 escalation events: ICU 84.5%,
**procedural 8.0%**, **ED 6.9%**, ward 0.3%, uncovered 0.3%. The **ED is
attributed** — it is a ward-type setting where NVPS and RASS are meant to be
charted, and it behaves like one with a different culture: documented
indication 14.5% against the ICU's 32.2%, while bolus co-administration is
nearly identical (11.5% vs 13.4%). A charting difference, not a dosing one.
**Procedural is deliberately excluded**: under general anaesthesia nobody charts
an NVPS or a RASS, so its 7.8% indication rate measures the absence of the
instrument, not the absence of an indication, and a funnel would flag it as a
dramatic outlier on an artifact. Its co-administration of 21.8% — the highest
anywhere — is a real fact about anaesthesia practice and is left for a separate
question. Ward and l&d fall below the suppression floor on volume, not by rule.

**The unit label is `care_setting`** = `location_type` where CLIF populates it,
else `location_category`. `location_type` is non-null only for ICU rows, so the
coalesce yields the six mCIDE ICU types plus `ed`. Both inputs are required,
controlled CLIF columns, so the label set is reproducible elsewhere. Phase 0
emits the raw columns alongside the derived one rather than overwriting
`location_type` — silently redefining a CLIF column name is how the next reader
is misled.

**Attribution is per event, not per episode.** A titration event is a moment in
time and episodes move between units — 6.4% of attributed episodes do so inside
the 72h window (11.0% of whole ICU hospitalizations do, over their longer span;
the two denominators are not interchangeable) — so each rate change is credited
to the unit the patient was physically in when it happened.
Phase 0 emits `icu_intervals.parquet` on the same relative clock the events use,
which keeps the rule changeable for the cost of an R re-run. Events charted while
the patient was ventilated but not in an ICU location are an **unattributed row**,
never a silent drop — and that row is not small: **14.3% of increase events
(3,445 of 24,061) at UCMC** were charted off the unit, in an operating,
procedural or radiology location. The unit views run on the remaining 20,616;
`attribution_funnel.csv` carries the accounting and the script asserts it closes.

**Two unit keys, because they travel differently.** `location_type` is a
required, mCIDE-controlled CLIF column and is the only one another site can
reproduce or pool on; `location_name` is optional site-local free text. Both are
drawn; anything audience-facing uses the type view or anonymised labels.

**It is not a quality measure and must not be called one.** No evidence
establishes that pairing a bolus with an up-titration is better care, so a unit
with a lower rate is not worse. `tests/test_covariates.py` bans *quality*,
*performance* and *benchmark* from the script that builds the split — the same
guard, for the same reason, as §5c's.

**The intervals describe; they do not adjust.** A cluster bootstrap resamples
**episodes**, because §5b records that event-level percentages "carry no valid
confidence interval without accounting for clustering within episode." It says
how much a unit's observed rate would wobble if the care process re-ran. It does
**not** adjust for case mix, and a unit's rate reflects who it admits as much as
how it practises. The association between unit and the metric is deliberately
**not modelled**; `covariates.json` → `unit_variation._STATUS_model` records what
the deferred method is (Myers et al., PMID 39018285: hierarchical logistic
regression, hospital random intercept, conditional modes on the log-odds scale —
a shrunken, risk-adjusted estimand) and how its question differs from this one.
`docs/unit_adjustment_design.md` carries the DAG, the confounder/mediator split
and the model specification for whenever that analysis is taken up.

Both a caterpillar and a funnel view ship. They answer different questions — how
far apart the units are, and whether any sits outside what its volume would
predict — and the funnel is the one the methods literature prefers for
institutional comparison, because it puts volume on an axis instead of hiding it
inside an ordering.

**The caterpillar carries three series per unit** (added 2026-09-30): an
*initiation*, an *uptitration*, and *any increase*. The third is the two pooled,
not a third kind of event, so it necessarily lies between its components at
their event-weighted mean — say so in the caption, or a reader will look for a
finding in the fact that it does. Units are ordered by the combined estimate so
one ranking governs all three rows. The funnel and the year series stay on the
combined figure alone: the funnel's envelope reading depends on one point per
unit, and three overlapping funnels defeat it.

## 5e. Was the dose change indicated?

Added 2026-09-30 — the other half of §5b. Co-administration asks whether a bolus
**accompanied** an increase; this asks whether anything documented **prompted**
it. A qualifying score near the event is the trigger; the bolus is the response,
and neither alone describes the practice.

Dichotomised at each instrument's own threshold (SG, 2026-09-30): **NVPS ≥ 4**
and **CPOT ≥ 3** for pain, **RASS ≥ +1** for agitation. The scales are not
comparable — NVPS 0–10, CPOT 0–8, RASS −5..+4 — so raw values are never pooled;
cutting each at its own threshold is exactly what makes them combinable, because
after the cut every instrument contributes the same yes/no.

**It measures a LEVEL, not a change**, and the distinction is not cosmetic. A
score held above threshold across three readings qualifies at every one of them;
a rise from 0 to 3 qualifies at none. A level is the right reading of "was there
an indication" — a patient in pain is in pain whether or not it is new — but the
word *change* must not appear in the caption.

**The window is hours, not minutes.** Assessments are charted hourly to q4h
(measured at UCMC: NVPS median gap 1.0 h, RASS median 1.0 h with p90 12.1 h)
while titrations are charted to roughly the nearest 5 minutes. The ±30 min
window §5b uses would find almost nothing here, and the nothing would be an
artifact of cadence rather than a fact about practice. Primary ±1 h, swept to
4 h for the same reason §5b sweeps 5–60 min: a result that only appears at the
widest window is measuring assessment frequency.

**Two denominators, the same distinction §4 draws between F3 and F5.** Against
*all* events, the indicated fraction conflates "no score was taken" with "a
score was taken and did not justify it." Against events that *had* an
assessment in the window, it separates documentation from indication. Both
ship. A backward-only variant, restricted to scores at or before the event, is
the literal causal reading and ships as a third column; the gap between it and
the symmetric window says how much charting happens after the fact.

**It is not a quality measure.** The literature on bolus-with-titration is
pharmacokinetic rather than clinical (SG, 2026-09-30), so no evidence
establishes that an unindicated increase is wrong care. A low indicated fraction
may mean scores are not charted near the decision — which is why both
denominators are reported rather than one.

**An absent instrument is reported, not fatal.** CPOT is not charted at UCMC.
The instruments are read through a tolerant path that prints a zero-record line
and continues, deliberately *not* through `time_varying`, whose
`assessment_covariates()` raises on a category matching no rows. That behaviour
is right for a required covariate and would kill the build for an optional
site-specific one.

## 6. Open decisions

- **The tracheostomy rule.** Identification is specified; the rule is not. It
  blocks any competing-risks outcome work.
- **The RASS/NVPS LOCF cap**, currently 4h (one window). Looser bridges more
  windows but asserts a stale score still describes the patient.
- **Fentanyl or analgosedation?** Midazolam is 98% zero at the window level and
  dexmedetomidine 83%, so only propofol is charted often enough to carry a
  "sedation was swapped, not withdrawn" story. Hydromorphone, ketamine,
  remifentanil and morphine are not collected. A scope decision.
- **The weight lag cap.** 11% of blocks take the dose denominator from a weight
  charted after the anchor; median 7h, max ~804h.
- **Window-width sensitivity** for the state definitions.

## 7. Where things live

| | |
|---|---|
| Definitions, thresholds, aggregation rules, and their rationale | `config/covariates.json`, `config/outlier_config.json` |
| Site paths and local settings (gitignored) | `config/config.json` |
| Mechanics | `code/`, read standalone |
| The state definitions, once | `code/utils/states.R` |
| Numbers and figures | `output/final_no_phi/` |
| Literature | [`references.md`](references.md) |
| How to run it, and what a new site needs | [`../README.md`](../README.md) |
| Session log, todo, lessons | `.claude/` (gitignored) |
