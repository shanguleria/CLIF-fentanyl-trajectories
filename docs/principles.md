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
complicated to tackle now. `code/04_gbmt_classes.R` and
`code/05_lcmm_classes.R` remain in the repo and still run, but they are tabled
and are not in either runner. Their outputs on disk predate their own scripts —
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

**All five figures are drawn journal-style** — no title or subtitle on the
panel. Each run writes its captions to `captions.md` beside the figures, with a
guard that no figure may ship without one.

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
outcome. `code/04_landmark_cohort.R` implements it, and T-sensitivity is reported
alongside the primary because the choice of T is a judgement.

The landmark is **available, not current**: the descriptive account comes first,
and returning to it is deliberate future work.

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
