# References

Literature this project's decisions rest on, grouped by what it is used for.

`references/` holds local reading copies and is **gitignored in full** — publisher
PDFs must not enter the history of a repo that ships to other sites. So every
entry carries a **DOI**, which is what a reader without the folder has; the local
filename is given where one exists, as a convenience, and nothing in the pipeline
may depend on it.

**Verification levels.** ✓✓ = read from the PDF itself. ✓ = verified against
PubMed/Crossref. "as cited in the design notes" = carried forward from an inline
citation in the retired `design_notes.md` and **not** re-verified — treat those as
leads, not as checked citations.

*A caution this file exists because of.* On 2026-09-23 this section was first
assembled from web searches while the PDFs sat unread in `references/`; it
declared two of those papers unfindable and a third unreachable, and 9 of 10
local PDFs were cited nowhere. Read the local copy before searching, and never
write a negative literature finding on the strength of a web search.

## Trajectory and latent-class methods

- **Proust-Lima C, Philipps V, Liquet B.** Estimation of extended mixed models
  using latent classes and latent processes: the R package `lcmm`. *J Stat Softw*
  2017;78(2):1–56. [doi:10.18637/jss.v078.i02](https://doi.org/10.18637/jss.v078.i02) ✓
  — the method citation for the latent-class arm. Cite **this**, not the manual, for the model.
- **Proust-Lima C, Philipps V, Diakite A, Liquet B.** *lcmm: Extended Mixed
  Models Using Latent Classes and Latent Processes.* R package v**2.2.2**
  (2025-11-20); manual generated 2026-05-08. `references/lcmm_R.pdf` ✓✓
  — the version pinned in `renv.lock` **and** the version actually loaded by the
  run: `logs/05_lcmm_classes_sessioninfo.txt` reads `lcmm_2.2.2`. The manual is
  the right edition for the numbers on disk. Source of the `gridsearch()` /
  multiple-start remedy that this project's `nstart >= 50` decision rests on.
- **`gbmt` package documentation** — `references/gbmt_R.pdf`; `scaling` formulas
  and the BIC definition `-2*logLik + npar*log(ss)` (`gbmt:::icCalc`). Used
  in the tabled `gbmt` arm. `renv.lock` pins **gbmt 0.1.4**.
- **Jones BL, Nagin DS.** A note on a Stata plugin for estimating group-based
  trajectory models. *Sociol Methods Res* 2013;42(4):608–613.
  [doi:10.1177/0049124113503141](https://doi.org/10.1177/0049124113503141) ✓
  — the `traj` plug-in: Nagin-style LCGA, the family `gbmt` is the R analogue of
  and the one Yang et al. used.

## Why our latent-class fits did not converge — the two papers SG raised *(2026-09-23)*

Both were raised to ask why their latent-class fits converged where our
`hlme` returned `conv = 4` at ng = 3, 4 and 5. **Neither supports reading our
non-convergence as a property of `lcmm`.** One of them used a different model
family entirely; the other used *the same package* and converged at every class
count, because of two design choices we did not make.

- **Xiao S, Zhuang Q, Li Y, Xue Z.** Longitudinal vasoactive inotrope score
  trajectories and their prognostic significance in critically ill sepsis
  patients: a retrospective cohort analysis. *Clin Ther* 2024;46(9):711–716.
  [doi:10.1016/j.clinthera.2024.07.006](https://doi.org/10.1016/j.clinthera.2024.07.006);
  PMID 39153910. `references/Xiao_lcmm_vasoactives.pdf` ✓✓

  MIMIC-IV v2.2, Sepsis-3, **n = 6,802**, VIS computed **bi-hourly over 72 h**
  (36 occasions/subject ≈ 244,872 observations), **R `lcmm`** under R 4.1.3,
  **4 classes** (52.1 / 23.0 / 13.2 / 11.6%). Same package as ours. `conv = 1`
  at **every** G from 1 to 5, and **no starting-value strategy is reported** — no
  `gridsearch`, no random starts. So the package is not the difference and
  neither is `nstart`. Two things are:

  **(1) The zero mass was removed by cohort construction.** Receipt of vasoactive
  therapy within 24 h *lasting at least six hours* was an **inclusion criterion**:
  of 33,411 screened, **4,842 were excluded for no vasoactive in the first 24 h**
  and 104 more for < 6 h. Minimum baseline VIS in the analytic cohort is
  **0.06** — nobody starts at zero. Our panel is **44.8% exact zeros** with no
  exposure-based entry filter.

  **(2) The outcome was banded, and they say why.** Verbatim: *"Due to the skewed
  distribution and presence of extreme values in VIS, it was transformed into a
  scale with six ranges based on five clinically meaningful thresholds: 0, 1–5,
  6–10, 11–20, 21–40, and greater than 40."* That 1–6 integer code was then
  modelled as continuous. VIS spanned 0.06–275, skew comparable to our 0–1383.
  Banding bounds the range, flattens the tail, and — the part that matters —
  **closes the gap between the zero spike and the rest of the distribution**:
  "0" becomes a band adjacent to "1–5" rather than a point mass separated by
  empty space. Their lowest class asymptotes at ≈ 2.40 on that 1–6 scale, 1.4
  units clear of the floor; **no class sits against a boundary.**

  *Specification, reconstructed from `npm` (the paper reports none of it):* `npm`
  runs 4, 8, 12, 16, 20 for G = 1…5 — exactly +4 per class with nothing left over,
  which gives class-specific **quadratic** time, one **shared** residual SD, and
  **no random effects at all**. Their Methods says "LGMM" while the Discussion
  says "Group-Based Trajectory Modeling"; the parameter count sides with GBTM. So
  they also did not carry the random intercept we did — a second, independent
  lever.

  **Before citing their model selection, know that it does not say what it
  claims.** BIC falls monotonically through G = 5 (403,998.6 → 342,476.2 →
  320,790.6 → 310,763.4 → **301,246.7**) and entropy falls monotonically
  (1.000 → 0.954 → 0.928 → 0.903 → 0.884), yet the text reads *"the entropy value
  showed a turning point at the four-class model"*. There is no turning point in
  the printed table; the 4-class choice rests on a BIC elbow plus clinical
  interpretability. **No APPA, no odds of correct classification, and no
  minimum-class-size rule are reported.** Contrast this project's own class-enumeration criteria, which required all three.

  Their acknowledgments thank Zhongheng Zhang for source code deposited with
  `doi:10.1186/s13054-020-2768-z` — i.e. **the model spec is inherited from Zhang
  et al. below**, not designed for this paper.

- **Yang J, Liujiao Y, Zhang X, Xiong J, Wang F, Shen F.** High NE dose
  trajectory is associated with new onset of acute kidney injury patients: a
  group-based trajectory modeling analysis. *PLoS One* 2025;20(5):e0323431.
  [doi:10.1371/journal.pone.0323431](https://doi.org/10.1371/journal.pone.0323431);
  PMC12074548 ✓ — **no local PDF; sourced from PMC.**

  MIMIC-IV, n = 3,462 septic shock, NE-equivalent dose at **8 points on a fixed
  12-hourly grid over 96 h**, complete by construction (≥ 96 h LOS required for
  entry), **3 classes** (47.3 / 41.5 / 11.2%). **Not an LCMM.** Verbatim: *"We
  utilized the traj plug-in in STATA to perform GBTM for estimating NE dose
  trajectories."* Stata `traj` is Nagin-style LCGA with **no random effects**, and
  their lowest class sits at **exactly 0.000 µg/kg/min from 48 to 96 h** — a
  heavily zero-inflated outcome fitted successfully *because* the model has no
  within-class variance component to collapse. Reporting gaps before using it as
  a template: the `traj` distributional model (`cnorm`/`zip`/`logit`) and
  censoring bound are never stated, the winning polynomial order is never
  reported, and **no convergence or starting-value information appears at all.**

- **Zhang Z, Ho KM, Gu H, Hong Y, Yu Y.** Defining persistent critical illness
  based on growth trajectories in patients with sepsis. *Crit Care* 2020;24:57.
  [doi:10.1186/s13054-020-2768-z](https://doi.org/10.1186/s13054-020-2768-z).
  `references/Zhang_lcmm_sepsis.pdf` ✓✓
  — eICU, **n = 22,868**, daily SOFA, latent growth mixture modelling, **5
  classes**; persistent critical illness transition at day 15; 643 patients (2.8%)
  developed PCI and consumed 19% of ICU bed-days. **The methodological parent of
  Xiao et al.** — its published code is what Xiao ran. Its residual value after
  the pivot is the transition-time framing: it turns trajectory classes into a
  statement about *when* patients change state, which is adjacent to the
  delivery-state approach.

**What this means if the latent-class arm is ever untabled.** The fix is already
in this repo. Xiao's move — recode a skewed, zero-heavy dose into clinically
defined **absolute** bands and model the band index — is exactly
`exposure.dose_states` (`window_mcg` cut at 200 and 400 mcg per 4 h window;
shares 50.8 / 21.3 / 12.5 / 15.5%), added in `c1c4d29` for the intensity-band
state definition. Feed **that**, not raw `total_dose`. Four levels over 18
occasions is coarser than their six over 36, so it would likely want more cuts —
and those cuts are federation-critical, so widening them is a consortium
decision (see `covariates.json exposure.dose_states._cuts_MUST_BE_ABSOLUTE`).
Dropping the random intercept is the independent second lever.

## Descriptive-figure sources

Both figures SG named by shorthand are here, and **both are from outside the
sedation literature** — which is why searching on topic found neither. The
lesson is recorded in `.claude/lessons.md`: a shorthand names a figure's *form*,
not its subject.

- **Baker L, Maley JH, Arévalo A, DeMichele F 3rd, Mateo-Collado R, Finkelstein S,
  Celi LA.** Real-world characterization of blood glucose control and insulin use
  in the intensive care unit. *Sci Rep* 2020;10:10718.
  [doi:10.1038/s41598-020-67864-z](https://doi.org/10.1038/s41598-020-67864-z).
  `references/Baker_insulin_single_trajectory.pdf` ✓✓
  — **the model for F1, the single-patient exemplar.** MIMIC-III v1.4, n = 19,694
  ICU admissions. Figure 1, captioned only *"Representative ICU admission"*, is
  the paper's sole patient-level figure and exists as a proof-of-granularity
  exhibit before any aggregate.

  *Construction:* one panel, **three y-axes** — glucose (left), insulin
  "units or units/hour" (right), D50 mL (a second right spine offset ≈32 pt
  outward). x = **`Time since entering ICU (days)`**, relative, 0–27 d.
  Glucose is a red straight-line-joined polyline with a filled circle at every
  measurement. The three insulin classes are **filled diamonds, all exactly the
  same size**, so dose is encoded by **vertical position only** and class by
  colour. The infusion is a true **step function with square risers, unfilled,
  closed to zero at both ends** — so "off drug" is the absence of the trace at
  baseline. D50 boluses are pentagons on the third axis. Legend sits
  **inside** the panel, semi-transparent, and visibly occludes data at days 5–8.
  Built in **matplotlib** (Python 3.7), notebook in the MIMIC Code Repository;
  a real de-identified MIMIC patient, with relative time partly because MIMIC
  date-shifts.

  *What F1 takes:* the relative clock anchored to the exposure window (for us,
  hours since intubation); the step-with-square-risers for mcg/hr closed to zero
  at both ends; instantaneous boluses as a reserved marker shape with dose read
  off position; and the overplotted "staircase" of near-simultaneous boluses,
  which is real signal about how hard a patient was being chased — do not jitter
  it away.

  *What F1 must not take:* **the three overlaid y-axes.** Baker got away with it
  because glucose, insulin units and dextrose mL are all 0–50 numbers on a 0–400
  backdrop; our quantities do not cohere (rate 0–400 mcg/hr, boluses 25–100 mcg,
  RASS **signed** −5…+4, NVPS 0–10). F1 uses **stacked facets on one shared x
  axis** instead. Also: the shared `"units or units/hour"` axis is a unit error
  papered over with an "or" — a rate and an amount must not share a numeric axis,
  so bolus dose gets its own spike panel or a size/label encoding. And Baker
  draws RASS-equivalents the wrong way for us: **straight-line interpolation
  between ordinal assessments asserts the patient passed through intermediate
  values at specific times.** RASS and NVPS get points plus
  `geom_step(direction = "hv")` — last documented value carried forward, the same
  semantics as the infusion step — on integer breaks, with a reference rule at
  RASS 0. Because Baker shows missingness as nothing at all, F1 must **break the
  step across gaps** beyond a stated threshold rather than assert a score
  persisted; the 4 h LOCF cap on `rass`/`nvps` is the same decision in the panel.
  Finally, *"Representative ICU admission"* with **no stated selection rule** is
  the figure's weakest point and the one thing F1 must improve on: state how the
  exemplar was chosen, and state the de-identification.

- **Iyer S, Kennedy JN, Nauka PC, Senussi MH, Seymour CW.** Epidemiology of
  β-blocker use among critically ill patients during and after septic shock.
  *Crit Care* 2024;28:364.
  [doi:10.1186/s13054-024-05145-1](https://doi.org/10.1186/s13054-024-05145-1).
  `references/Iyer_Bblocker_trajectory.pdf` ✓✓
  — **the model for F2, the per-patient state raster.** A 3-page Correspondence:
  12 UPMC hospitals 2010–2014, septic shock among patients on chronic outpatient
  β-blockers, n = 3,748 of 22,208; **Stata 18.0**.

  **Figure 1A is a table, not a plot, and there is no alluvial anywhere in the
  paper.** The alluvial-plus-raster pairing is ours to invent; only the raster
  half has a precedent here. The three prescribing strata that 1A tabulates are
  **not encoded in 1B at all.**

  *Panel 1B:* one row per patient, **"100 randomly selected patients"** (stated in
  the caption — the only methodological sentence about the figure). x = **Hospital
  Day 1…14**, every integer ticked and labelled, one tile per calendar day.
  **No y axis whatever** — no title, no ticks, no labels, no identifiers: the
  de-identification is structural. A continuous raster with **no tile borders**.
  Four qualitative colours — β-blockers `#03349B`, vasopressors `#30A6FE`, both
  `#FA0405`, neither `#FCD647` — with **white for out-of-observation time, absent
  from the legend.** The legend is a single outlined four-block colour bar
  outside the panel, no title.

  *Row ordering — not stated in the paper, reverse-engineered from the figure:*
  **non-increasing total number of days on vasopressors** (vasopressor-only plus
  "both" days), descending, **with ties left unbroken.** Not length of stay, not
  first state, not leading run length, and no clustering. The unbroken ties are
  why the panel looks organic rather than combed.

  *What F2 takes:* the form itself; no y axis at all; the four-sided panel border
  with no grid; `geom_raster` over `geom_tile` with no border colour; and stating
  the sample and its draw in the caption.

  *What F2 must change.* **Row ordering:** a single count-sort cannot carry seven
  states, because ours are two different families — a 4-level fentanyl exposure
  set and a 3-level disposition set — and sorting on "any fentanyl windows" would
  scatter deaths and discharges down the panel. F2 uses a **hierarchical key**:
  terminal state at 72 h, then time to first terminal event (descending), then
  total fentanyl-exposed windows as the tie-breaker — Iyer's key, demoted to
  where it belongs — and ties are broken rather than left arbitrary.
  **Terminal states:** Iyer lets rows stop and the page show through, so white
  does double duty for discharged, died *and* censored at day 14. We must not
  copy that: our three exits are **named states and are the point of the
  figure**, `extubated` is explicitly **not** absorbing (patients are
  reintubated), and a ragged right edge would make the raster and the alluvial
  disagree visually even where they agree numerically. Each terminal state gets
  its own colour, **carried forward to the 72 h horizon** — which is what
  `trajectory_long` already does — leaving the panel a full rectangle that can be
  read column-wise, and white reserved for genuinely unobserved time (of which
  there should then be none). **Palette:** four colours do not extend to seven;
  take the *logic* (family in hue, level in lightness) and check it in greyscale
  and for deuteranopia. The raster and the alluvial must share one palette and
  one factor-level order.

  *One caveat worth inheriting as a warning:* the cohort requires vasopressors
  within 24 h, yet 10 of the 100 rows begin in "neither" — so "hospital day 1" is
  a calendar-day bin anchored to admission, not to shock onset. Our 4 h windows
  are anchored to first IMV, and the anchor must be stated explicitly.

- **Swihart BJ, Caffo B, James BD, Strand M, Schwartz BS, Punjabi NM.** Lasagna
  plots: a saucy alternative to spaghetti plots. *Epidemiology*
  2010;21(5):621–625.
  [doi:10.1097/EDE.0b013e3181e5b06a](https://doi.org/10.1097/EDE.0b013e3181e5b06a);
  PMC2937254 ✓ — **no local PDF; sourced from PMC.**
  — the **methods** citation for F2's construction, complementing Iyer's worked
  example: the `m × n` history matrix, five named row sortings (**entire-row**
  sorting is the one preserving both patient identity and time order), the rule
  that absent time takes the *background* colour rather than a palette colour,
  and HCL-space palette guidance. It also states the precondition our 4 h grid
  already satisfies: a lasagna plot needs a discretised, shared time axis.
- **Barker AK, Valley TS, Kenes MT, Sjoding MW.** Early deep sedation practices
  worsened during the pandemic among adult patients without COVID-19. *Chest*
  2024;165(6):1429–1438.
  [doi:10.1016/j.chest.2024.01.019](https://doi.org/10.1016/j.chest.2024.01.019);
  PMC11317814 ✓ — **no local PDF.**
  — the closest published precedent for **a windowed ordinal sedation score as a
  state**: an alluvial over the first 48 h whose state is a 12 h weighted-average
  RASS binned deep (−3/−4/−5) vs. target (−2 to +1). The paper
  `covariates.json time_varying.rass._MEAN_DELIBERATELY_NOT_COLLECTED` points at
  as the likely eventual need for a central-tendency RASS.

## ICU sedation and analgesia — clinical comparators

- **Aragón RE, Proaño A, Mongilardi N, et al.** Sedation practices and clinical
  outcomes in mechanically ventilated patients in a prospective multicenter
  cohort. *Crit Care* 2019;23:130.
  [doi:10.1186/s13054-019-2394-9](https://doi.org/10.1186/s13054-019-2394-9).
  `references/Aragon_sedation.pdf` ✓✓
  — **the closest comparator for the descriptive paper.** Prospective, 5 ICUs in
  4 public hospitals in Lima, n = 1,657 ventilated adults, **1,338 (81%) with
  RASS scores over 18,645 ICU days**, followed to day 90. Reports **deep sedation
  in 98% of participants at some point**, cumulative benzodiazepine 774.5 mg and
  opioid 16.8 g, and 41% higher mortality at the 75th vs 25th percentile of
  benzodiazepine dose. Regression-based, not latent-class. **This is the
  benchmark for our deep-sedation figure** (41.3% of ventilated windows
  carrying a RASS are ≤ −3 at UCMC, 2026-09-23) — note their 98% is a
  *patient-ever* denominator and ours is a *window* denominator; the two are not
  comparable without restating one of them.
- **Su L, Liu C, Chang F, et al.** Selection strategy for sedation depth in
  critically ill patients on mechanical ventilation. *BMC Med Inform Decis Mak*
  2021;21(Suppl 2):79.
  [doi:10.1186/s12911-021-01452-7](https://doi.org/10.1186/s12911-021-01452-7).
  `references/Sedation LCA Paper.pdf` ✓✓
  — **despite the filename this is latent *profile* analysis**, cross-sectional
  over 36 characteristic variables in MIMIC-III, then PCA to 9 — **not** a
  longitudinal trajectory model, so it is a weak methods precedent for
  `gbmt`/`lcmm`. Two phenotypes. Its real value is as a second instance of the
  critique our own latent-class work ran into: two "phenotypes" that largely track a continuous
  severity and ventilator-load gradient.
- **Myers LC, Bosch NA, … Walkey AJ.** Opioid administration practice patterns in
  patients with acute respiratory failure who undergo invasive mechanical
  ventilation. *Crit Care Explor* 2024;6(7):e1123; PMC11257673 ✓ — **no local
  PDF.**
  — nearest epidemiology comparator: 21 ICUs, infusion *and* bolus opioid in
  morphine-milligram equivalents. Cite for context, not figure design (its
  figures are hospital-level caterpillar plots).

## Multistate and transition models

- **Lyons PG, Mody A, Bewley AF, et al.** Multistate modeling of clinical
  trajectories and outcomes in the ICU: a proof-of-concept evaluation of AKI
  among critically ill patients with COVID-19. *Crit Care Explor* 2022;4(12):e0784.
  [doi:10.1097/CCE.0000000000000784](https://doi.org/10.1097/CCE.0000000000000784).
  `references/lyons_multistate.pdf`, `references/lyons_multistate_supplement.pdf`
  — two multistate models on one cohort (AKI stage, then AKI × IMV) as
  triangulation: the pattern `code/utils/states.R` follows for the
  route/intensity pair, and the source of the transition-hazard curve style.
  *(As cited in the design notes.)*
- **Lyons PG, Bhavani SV, … Sinha P.** Hospital trajectories and early predictors
  of clinical outcomes differ between SARS-CoV-2 and influenza pneumonia.
  *EBioMedicine* 2022;85:104295; PMC9527494 ✓ — **no local PDF.**
  — carries the device F2 adopts: **every patient is given the full observation
  horizon inclusive of time after discharge or death**, which removes ragged
  right edges from a per-patient display.

## Measurement, conversion and reporting

- **Goradia S, et al.** Vasopressor dose equivalence: a scoping review and
  suggested formula. *J Crit Care* 2021;61:233–240.
  [doi:10.1016/j.jcrc.2020.11.002](https://doi.org/10.1016/j.jcrc.2020.11.002);
  PMID 33220576 — the NEE coefficients in `covariates.json`. *(Supplied by SG,
  2026-09-05; full citation already in the design notes.)*
- **Brown SM, et al.** *Chest* 2016; PMID 26836924, validated PMID 28538439 — the
  Severinghaus-based S/F → P/F mapping used instead of Rice. *(As cited in the design notes.)*
- **Rice TW, et al.** *Chest* 2007 — the linear `S/F = 64 + 0.84 × P/F` this
  project deliberately does **not** use. *(As cited in the design notes.)*
- **Chong SL, et al.** PMID 33780397 — arterial/venous pH offset. *(As cited in the design notes.)*
- **Yehya N, Harhay MO, et al.** *AJRCCM* 2019 — ventilator-free-day reporting;
  report median (IQR), rank-based tests. *(As cited in the design notes.)*
- **Anderson JR, Cain KC, Gelber RD.** *J Clin Oncol* 1983; and **Dafni U.**
  *Circ Cardiovasc Qual Outcomes* 2011 — landmark analysis and immortal-time
  bias, the basis of the landmark analysis. *(As cited in the design notes.)*
