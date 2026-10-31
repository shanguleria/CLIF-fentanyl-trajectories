# ==============================================================================
# 02_descriptive_cohort.R  --  cohort description: who, how long, how much fentanyl
#
# Purpose : Who the cohort is and how it is ventilated -- the retention table
#           that supplies the at-risk denominator, the dose trajectory over the
#           granular grid, per-admission ventilation and sedation over the whole
#           course, the 28-day liberation picture, and the pooling exports.
# Author  : Shan Guleria
# Created : 2026-09-05
# Inputs  : output/intermediate_phi/trajectory_long.parquet, time_to_event.parquet
# Outputs : output/final_no_phi/02_descriptive/ : the CSVs and figures listed in OWNED
# ==============================================================================

# Run this in a FRESH R session; rm(list = ls()) does not unload packages.

# ---- 1. Packages -------------------------------------------------------------

pkgs <- c("here", "jsonlite", "ggplot2", "arrow", "ggalluvial")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

source(here("code", "utils", "paths.R"))
source(here("code", "utils", "pooling.R"))
source(here("code", "utils", "figures.R"))


# ---- 2. Config (never setwd(); here() anchors to the .Rproj) -----------------

config <- fromJSON(here("config", "config.json"), simplifyVector = FALSE)
COV       <- fromJSON(here("config", "protocol.json"), simplifyVector = FALSE)
set.seed(COV$model$seed)

WINDOW_H  <- COV$cohort$window_hours
EXTENT_H  <- COV$cohort$granular_extent_hours
LANDMARK  <- COV$cohort$landmark_hours

SED_UNITS <- COV$exposure$sedatives$units

# Drug list comes from the config; a hardcoded one makes a declared sedative a no-op.
SED_COLS <- unlist(COV$exposure$sedatives$columns)
DRUGS <- c(fentanyl = "total_dose",
           setNames(SED_COLS, sub("_dose$", "", SED_COLS)))
UNITS <- c(fentanyl = COV$exposure$units,
           setNames(vapply(SED_COLS, function(c) SED_UNITS[[c]], character(1)),
                    sub("_dose$", "", SED_COLS)))

MIN_CELL <- COV$reporting$small_cell_min_den

# Where a categorical ships in both a full and a collapsed form, the full form is
# published at `overall` only: both views over the same strata subtract.
FULL_SCOPE <- COV$disclosure$full_levels_strata_scope
stopifnot("protocol.json must declare disclosure.full_levels_strata_scope" =
            !is.null(FULL_SCOPE))

# Absolute config-declared bins, never local quantiles -- those cannot be pooled.
HSPEC <- COV$pooling$histograms$variables
stopifnot("protocol.json must declare pooling.histograms.variables" =
            !is.null(HSPEC))
HIST_MINE <- names(HSPEC)[vapply(HSPEC, function(s)
  identical(as.character(s$owner), "02"), logical(1))]

POOL <- list()



# ---- 3. Paths and provenance -------------------------------------------------

dirs <- site_dirs()
dirs$phase <- phase_dir(dirs, "02_descriptive")   # shareable outputs, subdivided by script
prov <- provenance(config)

message(sprintf("[02_descriptive_cohort] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))


# ---- 4. Guards ---------------------------------------------------------------

manifest <- require_manifest(dirs, here())
message(sprintf("  reading cohort-build outputs from code %s, generated %s",
                manifest$code_version, manifest$generated))

OWNED <- list(phase = c(
  "retention.csv",
  "dose_summary.csv", "course_summary.csv",
  "imv_encounter_duration.csv",
  "liberation_survival.csv",
  "pooling_continuous.csv",
  "pooling_histograms.csv",
  "provenance.json",
  "fentanyl_curves.png", "liberation_survival.png", "captions.md"))

# Renames leave stale twins that clear_owned_outputs no longer names; retire them here.
RETIRED <- c(
  # Retired 2026-10-30 (SG): the balanced panels and the sedative curves were
  # trajectory-era figures that were hard to read and that nothing downstream
  # consumed. Replaced by course_summary.csv, which describes the same drugs as
  # per-admission distributions over the WHOLE ventilation course rather than as
  # overlaid window curves. The pooling_spec entry went with them.
  file.path("output", "final_no_phi", "02_descriptive",
            c("balanced_panels.csv", "fentanyl_balanced_panels.png",
              "sedative_curves.png")),
  # Moved to 07_outcome_model/ on 2026-10-28 (protocol 0.17.0). Declared so a
  # copy left by an earlier run is cleared rather than read back as current.
  file.path("output", "final_no_phi", "02_descriptive",
            c("baseline_characteristics.csv", "pooling_categorical.csv")),
  file.path("output", "final_no_phi",
            paste0("phase1_", c("dose_curves.png", "pct_receiving.png",
                                "dose_distribution.png", "balanced_panels.png"))),
  file.path("output", "final_no_phi", "02_descriptive",
            paste0("phase1_", c(
              "baseline_characteristics.csv", "retention.csv",
              "dose_summary.csv", "dose_distribution.csv",
              "balanced_panels.csv", "zero_fraction.csv",
              "imv_episodes.csv", "choosing_T.csv",
              "pooling_continuous.csv", "pooling_categorical.csv",
              "provenance.json",
              "fentanyl_curves.png", "fentanyl_balanced_panels.png",
              "fentanyl_distribution.png", "sedative_curves.png",
              "state_prevalence.csv", "state_transitions.csv",
              "state_alluvial.png", "state_prevalence.png",
              "dose_state_prevalence.csv", "dose_state_transitions.csv",
              "dose_state_alluvial.png", "dose_state_prevalence.png"))),
  file.path("output", "final_no_phi", "02_descriptive",
            c("state_prevalence.csv", "state_transitions.csv",
              "dose_state_prevalence.csv", "dose_state_transitions.csv",
              "state_prevalence.png", "state_alluvial.png",
              "dose_state_prevalence.png", "dose_state_alluvial.png")))
n_cleared <- clear_owned_outputs(dirs, OWNED, retired = RETIRED)
if (n_cleared) message(sprintf("  cleared %d output(s) from a previous run", n_cleared))


# ---- 5. Read -----------------------------------------------------------------

long <- as.data.frame(read_parquet(
  file.path(phi_dir(dirs, "01_cohort"), "trajectory_long.parquet"),
  col_select = c("encounter_block", "patient_id", "window_idx", "window_start_hr",
                 "alive_admitted", "imv_status", "inf_dose", "bolus_dose",
                 "total_dose", "window_mcg", "propofol_dose", "midazolam_dose",
                 "dexmedetomidine_dose",
                 "age", "sex", "race", "cci", "bmi_admission", "weight_kg",
                 "sofa_total", "nee", "oxygenation", "lactate",
                 "first_imv_episode_hours", "n_imv_episodes",
                 "died")))          # terminal states, section 12b below

tte <- as.data.frame(read_parquet(
  file.path(phi_dir(dirs, "01_cohort"), "time_to_event.parquet"),
  col_select = c("encounter_block", "landmark_eligible",
                 # Section 11b. Everything the 28-day liberation figure needs
                 # about a block's fate; the WITHIN-block timing comes from
                 # imv_intervals.parquet, read there.
                 "died", "death_hr", "hours_to_discharge", "ext_hr",
                 "vfd_28", "vfd_hours", "vfd_zero_lost_on_vent",
                 "vfd_zero_for_death")))

n_windows_expected <- EXTENT_H / WINDOW_H
stopifnot(
  "window count on disk disagrees with the config grid" =
    length(unique(long$window_idx)) == n_windows_expected,
  "the grid does not span granular_extent_hours" =
    max(long$window_start_hr) + WINDOW_H == EXTENT_H
)

# This script's denominator is ventilated, not alive_admitted.
long$ventilated <- !is.na(long$imv_status) & long$imv_status == 1
stopifnot(
  "a window cannot be ventilated without being alive and admitted" =
    all(long$ventilated <= long$alive_admitted)
)

missing_cols <- setdiff(unname(DRUGS), names(long))
if (length(missing_cols)) {
  stop("protocol.json declares dose columns that trajectory_long does not carry: ",
       paste(missing_cols, collapse = ", "),
       ". Re-run code/01_build_cohort.py rather than dropping them silently.",
       call. = FALSE)
}

anchor_n <- length(unique(long$encounter_block))
elig     <- tte$encounter_block[tte$landmark_eligible]

cat(sprintf("\nCohort: %s ventilation episodes; %s landmark-eligible at T=%dh\n",
            format(anchor_n, big.mark = ","), format(length(elig), big.mark = ","),
            LANDMARK))


# ---- 6. Retention ------------------------------------------------------------

retention <- do.call(rbind, lapply(sort(unique(long$window_idx)), function(w) {
  x <- long[long$window_idx == w, ]
  data.frame(window_idx = w,
             window_start_hr = x$window_start_hr[1],
             n_ventilated = sum(x$ventilated),
             n_alive_admitted = sum(x$alive_admitted),
             pct_of_anchor_cohort = round(100 * sum(x$ventilated) / anchor_n, 1))
}))

cat("\nRetention (at-risk = still ventilated)\n")
print(retention, row.names = FALSE)


# ---- 7. Dose summaries -------------------------------------------------------

summarise_dose <- function(v, drug, denom, w, hr) {
  v <- v[!is.na(v)]
  if (!length(v)) return(NULL)
  q <- unname(quantile(v, c(0.25, 0.5, 0.75)))
  # n_receiving is the numerator behind pct_receiving_any, so the proportion
  # pools exactly.
  k <- sum(v > 0)
  # receivers_only is already filtered to v > 0, so the percentage is structurally 100.
  receivers <- identical(denom, "receivers_only")
  data.frame(drug = drug, unit = UNITS[[drug]], denominator = denom,
             window_idx = w, window_start_hr = hr, n = length(v),
             n_receiving = k,
             median = round(q[2], 3), q1 = round(q[1], 3), q3 = round(q[3], 3),
             mean = round(mean(v), 3), sd = round(stats::sd(v), 3),
             pct_receiving_any = if (receivers) NA_real_
                                 else round(100 * k / length(v), 1))
}

dose_summary <- do.call(rbind, lapply(names(DRUGS), function(drug) {
  col <- DRUGS[[drug]]
  do.call(rbind, lapply(sort(unique(long$window_idx)), function(w) {
    x <- long[long$window_idx == w & long$ventilated, ]
    rbind(summarise_dose(x[[col]], drug, "all_ventilated", w, x$window_start_hr[1]),
          summarise_dose(x[[col]][x[[col]] > 0], drug, "receivers_only", w,
                         x$window_start_hr[1]))
  }))
}))

TV_POOL <- c(sofa_total = "SOFA", nee = "Norepinephrine equivalent (mcg/kg/min)",
             oxygenation = "P/F ratio", lactate = "Lactate (mmol/L)")

for (drug in names(DRUGS)) {
  col <- DRUGS[[drug]]
  for (w in sort(unique(long$window_idx))) {
    x <- long[long$window_idx == w & long$ventilated, ]
    hr <- x$window_start_hr[1]
    POOL[[length(POOL) + 1]] <- pool_row(
      "by_window", sprintf("%s dose", drug), UNITS[[drug]], "all_ventilated",
      w, hr, x[[col]], MIN_CELL)
    POOL[[length(POOL) + 1]] <- pool_row(
      "by_window", sprintf("%s dose", drug), UNITS[[drug]], "receivers_only",
      w, hr, x[[col]][!is.na(x[[col]]) & x[[col]] > 0], MIN_CELL)
  }
}
for (v in names(TV_POOL)) {
  for (w in sort(unique(long$window_idx))) {
    x <- long[long$window_idx == w & long$ventilated, ]
    POOL[[length(POOL) + 1]] <- pool_row(
      "by_window", TV_POOL[[v]], NA_character_, "ventilated",
      w, x$window_start_hr[1], x[[v]], MIN_CELL)
  }
}

cat("\nFentanyl dose by window (all ventilated, zeros included)\n")
print(dose_summary[dose_summary$drug == "fentanyl" &
                     dose_summary$denominator == "all_ventilated",
                   c("window_start_hr", "n", "median", "q1", "q3",
                     "pct_receiving_any")], row.names = FALSE)


# ---- 9. The whole ventilation course, per admission --------------------------
# Replaces the balanced panels and the sedative curves, retired 2026-10-30.
# Those
# described dose as overlaid WINDOW curves inside the first 72h; this describes
# the ADMISSION over the whole course, which is what the cohort's shape actually
# turns on.
#
# FOUR NESTED COHORTS, NOT A PARTITION: all > landmark_24h > landmark_48h >
# landmark_72h. Reported together so a reader can see what the landmark
# restriction selects on EXPOSURE, not only on outcome. Never sum across them.
#
# The numbers come from course_exposure.parquet, which rebuilds the SAME hourly
# exposure grid at the 28-day horizon rather than the 72h one -- see
# course_exposure() in 01_build_cohort.py, which asserts that restricting it to
# the first 72h reproduces trajectory_long's window_mcg exactly.

course <- as.data.frame(read_parquet(
  file.path(phi_dir(dirs, "01_cohort"), "course_exposure.parquet")))
stopifnot(
  "course_exposure must cover the analytic cohort -- re-run code/01_build_cohort.py" =
    setequal(course$encounter_block, unique(long$encounter_block)))

# Landmark membership. THE ELIGIBILITY RULE IS NOT REIMPLEMENTED HERE: 72h comes
# from time_to_event.parquet, and 24h/48h from the per-block membership
# code/01_build_cohort.py already wrote from _landmark_eligible(). A second rule
# in R would be a second definition.
lm_mem <- as.data.frame(read_parquet(
  file.path(phi_dir(dirs, "01_cohort"), "landmark_membership.parquet")))
COHORTS <- c("all", sprintf("landmark_%dh", sort(unique(lm_mem$T_hours))))

course_members <- function(tag) {
  if (identical(tag, "all")) return(course$encounter_block)
  th <- as.numeric(sub("^landmark_(\\d+)h$", "\\1", tag))
  lm_mem$encounter_block[lm_mem$T_hours == th]
}

SED_CATS <- unlist(COV$exposure$sedatives$course_summary_categories)
SED_RATE_UNITS <- vapply(SED_CATS,
                         function(c) SED_UNITS[[paste0(c, "_dose")]],
                         character(1))
MEASURES <- c(
  imv_hours                 = "hours",
  fentanyl_infusion_hours   = "hours",
  fentanyl_mcg              = "mcg",
  fentanyl_mcg_per_imv_hour = "mcg/hour of IMV",
  setNames(rep("hours", length(SED_RATE_UNITS)),
           paste0(names(SED_RATE_UNITS), "_hours")),
  setNames(unname(SED_RATE_UNITS), paste0(names(SED_RATE_UNITS), "_rate_mean")))
stopifnot("course_exposure is missing a declared measure" =
            all(names(MEASURES) %in% names(course)))

course_summary <- do.call(rbind, lapply(COHORTS, function(tag) {
  x <- course[course$encounter_block %in% course_members(tag), ]
  do.call(rbind, lapply(names(MEASURES), function(m) {
    v <- x[[m]]
    v <- v[!is.na(v)]            # per-IMV-hour is undefined without IMV time
    if (length(v) < MIN_CELL) return(NULL)
    q <- unname(quantile(v, c(0.25, 0.5, 0.75, 0.90, 0.95)))
    data.frame(
      cohort = tag, measure = m, unit = unname(MEASURES[[m]]),
      n = length(v), n_nonzero = sum(v > 0),
      median = round(q[2], 3), q1 = round(q[1], 3), q3 = round(q[3], 3),
      p90 = round(q[4], 3), p95 = round(q[5], 3),
      mean = round(mean(v), 3), sd = round(stats::sd(v), 3),
      # The pair that makes this poolable: an exact pooled mean and SD need no
      # assumption about equal variances. The quantiles beside them do not pool.
      sum = round(sum(v), 3), sum_sq = round(sum(v^2), 3),
      min = round(min(v), 3), max = round(max(v), 3))
  }))
}))

cat("\nVentilation and sedation over the whole course, by cohort\n")
print(course_summary[course_summary$measure %in%
                       c("imv_hours", "fentanyl_mcg",
                         "fentanyl_mcg_per_imv_hour"),
                     c("cohort", "measure", "unit", "n", "median", "q1", "q3")],
      row.names = FALSE)


# ---- 11. IMV episodes --------------------------------------------------------

blocks <- long[!duplicated(long$encounter_block),
               c("encounter_block", "patient_id", "first_imv_episode_hours",
                 "n_imv_episodes")]
blocks$landmark_eligible <- blocks$encounter_block %in% elig

q_row <- function(v, name, unit) {
  v <- v[!is.na(v)]
  q <- unname(quantile(v, c(0.25, 0.5, 0.75, 0.9, 0.95)))
  data.frame(measure = name, unit = unit, n = length(v),
             median = round(q[2], 2), q1 = round(q[1], 2), q3 = round(q[3], 2),
             p90 = round(q[4], 2), p95 = round(q[5], 2),
             min = round(min(v), 2), max = round(max(v), 2))
}

per_patient <- table(blocks$patient_id)
imv_encounter_duration <- rbind(
  q_row(blocks$first_imv_episode_hours, "first_imv_episode_hours", "hours"),
  q_row(blocks$n_imv_episodes, "n_imv_episodes_per_block", "count"),
  q_row(as.numeric(per_patient), "blocks_per_patient", "count")
)

cat("\nVentilation episode duration and encounter structure\n")
print(imv_encounter_duration, row.names = FALSE)
cat(sprintf("  patients contributing more than one block: %s of %s (%.1f%%)\n",
            format(sum(per_patient > 1), big.mark = ","),
            format(length(per_patient), big.mark = ","),
            100 * mean(per_patient > 1)))


# ---- 11b. Liberation and survival to day 28 ----------------------------------
# After MODE (Semler et al., NEJM Evid 2024) Figure 3. The question it answers
# is not "what is the median VFD" but "where did the VFD go": a VFD-28 of 0 can
# mean the patient died, or stayed ventilated, or was discharged alive while
# respiratory_support was still charting IMV -- and those are different facts
# with the same number. Each is a band here.
#
# DENOMINATOR IS THE WHOLE LANDMARK COHORT, not 07's complete-case modelling
# frame. The figure describes the cohort's outcome, not the model's.
#
# DISCHARGE ALIVE IS CARRIED FORWARD TO DAY 28, in whatever ventilation state
# the patient left in. That is not a new assumption: it is exactly what
# ventilator_free_days() already does, so this decomposes VFD rather than
# estimating it a second way. It also means EVERY block has a known state at
# every point: there is no censoring, so the Aalen-Johansen cumulative
# incidence collapses to an EXACT PROPORTION, which is what is computed below.
# No survival package is involved -- calling survfit() on uncensored data to
# borrow a method's name would be decoration. See code/utils/dependencies.R.

iv <- as.data.frame(read_parquet(file.path(phi_dir(dirs, "01_cohort"), "imv_intervals.parquet")))
stopifnot(
  "imv_intervals must carry the qualifying flag -- re-run code/01_build_cohort.py" =
    all(c("t_start_hr", "t_end_hr", "is_imv", "qualifies", "accrued_h") %in% names(iv)),
  # Containment, not equality: a landmark block whose every interval fell
  # outside [0, horizon] is legitimately absent from iv. An interval belonging
  # to no landmark block is not -- that is a join that went wrong.
  "imv_intervals carries a block the landmark cohort does not" =
    all(iv$encounter_block %in% elig))

HORIZON_H  <- as.numeric(COV$outcomes$vfd$horizon_days) * 24
GRID_STEP_H <- 6            # smooth enough to draw, small enough to pool
grid_h <- seq(0, HORIZON_H, by = GRID_STEP_H)

# Read, not restated: the caption quotes this site's carry-forward share, and a
# number typed into a caption is a number that goes stale on the next rerun.
wf_f <- file.path(dirs$out_final, "01_cohort", "diagnostics",
                  "waterfall_carryforward.csv")
stopifnot("no waterfall diagnostic -- re-run code/01_build_cohort.py" =
            file.exists(wf_f))
wf <- read.csv(wf_f, stringsAsFactors = FALSE)
WF_CARRIED_PCT <- wf$pct[wf$metric == "cohort IMV-hours that are carry-forward"]
stopifnot("the waterfall diagnostic has no carry-forward share" =
            length(WF_CARRIED_PCT) == 1 && is.finite(WF_CARRIED_PCT))

lib <- tte[tte$landmark_eligible, ]
lib <- lib[order(lib$encounter_block), ]
n_lib <- nrow(lib)

free_iv <- iv[iv$qualifies, c("encounter_block", "t_start_hr", "t_end_hr")]
# Block x grid-point occupancy, filled by walking the ~2.5k qualifying SPELLS
# rather than the 5.8k blocks x 113 points: the same answer, two orders of
# magnitude fewer R-level calls.
blk <- match(free_iv$encounter_block, lib$encounter_block)
stopifnot("a qualifying spell belongs to no landmark block" = !anyNA(blk))
in_free <- matrix(FALSE, nrow = n_lib, ncol = length(grid_h))
# AT THE HORIZON, ASK ABOUT THE INSTANT BEFORE IT. Intervals are half-open and
# clipped to the horizon, so a literal `t < t_end_hr` test at t = HORIZON_H is
# false for everyone and the last column of the figure empties out -- the
# day-28 bar would show every survivor back on a ventilator.
g_eval <- pmin(grid_h, HORIZON_H - 1e-9)
for (k in seq_len(nrow(free_iv))) {
  hit <- g_eval >= free_iv$t_start_hr[k] & g_eval < free_iv$t_end_hr[k]
  if (any(hit)) in_free[blk[k], hit] <- TRUE
}

# First qualifying liberation per block: the cumulative-incidence event time.
first_lib <- rep(Inf, n_lib)
starts <- tapply(free_iv$t_start_hr, free_iv$encounter_block, min)
first_lib[match(names(starts), lib$encounter_block)] <- as.numeric(starts)

# Occupancy, in PRIORITY order -- each test excludes the ones before it, so the
# four states are mutually exclusive by construction rather than by hoping they
# are. The partition is then asserted below anyway.
occ <- do.call(rbind, lapply(seq_along(grid_h), function(j) {
  t <- grid_h[j]
  dead  <- !is.na(lib$death_hr) & lib$death_hr <= t
  gone  <- !dead & (t >= lib$hours_to_discharge)
  lost  <- gone & lib$vfd_zero_lost_on_vent
  # Discharged alive OFF the ventilator: still alive and still free, by the
  # same extension ventilator_free_days() applies.
  freed <- !dead & !lost & (in_free[, j] | gone)
  data.frame(
    hour = t, day = t / 24,
    measure = c("alive, on IMV", "alive, free of IMV",
                "discharged alive, charted on IMV", "died", "cum_liberated"),
    n_admissions = c(sum(!dead & !lost & !freed), sum(freed), sum(lost),
                     sum(dead), sum(first_lib <= t)),
    denominator = n_lib, stringsAsFactors = FALSE)
}))
occ$pct <- round(100 * occ$n_admissions / occ$denominator, 2)

states4 <- setdiff(unique(occ$measure), "cum_liberated")
sums <- tapply(occ$n_admissions[occ$measure %in% states4],
               occ$hour[occ$measure %in% states4], sum)
stopifnot(
  "the four occupancy states must partition the landmark cohort at every point" =
    all(sums == n_lib),
  "the cumulative incidence of liberation must be monotone" =
    !is.unsorted(occ$n_admissions[occ$measure == "cum_liberated"]))

# THE AREA IDENTITY, verified rather than asserted in a caption. The free-of-IMV
# band integrated from the landmark to the horizon is mean ALIVE-AND-IMV-FREE
# DAYS. That is NOT mean VFD-28, and the gap between them is the whole point:
# VFD retroactively zeroes anyone who dies inside the horizon, so the lived
# ventilator-free time of everyone who died is real and is scored 0. The counts
# are printed below rather than written here, where they would go stale.
free_n <- occ$n_admissions[occ$measure == "alive, free of IMV"]
# A grid point stands for the step that STARTS there, so the horizon itself is
# not a step and must not be summed -- it is reported, not integrated.
in_win <- grid_h >= LANDMARK & grid_h < HORIZON_H
area_d <- sum(free_n[in_win]) * GRID_STEP_H / 24 / n_lib   # left Riemann sum
exact_d <- sum(iv$accrued_h) / 24 / n_lib
cat(sprintf(
  paste0("\nLiberation and survival to day %g\n",
         "  alive and IMV-free days per admission : %.2f  (grid) vs %.2f (exact)\n",
         "  mean VFD-28 as scored                 : %.2f\n",
         "    the difference is the %s admission(s) whose accrued free time\n",
         "    VFD zeroes for death inside the horizon, plus %s zeroed as\n",
         "    discharged alive while still charted on IMV\n"),
  HORIZON_H / 24, area_d, exact_d, mean(lib$vfd_28),
  format(sum(lib$vfd_zero_for_death), big.mark = ","),
  format(sum(lib$vfd_zero_lost_on_vent), big.mark = ",")))
stopifnot(
  "the drawn band and the exported intervals disagree about free time" =
    abs(area_d - exact_d) < 0.25)

liberation_survival <- occ[, c("hour", "day", "measure", "n_admissions",
                               "denominator", "pct")]

register_caption("liberation_survival.png",
  "Liberation from invasive ventilation and survival through day 28",
  sprintf(paste0(
    "All %s admissions ventilated at the %gh landmark, followed to day %g from ",
    "the start of invasive ventilation; the dashed rule marks the landmark. ",
    "UPPER PANEL, after MODE (Semler et al., NEJM Evid 2024) Figure 3: the ",
    "share alive, and the cumulative incidence of a first SUCCESSFUL ",
    "liberation -- at least %gh free of invasive ventilation -- with death as ",
    "a competing risk. LOWER PANEL: where the cohort actually is, hour by ",
    "hour. Discharge alive is carried forward to day %g in the ventilation ",
    "state the patient left in, which is the rule ventilator_free_days() ",
    "applies, so there is no censoring and the cumulative incidence is an ",
    "exact proportion. The orange band is the %s admission(s) discharged ",
    "alive while respiratory_support was still charting invasive ventilation; ",
    "all of them are scored VFD-28 = 0 by rule, and %s%% of the cohort's ",
    "invasive-ventilation hours are waterfall carry-forward rather than ",
    "charted device records (01_cohort/diagnostics/waterfall_carryforward.csv). ",
    "The free-of-IMV band integrates to %.1f alive and IMV-free days per ",
    "admission between the landmark and day %g; mean VFD-28 is %.1f days, ",
    "lower because VFD retroactively scores 0 for the %s admission(s) who ",
    "died inside the horizon. Mortality is IN-HOSPITAL and so a lower bound ",
    "on %g-day mortality. Grid %gh; denominator %s throughout."),
    format(n_lib, big.mark = ","), LANDMARK, HORIZON_H / 24,
    as.numeric(COV$outcomes$successful_extubation_hours), HORIZON_H / 24,
    format(sum(lib$vfd_zero_lost_on_vent), big.mark = ","),
    WF_CARRIED_PCT, area_d, HORIZON_H / 24, mean(lib$vfd_28),
    format(sum(lib$vfd_zero_for_death), big.mark = ","),
    HORIZON_H / 24, GRID_STEP_H, format(n_lib, big.mark = ",")))

p_lib <- liberation_panels(liberation_survival, denom = n_lib,
                           landmark_d = LANDMARK / 24,
                           horizon_d = HORIZON_H / 24)
ggsave(file.path(dirs$phase, "liberation_survival.png"), p_lib,
       width = 7.2, height = 8.4, dpi = 300, bg = PAPER)


# ---- 12. Baseline characteristics -- MOVED TO 07 --------------------------
# Table 1 is generated by code/07_outcome_model.R as of protocol 0.17.0: it is
# the only script holding the exposures and the fitted model frame, and it
# reads trajectory_long window 0 unrestricted so it can describe the whole
# cohort too. One generator, so the table and the model cannot describe
# different populations. Row builders live in code/utils/table1.R.


# ---- 12b. Pooling histograms -------------------------------------------------
# Counts on fixed absolute bins pool where a median does not. See protocol.json
# pooling.histograms.
pool_histograms <- do.call(rbind, lapply(HIST_MINE, function(v) {
  spec <- HSPEC[[v]]
  col <- sub("^trajectory_long\\.", "", as.character(spec$source))
  stopifnot("a declared histogram names a column trajectory_long does not carry" =
              col %in% names(long))
  x <- if (identical(as.character(spec$restrict), "one row per episode")) {
    # One value per episode: off the window grid unfiltered it would weight long stays.
    long[[col]][!duplicated(long$encounter_block)]
  } else {
    long[[col]][long$ventilated]
  }
  pool_hist(v, x, spec)
}))
stopifnot(
  "02 must emit every histogram the config assigns it" =
    setequal(unique(pool_histograms$variable), HIST_MINE))
cat(sprintf("\nPooling histograms: %d variables, %d bins total\n",
            length(HIST_MINE), nrow(pool_histograms)))


# ---- 13. Figures -------------------------------------------------------------

FENT_U <- UNITS[["fentanyl"]]
DOSE_COLS <- c("All ventilated, zeros included" = "#14427e",
               "Receivers only" = "#4a8bd8")

curves <- dose_summary
curves$series <- ifelse(curves$denominator == "all_ventilated",
                        "All ventilated, zeros included", "Receivers only")

# --- Primary figure: fentanyl, all three dose curves -------------------------
# Stacked facets, not sec_axis: sec_axis applies one transform to every facet and
# would draw the proportion against the wrong scale under free scales.
f <- curves[curves$drug == "fentanyl", ]
f_pct <- f[f$denominator == "all_ventilated", ]

SERIES <- c("Median dose, all ventilated (zeros included)" = "#14427e",
            "Median dose, receivers only"                  = "#4a8bd8",
            "% receiving any fentanyl"                     = "#eb6834")

fent_df <- rbind(
  data.frame(window_start_hr = f$window_start_hr,
             value = f$median, lo = f$q1, hi = f$q3,
             series = ifelse(f$denominator == "all_ventilated",
                             names(SERIES)[1], names(SERIES)[2]),
             quantity = sprintf("Dose (%s), median and IQR", FENT_U)),
  data.frame(window_start_hr = f_pct$window_start_hr,
             value = f_pct$pct_receiving_any, lo = NA_real_, hi = NA_real_,
             series = names(SERIES)[3],
             quantity = "% of ventilated episodes receiving any")
)
fent_df$series <- factor(fent_df$series, levels = names(SERIES))
fent_df$quantity <- factor(fent_df$quantity, levels = unique(fent_df$quantity))

p_fent <- house(
  ggplot(fent_df, aes(window_start_hr, value, colour = series, fill = series)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.13, colour = NA,
                na.rm = TRUE, show.legend = FALSE) +
    geom_line(linewidth = 1.0) + geom_point(size = 1.6) +
    facet_wrap(~ quantity, scales = "free_y", ncol = 1) +
    scale_colour_manual(values = SERIES) +
    scale_fill_manual(values = SERIES) +
    guides(fill = "none",
           colour = guide_legend(nrow = 2, override.aes = list(fill = NA))) +
    labs(x = "Hours since first IMV episode", y = NULL, colour = NULL))

register_caption("fentanyl_curves.png",
  "Fentanyl dose over the first 72h of ventilation",
  paste0("Denominator is episodes STILL VENTILATED in each window. The two ",
         "medians diverge: exposure narrows to fewer episodes rather than ",
         "falling within them."))
ggsave(file.path(dirs$phase, "fentanyl_curves.png"), p_fent,
       width = 7.5, height = 7.6, dpi = 200)

# ---- Write -------------------------------------------------------------------

write_out <- function(x, name) {
  f <- file.path(dirs$phase, name)
  write.csv(x, f, row.names = FALSE)
  cat(sprintf("written: %s\n", name))
}

cat("\n")
write_out(retention, "retention.csv")
write_out(dose_summary, "dose_summary.csv")
write_out(imv_encounter_duration, "imv_encounter_duration.csv")
write_out(course_summary, "course_summary.csv")
write_out(liberation_survival, "liberation_survival.csv")

pooling_continuous <- do.call(rbind, POOL)
write_out(pooling_continuous, "pooling_continuous.csv")
write_out(pool_histograms, "pooling_histograms.csv")

write_captions(file.path(dirs$phase, "captions.md"), "02_descriptive_cohort.R",
               grep("\\.png$", OWNED$phase, value = TRUE), prov)
cat("written: captions.md\n")

write_json(prov, file.path(dirs$phase, "provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)
cat("written: provenance.json\n")


# ---- Provenance --------------------------------------------------------------

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/02_descriptive_cohort.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "02_descriptive_cohort_sessioninfo.txt")
)
