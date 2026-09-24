# ==============================================================================
# 05_landmark_cohort.R  --  Phase 2 -- apply landmark T, retention reporting
#
# Purpose : Restrict to episodes alive and ventilated at landmark T; report the landmark flow, failed-extubation counts, the [0,T] dose curve and the T sensitivity sweep; write the analytic table Phases 3-6 read.
# Author  : Shan Guleria
# Created : 2026-09-05
# Inputs  : output/intermediate_phi/trajectory_long.parquet, time_to_event.parquet
# Outputs : output/intermediate_phi/landmark_cohort.parquet; output/final_no_phi/04_landmark/
#
# ==============================================================================

# Run this in a FRESH R session (RStudio: Cmd+Shift+F10).
# Do not use rm(list = ls()) -- it does not unload packages or reset options,
# so it only gives the appearance of a clean slate.


# ---- 1. Packages -------------------------------------------------------------

pkgs <- c("here", "jsonlite", "ggplot2", "arrow")
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
set.seed(config$model$seed)

WINDOW_H <- config$cohort$window_hours
EXTENT_H <- config$cohort$granular_extent_hours
LANDMARK <- config$cohort$landmark_hours
GAP_H    <- config$cohort$imv_episode_gap_hours
MIN_CELL <- config$reporting$small_cell_min_den
SUCC_H   <- config$outcomes$successful_extubation_hours

COV       <- fromJSON(here("config", "covariates.json"), simplifyVector = FALSE)
SED_COLS  <- unlist(COV$exposure$sedatives$columns)
SED_UNITS <- COV$exposure$sedatives$units
DRUGS <- c(fentanyl = "total_dose",
           setNames(SED_COLS, sub("_dose$", "", SED_COLS)))
UNITS <- c(fentanyl = COV$exposure$units,
           setNames(vapply(SED_COLS, function(c) SED_UNITS[[c]], character(1)),
                    sub("_dose$", "", SED_COLS)))
RACE_COLLAPSE <- COV$time_invariant$race$reporting_collapse

# Section 8 requires the sensitivity reported together with the primary, so T is
# swept even though config.landmark_hours fixes the analysis at 72h (SG).
T_SWEEP <- sort(unique(c(24, 48, LANDMARK)))


# ---- 3. Paths and provenance -------------------------------------------------
# One site, one output tree. site_dirs() creates them and labels the PHI ones.

dirs <- site_dirs()
dirs$phase <- phase_dir(dirs, "05_landmark")   # shareable outputs, subdivided by script
prov <- provenance(config)

message(sprintf("[05_landmark_cohort] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))


# ---- 4. Guards ---------------------------------------------------------------

manifest <- require_manifest(dirs, here())
message(sprintf("  reading Phase 0 outputs from code %s, generated %s",
                manifest$code_version, manifest$generated))

OWNED <- list(
  # landmark_cohort.csv was dropped from the writer but left here until
  # 2026-09-23. A declared output nothing writes is the same defect as a
  # declared config key nothing reads.
  out_phi = c("landmark_cohort.parquet"),
  phase = c("landmark_flow.csv", "landmark_flow.txt",
                "T_sensitivity.csv", "failed_extubation.csv",
                "dose_curve.csv", "dose_curve.png",
                "dependence.csv",
                "pooling_continuous.csv", "pooling_categorical.csv",
                "provenance.json", "captions.md"))
# This script was 03_landmark_cohort.R writing phase2_* into 03_landmark/ until
# 2026-09-23, when 02 split in two and took the 03 slot. Both the directory and
# the prefix changed, so every old path is retired here -- otherwise the whole of
# 03_landmark/ sits in the shareable tree forever, looking current.
# Renumbered again 04 -> 05 on 2026-09-24, when the exemplar took a slot ahead
# of the states script. Both old folders are retired, each built with paste0 so
# a name-level find-and-replace cannot reach into the list of old names -- the
# exact failure recorded as lessons.md #13.
RETIRED <- c(
  file.path("output", "final_no_phi", paste0("03_", "landmark"),
            paste0("phase2_", c(
              "landmark_flow.csv", "landmark_flow.txt",
              "T_sensitivity.csv", "failed_extubation.csv",
              "dose_curve.csv", "dose_curve.png",
              "dependence.csv",
              "pooling_continuous.csv", "pooling_categorical.csv",
              "provenance.json"))),
  file.path("output", "final_no_phi", paste0("04_", "landmark"), OWNED$phase))
n_cleared <- clear_owned_outputs(dirs, OWNED, retired = RETIRED)
if (n_cleared) message(sprintf("  cleared %d output(s) from a previous run", n_cleared))


# ---- 5. Read -----------------------------------------------------------------

long <- as.data.frame(read_parquet(file.path(dirs$out_phi, "trajectory_long.parquet")))
tte  <- as.data.frame(read_parquet(file.path(dirs$out_phi, "time_to_event.parquet")))

stopifnot(
  "window count on disk disagrees with the config grid" =
    length(unique(long$window_idx)) == EXTENT_H / WINDOW_H,
  "landmark_hours is not on the window grid" = LANDMARK %% WINDOW_H == 0,
  "landmark_hours exceeds the grid this table covers" = LANDMARK <= EXTENT_H
)

long$ventilated <- !is.na(long$imv_status) & long$imv_status == 1
anchor_n <- length(unique(long$encounter_block))


# ---- 6. Landmark eligibility, at any T ---------------------------------------
# Eligible = alive, admitted and STILL VENTILATED in the window that ends at T.
# Same rule Phase 0 uses to build time_to_event; recomputed here so the sweep can
# reach T values that table does not cover.

eligible_at <- function(t_hours) {
  w <- t_hours / WINDOW_H - 1
  x <- long[long$window_idx == w & long$alive_admitted & long$ventilated, ]
  unique(x$encounter_block)
}

elig <- eligible_at(LANDMARK)
stopifnot(
  "recomputed eligibility disagrees with Phase 0's time_to_event table" =
    setequal(elig, tte$encounter_block[tte$landmark_eligible])
)
cat(sprintf("\nLandmark T = %dh: %s of %s episodes eligible (%.1f%%)\n",
            LANDMARK, format(length(elig), big.mark = ","),
            format(anchor_n, big.mark = ","), 100 * length(elig) / anchor_n))


# ---- 7. Failed extubation inside [0, T] --------------------------------------
# Two definitions, reported together, because they bound the same quantity from
# opposite sides. A single window of imv_status == 0 may be a charting gap rather
# than an extubation; the episode builder treats a break as real only past
# imv_episode_gap_hours, so the gap-aligned count requires that many hours of
# consecutive zero windows.

MIN_RUN <- ceiling(GAP_H / WINDOW_H)

max_zero_run <- function(z) {
  if (!length(z) || !any(z)) return(0L)
  r <- rle(z)
  max(c(0L, r$lengths[r$values]))
}

pre_T <- long[long$encounter_block %in% elig & long$window_start_hr < LANDMARK, ]
pre_T <- pre_T[order(pre_T$encounter_block, pre_T$window_idx), ]
runs <- do.call(rbind, lapply(split(pre_T, pre_T$encounter_block), function(g) {
  z <- !g$ventilated
  data.frame(encounter_block = g$encounter_block[1],
             n_zero_windows = sum(z),
             max_zero_run = max_zero_run(z))
}))
runs$failed_liberal <- runs$n_zero_windows > 0
runs$failed_gap_aligned <- runs$max_zero_run >= MIN_RUN

failed <- data.frame(
  definition = c("any non-ventilated window in [0,T)",
                 sprintf("run of >= %d windows (>= %dh, the episode gap rule)",
                         MIN_RUN, MIN_RUN * WINDOW_H)),
  n_episodes = c(sum(runs$failed_liberal), sum(runs$failed_gap_aligned)),
  pct_of_landmark_cohort = round(
    100 * c(mean(runs$failed_liberal), mean(runs$failed_gap_aligned)), 1),
  n_remaining_if_excluded = length(elig) -
    c(sum(runs$failed_liberal), sum(runs$failed_gap_aligned)),
  stringsAsFactors = FALSE)

cat("\nFailed extubations inside [0, T] -- % lost if they were excluded\n")
print(failed, row.names = FALSE)
cat("  Section 10a keeps them: cohort membership is 'ventilated at T', which\n")
cat("  needs no look-ahead. The counts above are the required disclosure.\n")


# ---- 8. T sensitivity (reported alongside the primary) -----------------------

sens <- do.call(rbind, lapply(T_SWEEP, function(t_hours) {
  ids <- eligible_at(t_hours)
  x <- long[long$encounter_block %in% ids & long$window_start_hr < t_hours, ]
  b <- x[x$window_idx == 0, ]
  d <- x$total_dose[x$ventilated & !is.na(x$total_dose)]
  pre <- x[order(x$encounter_block, x$window_idx), ]
  fr <- vapply(split(pre$ventilated, pre$encounter_block),
               function(v) as.integer(max_zero_run(!v) >= MIN_RUN), integer(1))
  data.frame(
    T_hours = t_hours,
    is_primary = t_hours == LANDMARK,
    n_eligible = length(ids),
    pct_of_intubated = round(100 * length(ids) / anchor_n, 1),
    windows_available = t_hours / WINDOW_H,
    n_failed_extubation = sum(fr),
    pct_failed_extubation = round(100 * mean(fr), 1),
    age_median = round(median(b$age, na.rm = TRUE), 1),
    pct_female = round(100 * mean(b$sex == "Female", na.rm = TRUE), 1),
    sofa_median = round(median(b$sofa_total, na.rm = TRUE), 1),
    dose_mean = round(mean(d), 2), dose_median = round(median(d), 2),
    pct_windows_zero_dose = round(100 * mean(d == 0), 1))
}))

cat("\nT sensitivity\n")
print(sens, row.names = FALSE)


# ---- 9. The [0,T] dose curve, and the Phase 1 consistency check ---------------
# Section 10 says this is the same cohort as Phase 1's >=T balanced panel, so the
# curve should reproduce it. That is the check, not a new finding.

lm_long <- long[long$encounter_block %in% elig & long$window_start_hr < LANDMARK, ]

# NO ventilated filter here, deliberately. This is a BALANCED PANEL: the same
# episodes contribute every window, which is (a) what makes the Phase 1 panel a
# panel at all -- freezing the denominator is its whole point -- and (b) what
# gbmt requires, since an unbalanced panel silently caps the polynomial degree at
# (shortest unit's windows - 1). A transiently extubated window carries dose 0 by
# the extubated-gap rule, which is the settled treatment, not missingness. The
# ventilated-only counts ride alongside so the difference stays visible.
curve <- do.call(rbind, lapply(names(DRUGS), function(drug) {
  col <- DRUGS[[drug]]
  do.call(rbind, lapply(sort(unique(lm_long$window_idx)), function(w) {
    x <- lm_long[lm_long$window_idx == w, ]
    v <- x[[col]][!is.na(x[[col]])]
    if (!length(v)) return(NULL)
    q <- unname(quantile(v, c(0.25, 0.5, 0.75)))
    vent <- x[[col]][x$ventilated & !is.na(x[[col]])]
    data.frame(drug = drug, unit = UNITS[[drug]], window_idx = w,
               window_start_hr = x$window_start_hr[1], n = length(v),
               mean = round(mean(v), 3), sd = round(stats::sd(v), 3),
               median = round(q[2], 3), q1 = round(q[1], 3), q3 = round(q[3], 3),
               pct_receiving_any = round(100 * mean(v > 0), 1),
               n_ventilated = length(vent),
               mean_ventilated_only = round(mean(vent), 3))
  }))
}))

# Phase 1's output, not this script's -- read it from ITS folder.
p1 <- file.path(phase_dir(dirs, "02_descriptive"), "balanced_panels.csv")
if (file.exists(p1)) {
  bp <- read.csv(p1)
  bp <- bp[bp$panel_hours == LANDMARK, ]
  f <- curve[curve$drug == "fentanyl", ]
  m <- merge(f[, c("window_start_hr", "mean", "n")],
             bp[, c("window_start_hr", "mean", "n")],
             by = "window_start_hr", suffixes = c("_phase2", "_phase1"))
  d_mean <- max(abs(m$mean_phase2 - m$mean_phase1))
  d_n <- max(abs(m$n_phase2 - m$n_phase1))
  cat(sprintf(
    "\nConsistency with the Phase 1 >=%dh balanced panel: max |mean diff| %.4f %s, max |n diff| %d\n",
    LANDMARK, d_mean, UNITS[["fentanyl"]], d_n))
  if (d_mean > 1e-6 || d_n > 0) {
    stop(sprintf(paste0(
      "Phase 2 disagrees with the Phase 1 >=%dh balanced panel (max mean diff ",
      "%.4f, max n diff %d). They are the same episodes over the same windows, ",
      "so a difference is a bug in one of the two -- not a finding."),
      LANDMARK, d_mean, d_n), call. = FALSE)
  }
  cat("  identical, as they must be.\n")
} else {
  cat("\nPhase 1 balanced panels absent; consistency check skipped.\n")
}


# ---- 10. Repeat-episode dependence, within the landmark cohort ---------------
# Section 11 requires this: keeping every encounter block is defensible only if
# the dependence it admits is measured. Two of the three diagnostics land here;
# the third (classes holding two episodes from one patient) needs Phase 3.

lm_blocks <- lm_long[!duplicated(lm_long$encounter_block),
                     c("encounter_block", "patient_id")]
per_pt <- table(lm_blocks$patient_id)
dep <- data.frame(
  measure = c("episodes", "patients", "patients with >1 episode",
              "episodes from a repeating patient", "max episodes per patient"),
  value = c(nrow(lm_blocks), length(per_pt), sum(per_pt > 1),
            sum(per_pt[per_pt > 1]), max(per_pt)),
  pct = c(NA, NA, round(100 * mean(per_pt > 1), 1),
          round(100 * sum(per_pt[per_pt > 1]) / nrow(lm_blocks), 1), NA))

cat("\nRepeat-episode dependence in the landmark cohort\n")
print(dep, row.names = FALSE)


# ---- 11. Landmark flow -------------------------------------------------------

flow <- data.frame(
  step = c("Analytic cohort (Phase 0)",
           sprintf("Alive and ventilated at T = %dh", LANDMARK),
           "Analysed (Phases 3-6)"),
  n = c(anchor_n, length(elig), length(elig)),
  excluded = c(NA, anchor_n - length(elig), 0),
  reason = c(NA, sprintf("extubated, died or discharged before T = %dh", LANDMARK),
             "none -- failed extubations are retained"),
  stringsAsFactors = FALSE)

cat("\nLandmark flow\n")
print(flow, row.names = FALSE)

writeLines(c(
  sprintf("LANDMARK FLOW  (T = %dh)", LANDMARK),
  strrep("=", 60),
  sprintf("%-46s %9s", "Analytic cohort (Phase 0)", format(anchor_n, big.mark = ",")),
  sprintf("%-46s %9s", "  excluded before T",
          format(anchor_n - length(elig), big.mark = ",")),
  sprintf("%-46s %9s", sprintf("Alive and ventilated at T = %dh", LANDMARK),
          format(length(elig), big.mark = ",")),
  sprintf("%-46s %9s", "Analysed (Phases 3-6)", format(length(elig), big.mark = ",")),
  "",
  sprintf("Failed extubations retained: %s (%.1f%%) by the %dh gap rule",
          format(failed$n_episodes[2], big.mark = ","),
          failed$pct_of_landmark_cohort[2], MIN_RUN * WINDOW_H),
  sprintf("Estimand: conditional on being alive and mechanically ventilated at T = %dh.",
          LANDMARK),
  "The unit is the ventilation episode, not the patient."),
  file.path(dirs$phase, "landmark_flow.txt"))


# ---- 12. Figure --------------------------------------------------------------

# house() and the palette come from code/utils/figures.R, sourced at the top.
# This script carried its OWN copy until 2026-09-24 -- a verbatim duplicate that
# had already stopped matching, since it never got the legend changes the shared
# theme did. figures.R exists precisely so a figure drawn here looks like one
# drawn in 02 or 03; a private copy is how that quietly stops being true.

f <- curve[curve$drug == "fentanyl", ]
fig <- rbind(
  data.frame(window_start_hr = f$window_start_hr, value = f$mean,
             quantity = sprintf("MEAN dose (%s)", UNITS[["fentanyl"]])),
  data.frame(window_start_hr = f$window_start_hr, value = f$median,
             quantity = sprintf("MEDIAN dose (%s)", UNITS[["fentanyl"]])),
  data.frame(window_start_hr = f$window_start_hr, value = f$pct_receiving_any,
             quantity = "% receiving any fentanyl"))
fig$quantity <- factor(fig$quantity, levels = unique(fig$quantity))

p_curve <- house(
  ggplot(fig, aes(window_start_hr, value)) +
    geom_line(linewidth = 1.0, colour = "#14427e") +
    geom_point(size = 1.5, colour = "#14427e") +
    facet_wrap(~ quantity, scales = "free_y", ncol = 1) +
    labs(x = "Hours since first IMV episode", y = NULL))

register_caption("dose_curve.png",
  sprintf("Fentanyl over [0, %dh] in the landmark cohort", LANDMARK),
  sprintf(paste0("%s ventilation episodes alive and ventilated at T = %dh ",
                 "(%.1f%% of the intubated cohort). This is the Phase 1 >=%dh ",
                 "balanced panel promoted to the analytic set, not new analysis."),
          format(length(elig), big.mark = ","), LANDMARK,
          100 * length(elig) / anchor_n, LANDMARK))
ggsave(file.path(dirs$phase, "dose_curve.png"), p_curve,
       width = 7.5, height = 7.6, dpi = 200)


# ---- 13. Pooling exports for the landmark cohort -----------------------------

base <- lm_long[lm_long$window_idx == 0, ]
one <- rep("landmark", nrow(base))

POOL <- list()
for (v in list(c("Age, years", "age", "years"),
               c("Charlson Comorbidity Index", "cci", "index"),
               c("BMI at admission, kg/m2", "bmi_admission", "kg/m2"),
               c("Weight, kg", "weight_kg", "kg"),
               c("SOFA, first window", "sofa_total", "points"),
               c("Norepinephrine equivalent", "nee", "mcg/kg/min"),
               c("P/F ratio, first window", "oxygenation", "mmHg"),
               c("Lactate", "lactate", "mmol/L"),
               c("First IMV episode, hours", "first_imv_episode_hours", "hours"))) {
  POOL[[length(POOL) + 1]] <- pool_row("baseline", v[1], v[3], "landmark",
                                       NA, NA, base[[v[2]]], MIN_CELL)
}
for (drug in names(DRUGS)) {
  col <- DRUGS[[drug]]
  for (w in sort(unique(lm_long$window_idx))) {
    x <- lm_long[lm_long$window_idx == w, ]      # balanced panel, as above
    hr <- x$window_start_hr[1]
    POOL[[length(POOL) + 1]] <- pool_row(
      "by_window", sprintf("%s dose", drug), UNITS[[drug]], "balanced_panel",
      w, hr, x[[col]], MIN_CELL)
    POOL[[length(POOL) + 1]] <- pool_row(
      "by_window", sprintf("%s dose", drug), UNITS[[drug]], "receivers_only",
      w, hr, x[[col]][!is.na(x[[col]]) & x[[col]] > 0], MIN_CELL)
  }
}
pooling_continuous <- do.call(rbind, POOL)

pooling_categorical <- rbind(
  pool_cat("Sex", ifelse(is.na(base$sex), "Missing", base$sex), one, MIN_CELL),
  pool_cat("Race", ifelse(is.na(base$race), "Missing", base$race), one, MIN_CELL),
  pool_cat("Race (collapsed)",
           collapse_levels(ifelse(is.na(base$race), "Missing", base$race),
                           RACE_COLLAPSE), one, MIN_CELL))
pooling_categorical <- pooling_categorical[pooling_categorical$stratum != "overall", ]


# ---- 14. The analytic table Phases 3-6 read ----------------------------------
# id_num is a dense integer rank of encounter_block: gbmt and lcmm want a numeric
# unit. Both id columns ride along so the blocks-per-patient choice stays open.

lm_out <- lm_long[order(lm_long$encounter_block, lm_long$window_idx), ]
lm_out$id_num <- as.integer(factor(lm_out$encounter_block))
stopifnot(
  "every landmark episode must contribute every window in [0,T)" =
    nrow(lm_out) == length(elig) * (LANDMARK / WINDOW_H),
  "id_num must be a dense rank of encounter_block" =
    max(lm_out$id_num) == length(elig)
)
write_parquet(lm_out, file.path(dirs$out_phi, "landmark_cohort.parquet"))
cat(sprintf("\nwritten: landmark_cohort.parquet  %s rows x %d cols (PHI)\n",
            format(nrow(lm_out), big.mark = ","), ncol(lm_out)))


# ---- 15. Write ---------------------------------------------------------------

write_out <- function(x, name) {
  f <- file.path(dirs$phase, name)
  write.csv(x, f, row.names = FALSE)
  cat(sprintf("written: %s\n", name))
}

write_out(flow, "landmark_flow.csv")
write_out(sens, "T_sensitivity.csv")
write_out(failed, "failed_extubation.csv")
write_out(curve, "dose_curve.csv")
write_out(dep, "dependence.csv")
write_out(pooling_continuous, "pooling_continuous.csv")
write_out(pooling_categorical, "pooling_categorical.csv")

write_captions(file.path(dirs$phase, "captions.md"), "05_landmark_cohort.R",
               grep("\\.png$", OWNED$phase, value = TRUE), prov)
cat("written: captions.md\n")

write_json(prov, file.path(dirs$phase, "provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)
cat("written: provenance.json\n")


# ---- 16. Provenance ----------------------------------------------------------
# Which package versions produced these numbers?

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/05_landmark_cohort.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "05_landmark_cohort_sessioninfo.txt")
)
