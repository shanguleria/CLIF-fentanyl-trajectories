# ==============================================================================
# 07_outcome_model.R  --  per-admission titration exposures vs ventilator-free days
#
# Purpose : Two per-admission exposures over the first 72h -- the share of fentanyl
#           increases paired with a bolus (Pbolus) and the share with a documented
#           pain or agitation score (Pindicated) -- against VFD-28, by
#           proportional-odds ordinal regression with ICU fixed effects.
# Author  : Shan Guleria
# Created : 2026-10-27
# Inputs  : output/intermediate_phi/ : titration_adherence.parquet,
#           indication_events.parquet, time_to_event.parquet,
#           trajectory_long.parquet, care_intervals.parquet
# Outputs : output/final_no_phi/07_outcome_model/ : the files listed in OWNED
#           output/intermediate_phi/outcome_model_frame.parquet (PHI, per admission)
# ==============================================================================

# Run this in a FRESH R session (RStudio: Cmd+Shift+F10).

# ---- 1. Packages -------------------------------------------------------------

# No metafor: pooling is coordinating-centre code. See protocol.json _pooling_why.
pkgs <- c("here", "jsonlite", "ggplot2", "arrow", "ordinal")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

source(here("code", "utils", "paths.R"))
source(here("code", "utils", "figures.R"))
source(here("code", "utils", "pooling.R"))
source(here("code", "utils", "table1.R"))


# ---- 2. Config ---------------------------------------------------------------

config <- fromJSON(here("config", "config.json"), simplifyVector = FALSE)
COV <- fromJSON(here("config", "protocol.json"), simplifyVector = FALSE)
SEED <- COV$model$seed
set.seed(SEED)

SPEC <- COV$model$outcome_model
stopifnot(
  "protocol.json must declare model.outcome_model" = !is.null(SPEC),
  "this script implements proportional-odds ordinal logistic regression only" =
    identical(SPEC$family, "ordinal_logistic_proportional_odds"),
  "the ICU attribution rule must be the interval containing hour 0" =
    identical(SPEC$icu_attribution, "interval_containing_hour_0"),
  "only complete-case is implemented" = identical(SPEC$missing_data, "complete_case"))

EXPOSURES  <- unlist(SPEC$exposures)
COVARIATES <- unlist(SPEC$covariates)
MIN_UNIT_N <- as.integer(SPEC$min_admissions_per_unit)
SENS_MIN   <- as.integer(SPEC$sensitivity_min_increase_events)

# Reporting unit: MODEL terms are the share / EXP_SCALE and say so in their names;
# tables and figures keep the [0, 1] share. See _exposure_scale_why.
EXP_SCALE <- as.numeric(SPEC$exposure_scale_per)
stopifnot(
  "protocol.json must declare model.outcome_model.exposure_scale_per" =
    length(EXP_SCALE) == 1 && is.finite(EXP_SCALE),
  "the exposure scale must be a positive fraction of the [0,1] share" =
    EXP_SCALE > 0 && EXP_SCALE <= 1)
EXP_SCALE_PP <- round(EXP_SCALE * 100)
MODEL_EXPOSURES <- sprintf("%s_per%dpp", EXPOSURES, EXP_SCALE_PP)
LANDMARK_H <- COV$cohort$landmark_hours
MIN_CELL   <- COV$reporting$small_cell_min_den

# The ICU-only sensitivity restricts on the mCIDE location_category.
ICU_CAT <- tolower(as.character(SPEC$sensitivity_icu_only_location_category))
ATTR_CATS <- tolower(unlist(COV$unit_variation$attributed_location_categories))
stopifnot(
  "protocol.json must declare sensitivity_icu_only_location_category" =
    length(ICU_CAT) == 1 && nzchar(ICU_CAT),
  "the ICU-only category must be one the study attributes" = ICU_CAT %in% ATTR_CATS)

# Descriptive table and the four_band sensitivity only, not the model outcome.
CUTS   <- as.numeric(unlist(SPEC$vfd_ordinal_breaks))
LABELS <- unlist(SPEC$vfd_ordinal_labels)
stopifnot(
  "protocol.json must declare vfd_ordinal_breaks and vfd_ordinal_labels" =
    length(CUTS) > 1 && length(LABELS) == length(CUTS) - 1,
  "the breaks must be increasing" = !is.unsorted(CUTS, strictly = TRUE))

# One ordinal level per whole day; see _outcome_scale_why and _day_rounding_why.
OUTCOME_SCALE <- as.character(SPEC$outcome_scale)
DAY_ROUNDING  <- as.character(SPEC$day_rounding)
BAND_VALUES   <- as.character(SPEC$band_values)
CUTS4   <- as.numeric(unlist(SPEC$four_band_breaks))
LABELS4 <- unlist(SPEC$four_band_labels)
stopifnot(
  "this script implements outcome_scale = declared_bands only" =
    identical(OUTCOME_SCALE, "declared_bands"),
  "this script implements day_rounding = floor only" =
    identical(DAY_ROUNDING, "floor"),
  "only band_values = observed_within_band_mean is implemented" =
    identical(BAND_VALUES, "observed_within_band_mean"),
  "protocol.json must declare four_band_breaks and four_band_labels" =
    length(CUTS4) > 1 && length(LABELS4) == length(CUTS4) - 1,
  "the four-band breaks must be increasing" =
    !is.unsorted(CUTS4, strictly = TRUE))

# Declared with the UNSCALED exposure names and matched to the model terms below,
# so the gate cannot drift from SPEC$exposures. Enforced at the end of the script.
GATE       <- SPEC$proportional_odds_gate
GATE_TERMS <- unlist(GATE$terms)
GATE_ALPHA <- as.numeric(GATE$alpha)
stopifnot(
  "protocol.json must declare model.outcome_model.proportional_odds_gate" =
    !is.null(GATE),
  "the gate must name terms that are declared exposures" =
    length(GATE_TERMS) > 0 && all(GATE_TERMS %in% EXPOSURES),
  "the gate alpha must be a probability" =
    length(GATE_ALPHA) == 1 && GATE_ALPHA > 0 && GATE_ALPHA < 1,
  "only proportional_odds_gate.on_failure = fail is implemented" =
    identical(as.character(GATE$on_failure), "fail"))
GATE_MODEL_TERMS <- MODEL_EXPOSURES[match(GATE_TERMS, EXPOSURES)]

PRED          <- SPEC$prediction
PRED_EXPOSURE <- as.character(PRED$exposure)
PRED_GRID_PCT <- as.numeric(unlist(PRED$exposure_grid_pct))
PRED_DRAWS    <- as.integer(PRED$bootstrap_draws)
BIN_BREAKS    <- as.numeric(unlist(PRED$unadjusted_bin_breaks))
BIN_LABELS    <- unlist(PRED$unadjusted_bin_labels)
stopifnot(
  "protocol.json must declare model.outcome_model.prediction" = !is.null(PRED),
  "the prediction exposure must be one of the declared exposures" =
    length(PRED_EXPOSURE) == 1 && PRED_EXPOSURE %in% EXPOSURES,
  "the prediction grid must be percentages in [0, 100]" =
    length(PRED_GRID_PCT) > 1 && all(PRED_GRID_PCT >= 0 & PRED_GRID_PCT <= 100),
  "the prediction grid must be increasing" =
    !is.unsorted(PRED_GRID_PCT, strictly = TRUE),
  "only prediction.method = standardization is implemented" =
    identical(as.character(PRED$method), "standardization"),
  "bootstrap_draws must be a positive count" =
    length(PRED_DRAWS) == 1 && is.finite(PRED_DRAWS) && PRED_DRAWS > 0,
  "one unadjusted bin label per interval" =
    length(BIN_LABELS) == length(BIN_BREAKS) - 1,
  "the unadjusted bin breaks must be increasing" =
    !is.unsorted(BIN_BREAKS, strictly = TRUE),
  "the first bin break must sit below zero so exactly 0 gets its own bin" =
    BIN_BREAKS[1] < 0)
PRED_MODEL_EXPOSURE <- MODEL_EXPOSURES[match(PRED_EXPOSURE, EXPOSURES)]


# ---- 3. Paths and provenance -------------------------------------------------

dirs <- site_dirs()
dirs$phase <- phase_dir(dirs, "07_outcome_model")
dirs$phi <- phi_dir(dirs, "07_outcome_model")   # patient-level, same subdivision
dirs$diagnostics <- phase_dir(dirs, file.path("07_outcome_model", "diagnostics"))
prov <- provenance(config)


# ---- 4. Guards ---------------------------------------------------------------

manifest <- require_manifest(dirs, here())
message(sprintf("  reading cohort-build outputs from code %s, generated %s",
                manifest$code_version, manifest$generated))

OWNED <- list(
  phi = c("outcome_model_frame.parquet"),
  # Reported directly.
  phase = c("exposure_distribution.csv", "exposure_by_unit.csv",
            "exposure_correlation.csv", "vfd_ordinal.csv",
            "table1.csv", "table1_continuous.csv", "table1_categorical.csv",
            "vfd_day_distribution.csv", "adjusted_prediction.csv",
            "unadjusted_bins.csv",
            "model_coefficients.csv", "model_covariance.csv", "strobe.csv",
            "exposure_histogram.png", "exposure_by_unit.png",
            "exposure_scatter.png", "vfd_vs_pindicated.png",
            "provenance.json", "captions.md"),
  # Method checks: evidence about whether the model is admissible, not results.
  diagnostics = c("model_diagnostics.csv", "proportional_odds_test.csv",
                  "proportional_odds_gate.csv"))
# Flat intermediate_phi locations, before it was subdivided by writing script
# on 2026-10-30. Declared so a stale twin is deleted, not read as current.
RETIRED_PHI <- file.path("output", "intermediate_phi",
                         c("outcome_model_frame.parquet"))

RETIRED <- character(0)
n_cleared <- clear_owned_outputs(dirs, OWNED, retired = c(RETIRED, RETIRED_PHI))
if (n_cleared) message(sprintf("  cleared %d output(s) from a previous run", n_cleared))

FLOW <- list()
flow <- function(label, n, dropped_because = NA_character_) {
  FLOW[[length(FLOW) + 1]] <<- data.frame(
    step = label, n = n, dropped_because = dropped_because,
    stringsAsFactors = FALSE)
  cat(sprintf("  %-58s %7s\n", label, format(n, big.mark = ",")))
}


# ---- 5. Read -----------------------------------------------------------------

# 05's `overall` scope IS the first-72h exposure, so no 72h filter is applied here.
adh <- as.data.frame(read_parquet(file.path(phi_dir(dirs, "05_titration"), "titration_adherence.parquet")))
stopifnot("titration_adherence must carry the overall scope" =
            "overall" %in% adh$scope)
pb <- adh[adh$scope == "overall", c("encounter_block", "n_events", "n_paired")]
names(pb)[names(pb) == "n_events"] <- "n_increase"

ind <- as.data.frame(read_parquet(file.path(phi_dir(dirs, "05_titration"), "indication_events.parquet")))
stopifnot(
  "indication events must be bounded by the granular extent" =
    max(ind$t_hr) <= COV$cohort$granular_extent_hours + 1e-9,
  "indication events must start at or after the anchor" = min(ind$t_hr) >= 0)

tte <- as.data.frame(read_parquet(
  file.path(phi_dir(dirs, "01_cohort"), "time_to_event.parquet"),
  col_select = c("encounter_block", "vfd_28", "landmark_eligible",
                 "extubation_event", "died")))

base <- as.data.frame(read_parquet(
  file.path(phi_dir(dirs, "01_cohort"), "trajectory_long.parquet"),
  col_select = c("encounter_block", "window_idx", "age", "sex", "cci",
                 "sofa_baseline", "sofa_baseline_n_components",
                 # Table 1 only -- 07 is its sole generator as of 0.17.0.
                 "race", "first_imv_episode_hours")))
base <- base[base$window_idx == 0, setdiff(names(base), "window_idx")]
stopifnot("trajectory_long must carry sofa_baseline -- re-run code/01_build_cohort.py" =
            "sofa_baseline" %in% names(base),
          "window 0 must hold each block once" = !anyDuplicated(base$encounter_block))

icu <- as.data.frame(read_parquet(file.path(phi_dir(dirs, "01_cohort"), "care_intervals.parquet")))


# ---- 6. The two exposures ----------------------------------------------------
# A share of INCREASE events; the shared denominator is asserted, not assumed.

inc <- ind[ind$event_type == "increase", ]
pi_ <- aggregate(cbind(n_inc = rep(1, nrow(inc)),
                       n_indicated = as.integer(inc$any_indicated)) ~ encounter_block,
                 data = inc, FUN = sum)
ups <- aggregate(list(n_uptitration = as.integer(inc$kind == "uptitration")),
                 by = list(encounter_block = inc$encounter_block), FUN = sum)

ex <- merge(pb, pi_, by = "encounter_block")
ex <- merge(ex, ups, by = "encounter_block")
stopifnot(
  "the two handoffs disagree about the number of increase events" =
    all(ex$n_increase == ex$n_inc),
  "a numerator exceeds its denominator" =
    all(ex$n_paired <= ex$n_increase) && all(ex$n_indicated <= ex$n_increase))
ex$Pbolus     <- ex$n_paired / ex$n_increase
ex$Pindicated <- ex$n_indicated / ex$n_increase
ex$n_inc <- NULL

cat(sprintf("\nAdmissions with a fentanyl infusion and at least one increase: %s\n",
            format(nrow(ex), big.mark = ",")))


# ---- 7. The analysis set -----------------------------------------------------
# Three restrictions, each a STROBE row. The landmark one is not optional: VFD-28
# is defined only on IMV at T, and it puts the exposure window before the outcome.

flow("Admissions with >= 1 fentanyl increase event", nrow(ex))

d <- merge(ex, tte[tte$landmark_eligible %in% TRUE, ], by = "encounter_block")
flow(sprintf("Landmark cohort (VFD-28 defined, on IMV at %dh)", LANDMARK_H),
     nrow(d), "not in the landmark cohort")

# Table 1's middle scope, captured BEFORE the up-titration restriction: this is
# "on IMV at 72h with at least one fentanyl start or increase", which is NOT the
# frame the regression fits. Covariates are merged on here rather than held by
# reference, because the columns 07 adds later would not reach a stale copy.
d_landmark <- merge(d, base, by = "encounter_block", all.x = TRUE)
d <- d[d$n_uptitration >= 1, ]
flow("With >= 1 UPTITRATION (not initiation alone)", nrow(d),
     "initiation only -- no titration decision to describe")

# ICU at hour 0, not the modal unit over the window: see _icu_attribution_why.
icu <- icu[icu$t_start_hr <= 0 & icu$t_end_hr > 0, ]
# location_category is the mCIDE key the ICU-only sensitivity restricts on.
# care_setting cannot serve: it gives the subtype ("which ICU"), not "is it one".
icu <- icu[!duplicated(icu$encounter_block),
           c("encounter_block", "location_name", "location_category")]
d <- merge(d, icu, by = "encounter_block", all.x = TRUE)

# `unattributed` is a real group, kept as its own level and labelled on BOTH keys
# so the ICU-only restriction excludes it explicitly, not by an NA comparison.
d$location_name[is.na(d$location_name)] <- "unattributed"
d$location_category[is.na(d$location_category)] <- "unattributed"
d$icu <- factor(d$location_name)
cat(sprintf("  ICU levels before collapsing: %d\n", nlevels(d$icu)))

# Care-category composition: the fixed effects are named "ICU" but the largest
# level need not be one, and this is the ICU-only sensitivity's denominator.
comp <- as.data.frame(table(location_category = tolower(d$location_category)),
                      stringsAsFactors = FALSE)
comp$pct <- round(100 * comp$Freq / nrow(d), 1)
cat("  Care category at hour 0:\n")
for (i in seq_len(nrow(comp))) {
  cat(sprintf("    %-14s %6s  %5.1f%%%s\n", comp$location_category[i],
              format(comp$Freq[i], big.mark = ","), comp$pct[i],
              if (identical(comp$location_category[i], ICU_CAT))
                "   <- the ICU-only sensitivity" else ""))
}
# The unit collapse happens PER MODEL FRAME: each model drops a different set of
# incomplete cases, so it must be sized on the frame fitted. See collapse_units().

# Baseline covariates, at or before intubation. sofa_baseline is a cohort-build product
# (protocol.json time_invariant.sofa_baseline), not derived here.
d <- merge(d, base, by = "encounter_block", all.x = TRUE)
stopifnot("the covariate merge changed the admission count" =
            nrow(d) == nrow(unique(d["encounter_block"])),
          "every model covariate must be present after the merge" =
            all(c("age", "sex", "cci", "sofa_baseline") %in% names(d)))
d$sex <- factor(d$sex)

# Three outcome columns. vfd_ord is the PRIMARY: the declared bands. vfd_band4
# and vfd_day are the two scale sensitivities. cut() is right-closed, so the first
# break sits BELOW zero to give exactly 0 its own level. Ceiling derived, not typed.
VFD_CEILING_D <- as.numeric(COV$outcomes$vfd$horizon_days) - LANDMARK_H / 24
stopifnot("the VFD ceiling must be a positive whole number of days" =
            VFD_CEILING_D > 0 && abs(VFD_CEILING_D - round(VFD_CEILING_D)) < 1e-9)

d$vfd_ord   <- cut(d$vfd_28, breaks = CUTS,  labels = LABELS,  ordered_result = TRUE)
d$vfd_band4 <- cut(d$vfd_28, breaks = CUTS4, labels = LABELS4, ordered_result = TRUE)
day <- floor(d$vfd_28 + 1e-9)
d$vfd_day <- factor(day, levels = sort(unique(day)), ordered = TRUE)
stopifnot(
  "a VFD-28 value fell outside the declared breaks" = !any(is.na(d$vfd_ord)),
  "a VFD-28 value fell outside the four-band breaks" = !any(is.na(d$vfd_band4)),
  "the model outcome must be ordered" = is.ordered(d$vfd_ord),
  "the four-band sensitivity outcome must be ordered" = is.ordered(d$vfd_band4),
  "a VFD-28 value produced no integer day" = !any(is.na(d$vfd_day)),
  "flooring must not invent a day above the VFD ceiling" =
    max(day) <= VFD_CEILING_D)

vfd_tab <- as.data.frame(table(level = d$vfd_ord), stringsAsFactors = FALSE)
names(vfd_tab)[2] <- "n_admissions"
vfd_tab$denominator <- nrow(d)
vfd_tab$pct <- round(100 * vfd_tab$n_admissions / nrow(d), 1)
# The band -> days mapping for the adjusted expectation, measured here and
# SHIPPED, because a band carries no days of its own. See protocol.json
# band_values. A band with no admissions has no observed mean; the model cannot
# fit an empty level either, so both are reported rather than imputed.
vfd_tab$mean_vfd_28 <- round(as.numeric(tapply(d$vfd_28, d$vfd_ord, mean)), 4)
cat("\nVFD-28 as the PRIMARY ordinal outcome\n")
print(vfd_tab, row.names = FALSE)
if (any(vfd_tab$n_admissions < MIN_CELL)) {
  cat(sprintf("  NOTE %d band(s) below %d admissions. A thin band does not ",
              sum(vfd_tab$n_admissions < MIN_CELL), MIN_CELL),
      "destabilise the fit -- proportional-odds thresholds are estimated from\n",
      "  CUMULATIVE splits, not cell counts -- but it is thin to report.\n", sep = "")
}

# One row per attainable day. The empty levels are information: a gap at day 1 is
# the liberation threshold showing through, not missing data.
day_tab <- data.frame(day = seq.int(0L, as.integer(round(VFD_CEILING_D))),
                      stringsAsFactors = FALSE)
day_tab$n_admissions <- as.integer(table(factor(day, levels = day_tab$day)))
day_tab$denominator <- nrow(d)
day_tab$pct <- round(100 * day_tab$n_admissions / nrow(d), 2)
day_tab$is_model_level <- day_tab$day %in% as.numeric(levels(d$vfd_day))
cat(sprintf("VFD-28 by whole day: %d of %d day levels observed (absent: %s)\n",
            nlevels(d$vfd_day), nrow(day_tab),
            paste(day_tab$day[!day_tab$is_model_level], collapse = ", ")))

# Scaling the COLUMN, not the coefficient afterwards, keeps estimate, SE, CI and
# vcov on one scale -- a variance scales by the SQUARE of the factor.
for (i in seq_along(EXPOSURES)) {
  d[[MODEL_EXPOSURES[i]]] <- d[[EXPOSURES[i]]] / EXP_SCALE
}
stopifnot("the scaled exposures must reach the model frame" =
            all(MODEL_EXPOSURES %in% names(d)))
cat(sprintf("\nExposure reporting unit: one OR per %d percentage points (%s)\n",
            EXP_SCALE_PP, paste(MODEL_EXPOSURES, collapse = ", ")))

# Complete case, derived from the config so a new covariate is not a silent no-op.
NEED <- c(EXPOSURES, COVARIATES)
d$.complete <- stats::complete.cases(d[, NEED])
flow("Complete on every model covariate", sum(d$.complete),
     "missing a covariate -- sofa_baseline is the binding one")
for (v in NEED) {
  miss <- sum(is.na(d[[v]]))
  if (miss) cat(sprintf("    %-28s missing %s (%.1f%%)\n", v,
                        format(miss, big.mark = ","), 100 * miss / nrow(d)))
}


# ---- 8. Exposure description -------------------------------------------------
# BINNED COUNTS, not a density: counts on fixed bins pool across sites by addition.
BINS <- seq(0, 1, by = 0.05)
dist_rows <- do.call(rbind, lapply(EXPOSURES, function(v) {
  b <- cut(d[[v]], breaks = BINS, include.lowest = TRUE, right = TRUE)
  tb <- as.data.frame(table(bin = b), stringsAsFactors = FALSE)
  data.frame(exposure = v, bin = tb$bin,
             bin_lo = BINS[-length(BINS)], bin_hi = BINS[-1],
             n = tb$Freq, n_total = nrow(d),
             pct = round(100 * tb$Freq / nrow(d), 2), stringsAsFactors = FALSE)
}))

summ <- do.call(rbind, lapply(EXPOSURES, function(v) {
  x <- d[[v]]
  data.frame(exposure = v, n_admissions = length(x),
             median = round(median(x), 3),
             q1 = round(quantile(x, .25, names = FALSE), 3),
             q3 = round(quantile(x, .75, names = FALSE), 3),
             n_at_zero = sum(x == 0), n_at_one = sum(x == 1),
             pct_at_zero = round(100 * mean(x == 0), 1),
             pct_at_one = round(100 * mean(x == 1), 1),
             stringsAsFactors = FALSE)
}))
cat("\nExposure distributions (per admission)\n")
print(summ, row.names = FALSE)
cat(sprintf("  shared denominator (increase events): median %.0f, IQR %.0f-%.0f, max %.0f\n",
            median(d$n_increase), quantile(d$n_increase, .25),
            quantile(d$n_increase, .75), max(d$n_increase)))

# Per unit, for the violin figure. Quantiles do NOT pool; the counts beside them
# do. Units below the reporting floor are suppressed, not published.
by_unit <- do.call(rbind, lapply(EXPOSURES, function(v) {
  do.call(rbind, lapply(levels(d$icu), function(u) {
    x <- d[[v]][d$icu == u]
    if (!length(x)) return(NULL)
    # location_category SHIPS BESIDE THE NAME: location_name is site-local free
    # text, so another site cannot tell an ICU from an ED from a ward.
    cat_u <- unique(tolower(d$location_category[d$icu == u]))
    stopifnot("a location_name maps to more than one location_category" =
                length(cat_u) == 1)
    data.frame(exposure = v, location_name = u, location_category = cat_u,
               n_admissions = length(x),
               # The POINT MASS, as a count so it pools and the label is checkable.
               n_at_zero = sum(x == 0),
               median = round(median(x), 3),
               q1 = round(quantile(x, .25, names = FALSE), 3),
               q3 = round(quantile(x, .75, names = FALSE), 3),
               mean = round(mean(x), 3), sd = round(sd(x), 3),
               stringsAsFactors = FALSE)
  }))
}))
suppressed <- by_unit$n_admissions < MIN_CELL
if (any(suppressed)) {
  cat(sprintf("  suppressing %d unit-exposure cell(s) below reporting floor %d\n",
              sum(suppressed), MIN_CELL))
  by_unit[suppressed, c("n_at_zero", "median", "q1", "q3", "mean", "sd")] <- NA
}

# Spearman as well as Pearson: both exposures are proportions massed at zero.
ct_p <- stats::cor.test(d$Pbolus, d$Pindicated, method = "pearson")
ct_s <- suppressWarnings(stats::cor.test(d$Pbolus, d$Pindicated, method = "spearman",
                                         exact = FALSE))
corr <- data.frame(
  n_admissions = nrow(d),
  pearson_r = round(unname(ct_p$estimate), 4),
  pearson_lo = round(ct_p$conf.int[1], 4), pearson_hi = round(ct_p$conf.int[2], 4),
  pearson_p = signif(ct_p$p.value, 3),
  spearman_rho = round(unname(ct_s$estimate), 4),
  spearman_p = signif(ct_s$p.value, 3),
  vif_two_predictor = round(1 / (1 - unname(ct_p$estimate)^2), 3),
  stringsAsFactors = FALSE)
cat(sprintf("\nCollinearity: Pearson r = %+.3f (%.3f to %.3f), Spearman rho = %+.3f, VIF = %.2f\n",
            corr$pearson_r, corr$pearson_lo, corr$pearson_hi,
            corr$spearman_rho, corr$vif_two_predictor))
cat("  VIF below 5 means both exposures can enter one model; above it, split them.\n")


# ---- 9. Figures --------------------------------------------------------------

EXP_LAB <- c(Pbolus = "P(bolus) -- share of increases paired with a bolus",
             Pindicated = "P(indicated) -- share of increases with a qualifying score")
EXP_COL <- c(Pbolus = STATE_COLOURS[["bolus only"]],
             Pindicated = STATE_COLOURS[["continuous only"]])

long_ex <- do.call(rbind, lapply(EXPOSURES, function(v)
  data.frame(exposure = v, value = d[[v]], icu = d$icu, stringsAsFactors = FALSE)))
long_ex$exposure <- factor(long_ex$exposure, levels = EXPOSURES,
                           labels = unname(EXP_LAB[EXPOSURES]))

register_caption("exposure_histogram.png",
  "How often is a fentanyl increase paired with a bolus, or justified by a score?",
  sprintf(paste0("Distribution across %s admissions of the two per-admission ",
                 "exposures, each the share of INCREASE events (initiation union ",
                 "up-titration) in the first %gh meeting its condition. Bin width ",
                 "0.05. BOTH DISTRIBUTIONS ARE DOMINATED BY ZERO -- %.1f%% of ",
                 "admissions never pair a bolus with an increase and %.1f%% never ",
                 "have a qualifying score at one. The shared denominator is small ",
                 "(median %.0f increases, IQR %.0f-%.0f), so a proportion on two or ",
                 "three events can only take a few values; the sensitivity analysis ",
                 "restricted to >= %d events is the check on that."),
          format(nrow(d), big.mark = ","), LANDMARK_H,
          summ$pct_at_zero[summ$exposure == "Pbolus"],
          summ$pct_at_zero[summ$exposure == "Pindicated"],
          median(d$n_increase), quantile(d$n_increase, .25),
          quantile(d$n_increase, .75), SENS_MIN))

p_hist <- house(
  ggplot(long_ex, aes(value, fill = exposure)) +
    geom_histogram(breaks = BINS, colour = PAPER, linewidth = 0.2) +
    facet_wrap(~exposure, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = unname(EXP_COL[EXPOSURES]), guide = "none") +
    scale_x_continuous(limits = c(-0.025, 1.025),
                       breaks = seq(0, 1, 0.25), labels = scales::percent_format(1)) +
    labs(x = "Share of increase events", y = "Admissions"))
ggsave(file.path(dirs$phase, "exposure_histogram.png"), p_hist,
       width = 7.0, height = 5.2, dpi = 200)

# Violin by unit, with the n under each: a wide violin on twelve admissions is not
# the same claim as one on four hundred.
unit_n <- as.data.frame(table(icu = d$icu), stringsAsFactors = FALSE)
names(unit_n)[2] <- "n"
long_ex$icu_lab <- factor(long_ex$icu, levels = unit_n$icu,
                          labels = sprintf("%s\n(n = %s)", unit_n$icu,
                                           format(unit_n$n, big.mark = ",")))

register_caption("exposure_by_unit.png",
  "Does titration practice differ between ICUs, or only between patients?",
  sprintf(paste0("Per-admission exposure by the ICU the patient was INTUBATED IN ",
                 "(the attributed care location containing hour 0, not the modal ",
                 "unit over the exposure window, which would be partly determined ",
                 "by the exposure). ONE POINT PER ADMISSION, jittered ",
                 "HORIZONTALLY ONLY so that a value of exactly zero is drawn AT ",
                 "zero -- both exposures have a point mass on that bound, and ",
                 "vertical jitter would read as documentation that did not ",
                 "happen. Diamond is the median; the percentage above each unit ",
                 "is the share of its admissions at exactly zero. POINTS ARE ",
                 "SINA-POSITIONED, offset by the LOCAL density on the same ",
                 "width normalisation as the violin outline, so every admission ",
                 "sits inside the envelope; uniform jitter ignores the envelope ",
                 "and sprays points past the waist above the zero spike. The ",
                 "density is drawn at the DEFAULT bandwidth -- smoothing it ",
                 "further made the panel tidier and misleading, rendering a unit ",
                 "with 61%% of admissions at exactly zero as an almost uniform ",
                 "distribution. NO BOX: with half of this variable at one value ",
                 "the median and Q1 coincide on the floor for most units. The ",
                 "zero share is LABELLED rather than drawn because neither a ",
                 "violin nor overlapping dots distinguish 50%% from 80%%, and ",
                 "it ships as n_at_zero beside n_admissions so the label can be ",
                 "recomputed. A UNIT NEEDS AT LEAST %d ADMISSIONS TO APPEAR ",
                 "HERE AT ALL and is otherwise absent from the panel: below ",
                 "that the point, the median and the zero share are all the ",
                 "same one or two patients' exposure against a named unit. Such ",
                 "units still ship in exposure_by_unit.csv with their n and ",
                 "their statistics blanked, so the count is auditable even ",
                 "though the unit is not drawn. HORIZONTAL ",
                 "BANDS ARE REAL: the median admission contributes 4 increase ",
                 "events, so attainable shares form a sparse rational grid and ",
                 "the banding is the small-denominator limitation, not a ",
                 "rendering artifact. n under each label. ",
                 "UNITS ARE SHOWN AS CHARTED HERE -- the collapsing of units ",
                 "below %d admissions into 'other' happens inside each MODEL, on ",
                 "the frame actually fitted, because complete-case filtering ",
                 "changes which units are too small. READ THE WITHIN-UNIT SPREAD ",
                 "FIRST: if it dwarfs the between-unit differences, unit is not ",
                 "where the variation lives. Units carry their real names here; ",
                 "audience-facing versions anonymise them."),
          MIN_CELL, MIN_UNIT_N))

# One point per admission, not a box: both exposures are massed at exactly zero, so
# median and Q1 coincide on the floor for most units. The offset is HORIZONTAL ONLY
# -- a true zero moved vertically reads as documentation that did not happen. Points
# are SINA-POSITIONED, not jittered: a scale = "width" violin narrows to a waist
# above the zero spike, so each is offset by the LOCAL density. adjust stays at 1.
long_ex$n_unit <- unit_n$n[match(long_ex$icu, unit_n$icu)]

# A unit needs MIN_CELL admissions to appear at all, and droplevels is what removes
# the axis slot. It still ships in exposure_by_unit.csv with its statistics blanked.
pts <- droplevels(long_ex[long_ex$n_unit >= MIN_CELL, ])

VIOLIN_W <- 0.90
sina_offset <- function(v) {
  if (length(v) < 2) return(rep(0, length(v)))
  de <- stats::density(v, adjust = 1, from = min(v), to = max(v))
  w <- stats::approx(de$x, de$y, xout = v, rule = 2)$y
  w <- w / max(w, na.rm = TRUE)          # the scale = "width" normalisation
  stats::runif(length(v), -1, 1) * w * (VIOLIN_W / 2) * 0.95
}
pts$xoff <- stats::ave(pts$value, interaction(pts$exposure, pts$icu_lab, drop = TRUE),
                       FUN = sina_offset)
pts$xnum <- as.numeric(pts$icu_lab) + pts$xoff

# Derived FROM by_unit, not recomputed: the panel and the CSV cannot then disagree.
agg <- by_unit[by_unit$n_admissions >= MIN_CELL, ]
agg$icu_lab <- factor(agg$location_name, levels = unit_n$icu,
                      labels = levels(long_ex$icu_lab))
# Re-level onto the units the panel draws, or the diamonds land on a stale level set.
agg$icu_lab <- factor(as.character(agg$icu_lab), levels = levels(pts$icu_lab))
agg$exposure <- factor(agg$exposure, levels = EXPOSURES,
                       labels = unname(EXP_LAB[EXPOSURES]))
agg$med <- agg$median
# Integer-first: 100 * (n/d) can print a percent below (100*n)/d recomputed from
# the shipped counts, and the panel number has to match the CSV.
agg$p0  <- (100 * agg$n_at_zero) / agg$n_admissions

p_unit <- house(
  ggplot(pts, aes(icu_lab, value)) +
    geom_violin(data = pts, aes(fill = exposure), scale = "width",
                width = VIOLIN_W, adjust = 1, trim = TRUE,
                colour = INK, linewidth = 0.35, alpha = 0.70) +
    geom_point(data = pts, aes(x = xnum), size = 0.45, alpha = 0.22,
               shape = 16, colour = INK) +
    geom_point(data = agg, aes(y = med, fill = exposure), shape = 23,
               size = 2.0, colour = INK, stroke = 0.4) +
    geom_text(data = agg, aes(y = 1.10, label = sprintf("%.0f%% at 0", p0)),
              size = 2.3, colour = INK) +
    facet_wrap(~exposure, ncol = 1) +
    scale_fill_manual(values = unname(EXP_COL[EXPOSURES]), guide = "none") +
    # Headroom for the "% at 0" row: ggplot DROPS data outside the scale limits.
    scale_y_continuous(limits = c(0, 1.14), breaks = seq(0, 1, 0.25),
                       labels = scales::percent_format(1)) +
    labs(x = NULL, y = "Share of increase events") +
    theme(axis.text.x = element_text(size = 8, colour = INK)))
ggsave(file.path(dirs$phase, "exposure_by_unit.png"), p_unit,
       width = 9.0, height = 6.2, dpi = 200)

# Jittered because both axes are proportions on a small denominator.
register_caption("exposure_scatter.png",
  "Are the two exposures measuring the same thing?",
  sprintf(paste0("P(bolus) against P(indicated), one point per admission, ",
                 "jittered because both are proportions on a median denominator ",
                 "of %.0f and many admissions share exact coordinates. Pearson ",
                 "r = %+.3f (%.3f to %.3f), Spearman rho = %+.3f, two-predictor ",
                 "VIF %.2f -- WEAKLY correlated, so both enter one model rather ",
                 "than being fitted separately. The solid line is an ordinary ",
                 "least-squares fit shown for orientation only; neither variable ",
                 "is an outcome of the other."),
          median(d$n_increase), corr$pearson_r, corr$pearson_lo,
          corr$pearson_hi, corr$spearman_rho, corr$vif_two_predictor))

p_sc <- house(
  ggplot(d, aes(Pindicated, Pbolus)) +
    geom_jitter(width = 0.018, height = 0.018, size = 0.7, alpha = 0.22,
                colour = INK) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = STATE_COLOURS[["bolus only"]], linewidth = 0.7) +
    scale_x_continuous(limits = c(-0.04, 1.04), breaks = seq(0, 1, 0.25),
                       labels = scales::percent_format(1)) +
    scale_y_continuous(limits = c(-0.04, 1.04), breaks = seq(0, 1, 0.25),
                       labels = scales::percent_format(1)) +
    labs(x = "P(indicated)", y = "P(bolus)"))
ggsave(file.path(dirs$phase, "exposure_scatter.png"), p_sc,
       width = 5.6, height = 5.2, dpi = 200)


# ---- 10. The models ----------------------------------------------------------
# Proportional-odds ordinal logistic. clm() models cumulative odds of being at or
# below a level and nothing is negated: positive means higher odds of a BETTER band.

# The exposure terms are the SCALED ones, so coefficient names carry their unit.
FORM_PRIMARY <- stats::as.formula(paste(
  "vfd_ord ~", paste(c(MODEL_EXPOSURES, COVARIATES), collapse = " + ")))
FORM_COVERAGE <- stats::update(FORM_PRIMARY, . ~ . + sofa_baseline_n_components)

# Rare units collapse to `other`; if `other` is ITSELF below the threshold those
# admissions are DROPPED rather than relabelled. The count is reported.
collapse_units <- function(frame) {
  n_by <- table(as.character(frame$location_name))
  small <- names(n_by)[n_by < MIN_UNIT_N]
  lab <- as.character(frame$location_name)
  if (length(small)) lab[lab %in% small] <- "other"
  n_by2 <- table(lab)
  still <- names(n_by2)[n_by2 < MIN_UNIT_N]
  dropped <- lab %in% still
  frame$icu <- factor(lab)
  attr(frame, "collapsed") <- small
  attr(frame, "dropped_n") <- sum(dropped)
  if (any(dropped)) frame <- frame[!dropped, ]
  frame$icu <- droplevels(factor(as.character(frame$icu)))
  frame
}

fit_one <- function(frame, form, label, why, nominal = NULL) {
  # The nominal terms are NOT in `form`, so name them here too or the complete-case
  # filter misses them and clm drops rows itself, after the collapse was sized.
  vars <- setdiff(c(all.vars(form), all.vars(nominal)), c("vfd_ord", "icu"))
  frame <- frame[stats::complete.cases(frame[, vars, drop = FALSE]), ]
  frame <- collapse_units(frame)
  n_collapsed <- length(attr(frame, "collapsed"))
  n_dropped <- attr(frame, "dropped_n")
  if (n_collapsed || n_dropped) {
    cat(sprintf("  %s: collapsed %d unit(s) into 'other'%s\n", label, n_collapsed,
                if (n_dropped) sprintf(", dropped %d admission(s) still below %d",
                                       n_dropped, MIN_UNIT_N) else ""))
  }
  frame <- droplevels(frame)
  if (nrow(frame) < 50 || nlevels(frame$vfd_ord) < 3) {
    cat(sprintf("  SKIP %s: %d admissions, %d outcome levels -- too few to fit\n",
                label, nrow(frame), nlevels(frame$vfd_ord)))
    return(NULL)
  }
  # do.call, NOT clm(form, data = frame): nominal_test() re-evaluates clm's verbatim
  # call, where `frame` resolves to graphics::frame and every refit silently NAs.
  args <- list(formula = form, data = frame, link = "logit")
  if (!is.null(nominal)) args$nominal <- nominal
  fit <- tryCatch(do.call(ordinal::clm, args),
                  error = function(e) {cat(sprintf("  FAILED %s: %s\n", label,
                                                   conditionMessage(e))); NULL})
  if (is.null(fit)) return(NULL)
  list(label = label, why = why, fit = fit, frame = frame, n = nrow(frame),
       n_units = nlevels(frame$icu), form = paste(deparse(form), collapse = " "))
}

MODELS <- list(
  fit_one(d, FORM_PRIMARY, "primary",
          "all admissions with >= 1 uptitration"),
  fit_one(d[d$n_increase >= SENS_MIN, ], FORM_PRIMARY,
          sprintf("sensitivity_min_%d_events", SENS_MIN),
          sprintf(paste0("restricted to >= %d increase events -- SG's check on ",
                         "the small denominator"), SENS_MIN)),
  # ICU-only. The ED documents indication far more sparsely than it doses, so part
  # of the between-level contrast in Pindicated is CHARTING. A sensitivity only.
  fit_one(d[tolower(d$location_category) == ICU_CAT, ], FORM_PRIMARY,
          "sensitivity_icu_only",
          sprintf(paste0("restricted to admissions whose care location at hour ",
                         "0 has location_category == '%s' -- excludes the ED ",
                         "and unattributed, where documented indication is ",
                         "measurably sparser than dosing"), ICU_CAT)),
  # The previously shipped primary on the four-band outcome, so the resolution
  # change is visible: same frame and formula, only the response differs.
  fit_one(local({ dd <- d; dd$vfd_ord <- dd$vfd_day; dd }), FORM_PRIMARY,
          "sensitivity_day_scale",
          paste0("VFD-28 at whole-day resolution, the primary until 0.16.0 -- ",
                 "expect not_testable PO terms, which is a refit that will not ",
                 "fit 22 cut-points, not a violation")),
  fit_one(local({ dd <- d; dd$vfd_ord <- dd$vfd_band4; dd }), FORM_PRIMARY,
          "sensitivity_four_band",
          paste0("the primary refitted on the four-band cut declared in ",
                 "four_band_breaks -- the primary from 0.14.0 to 0.15.0 -- so ",
                 "the resolution change is visible, not silent")),
  fit_one(d, FORM_COVERAGE, "sensitivity_sofa_coverage",
          paste0("primary plus sofa_baseline_n_components; NOT a primary ",
                 "covariate, see protocol.json ",
                 "model.outcome_model.sensitivity_models.sofa_coverage")))
MODELS <- Filter(Negate(is.null), MODELS)
stopifnot("no model converged -- nothing to report" = length(MODELS) > 0)

tidy_one <- function(m) {
  s <- summary(m$fit)$coefficients
  # clm puts the thresholds first; they are nuisance and MUST NOT BE SHIPPED.
  # Matched by name via names(fit$alpha) -- there is no `fit$alpha.names`.
  stopifnot("clm did not expose its thresholds as names(fit$alpha)" =
              length(m$fit$alpha) > 0 && !is.null(names(m$fit$alpha)))
  thr <- rownames(s) %in% names(m$fit$alpha)
  stopifnot("the threshold rows were not found in the coefficient table" =
              sum(thr) == length(m$fit$alpha))
  co <- s[!thr, , drop = FALSE]
  ci <- suppressMessages(tryCatch(stats::confint(m$fit, type = "Wald"),
                                  error = function(e) NULL))
  out <- data.frame(
    model = m$label, term = rownames(co),
    estimate = round(co[, 1], 6), std_error = round(co[, 2], 6),
    z = round(co[, 3], 4), p_value = signif(co[, 4], 4),
    odds_ratio = round(exp(co[, 1]), 4),
    n_admissions = m$n, n_units = m$n_units,
    is_exposure = rownames(co) %in% MODEL_EXPOSURES,
    stringsAsFactors = FALSE)
  # TWO SCALES. estimate/std_error/ci_* are LOG ODDS, which is what pooling needs;
  # odds_ratio/or_ci_* are exp() of them.
  if (!is.null(ci)) {
    out$ci_lo <- round(ci[match(out$term, rownames(ci)), 1], 6)
    out$ci_hi <- round(ci[match(out$term, rownames(ci)), 2], 6)
  } else {
    out$ci_lo <- NA_real_; out$ci_hi <- NA_real_
  }
  out$or_ci_lo <- round(exp(out$ci_lo), 4)
  out$or_ci_hi <- round(exp(out$ci_hi), 4)
  out$scale_note <- sprintf(paste0("estimate/std_error/ci_* are log odds; ",
    "odds_ratio/or_ci_* are exp(). EXPOSURE TERMS ARE PER %d PERCENTAGE ",
    "POINTS of the share, not per unit share -- see protocol.json ",
    "model.outcome_model._exposure_scale_why"), EXP_SCALE_PP)
  rownames(out) <- NULL
  out
}
coefs <- do.call(rbind, lapply(MODELS, tidy_one))

# ICU fixed effects are site-local nuisance -- `location_name` is free text. BLANKED
# rather than dropped, so the pooled table can see they were fitted.
is_icu <- grepl("^icu", coefs$term)
coefs$poolable <- !is_icu
coefs$estimate[is_icu] <- NA; coefs$std_error[is_icu] <- NA
coefs$z[is_icu] <- NA; coefs$p_value[is_icu] <- NA
coefs$odds_ratio[is_icu] <- NA
coefs$ci_lo[is_icu] <- NA; coefs$ci_hi[is_icu] <- NA
coefs$or_ci_lo[is_icu] <- NA; coefs$or_ci_hi[is_icu] <- NA

cat("\nModel coefficients (ICU fixed effects fitted but not shipped)\n")
show <- coefs[coefs$poolable, c("model", "term", "estimate", "std_error",
                               "odds_ratio", "or_ci_lo", "or_ci_hi", "p_value",
                               "n_admissions")]
print(show, row.names = FALSE)

# The 2x2 between the exposures, so the centre can test them jointly.
vcov_rows <- do.call(rbind, lapply(MODELS, function(m) {
  V <- tryCatch(stats::vcov(m$fit), error = function(e) NULL)
  if (is.null(V)) return(NULL)
  present <- intersect(MODEL_EXPOSURES, rownames(V))
  if (length(present) < 2) return(NULL)
  g <- expand.grid(term_row = present, term_col = present,
                   stringsAsFactors = FALSE)
  data.frame(model = m$label, g,
             covariance = round(mapply(function(a, b) V[a, b],
                                       g$term_row, g$term_col), 10),
             n_admissions = m$n, stringsAsFactors = FALSE)
}))

# Proportional odds is tested, not assumed, ONE TERM AT A TIME via
# nominal_test(scope = ): the all-terms call dies under R >= 4.3, since clm can
# return a length-2 convergence code and ordinal guards its refit with `&&`.
po_rows <- do.call(rbind, lapply(MODELS, function(m) {
  terms_in <- attr(stats::terms(m$fit), "term.labels")
  base <- data.frame(
    model = m$label, term = "<none>", df = NA_integer_,
    logLik = round(as.numeric(stats::logLik(m$fit)), 3),
    AIC = round(stats::AIC(m$fit), 3), LRT = NA_real_, p_value = NA_real_,
    status = "reference", note = NA_character_, stringsAsFactors = FALSE)
  if (!length(terms_in)) {
    base$status <- "errored"
    base$note <- "the fit exposes no model terms to test"
    return(base)
  }
  rows <- lapply(terms_in, function(tm) {
    # nominal_test swallows the refit's warnings into an NA row; collected here.
    warns <- character(0)
    res <- withCallingHandlers(
      tryCatch({
        nt <- as.data.frame(ordinal::nominal_test(
          m$fit, scope = stats::as.formula(paste("~", tm))))
        r <- nt[rownames(nt) == tm, , drop = FALSE]
        if (!nrow(r)) stop("the scoped table carries no row for this term")
        list(df = r[["Df"]], logLik = r[["logLik"]], AIC = r[["AIC"]],
             LRT = r[["LRT"]], p = r[["Pr(>Chi)"]], err = NA_character_)
      }, error = function(e) list(df = NA, logLik = NA, AIC = NA, LRT = NA,
                                  p = NA, err = conditionMessage(e))),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
      })
    p <- suppressWarnings(as.numeric(res$p))
    note <- if (!is.na(res$err)) res$err
            else if (is.na(p) && length(warns))
              paste(unique(trimws(warns)), collapse = " | ")
            else NA_character_
    data.frame(
      model = m$label, term = tm,
      df = suppressWarnings(as.integer(res$df)),
      logLik = round(suppressWarnings(as.numeric(res$logLik)), 3),
      AIC = round(suppressWarnings(as.numeric(res$AIC)), 3),
      LRT = round(suppressWarnings(as.numeric(res$LRT)), 4),
      p_value = signif(p, 4),
      status = if (is.na(p)) "not_testable" else "tested",
      note = note, stringsAsFactors = FALSE)
  })
  do.call(rbind, c(list(base), rows))
}))

GATE_PASSED <- TRUE
gate_rows <- NULL
if (!is.null(po_rows)) {
  tested <- po_rows[po_rows$status == "tested", ]
  untested <- po_rows[po_rows$status == "not_testable", ]
  bad <- tested[tested$p_value < GATE_ALPHA, ]
  cat(sprintf(paste0("\nProportional-odds test: %d term(s) tested, %d not ",
                     "testable, %d with p < %.2f\n"),
              nrow(tested), nrow(untested), nrow(bad), GATE_ALPHA))
  if (nrow(untested)) {
    cat(sprintf("  NOT TESTABLE: %s\n",
                paste(sprintf("%s/%s", untested$model, untested$term),
                      collapse = ", ")),
        "  The per-term refit returned nothing, so proportional odds is\n",
        "  UNTESTED for those terms -- which is not the same as satisfied.\n",
        "  A zero violation count is not a pass.\n", sep = "")
  }
  if (nrow(bad)) {
    print(bad[, c("model", "term", "LRT", "p_value")], row.names = FALSE)
    cat("  A violated PO assumption means ONE odds ratio does not describe that\n",
        "  term across all cut points. Report it; do not quietly keep the OR.\n")
  }

  # ---- The proportional-odds gate --------------------------------------------
  # The partial-PO sensitivity was dropped conditional on both exposures passing,
  # so this is enforced at the VERY END, after every output is on disk. Scope is
  # the primary -- see protocol.json proportional_odds_gate._scope_why.
  GATE_SCOPE <- as.character(GATE$scope)
  stopifnot("only proportional_odds_gate.scope = primary is implemented" =
              identical(GATE_SCOPE, "primary"))
  gr <- po_rows[po_rows$model == GATE_SCOPE & po_rows$term %in% GATE_MODEL_TERMS, ]
  absent <- setdiff(GATE_MODEL_TERMS, gr$term)
  gate_rows <- rbind(
    if (nrow(gr)) data.frame(
      model = GATE_SCOPE, term = gr$term, status = gr$status,
      p_value = gr$p_value, alpha = GATE_ALPHA,
      passed = gr$status == "tested" & !is.na(gr$p_value) &
        gr$p_value >= GATE_ALPHA, stringsAsFactors = FALSE),
    if (length(absent)) data.frame(
      model = GATE_SCOPE, term = absent, status = "absent",
      p_value = NA_real_, alpha = GATE_ALPHA, passed = FALSE,
      stringsAsFactors = FALSE))
  GATE_PASSED <- !is.null(gate_rows) && nrow(gate_rows) > 0 &&
    all(gate_rows$passed)
  cat(sprintf("\nProportional-odds gate on the exposures, %s model (alpha %.2f): %s\n",
              GATE_SCOPE, GATE_ALPHA,
              if (GATE_PASSED) "PASSED" else "*** FAILED ***"))
  print(gate_rows[, c("term", "status", "p_value", "passed")], row.names = FALSE)

  other_bad <- bad[bad$model != GATE_SCOPE & bad$term %in% GATE_MODEL_TERMS, ]
  if (nrow(other_bad)) {
    cat("  WARNING, not a halt: an exposure violates PO in a SENSITIVITY model.\n")
    print(other_bad[, c("model", "term", "p_value")], row.names = FALSE)
  }
}

diagnostics <- do.call(rbind, lapply(MODELS, function(m) data.frame(
  model = m$label, why = m$why, formula = m$form,
  n_admissions = m$n, n_units = m$n_units,
  n_outcome_levels = nlevels(m$fit$model[[1]]),
  converged = isTRUE(m$fit$convergence$code == 0),
  convergence_code = m$fit$convergence$code,
  convergence_message = paste(m$fit$convergence$messages, collapse = "; "),
  cond_H = signif(m$fit$cond.H, 4),
  max_abs_gradient = signif(max(abs(m$fit$gradient)), 4),
  n_collapsed_units = length(attr(m$frame, "collapsed")),
  n_dropped_small_unit = attr(m$frame, "dropped_n"),
  logLik = round(as.numeric(stats::logLik(m$fit)), 3),
  AIC = round(stats::AIC(m$fit), 3),
  link = m$fit$link, stringsAsFactors = FALSE)))
cat("\nModel diagnostics\n")
print(diagnostics[, c("model", "n_admissions", "n_units", "converged",
                      "cond_H", "AIC")], row.names = FALSE)
if (any(!diagnostics$converged)) {
  cat("\n  A MODEL THAT DID NOT CONVERGE IS NOT REPORTABLE. cond_H above ~1e10\n")
  cat("  means the parameters are not uniquely determined -- usually a fixed\n")
  cat("  effect level with too few admissions, or one that is separated on the\n")
  cat("  outcome. Check n_dropped_small_unit and the convergence_message.\n")
}


# ---- 10b. The adjusted prediction, and the unadjusted comparison -------------
# Adjusted mean VFD-28 across plausible exposure levels, with the unadjusted bins
# beside it. The day-scale expectation works only because a level IS a day.

pred_model <- Filter(function(m) identical(m$label, "primary"), MODELS)
stopifnot("the primary model is required for the adjusted prediction" =
            length(pred_model) == 1)
pred_model <- pred_model[[1]]
pfit   <- pred_model$fit
pframe <- pred_model$frame
# A BAND CARRIES NO DAYS. The expectation needs a days value per level, and it
# is the observed within-band mean ON THE MODEL FRAME, so the mapping matches the
# admissions actually fitted. Shipped in adjusted_prediction.csv.
DVALS  <- as.numeric(tapply(pframe$vfd_28, pframe$vfd_ord, mean))
N_ALPHA <- length(pfit$alpha)
stopifnot(
  "every fitted band must have an observed mean to map to" = !any(is.na(DVALS)),
  "the band values must increase with the band" = !is.unsorted(DVALS, strictly = TRUE),
  "the lowest level must be day 0, or p_zero does not mean VFD = 0" =
    isTRUE(all.equal(DVALS[1], 0)),
  "a band value fell outside the attainable VFD range" =
    DVALS[length(DVALS)] <= VFD_CEILING_D,
  "there must be one threshold per gap between levels" =
    N_ALPHA == length(DVALS) - 1)

# Design matrix rebuilt from the fit's own terms and xlevels: the bootstrap
# evaluates at a DRAWN coefficient vector, which predict.clm cannot. Column names
# are ASSERTED, or a reordered contrast multiplies the wrong column by a beta.
design_of <- function(fit, nd) {
  tt <- stats::delete.response(stats::terms(fit))
  X <- stats::model.matrix(tt, data = nd, xlev = fit$xlevels)
  X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
  stopifnot("the rebuilt design matrix does not match the fitted coefficients" =
              identical(colnames(X), names(fit$beta)))
  X
}

# Full band probabilities for a proportional-odds clm: P(Y <= k) =
# plogis(alpha_k - X beta), differenced across the thresholds.
po_probs <- function(X, alpha, beta) {
  eta <- as.vector(X %*% beta)
  cum <- stats::plogis(outer(eta, alpha, function(e, a) a - e))
  cbind(cum, 1) - cbind(0, cum)
}

# ONE plogis call over a (thresholds x admissions) outer difference -- a per-threshold
# loop is the entire runtime. Averaging cumulatively then differencing equals
# averaging the band probabilities, both steps being linear.
std_summary <- function(X, par) {
  alpha <- par[seq_len(N_ALPHA)]
  beta  <- par[-seq_len(N_ALPHA)]
  eta   <- as.vector(X %*% beta)
  cm    <- rowMeans(stats::plogis(outer(alpha, eta, "-")))
  c(mean_vfd = unname(DVALS[length(DVALS)] - sum(diff(DVALS) * cm)),
    p_zero = unname(cm[1]))
}

# Three routes to the same quantity at the fitted coefficients -- predict.clm, the
# explicit band matrix, and the cumulative shortcut -- must agree before any draw.
X0 <- design_of(pfit, pframe)
P_manual <- po_probs(X0, pfit$alpha, pfit$beta)
nd0 <- pframe[, setdiff(names(pframe), "vfd_ord"), drop = FALSE]
# predict.clm returns P(Y = the OBSERVED y) when the response is present in newdata
# and the full band matrix only when it is absent. Dropping the column is required.
P_ref <- stats::predict(pfit, newdata = nd0, type = "prob")$fit
if (!is.matrix(P_ref)) P_ref <- matrix(P_ref, nrow = nrow(X0))
stopifnot(
  "the hand-built band probabilities disagree with predict.clm" =
    max(abs(P_manual - P_ref)) < 1e-9,
  "the band probabilities must sum to 1 per admission" =
    max(abs(rowSums(P_manual) - 1)) < 1e-9,
  "the cumulative shortcut disagrees with the explicit band matrix" =
    max(abs(std_summary(X0, c(pfit$alpha, pfit$beta)) -
            c(mean(as.vector(P_manual %*% DVALS)), mean(P_manual[, 1])))) < 1e-9)
cat(sprintf("\n  prediction machinery verified against predict.clm (max |diff| %.2e)\n",
            max(abs(P_manual - P_ref))))

# STANDARDISATION, not at-means: each admission is set to the grid value with all
# other covariates observed, then averaged. An average sex or ICU does not exist.
PRED_GRID_SHARE <- PRED_GRID_PCT / 100
X_grid <- lapply(PRED_GRID_SHARE, function(sh) {
  nd <- nd0
  nd[[PRED_MODEL_EXPOSURE]] <- sh / EXP_SCALE
  design_of(pfit, nd)
})

set.seed(SEED)
cf0 <- c(pfit$alpha, pfit$beta)
V0  <- stats::vcov(pfit)
stopifnot(
  "vcov() is not ordered as (thresholds, coefficients)" =
    identical(rownames(V0), names(cf0)),
  "vcov() is not square over the whole parameter vector" =
    all(dim(V0) == length(cf0)))
Rchol <- tryCatch(chol(V0), error = function(e) NULL)
stopifnot("the covariance matrix is not positive definite -- cannot draw" =
            !is.null(Rchol))
DRAWS <- sweep(matrix(stats::rnorm(PRED_DRAWS * length(cf0)), nrow = PRED_DRAWS)
               %*% Rchol, 2, cf0, "+")
cat(sprintf("  %s bootstrap draws at seed %s\n",
            format(PRED_DRAWS, big.mark = ","), format(SEED)))

t_pred <- proc.time()[3]
adj_rows <- do.call(rbind, lapply(seq_along(PRED_GRID_SHARE), function(i) {
  X <- X_grid[[i]]
  point <- std_summary(X, cf0)
  reps <- vapply(seq_len(PRED_DRAWS), function(b) std_summary(X, DRAWS[b, ]),
                 numeric(2))
  q <- function(v) stats::quantile(v, c(0.025, 0.975), names = FALSE)
  qm <- q(reps["mean_vfd", ]); qp <- q(reps["p_zero", ])
  data.frame(
    model = pred_model$label,
    exposure = PRED_EXPOSURE,
    exposure_pct = PRED_GRID_PCT[i],
    exposure_share = PRED_GRID_SHARE[i],
    mean_vfd_28 = round(point[["mean_vfd"]], 4),
    mean_vfd_28_lo = round(qm[1], 4), mean_vfd_28_hi = round(qm[2], 4),
    pct_vfd_zero = round(100 * point[["p_zero"]], 2),
    pct_vfd_zero_lo = round(100 * qp[1], 2),
    pct_vfd_zero_hi = round(100 * qp[2], 2),
    band_values = paste(sprintf("%.4f", DVALS), collapse = "|"),
    n_admissions = nrow(X), n_draws = PRED_DRAWS, seed = SEED,
    stringsAsFactors = FALSE)
}))
cat(sprintf("  standardised over %s admissions in %.1fs\n",
            format(nrow(X0), big.mark = ","), proc.time()[3] - t_pred))
stopifnot(
  "a standardised interval does not bracket its own point estimate" =
    all(adj_rows$mean_vfd_28_lo <= adj_rows$mean_vfd_28 &
        adj_rows$mean_vfd_28 <= adj_rows$mean_vfd_28_hi),
  "a standardised mean fell outside the attainable VFD range" =
    all(adj_rows$mean_vfd_28 >= 0 & adj_rows$mean_vfd_28 <= VFD_CEILING_D))
cat("\nAdjusted predicted VFD-28, standardised over the covariate distribution\n")
print(adj_rows[, c("exposure_pct", "mean_vfd_28", "mean_vfd_28_lo",
                   "mean_vfd_28_hi", "pct_vfd_zero", "pct_vfd_zero_lo",
                   "pct_vfd_zero_hi")], row.names = FALSE)

# Exactly zero is its own bin -- a 0-20 bin would mix the point mass with a level of
# exposure. Each bin sits at its OBSERVED MEAN (_bins_plotted_at_observed_mean).
pframe$.bin <- cut(pframe[[PRED_EXPOSURE]], breaks = BIN_BREAKS,
                   labels = BIN_LABELS, ordered_result = TRUE)
stopifnot("an exposure value fell outside the declared bin breaks" =
            !any(is.na(pframe$.bin)))
# Wilson, not Wald: the zero share runs from near 0 to near 1 across these bins.
wilson <- function(k, n, z = 1.959964) {
  if (n == 0) return(c(NA_real_, NA_real_))
  p <- k / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(max(0, c - h), min(1, c + h))
}
bin_rows <- do.call(rbind, lapply(BIN_LABELS, function(lb) {
  s <- pframe[pframe$.bin == lb, ]
  n <- nrow(s); k <- sum(s$vfd_28 == 0)
  w <- wilson(k, n)
  se <- if (n > 1) stats::sd(s$vfd_28) / sqrt(n) else NA_real_
  data.frame(
    exposure = PRED_EXPOSURE, bin = lb, n_admissions = n,
    exposure_mean = if (n) round(mean(s[[PRED_EXPOSURE]]), 4) else NA_real_,
    exposure_pct_mean = if (n) round(100 * mean(s[[PRED_EXPOSURE]]), 2) else NA_real_,
    mean_vfd_28 = if (n) round(mean(s$vfd_28), 4) else NA_real_,
    mean_vfd_28_lo = if (!is.na(se)) round(mean(s$vfd_28) - 1.959964 * se, 4) else NA_real_,
    mean_vfd_28_hi = if (!is.na(se)) round(mean(s$vfd_28) + 1.959964 * se, 4) else NA_real_,
    median_vfd_28 = if (n) round(stats::median(s$vfd_28), 4) else NA_real_,
    n_vfd_zero = k,
    pct_vfd_zero = if (n) round(100 * k / n, 2) else NA_real_,
    pct_vfd_zero_lo = round(100 * w[1], 2), pct_vfd_zero_hi = round(100 * w[2], 2),
    stringsAsFactors = FALSE)
}))
# Small-cell suppression, this script's own rule: counts stay, statistics go.
thin <- bin_rows$n_admissions < MIN_CELL
if (any(thin)) {
  cat(sprintf("  %d bin(s) below %d admissions -- statistics suppressed, counts kept\n",
              sum(thin), MIN_CELL))
  for (cl in c("exposure_mean", "exposure_pct_mean", "mean_vfd_28",
               "mean_vfd_28_lo", "mean_vfd_28_hi", "median_vfd_28",
               "pct_vfd_zero", "pct_vfd_zero_lo", "pct_vfd_zero_hi")) {
    bin_rows[[cl]][thin] <- NA_real_
  }
}
bin_rows$suppressed <- thin
stopifnot("the bins must partition the model frame" =
            sum(bin_rows$n_admissions) == nrow(pframe))
cat("\nUnadjusted, by exposure bin\n")
print(bin_rows[, c("bin", "n_admissions", "exposure_pct_mean", "mean_vfd_28",
                   "n_vfd_zero", "pct_vfd_zero")], row.names = FALSE)


# ---- 10c. The prediction figure ----------------------------------------------
# Two panels on a shared exposure axis; only the top needs the day scale.

# ADJUSTED CURVE ONLY: a bin is an interval at its own observed mean and a grid
# value is a point, so drawing both reads as one series measured two ways; the bins
# ship as unadjusted_bins.csv. No legend and no colour AESTHETIC -- one series.
PRED_COL <- PREDICTION_COLOURS[["adjusted"]]
x_lab <- "% Fentanyl Increases with a Documented Indication"

panel <- function(y, ylo, yhi, y_title, show_x) {
  adj <- data.frame(x = adj_rows$exposure_pct, y = adj_rows[[y]],
                    lo = adj_rows[[ylo]], hi = adj_rows[[yhi]])
  p <- ggplot(adj, aes(x = x)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = PRED_COL, alpha = 0.16) +
    geom_line(aes(y = y), colour = PRED_COL, linewidth = 0.8) +
    geom_point(aes(y = y), colour = PRED_COL, size = 2.1) +
    # The axis ends where the grid does, now that no bin series needs the room.
    scale_x_continuous(limits = c(min(PRED_GRID_PCT) - 3, max(PRED_GRID_PCT) + 3),
                       breaks = PRED_GRID_PCT, expand = c(0, 0)) +
    labs(x = if (show_x) x_lab else NULL, y = y_title)
  p <- house(p) + theme(legend.position = "none")
  if (!show_x) p <- p + theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank())
  p
}

# KEEP THESE SHORT. A y axis title is set in the panel's own height, so a longer one
# overflows into the neighbouring panel's title. The denominator is in the caption.
p_top <- panel("mean_vfd_28", "mean_vfd_28_lo", "mean_vfd_28_hi",
               "Adjusted Mean VFD-28, Days", show_x = FALSE)
p_bot <- panel("pct_vfd_zero", "pct_vfd_zero_lo", "pct_vfd_zero_hi",
               "Adjusted % Admissions with VFD-28 = 0", show_x = TRUE)

# EQUAL PANEL HEIGHTS, SET EXPLICITLY: with heights = NULL the bottom grob carries
# the x axis and its title, so the two plotting regions differ in height.
g_pred <- stack_panels(list(p_top, p_bot), c(2.8, 2.8))
ggsave(file.path(dirs$phase, "vfd_vs_pindicated.png"), g_pred,
       width = 7.2, height = 6.8, dpi = 300, bg = PAPER)
cat("written: vfd_vs_pindicated.png\n")

register_caption(
  "vfd_vs_pindicated.png",
  "Does documenting an indication for a fentanyl increase track ventilator-free days?",
  sprintf(paste0(
    "Top: adjusted mean VFD-28 in days. Bottom: the adjusted share of ",
    "admissions with VFD-28 exactly zero. Both are the primary ",
    "proportional-odds model (n = %s admissions, %d ICU fixed effects) ",
    "STANDARDISED over the observed covariate distribution: every admission is ",
    "set to the exposure value on the x axis with age, sex, CCI, ",
    "pre-intubation SOFA, unit AND the other exposure (P(bolus)) left at their ",
    "observed values, and the predictions are averaged over admissions. NOT an ",
    "average patient -- sex and unit are factors and an average one does not ",
    "exist. NO LEGEND: there is one series, and \'adjusted\' is in the axis ",
    "titles because it describes the quantity rather than distinguishing a ",
    "mark; the bottom panel's denominator is the %s admissions in the model ",
    "frame, stated here because the axis title is set in the panel height and a ",
    "longer one overflows into the panel above. Points are the %d declared grid values, the line joins them, and the ",
    "band is a %s-draw percentile interval from MVN(coef, vcov) at seed %s. THE GRID STOPS AT ",
    "%d%% because the top of the exposure range is thinly populated and the ",
    "fixed effects would be extrapolated beyond it. THE UNADJUSTED ",
    "COMPARISON IS A TABLE, NOT A SERIES HERE -- see unadjusted_bins.csv, ",
    "which gives mean VFD-28, the share at zero and the admission count for ",
    "each of %d exposure bins (%s). It is deliberately not plotted on these ",
    "axes: a bin is an interval summarised at its own observed mean and a grid ",
    "value is a point, so drawing them together reads as one series measured ",
    "two ways. THE TOP PANEL RESTS ON A BAND-TO-DAYS MAPPING, because VFD-28 is ",
    "modelled as %d ordered bands and a band carries no days of its own. Each ",
    "band is assigned the OBSERVED MEAN VFD-28 of its own members on the model ",
    "frame (%s days), shipped in adjusted_prediction.csv so the mapping is ",
    "auditable. It moves the LEVEL of the curve by a fraction of a day, not its ",
    "shape. THE BOTTOM PANEL NEEDS NO MAPPING: it is the model's own first ",
    "level, and the lowest band holds exactly-zero and nothing else. AN ",
    "ASSOCIATION, NOT AN EFFECT -- confounding by indication is unaddressed, ",
    "and a frequently titrated patient is a patient someone was worried about."),
    format(nrow(X0), big.mark = ",", trim = TRUE), pred_model$n_units,
    format(nrow(X0), big.mark = ",", trim = TRUE),
    length(PRED_GRID_PCT), format(PRED_DRAWS, big.mark = ",", trim = TRUE),
    format(SEED), as.integer(max(PRED_GRID_PCT)),
    nrow(bin_rows), paste(bin_rows$bin, collapse = ", "),
    length(DVALS), paste(sprintf("%.2f", DVALS), collapse = " / ")))


# ---- 10d. Table 1 -------------------------------------------------------------
# 07 is the SOLE GENERATOR (protocol.json table1._role). It reads trajectory_long
# window 0 unrestricted, so it can describe the whole cohort, and it is the only
# script holding the exposures and the fitted model frame -- so the table and the
# model cannot end up describing different populations.
#
# THREE NESTED SCOPES ARE EXPORTED, two are displayed. Exporting all three makes
# the displayed pair a config edit rather than another re-run at every site.
T1 <- COV$table1
stopifnot(
  "protocol.json must declare table1" = !is.null(T1),
  "Table 1 must be stratified on cohort_scope" =
    identical(T1$stratify_on, "cohort_scope"),
  "table1 must declare no test -- every column is a subset of Overall" =
    is.null(T1$test_continuous) && is.null(T1$test_categorical))
T1_SCOPES <- names(T1$scopes)
T1_SHOW   <- unlist(T1$display_columns)
T1_MODEL_ONLY <- unlist(T1$model_scope_only)
T1_RACE   <- unlist(T1$race_levels)
RACE_COLLAPSE <- COV$time_invariant$race$reporting_collapse
stopifnot(
  "display_columns must name declared scopes" = all(T1_SHOW %in% T1_SCOPES),
  "protocol.json must declare race.reporting_collapse" = !is.null(RACE_COLLAPSE))

# The scope frames. `base` is window 0 for EVERY block, so it is the descriptive
# denominator including admissions that never received fentanyl; `d_landmark` is
# the cohort captured before the up-titration restriction; the model frame is the
# one the primary model actually fitted, after complete-case and the unit
# collapse. Their sizes differ by about a quarter and that is the point.
t1_frames <- list(
  overall = base,
  landmark_increase = d_landmark,
  model = pred_model$frame)
t1_n <- lapply(t1_frames, nrow)
cat("\nTable 1 scopes\n")
for (sc in T1_SCOPES)
  cat(sprintf("  %-20s %-46s %7s\n", sc, T1$scopes[[sc]]$label,
              format(t1_n[[sc]], big.mark = ",")))

# A column of values per scope. A variable absent from a frame yields NULL, which
# the builders render as an empty cell rather than a number over the wrong
# denominator -- that is how the model-only rows stay out of Overall.
t1_vals <- function(col, scopes = T1_SCOPES) {
  out <- lapply(scopes, function(sc) {
    fr <- t1_frames[[sc]]
    if (is.null(fr) || !(col %in% names(fr))) NULL else fr[[col]]
  })
  names(out) <- scopes
  Filter(Negate(is.null), out)
}
t1_scopes_for <- function(label)
  if (label %in% T1_MODEL_ONLY) "model" else T1_SCOPES

t1_cont <- function(label, col, unit, digits = 1) {
  t1_continuous(label, t1_vals(col, t1_scopes_for(label)), unit, digits,
                MIN_CELL, T1_SHOW)
}

T1_PARTS <- list(
  t1_cont("Age, years", "age", "years"),
  t1_categorical("Sex", t1_vals("sex"), min_cell = MIN_CELL,
                 scope_order = T1_SHOW, levels_shown = c("Female", "Male")),
  # Race reports the DECLARED levels; a level present in the data but not listed
  # still ships in the tidy export, so no pooled count is lost to a display choice.
  t1_categorical("Race", t1_vals("race"), collapse = RACE_COLLAPSE,
                 min_cell = MIN_CELL, scope_order = T1_SHOW,
                 levels_shown = T1_RACE),
  t1_cont("Charlson Comorbidity Index", "cci", "index", 0),
  t1_cont("SOFA at baseline (24h pre-intubation)", "sofa_baseline", "points", 0),
  t1_cont("First IMV episode, hours", "first_imv_episode_hours", "hours", 1),
  t1_cont("Ventilator-free days to day 28", "vfd_28", "days", 1),
  t1_cont("P(bolus)", "Pbolus", "share", 3),
  t1_cont("P(indicated)", "Pindicated", "share", 3))

t1_out <- t1_assemble(T1_PARTS, T1_SHOW,
                      setNames(lapply(T1_SHOW, function(sc) T1$scopes[[sc]]$label), T1_SHOW),
                      t1_n[T1_SHOW])
table1 <- t1_out$display
table1_continuous <- t1_out$continuous
table1_categorical <- t1_out$categorical
stopifnot(
  "Table 1 produced no rows" = nrow(table1) > 0,
  "the tidy continuous export must cover every displayed scope" =
    all(T1_SHOW %in% table1_continuous$stratum),
  "a model-only row must not carry an overall value" =
    !any(table1_continuous$stratum == "overall" &
         table1_continuous$variable %in% T1_MODEL_ONLY))
cat(sprintf("Table 1: %d display row(s), %d tidy continuous, %d tidy categorical\n",
            nrow(table1), nrow(table1_continuous), nrow(table1_categorical)))


# THE STROBE. 07 ships it because it is the only script that can: the funnel ends
# at the model frame, and the restrictions after the analytic cohort are 05's
# titration events and 07's own complete-case and unit collapse, none of which
# exist when 01 runs. 01's own funnel is cohort_construction.csv -- extract to
# analytic cohort -- and is read here rather than re-derived.
#
# TWO BRANCHES, NOT ONE CHAIN, because there are two denominators. The practice
# pattern results are per-escalation over every admission that escalated at all;
# the model is over the complete-case frame. Collapsing them into one column
# would imply both describe the same population.
cc_f <- file.path(dirs$out_final, "01_cohort", "cohort_construction.csv")
stopifnot("01_cohort/cohort_construction.csv is missing -- re-run code/01_build_cohort.py" =
            file.exists(cc_f))
cc <- utils::read.csv(cc_f, stringsAsFactors = FALSE)
stopifnot("cohort_construction.csv must carry the exclusion-style columns" =
            all(c("step", "n_before", "n_excluded", "reason", "n_after") %in% names(cc)))

# The trunk's last row is the analytic cohort, and it must be the denominator
# everything below descends from. A mismatch means 07 is reading a funnel from a
# different cohort-build run than the tables it analysed.
analytic <- suppressWarnings(as.numeric(
  cc$n_after[grepl("ANALYTIC COHORT", cc$step, fixed = TRUE)]))
stopifnot("cohort_construction.csv names no ANALYTIC COHORT row" =
            length(analytic) == 1 && is.finite(analytic))

trunk <- data.frame(branch = "trunk", cc[, c("step", "n_before", "n_excluded",
                                             "reason", "n_after")],
                    stringsAsFactors = FALSE)
# Drop 01's own landmark row from the trunk: it applies the landmark to the WHOLE
# cohort, whereas the model branch below applies it after the escalation
# restriction. Keeping both would show two different numbers for one step.
trunk <- trunk[!grepl("LANDMARK COHORT", trunk$step, fixed = TRUE), , drop = FALSE]

# FLOW holds 07's own steps as (step, n, dropped_because). Convert to the
# exclusion shape so one schema describes the whole funnel: consecutive retained
# counts give the exclusions.
f07 <- do.call(rbind, FLOW)
to_excl <- function(df, first_before, branch) {
  before <- c(first_before, df$n[-nrow(df)])
  data.frame(branch = branch, step = df$step, n_before = before,
             n_excluded = before - df$n, reason = df$dropped_because,
             n_after = df$n, stringsAsFactors = FALSE)
}
# The practice-pattern denominator is 07's first step: every admission that
# escalated at all. It is a BRANCH off the analytic cohort, not a step toward
# the model, which is why it is labelled and not chained.
descr <- to_excl(f07[1, , drop = FALSE], analytic, "practice_patterns")
descr$reason[1] <- "no fentanyl initiation or up-titration in the exposure window"
model <- to_excl(f07[-1, , drop = FALSE], f07$n[1], "model")
# The frame the model actually FITTED, after the unit collapse -- one more
# admission than the complete-case count at UCMC, and it is the number every
# coefficient is computed on.
model <- rbind(model, data.frame(
  branch = "model", step = "ANALYSED IN THE MODEL",
  n_before = f07$n[nrow(f07)], n_excluded = f07$n[nrow(f07)] - pred_model$n,
  reason = sprintf("in a unit below min_admissions_per_unit = %d", MIN_UNIT_N),
  n_after = pred_model$n, stringsAsFactors = FALSE))

strobe <- rbind(trunk, descr, model)
stopifnot(
  "the practice-pattern branch must descend from the analytic cohort" =
    descr$n_before[1] == analytic,
  "the model branch must descend from the escalation cohort" =
    model$n_before[1] == f07$n[1],
  "every step must retain no more than it started with" =
    all(strobe$n_after <= strobe$n_before | is.na(strobe$n_before)),
  "the fitted frame must be the last row" =
    strobe$n_after[nrow(strobe)] == pred_model$n)
cat(sprintf("\nSTROBE: %d trunk + %d practice-pattern + %d model row(s)\n",
            nrow(trunk), nrow(descr), nrow(model)))
print(strobe[, c("branch", "step", "n_after")], row.names = FALSE)

# ---- 11. Write ---------------------------------------------------------------

write_out <- function(x, name) {
  if (is.null(x) || !nrow(x)) {
    cat(sprintf("  SKIP %s: nothing to write\n", name)); return(invisible(NULL))
  }
  key <- if (name %in% OWNED$diagnostics) "diagnostics" else "phase"
  stopifnot("write_out called for a file OWNED does not declare" =
              name %in% c(OWNED$phase, OWNED$diagnostics))
  write.csv(x, file.path(dirs[[key]], name), row.names = FALSE)
  cat(sprintf("written: %s%s\n", if (key == "diagnostics") "diagnostics/" else "",
              name))
}

cat("\n")
write_out(dist_rows, "exposure_distribution.csv")
write_out(by_unit, "exposure_by_unit.csv")
write_out(corr, "exposure_correlation.csv")
write_out(vfd_tab, "vfd_ordinal.csv")
write_out(table1, "table1.csv")
write_out(table1_continuous, "table1_continuous.csv")
write_out(table1_categorical, "table1_categorical.csv")
write_out(day_tab, "vfd_day_distribution.csv")
write_out(adj_rows, "adjusted_prediction.csv")
write_out(bin_rows, "unadjusted_bins.csv")
write_out(coefs, "model_coefficients.csv")
write_out(vcov_rows, "model_covariance.csv")
write_out(diagnostics, "model_diagnostics.csv")
write_out(po_rows, "proportional_odds_test.csv")
write_out(gate_rows, "proportional_odds_gate.csv")

# PHI handoff, one row per admission, so the model can be refitted. Stays PHI-side.
keep <- c("encounter_block", "n_increase", "n_paired", "n_indicated",
          "n_uptitration", "Pbolus", "Pindicated", "vfd_28", "vfd_ord",
          "vfd_band4", "vfd_day",
          "age", "sex", "cci", "sofa_baseline", "sofa_baseline_n_components",
          "location_name", "extubation_event", "died", ".complete")
write_parquet(d[, intersect(keep, names(d))],
              file.path(phi_dir(dirs, "07_outcome_model"), "outcome_model_frame.parquet"))
cat(sprintf("written: outcome_model_frame.parquet (PHI -- %s admissions)\n",
            format(nrow(d), big.mark = ",")))

write.csv(strobe, file.path(dirs$phase, "strobe.csv"), row.names = FALSE)
cat("written: strobe.csv\n")

write_json(prov, file.path(dirs$phase, "provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)
cat("written: provenance.json\n")
write_captions(file.path(dirs$phase, "captions.md"), "07_outcome_model.R",
               grep("\\.png$", c(OWNED$phase, OWNED$diagnostics), value = TRUE),
               prov)
cat("written: captions.md\n")


# ---- 12. Provenance ----------------------------------------------------------

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/07_outcome_model.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "07_outcome_model_sessioninfo.txt")
)

cat(sprintf("\nDone. %s admissions analysed, %s complete on every covariate.\n",
            format(nrow(d), big.mark = ","), format(sum(d$.complete), big.mark = ",")))
cat("THIS IS AN ASSOCIATION. No randomisation, no instrument, and confounding by\n")
cat("indication is unaddressed -- a frequently titrated patient is a patient\n")
cat("someone was worried about. Say 'associated with', never 'improves'.\n")


# ---- 13. The gate, enforced --------------------------------------------------
# Last statement, deliberately: every output is already on disk when this fires.
if (!GATE_PASSED) {
  failed <- if (is.null(gate_rows)) "the exposures" else
    paste(gate_rows$term[!gate_rows$passed], collapse = " and ")
  stop(sprintf(paste0(
    "PROPORTIONAL-ODDS GATE FAILED for %s.\n",
    "  The partial-proportional-odds sensitivity was dropped on 2026-10-28 ",
    "CONDITIONAL on the exposures passing this test -- see protocol.json ",
    "model.outcome_model.proportional_odds_gate. That condition no longer ",
    "holds here, so a single odds ratio is not a valid summary of the ",
    "exposure and model_coefficients.csv does not mean what its column names ",
    "say.\n",
    "  Every output HAS been written. Read diagnostics/",
    "proportional_odds_test.csv and diagnostics/proportional_odds_gate.csv, ",
    "and DO NOT POOL this bundle. Reinstate a partial-PO model -- fit_one() ",
    "still takes `nominal` -- or raise it with the coordinating centre."),
    failed))
}
