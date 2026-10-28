# ==============================================================================
# 07_outcome_model.R  --  from events to admissions: does titration practice
#                         track ventilator-free days?
#
# Purpose : The first analysis here that leaves the EVENT and asks about the
#           ADMISSION. Two per-admission exposures measured over the first 72h --
#           the share of fentanyl increases paired with a bolus (Pbolus) and the
#           share with a documented pain or agitation score (Pindicated) -- against
#           VFD-28, by proportional-odds ordinal regression with ICU fixed effects.
# Author  : Shan Guleria
# Created : 2026-10-27
# Inputs  : output/intermediate_phi/ : titration_adherence.parquet,
#           indication_events.parquet, time_to_event.parquet,
#           trajectory_long.parquet, care_intervals.parquet
# Outputs : output/final_no_phi/07_outcome_model/ : the files listed in OWNED
#           output/intermediate_phi/outcome_model_frame.parquet (PHI, per admission)
#
# THE EXPOSURES ARE NOT REDEFINED HERE. 05_titration.R owns what counts as an
# increase, a pairing and an indication; this script aggregates its two
# per-event handoffs to the admission. A second implementation of "what is an
# uptitration" would be a drift hazard with no upside -- the same rule
# 06_unit_variation.R follows.
#
# THIS IS AN ASSOCIATION, NOT AN EFFECT. No randomisation, no instrument, and
# confounding by indication is the obvious threat: a patient who is titrated
# often is a patient someone is worried about. The ICU fixed effects absorb
# between-unit practice, not within-unit case mix. Say "associated with" and
# never "improves" or "causes".
#
# THE DENOMINATOR IS SMALL AND THAT IS THE MAIN THREAT TO THIS ANALYSIS.
# MEASURED at UCMC: median 4 increase events per admission (IQR 2-6), 9.4% with
# exactly one, and only 55 distinct Pbolus values across 2,897 admissions. A
# proportion built on a denominator of one to three is mostly denominator noise,
# which pulls the coefficient toward the null. That is what the sensitivity
# analysis is for, and why its result deserves as much space as the primary.
# Adjusting for the event count is NOT the fix and is deliberately not done --
# see protocol.json model.outcome_model._n_events_is_NOT_a_covariate.
# ==============================================================================

# Run this in a FRESH R session (RStudio: Cmd+Shift+F10).

# ---- 1. Packages -------------------------------------------------------------

# `ordinal` is the only new pinned dependency this study has taken on since the
# stack was fixed; clm() has no substitute in base R. metafor is deliberately NOT
# here -- the pooling is coordinating-centre code and that decision is deferred,
# so a site never installs it. See protocol.json model.outcome_model._pooling_why.
pkgs <- c("here", "jsonlite", "ggplot2", "arrow", "ordinal")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

source(here("code", "utils", "paths.R"))
source(here("code", "utils", "figures.R"))


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
LANDMARK_H <- COV$cohort$landmark_hours
MIN_CELL   <- COV$reporting$small_cell_min_den

# The ordinal cut points. Declared, never chosen here: two sites cutting VFD-28
# at different days are not estimating the same parameter.
CUTS   <- as.numeric(unlist(SPEC$vfd_ordinal_breaks))
LABELS <- unlist(SPEC$vfd_ordinal_labels)
stopifnot(
  "protocol.json must declare vfd_ordinal_breaks and vfd_ordinal_labels" =
    length(CUTS) > 1 && length(LABELS) == length(CUTS) - 1,
  "the breaks must be increasing" = !is.unsorted(CUTS, strictly = TRUE))


# ---- 3. Paths and provenance -------------------------------------------------

dirs <- site_dirs()
dirs$phase <- phase_dir(dirs, "07_outcome_model")
dirs$diagnostics <- phase_dir(dirs, file.path("07_outcome_model", "diagnostics"))
prov <- provenance(config)


# ---- 4. Guards ---------------------------------------------------------------

manifest <- require_manifest(dirs, here())
message(sprintf("  reading Phase 0 outputs from code %s, generated %s",
                manifest$code_version, manifest$generated))

OWNED <- list(
  out_phi = c("outcome_model_frame.parquet"),
  # REPORTED DIRECTLY.
  phase = c("exposure_distribution.csv", "exposure_by_unit.csv",
            "exposure_correlation.csv", "vfd_ordinal.csv",
            "model_coefficients.csv", "model_covariance.csv", "strobe.csv",
            "exposure_histogram.png", "exposure_by_unit.png",
            "exposure_scatter.png", "provenance.json", "captions.md"),
  # METHOD CHECKS. The proportional-odds test belongs here: it is evidence about
  # whether the model is admissible, not a result.
  diagnostics = c("model_diagnostics.csv", "proportional_odds_test.csv"))
RETIRED <- character(0)          # new script on 2026-10-27
n_cleared <- clear_owned_outputs(dirs, OWNED, retired = RETIRED)
if (n_cleared) message(sprintf("  cleared %d output(s) from a previous run", n_cleared))

FLOW <- list()
flow <- function(label, n, dropped_because = NA_character_) {
  FLOW[[length(FLOW) + 1]] <<- data.frame(
    step = label, n = n, dropped_because = dropped_because,
    stringsAsFactors = FALSE)
  cat(sprintf("  %-58s %7s\n", label, format(n, big.mark = ",")))
}


# ---- 5. Read -----------------------------------------------------------------

# Pbolus comes straight from 05's per-encounter handoff. Its `overall` scope IS
# the first-72h exposure: titration events are bounded at the granular extent and
# 05 asserts it, so there is no separate 72h filter to apply and no second place
# the window could drift.
adh <- as.data.frame(read_parquet(file.path(dirs$out_phi, "titration_adherence.parquet")))
stopifnot("titration_adherence must carry the overall scope" =
            "overall" %in% adh$scope)
pb <- adh[adh$scope == "overall", c("encounter_block", "n_events", "n_paired")]
names(pb)[names(pb) == "n_events"] <- "n_increase"

ind <- as.data.frame(read_parquet(file.path(dirs$out_phi, "indication_events.parquet")))
stopifnot(
  "indication events must be bounded by the granular extent" =
    max(ind$t_hr) <= COV$cohort$granular_extent_hours + 1e-9,
  "indication events must start at or after the anchor" = min(ind$t_hr) >= 0)

tte <- as.data.frame(read_parquet(
  file.path(dirs$out_phi, "time_to_event.parquet"),
  col_select = c("encounter_block", "vfd_28", "landmark_eligible",
                 "extubation_event", "died")))

base <- as.data.frame(read_parquet(
  file.path(dirs$out_phi, "trajectory_long.parquet"),
  col_select = c("encounter_block", "window_idx", "age", "sex", "cci",
                 "sofa_baseline", "sofa_baseline_n_components")))
base <- base[base$window_idx == 0, setdiff(names(base), "window_idx")]
stopifnot("trajectory_long must carry sofa_baseline -- re-run Phase 0" =
            "sofa_baseline" %in% names(base),
          "window 0 must hold each block once" = !anyDuplicated(base$encounter_block))

icu <- as.data.frame(read_parquet(file.path(dirs$out_phi, "care_intervals.parquet")))


# ---- 6. The two exposures ----------------------------------------------------
# Both are a share of INCREASE events: initiation union uptitration, excluding
# downtitrations. They share the identical denominator, which is asserted rather
# than assumed -- if the two handoffs ever disagree about what an increase is,
# the ratios stop being comparable and the collinearity check becomes meaningless.

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
# Three restrictions, each with a reason, each a STROBE row. The landmark one is
# not optional: VFD-28 is only defined for admissions still ventilated at T, and
# it is also what guarantees the exposure window closed before the outcome
# window opened.

flow("Admissions with >= 1 fentanyl increase event", nrow(ex))

d <- merge(ex, tte[tte$landmark_eligible %in% TRUE, ], by = "encounter_block")
flow(sprintf("Landmark cohort (VFD-28 defined, on IMV at %dh)", LANDMARK_H),
     nrow(d), "not in the landmark cohort")

d <- d[d$n_uptitration >= 1, ]
flow("With >= 1 UPTITRATION (not initiation alone)", nrow(d),
     "initiation only -- no titration decision to describe")

# ICU at hour 0. Not the modal unit over the exposure window: that would be
# partly determined by the exposure and would make a baseline covariate
# post-exposure. See protocol.json model.outcome_model._icu_attribution_why.
icu <- icu[icu$t_start_hr <= 0 & icu$t_end_hr > 0, ]
icu <- icu[!duplicated(icu$encounter_block), c("encounter_block", "location_name")]
d <- merge(d, icu, by = "encounter_block", all.x = TRUE)
d$location_name[is.na(d$location_name)] <- "unattributed"

# `unattributed` is a real group -- no attributed care location covers hour 0 --
# and is kept as its own level rather than dropped or merged, because merging it
# into a ward it was not in would be a fabrication.
d$location_name[is.na(d$location_name)] <- "unattributed"
d$icu <- factor(d$location_name)
cat(sprintf("  ICU levels before collapsing: %d\n", nlevels(d$icu)))
# THE COLLAPSE HAPPENS PER MODEL FRAME, not here. Each model drops a different
# set of incomplete cases, so a level that clears the threshold on 2,897
# admissions can fall to one on the 2,472 that are complete -- which is exactly
# what happened on the first run: a single collapsed admission sat alone in the
# top VFD band, the Hessian went singular (condition number 1.1e16) and the
# primary model did not converge. Sizing the collapse on the wrong frame is the
# bug; see collapse_units() in fit_one().

# Baseline covariates, all measured at or before intubation. age / sex / cci are
# time-invariant; sofa_baseline is the worst value in the 24h BEFORE the anchor
# (protocol.json time_invariant.sofa_baseline), which is why it is a Phase 0
# product and not derived here.
d <- merge(d, base, by = "encounter_block", all.x = TRUE)
stopifnot("the covariate merge changed the admission count" =
            nrow(d) == nrow(unique(d["encounter_block"])),
          "every model covariate must be present after the merge" =
            all(c("age", "sex", "cci", "sofa_baseline") %in% names(d)))
d$sex <- factor(d$sex)

# The ordinal outcome. cut() is right-closed, so the first break sits BELOW zero
# to put exactly 0 in its own bottom level -- half the cohort lives there and it
# must not leak into the next band.
d$vfd_ord <- cut(d$vfd_28, breaks = CUTS, labels = LABELS, ordered_result = TRUE)
stopifnot("a VFD-28 value fell outside the declared breaks" = !any(is.na(d$vfd_ord)),
          "the ordinal outcome must be ordered" = is.ordered(d$vfd_ord))

vfd_tab <- as.data.frame(table(level = d$vfd_ord), stringsAsFactors = FALSE)
names(vfd_tab)[2] <- "n_admissions"
vfd_tab$denominator <- nrow(d)
vfd_tab$pct <- round(100 * vfd_tab$n_admissions / nrow(d), 1)
cat("\nVFD-28 as an ordinal outcome\n")
print(vfd_tab, row.names = FALSE)

# Complete case, declared rather than left to clm()'s default.
NEED <- c("Pbolus", "Pindicated", "age", "sex", "cci", "sofa_baseline", "icu")
d$.complete <- stats::complete.cases(d[, NEED])
flow("Complete on every model covariate", sum(d$.complete),
     "missing a covariate -- sofa_baseline is the binding one")
for (v in NEED) {
  miss <- sum(is.na(d[[v]]))
  if (miss) cat(sprintf("    %-28s missing %s (%.1f%%)\n", v,
                        format(miss, big.mark = ","), 100 * miss / nrow(d)))
}


# ---- 8. Exposure description -------------------------------------------------
# Shipped as BINNED COUNTS, not as a density: counts on fixed bins pool across
# sites by addition, a kernel density does not.
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

# Per unit, for the violin figure and for between-site heterogeneity. Quantiles
# do NOT pool; the counts beside them do, and the RANGE across units is the thing
# to read. Units below the reporting floor are suppressed, not published.
by_unit <- do.call(rbind, lapply(EXPOSURES, function(v) {
  do.call(rbind, lapply(levels(d$icu), function(u) {
    x <- d[[v]][d$icu == u]
    if (!length(x)) return(NULL)
    data.frame(exposure = v, location_name = u, n_admissions = length(x),
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
  by_unit[suppressed, c("median", "q1", "q3", "mean", "sd")] <- NA
}

# Collinearity: the question SG asked before the model was written. Pearson for
# the linear-model concern, Spearman because both exposures are proportions with
# a point mass at zero and Pearson alone would be read too confidently.
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

# Violin by unit: within-ICU spread AND between-ICU differences in one panel,
# which is what SG asked these figures to show. A boxplot inside each violin
# because a violin alone hides the median, and the n under each unit because a
# wide violin on twelve admissions is not the same claim as one on four hundred.
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
                 "by the exposure). Violin shows the full within-unit ",
                 "distribution, box the median and IQR; n under each label. ",
                 "UNITS ARE SHOWN AS CHARTED HERE -- the collapsing of units ",
                 "below %d admissions into 'other' happens inside each MODEL, on ",
                 "the frame actually fitted, because complete-case filtering ",
                 "changes which units are too small. READ THE WITHIN-UNIT SPREAD ",
                 "FIRST: if it dwarfs the between-unit differences, unit is not ",
                 "where the variation lives. Units carry their real names here; ",
                 "audience-facing versions anonymise them."),
          MIN_UNIT_N))

p_unit <- house(
  ggplot(long_ex, aes(icu_lab, value, fill = exposure)) +
    geom_violin(scale = "width", width = 0.85, colour = NA, alpha = 0.55) +
    geom_boxplot(width = 0.16, outlier.shape = NA, colour = INK,
                 fill = PAPER, linewidth = 0.3) +
    facet_wrap(~exposure, ncol = 1) +
    scale_fill_manual(values = unname(EXP_COL[EXPOSURES]), guide = "none") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = scales::percent_format(1)) +
    labs(x = NULL, y = "Share of increase events") +
    theme(axis.text.x = element_text(size = 8, colour = INK)))
ggsave(file.path(dirs$phase, "exposure_by_unit.png"), p_unit,
       width = 9.0, height = 6.0, dpi = 200)

# Scatter: the collinearity check, drawn. Jittered because both axes are
# proportions on a small denominator and hundreds of admissions sit on the exact
# same coordinate -- an unjittered plot would show a dozen points for 2,897.
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
# Proportional-odds ordinal logistic. A positive coefficient means HIGHER odds of
# landing in a BETTER (higher) VFD band -- clm() models cumulative odds of being
# at or below a level, and we negate nothing, so read the sign against the
# ordering of vfd_ord and not against intuition.

FORM_PRIMARY <- stats::as.formula(
  "vfd_ord ~ Pbolus + Pindicated + age + sex + cci + sofa_baseline + icu")
FORM_COVERAGE <- stats::update(FORM_PRIMARY, . ~ . + sofa_baseline_n_components)

# A fixed effect needs enough admissions to be identified. Rare units collapse to
# `other`; if `other` is ITSELF below the threshold the remaining admissions are
# DROPPED rather than relabelled, because assigning a patient to a unit they were
# never in to make a model fit is not a modelling choice. The count is reported.
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
  # The nominal terms are NOT in `form`, so they must be named here as well or
  # the complete-case filter would miss them and clm would drop rows itself,
  # silently, after the unit collapse had already been sized.
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
  # do.call, NOT ordinal::clm(form, data = frame, ...). The difference is not
  # style: clm() stores its call verbatim, and nominal_test() REFITS the model
  # once per term by re-evaluating that call. A plain call records
  # `data = frame`, `frame` is a local that no longer exists once this function
  # returns -- and `frame` is also a BASE R FUNCTION (graphics::frame), so the
  # refit resolves it to that function instead of erroring on an unknown name.
  # Every refit then fails, and nominal_test returns a well-formed table with
  # every term's LRT and p-value NA. do.call inlines the data frame into the
  # recorded call, so the refits have the data they need. Verified: the
  # coefficients are identical to the plain-call fit at tolerance = 0.
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
  # PARTIAL PROPORTIONAL ODDS. Added after the PO test was repaired and turned
  # out to FAIL for age (p = 1.3e-08), cci (p = 9.9e-05) and sex (p = 0.028)
  # while HOLDING for both exposures and for sofa_baseline. A single odds ratio
  # is not a valid summary for those three, so here they get cut-point-specific
  # effects and the exposures keep one OR each.
  # A SENSITIVITY, not the primary: SG specified the model, the violations are in
  # adjustment terms rather than the estimands, and the exposure coefficient is
  # what travels to the pool. If the exposure ORs move materially between this
  # and the primary, the primary is the one to doubt.
  fit_one(d, stats::as.formula("vfd_ord ~ Pbolus + Pindicated + sofa_baseline + icu"),
          "sensitivity_partial_po",
          paste0("age, sex and cci given nominal (cut-point-specific) effects ",
                 "because they violate proportional odds; exposures keep one OR"),
          nominal = ~ age + sex + cci),
  fit_one(d, FORM_COVERAGE, "sensitivity_sofa_coverage",
          paste0("primary plus sofa_baseline_n_components; NOT a primary ",
                 "covariate, see protocol.json ",
                 "model.outcome_model.sensitivity_models.sofa_coverage")))
MODELS <- Filter(Negate(is.null), MODELS)
stopifnot("no model converged -- nothing to report" = length(MODELS) > 0)

tidy_one <- function(m) {
  s <- summary(m$fit)$coefficients
  # clm puts the thresholds first. They are NUISANCE AND MUST NOT BE SHIPPED: a
  # threshold depends on each site's covariate distribution, so pooling one
  # across sites would combine numbers that are not the same quantity. Matched by
  # NAME via names(fit$alpha) -- `fit$alpha.names` does not exist, which an
  # earlier version of this line assumed, and the thresholds leaked into the
  # coefficient table as a result.
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
    is_exposure = rownames(co) %in% EXPOSURES,
    stringsAsFactors = FALSE)
  # TWO SCALES, NAMED. estimate/std_error/ci_* are LOG ODDS -- that is the scale
  # inverse-variance pooling needs. odds_ratio/or_ci_* are the exponentiated
  # version for reading. An earlier draft shipped an OR beside a log-odds CI
  # under the bare names ci_lo/ci_hi, which invites pairing them.
  if (!is.null(ci)) {
    out$ci_lo <- round(ci[match(out$term, rownames(ci)), 1], 6)
    out$ci_hi <- round(ci[match(out$term, rownames(ci)), 2], 6)
  } else {
    out$ci_lo <- NA_real_; out$ci_hi <- NA_real_
  }
  out$or_ci_lo <- round(exp(out$ci_lo), 4)
  out$or_ci_hi <- round(exp(out$ci_hi), 4)
  out$scale_note <- "estimate/std_error/ci_* are log odds; odds_ratio/or_ci_* are exp()"
  rownames(out) <- NULL
  out
}
coefs <- do.call(rbind, lapply(MODELS, tidy_one))

# ICU fixed effects are site-local nuisance: `location_name` is free text and an
# `icuMICU` coefficient at one site names nothing at another. They are BLANKED
# rather than dropped so the pooled table can still see that they were fitted.
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

# The 2x2 between the two exposures, so the centre can test them jointly rather
# than only one at a time. Shipping SEs alone forecloses that.
vcov_rows <- do.call(rbind, lapply(MODELS, function(m) {
  V <- tryCatch(stats::vcov(m$fit), error = function(e) NULL)
  if (is.null(V)) return(NULL)
  present <- intersect(EXPOSURES, rownames(V))
  if (length(present) < 2) return(NULL)
  g <- expand.grid(term_row = present, term_col = present,
                   stringsAsFactors = FALSE)
  data.frame(model = m$label, g,
             covariance = round(mapply(function(a, b) V[a, b],
                                       g$term_row, g$term_col), 10),
             n_admissions = m$n, stringsAsFactors = FALSE)
}))

# PROPORTIONAL ODDS IS TESTED, NOT ASSUMED. With half the cohort at VFD = 0 the
# assumption is exactly where this model is most likely to fail, and a violated
# PO assumption makes a single odds ratio per exposure the wrong summary.
#
# The table CARRIES ITS OWN STATUS, because the way this diagnostic fails is by
# succeeding in appearance. nominal_test() returns a fully formed table -- the
# `<none>` row populated from the existing fit, one row per term -- and NA's the
# terms when the per-term refits fail. Counting "terms with p < 0.05" over that
# gives 0 and reads as "the assumption holds", which is the opposite of what
# happened. An untestable model must be distinguishable from a passing one in
# the shipped file, not just in the console.
po_rows <- do.call(rbind, lapply(MODELS, function(m) {
  nt <- tryCatch(ordinal::nominal_test(m$fit), error = function(e) NULL)
  if (is.null(nt)) {
    return(data.frame(model = m$label, term = NA_character_, df = NA_integer_,
                      logLik = NA_real_, AIC = NA_real_, LRT = NA_real_,
                      p_value = NA_real_, status = "errored",
                      stringsAsFactors = FALSE))
  }
  nt <- as.data.frame(nt)
  terms_only <- rownames(nt) != "<none>"
  n_tested <- sum(!is.na(nt[["Pr(>Chi)"]][terms_only]))
  data.frame(model = m$label, term = rownames(nt),
             df = nt[["Df"]], logLik = round(nt[["logLik"]], 3),
             AIC = round(nt[["AIC"]], 3),
             LRT = round(nt[["LRT"]], 4),
             p_value = signif(nt[["Pr(>Chi)"]], 4),
             status = if (n_tested == 0L) "not_testable" else "tested",
             stringsAsFactors = FALSE)
}))
if (!is.null(po_rows)) {
  untestable <- unique(po_rows$model[po_rows$status != "tested"])
  tested <- po_rows[po_rows$status == "tested" & po_rows$term != "<none>", ]
  bad <- tested[!is.na(tested$p_value) & tested$p_value < 0.05, ]
  cat(sprintf("\nProportional-odds test: %d of %d model(s) testable, %d term(s) with p < 0.05\n",
              length(MODELS) - length(untestable), length(MODELS), nrow(bad)))
  if (length(untestable)) {
    cat(sprintf("  NOT TESTABLE: %s\n", paste(untestable, collapse = ", ")),
        " The per-term refits returned nothing, so the proportional-odds\n",
        " assumption is UNTESTED for those models -- not confirmed. Do not read\n",
        " a zero count above as the assumption holding.\n", sep = "")
  }
  if (nrow(bad)) {
    print(bad[, c("model", "term", "LRT", "p_value")], row.names = FALSE)
    cat("  A violated PO assumption means ONE odds ratio does not describe that\n",
        "  term across all cut points. Report it; do not quietly keep the OR.\n")
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
write_out(coefs, "model_coefficients.csv")
write_out(vcov_rows, "model_covariance.csv")
write_out(diagnostics, "model_diagnostics.csv")
write_out(po_rows, "proportional_odds_test.csv")

# PHI handoff: one row per admission, so the model can be refitted or extended
# without re-deriving the exposures. Stays PHI-side -- it is per-patient.
keep <- c("encounter_block", "n_increase", "n_paired", "n_indicated",
          "n_uptitration", "Pbolus", "Pindicated", "vfd_28", "vfd_ord",
          "age", "sex", "cci", "sofa_baseline", "sofa_baseline_n_components",
          "location_name", "extubation_event", "died", ".complete")
write_parquet(d[, intersect(keep, names(d))],
              file.path(dirs$out_phi, "outcome_model_frame.parquet"))
cat(sprintf("written: outcome_model_frame.parquet (PHI -- %s admissions)\n",
            format(nrow(d), big.mark = ",")))

write.csv(do.call(rbind, FLOW), file.path(dirs$phase, "strobe.csv"),
          row.names = FALSE)
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
