# ==============================================================================
# 06_transition_model.R  --  Phase 5 -- discrete-time multinomial transition model
#
# Purpose : Model where a ventilation episode's fentanyl delivery state goes next, given where it is now, how long it has been there, and how sick the patient is. Whole analytic cohort, 4h windows.
# Author  : Shan Guleria
# Created : 2026-09-09
# Inputs  : output/intermediate_phi/trajectory_long.parquet
# Outputs : output/intermediate_phi/transition_model.rds; output/final_no_phi/06_transitions/
#
# WHY THIS AND NOT A TRAJECTORY CLASS. Phases 3 and 4 established that dose LEVEL
# carries 66% of the variance, is continuously distributed with no gaps, and that
# gbmt's four classes were an account of level (R^2 = 0.82 for predicting an
# episode's mean dose) rather than a discovery of subpopulations -- a random
# intercept absorbed them almost entirely. This script does not partition anyone.
# It models the moves between states, so there is no k to choose and no claim
# that groups exist.
#
# WHY NOT msm OR mstate. msm requires the Markov property: the next state depends
# only on the current one. That is not credible here -- an episode on hour 60 of
# an infusion is not exchangeable with one that started an hour ago. mstate
# relaxes it via a clock-reset (semi-Markov) time scale but needs exactly
# observed transition times, which a 4h window grid does not give. The states are
# DEFINED on the window, so the process is genuinely discrete-time: history goes
# in as covariates and no Markov assumption is made at all.
#
# ESTIMAND: none. This is descriptive. Coefficients are associations between the
# covariates and the next observed state, not causal effects of dosing.
# ==============================================================================

# Run this in a FRESH R session (RStudio: Cmd+Shift+F10).
# Do not use rm(list = ls()) -- it does not unload packages or reset options,
# so it only gives the appearance of a clean slate.


# ---- 1. Packages -------------------------------------------------------------
# nnet and MASS ship with R; nothing to install.

pkgs <- c("here", "jsonlite", "arrow", "ggplot2", "nnet")
library(splines)   # ships with R; ns() for the clock terms in section 6 below
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

source(here("code", "utils", "paths.R"))
source(here("code", "utils", "states.R"))


# ---- 2. Config (never setwd(); here() anchors to the .Rproj) -----------------

config <- fromJSON(here("config", "config.json"), simplifyVector = FALSE)
set.seed(config$model$seed)

WINDOW_H <- config$cohort$window_hours
EXTENT_H <- config$cohort$granular_extent_hours
MIN_CELL <- config$reporting$small_cell_min_den
COV      <- fromJSON(here("config", "covariates.json"), simplifyVector = FALSE)
FENT_U   <- COV$exposure$units
# Dose bands, read from the config and never defaulted here -- the cuts are
# federation-critical (covariates.json exposure.dose_states._cuts_MUST_BE_ABSOLUTE).
DS_SPEC   <- COV$exposure$dose_states
DOSE_CUTS <- as.numeric(unlist(DS_SPEC$cuts))
DOSE_LAB  <- unlist(DS_SPEC$labels)
stopifnot("covariates.json exposure.dose_states must declare cuts and labels" =
            length(DOSE_CUTS) > 0 && length(DOSE_LAB) == length(DOSE_CUTS) + 2,
          "dose states must be defined on window_mcg" =
            identical(DS_SPEC$variable, "window_mcg"))

# Cluster bootstrap replicates. 17 transitions from one episode are correlated,
# and some patients contribute several episodes, so naive multinom standard
# errors are too narrow. sandwich has no estfun method for multinom, so the
# clustering is handled by resampling PATIENTS with replacement.
# MEASURED 2026-09-09: one multinom fit over 178,126 transition rows takes ~29s,
# so 200 replicates is ~1.6h. Override for a quick pass:
#     TRANSITION_BOOT_REPS=25 Rscript code/06_transition_model.R
# Deliberately NOT in config.json: it does not change the estimand, only the
# precision of the interval, and putting it there would force every site through
# a Phase 0 re-run to alter a compute setting.
BOOT_REPS <- as.integer(Sys.getenv("TRANSITION_BOOT_REPS", "200"))
stopifnot("TRANSITION_BOOT_REPS must be a positive integer" =
            !is.na(BOOT_REPS) && BOOT_REPS >= 1)


# ---- 3. Paths and provenance -------------------------------------------------

dirs <- site_dirs()
dirs$phase <- phase_dir(dirs, "06_transitions")
prov <- provenance(config)

message(sprintf("[06_transition_model] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))


# ---- 4. Guards ---------------------------------------------------------------

manifest <- require_manifest(dirs, here())
message(sprintf("  reading Phase 0 outputs from code %s, generated %s",
                manifest$code_version, manifest$generated))

PHASE5_STEMS <- c("transition_matrix.csv", "transition_counts.csv",
                  "model_coefficients.csv", "predicted_transitions.csv",
                  "transition_hazards.csv", "model_fit.csv",
                  "measurement_by_origin.csv", "covariate_missingness.csv",
                  "predicted_transitions.png", "severity_gradient.png",
                  "transition_hazards.png")
# One list, both definitions. Deriving the owned set from the same stems the
# writer uses is the point: a stem added to one and not the other leaves a stale
# file that the next run silently keeps.
OWNED <- list(
  out_phi = c("transition_model.rds"),
  phase   = c(as.vector(outer(c("phase5_", "phase5_dose_"), PHASE5_STEMS, paste0)),
              "phase5_provenance.json"))
# A stem removed from the writer leaves a stale twin the owned list no longer
# names, so clear_owned_outputs walks past it forever. phase5_complete_case_loss
# was written while the missing-indicator decision was still open, and
# orphaned when it closed.
RETIRED <- file.path("output", "final_no_phi", "06_transitions",
                     "phase5_complete_case_loss.csv")
n_cleared <- clear_owned_outputs(dirs, OWNED, retired = RETIRED)
if (n_cleared) message(sprintf("  cleared %d output(s) from a previous run", n_cleared))


# ---- 5. Data -- the WHOLE analytic cohort ------------------------------------
# NOT the landmark set. The landmark conditions on being ventilated at T, so
# inside [0, T] nobody dies, is discharged or is permanently extubated: the three
# terminal states would be empty by construction and the liberation pathway --
# the thing worth modelling -- would be excluded. 8,169 of 14,897 episodes leave
# before 72h, median at 24h.

long <- as.data.frame(read_parquet(
  file.path(dirs$out_phi, "trajectory_long.parquet"),
  col_select = c("encounter_block", "patient_id", "window_idx", "window_start_hr",
                 "alive_admitted", "imv_status", "inf_dose", "bolus_dose",
                 "total_dose", "window_mcg", "died",
                 "age", "sex", "bmi_admission", "cci",
                 "nee", "oxygenation", "spo2_plateau", "crrt_status",
                 "bun", "bicarbonate", "lactate")))

stopifnot("window count on disk disagrees with the config grid" =
            length(unique(long$window_idx)) == EXTENT_H / WINDOW_H)

# States are derived INSIDE run_definition (section 11 below) -- each needs
# its own. add_history is definition-dependent too: hours_in_state counts a run
# of the CURRENT state, so it must be recomputed per definition.

n_ep <- length(unique(long$encounter_block))
cat(sprintf("\nAnalytic cohort: %s ventilation episodes x %d windows\n",
            format(n_ep, big.mark = ","), length(unique(long$window_idx))))


# ---- 6..10b as ONE function, run once per state definition -------------------
# Both definitions -- delivery ROUTE and intensity BAND -- go through identical
# machinery: per-origin multinomial, cluster bootstrap, predicted transitions,
# hazard curves. Writing it once is the only way the two stay comparable; a
# copied second script would drift on the first edit. Lyons et al. fit two
# multistate models on one cohort for the same reason.
#
# TAG prefixes the outputs ("" for route, "dose_" for intensity) and REF names
# the multinomial reference destination, which differs because the natural
# baseline differs: `no fentanyl` for routes, `zero` for bands.

run_definition <- function(long, TAG, LABEL, REF) {

LEVELS <- levels(long$state)
cat(sprintf("\n\n%s\n== %s states: %s\n%s\n",
            strrep("=", 78), LABEL, paste(LEVELS, collapse = " | "), strrep("=", 78)))

# ---- 6. Transition pairs -----------------------------------------------------

long <- add_history(long, window_h = WINDOW_H)
tp <- transition_pairs(long)
cat(sprintf("Transition pairs: %s rows (terminal origins dropped)\n",
            format(nrow(tp), big.mark = ",")))

tm <- transition_matrix(tp)
print(tm, row.names = FALSE)

counts <- as.data.frame(table(from = tp$state, to = tp$state_next))
names(counts)[3] <- "n"


# ---- 7. Model -- ONE FIT PER ORIGIN STATE ------------------------------------
# Five separate multinomial fits, one per origin, rather than a single pooled fit
# with `state` as a covariate. The two differ in whether covariate effects may
# vary by origin: pooling forces ONE coefficient for, say, SOFA across every
# origin, when clinically a rising SOFA should mean "keep sedating, do not
# extubate" from `continuous only` and "this patient is failing" from
# `extubated` -- plausibly opposite signs.
#
# NOT a matter of taste. MEASURED 2026-09-09 on 238,630 rows: the pooled fit has
# deviance 262,596.8 on 78 df, the fully interacted equivalent 258,200.3 on 270
# df -- a drop of 4,396.5 on 192 df, p < 1e-300, and AIC prefers per-origin by
# ~4,000 despite the extra parameters. The pooled specification was misspecified.
#
# Fitting them separately rather than as `state * (covariates)` is statistically
# identical but reads better: each model's coefficient table is directly "what
# drives moves out of THIS state".
#
# Reference destination is `no fentanyl` throughout -- observed from every origin,
# and the natural comparator ("relative to fentanyl being stopped").

# Covariate set fixed by SG 2026-09-09. Summary rules are Phase 0's and are
# declared in covariates.json -- NEE max-of-summed-step-function, oxygenation min
# (worst, with the Severinghaus S/F fallback), CRRT any-in-window, BUN max,
# bicarbonate min, lactate max.
#
# IMV STATUS IS DELIBERATELY ABSENT. It is not an oversight: ventilation status
# IS the state here, so imv_status is constant within every per-origin model
# (1 for all four fentanyl states, 0 for extubated) and carries no information a
# covariate could use. The information is retained by the state itself.
#
# SOFA is also absent. covariates.json records four verified defects in clifpy's
# SOFA, all silent and all biasing severity downward; the explicit markers below
# carry the severity signal instead.
TIME_INVARIANT <- c("age", "sex", "bmi_admission", "cci")
# spo2_plateau sits BESIDE oxygenation, not instead of it. Above the SpO2
# ceiling the dissociation curve is flat, so Severinghaus cannot invert a
# saturation into a PaO2 and Phase 0 leaves the P/F NA -- but NA is the wrong
# encoding, because the window says the patient is oxygenating WELL. Filling it
# with a median P/F imputes moderate hypoxaemia for the patients doing best.
# The flag carries that fact directly (SG, 2026-09-09).
TIME_VARYING   <- c("nee", "oxygenation", "spo2_plateau", "crrt_status",
                    "bun", "bicarbonate", "lactate")
HISTORY        <- c("hours_in_state", "cumulative_dose", "window_start_hr")
COVS <- c(HISTORY, TIME_INVARIANT, TIME_VARYING)
FORM <- as.formula(paste("state_next ~", paste(COVS, collapse = " + ")))
REF  <- "no fentanyl"

tp$state <- droplevels(tp$state)

# ---- Missing-indicator method (SG, 2026-09-09) --------------------------------
# Complete-case analysis would keep only 133,653 of 238,811 transitions (56.0%),
# and the loss is NOT uniform: labs are drawn far less often once a patient is
# extubated, so extubated-origin rows kept just 34.3% against 57-71% elsewhere.
# That preferentially discards reintubation, discharge and death -- the
# transitions the liberation story turns on. multinom drops them silently.
#
# So: every row is kept. For each covariate with any missingness, a binary
# `<var>_measured` flag is added and the value is filled with the cohort median
# (mode for a categorical). The flag is not a nuisance -- a drawn lactate means
# somebody was worried, so "was it measured" carries real information about the
# clinician's assessment, which is exactly what this model is describing.
#
# VALID HERE, NOT IN GENERAL. The missing-indicator method is known to be biased
# for causal estimation. This model is explicitly descriptive (no estimand), so
# the trade -- keeping every extubated-origin row against a bias that would
# matter only for a causal claim we are not making -- is the right way round.
# Fill values are computed ONCE from the full data and held fixed across
# bootstrap replicates; they are nuisance constants, not estimands.

missing_by_cov <- data.frame(
  covariate = COVS,
  pct_missing = round(100 * vapply(COVS, function(v) mean(is.na(tp[[v]])), numeric(1)), 2),
  row.names = NULL)
print(missing_by_cov, row.names = FALSE)

NEEDS_FLAG <- COVS[vapply(COVS, function(v) anyNA(tp[[v]]), logical(1))]
FILL <- lapply(NEEDS_FLAG, function(v) {
  x <- tp[[v]]
  if (is.character(x) || is.factor(x)) names(sort(table(x), decreasing = TRUE))[1]
  else stats::median(x, na.rm = TRUE)
})
names(FILL) <- NEEDS_FLAG

for (v in NEEDS_FLAG) {
  flag <- paste0(v, "_measured")
  tp[[flag]] <- as.integer(!is.na(tp[[v]]))
  tp[[v]][is.na(tp[[v]])] <- FILL[[v]]
}
FLAGS <- paste0(NEEDS_FLAG, "_measured")
COVS  <- c(COVS, FLAGS)
# The three clock terms enter as NATURAL SPLINES, not linear.
#
# WHY, measured 2026-09-09: with window_start_hr linear, the case-mix-adjusted
# hazard curve in section 10b below is a straight line in logit space and cannot bend.
# The observed hazard of `no fentanyl` -> `bolus only` FALLS from 0.133 at hour 0
# to 0.049 at hour 64; the linear-time adjusted curve ROSE across the same span.
# That contradiction was misspecification, not a case-mix story, and it would
# have been read as one. Time in an ICU course is not linear on the logit scale
# for anything -- early boluses are peri-intubation and late ones are not the
# same act -- so all three history terms get the same treatment.
#
# df = 4 is 3 interior knots at the quartiles: enough to turn once or twice
# across 72h, far short of the 18-window saturation that would just redraw the
# observed curve. Held in the config so every site splines identically.
NS_DF <- config$model$spline_df
stopifnot("model.spline_df must be an integer >= 3" =
            is.numeric(NS_DF) && NS_DF >= 3 && NS_DF == as.integer(NS_DF))
spl <- function(v) sprintf("ns(%s, df = %d)", v, NS_DF)
FORM <- as.formula(paste("state_next ~",
                         paste(c(vapply(HISTORY, spl, character(1)),
                                 setdiff(COVS, HISTORY)), collapse = " + ")))

cat(sprintf("\n  %d covariates carry missingness; a _measured flag was added for each\n",
            length(NEEDS_FLAG)))
cat(sprintf("  %s: %s\n", "flagged", paste(NEEDS_FLAG, collapse = ", ")))

# How often each origin has each covariate measured -- the reason for doing this.
cc <- do.call(rbind, lapply(levels(tp$state), function(s) {
  i <- tp$state == s
  d <- data.frame(origin = s, n_rows = sum(i))
  for (v in NEEDS_FLAG) d[[paste0("pct_", v)]] <-
    round(100 * mean(tp[[paste0(v, "_measured")]][i]), 1)
  d
}))
cat("\n% of rows with each covariate actually measured, by origin\n")
print(cc, row.names = FALSE)

stopifnot("the missing-indicator fill must leave no NA in the model frame" =
            !anyNA(tp[, c("state", "state_next", COVS)]),
          "no row may be lost" = nrow(tp) == sum(!is.na(tp$state_next)))
cat(sprintf("\n  ALL %s transition rows enter the models\n",
            format(nrow(tp), big.mark = ",")))

ORIGINS <- levels(droplevels(tp$state))

fit_origin <- function(d) {
  d$state_next <- droplevels(d$state_next)
  if (REF %in% levels(d$state_next)) d$state_next <- relevel(d$state_next, ref = REF)
  multinom(FORM, data = d, maxit = 800, trace = FALSE)
}

cat("\nFitting one multinomial model per origin state ...\n")
t0 <- Sys.time()
fits <- lapply(ORIGINS, function(s) fit_origin(tp[tp$state == s, ]))
names(fits) <- ORIGINS
cat(sprintf("  %d models in %.1f min\n", length(fits),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

model_fit <- do.call(rbind, lapply(ORIGINS, function(s) {
  d <- tp[tp$state == s, ]; f <- fits[[s]]
  data.frame(origin = s,
             n_transitions  = nrow(d),
             n_episodes     = length(unique(d$encounter_block)),
             n_destinations = nlevels(droplevels(d$state_next)),
             deviance       = round(f$deviance, 1),
             aic            = round(f$AIC, 1),
             n_coefficients = length(coef(f)),
             converged      = as.integer(f$convergence == 0))
}))
print(model_fit, row.names = FALSE)
if (any(model_fit$converged != 1)) {
  warning("these origin models did not converge: ",
          paste(model_fit$origin[model_fit$converged != 1], collapse = ", "))
}
cat("  NOTE: an origin with few at-risk windows has the least well estimated\n")
cat("  coefficients -- read n_transitions alongside them.\n")


# ---- 8. Cluster bootstrap ----------------------------------------------------
# Resample PATIENTS with replacement and refit ALL FIVE models per replicate.
# 17 transitions from one episode are correlated and 932 of 13,627 patients
# contribute more than one episode, so naive multinom standard
# errors are too narrow. sandwich has no estfun method for multinom.

flat <- function(fl) unlist(lapply(names(fl), function(s) {
  cf <- coef(fl[[s]])
  if (is.null(dim(cf))) cf <- matrix(cf, nrow = 1,
                                     dimnames = list(setdiff(levels(
                                       fl[[s]]$lev), REF)[1], names(cf)))
  setNames(as.numeric(cf),
           paste(s, rep(rownames(cf), times = ncol(cf)),
                 rep(colnames(cf), each = nrow(cf)), sep = "|"))
}))

est <- flat(fits)
pt  <- unique(tp$patient_id)
idx <- split(seq_len(nrow(tp)), tp$patient_id)

cat(sprintf("\nCluster bootstrap: %d replicates over %s patients ...\n",
            BOOT_REPS, format(length(pt), big.mark = ",")))
t0 <- Sys.time()
boot <- vapply(seq_len(BOOT_REPS), function(b) {
  take <- sample(pt, length(pt), replace = TRUE)
  s <- tp[unlist(idx[as.character(take)], use.names = FALSE), ]
  r <- try({
    fl <- lapply(ORIGINS, function(o) fit_origin(s[s$state == o, ]))
    names(fl) <- ORIGINS
    v <- flat(fl)
    v[names(est)]          # align; a destination unseen in this replicate is NA
  }, silent = TRUE)
  if (inherits(r, "try-error")) rep(NA_real_, length(est)) else as.numeric(r)
}, numeric(length(est)))

cat(sprintf("  %.1f min; %d of %d replicates converged\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            sum(!is.na(boot[1, ])), BOOT_REPS))

parts <- do.call(rbind, strsplit(names(est), "|", fixed = TRUE))
coefs <- data.frame(
  origin      = parts[, 1],
  destination = parts[, 2],
  term        = parts[, 3],
  estimate    = as.numeric(est),
  boot_se     = apply(boot, 1, stats::sd, na.rm = TRUE),
  lower       = apply(boot, 1, stats::quantile, 0.025, na.rm = TRUE),
  upper       = apply(boot, 1, stats::quantile, 0.975, na.rm = TRUE))
coefs$odds_ratio <- round(exp(coefs$estimate), 3)
coefs$or_lower   <- round(exp(coefs$lower), 3)
coefs$or_upper   <- round(exp(coefs$upper), 3)
coefs[c("estimate", "boot_se", "lower", "upper")] <-
  round(coefs[c("estimate", "boot_se", "lower", "upper")], 4)


# ---- 9. Predicted transition probabilities -----------------------------------
# More legible than 250-odd coefficients: hold the covariates fixed and read off
# where each state goes next. Predictions come from that origin's OWN model, so
# each row is a probability distribution over the destinations observed from it;
# a destination never seen from an origin is a structural zero, not a small
# estimate.

# Built from COVS so the profile can never drift from the model formula -- a
# hand-written list here is how a new covariate produces "object not found" at
# the prediction step, after the fits have already run.
profile_at <- function(sofa, nee_val, label) {
  vals <- lapply(COVS, function(v) {
    x <- tp[[v]]
    if (is.character(x) || is.factor(x)) names(sort(table(x), decreasing = TRUE))[1]
    else stats::median(x, na.rm = TRUE)
  })
  names(vals) <- COVS
  # Predict for a patient whose covariates WERE measured: the flags are set to 1,
  # so the profile is a real clinical scenario rather than an average over
  # measured and unmeasured windows.
  for (f in FLAGS) vals[[f]] <- 1L
  vals$nee <- nee_val
  d <- as.data.frame(vals, stringsAsFactors = FALSE)
  d$label <- label
  d
}

qq <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
# Severity is indexed by NEE and oxygenation, the two explicit markers, now that
# SOFA is out of the model.
sev <- function(p, label) {
  d <- profile_at(NULL, qq(tp$nee, p), label)
  d$oxygenation <- qq(tp$oxygenation, 1 - p)   # worse oxygenation = LOWER P/F
  d$lactate     <- qq(tp$lactate, p)
  d
}
PROFILES <- rbind(sev(.10, "least sick (p10)"),
                  sev(.50, "median (p50)"),
                  sev(.90, "sickest (p90)"))

predict_all <- function(prof) {
  do.call(rbind, lapply(ORIGINS, function(s) {
    pr <- predict(fits[[s]], newdata = prof, type = "probs")
    if (is.null(dim(pr))) pr <- setNames(as.numeric(pr), fits[[s]]$lev)
    p <- setNames(rep(0, length(LEVELS)), LEVELS)
    p[names(pr)] <- as.numeric(pr)
    data.frame(profile = prof$label, from = s, to = LEVELS,
               probability = round(as.numeric(p), 5), row.names = NULL)
  }))
}

predicted_all <- do.call(rbind, lapply(seq_len(nrow(PROFILES)),
                                       function(i) predict_all(PROFILES[i, ])))
predicted <- predicted_all[predicted_all$profile == "median (p50)", ]

# Each row must be a distribution. Tolerance is 1e-4, not 1e-6: the stored
# probabilities are rounded to 5 dp, and seven rounded values can drift further
# than 1e-6 from 1 without anything being wrong.
chk <- tapply(predicted_all$probability,
              paste(predicted_all$profile, predicted_all$from), sum)
stopifnot("each origin's predicted probabilities must sum to 1" =
            all(abs(chk - 1) < 1e-4))

cat("\nPredicted next-state probability, median covariate profile\n")
print(reshape(predicted[, c("from", "to", "probability")],
              idvar = "from", timevar = "to", direction = "wide"),
      row.names = FALSE)

cat("\nSeverity gradient: P(died next window) and P(extubated next window)\n")
g <- predicted_all[predicted_all$to %in% c("died", "extubated"), ]
print(reshape(g[, c("from", "profile", "to", "probability")],
              idvar = c("from", "to"), timevar = "profile", direction = "wide"),
      row.names = FALSE)


# ---- 10. Figure --------------------------------------------------------------

ink <- "#0b0b0b"; muted <- "#898781"; gridline <- "#e1e0d9"
predicted$to <- factor(predicted$to, levels = LEVELS)
predicted$from <- factor(as.character(predicted$from), levels = LEVELS)

p_heat <- ggplot(predicted, aes(to, from, fill = probability)) +
  geom_tile(colour = "#fcfcfb", linewidth = 0.6) +
  # "<0.01" rather than "0.00": several real cells sit near 1e-3 (you rarely leave
  # hospital directly from a ventilated state) and printing them as zero reads as
  # "impossible" when it means "small".
  geom_text(aes(label = ifelse(probability < 0.005, "<0.01",
                               sprintf("%.2f", probability)),
                colour = probability > 0.5), size = 3.1, show.legend = FALSE) +
  scale_fill_gradient(low = "#f2f6fb", high = "#14427e", limits = c(0, 1)) +
  # a discrete y axis puts level 1 at the BOTTOM; reverse it so the matrix reads
  # top-to-bottom in state order, matching the printed table
  scale_y_discrete(limits = rev) +
  scale_colour_manual(values = c(`FALSE` = ink, `TRUE` = "white")) +
  labs(title = sprintf("Probability of moving to each state in the next %d hours",
                       WINDOW_H),
       subtitle = paste0(
         "Median covariate profile. Rows are the current state, columns the next.\n",
         "A 4h probability of 0.005 compounds to roughly 9% over the 18-window course."),
       x = "State at t + 1", y = "State at t", fill = "P") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(colour = ink, face = "bold"),
        plot.subtitle = element_text(colour = muted, margin = margin(b = 10)),
        axis.title = element_text(colour = muted),
        axis.text = element_text(colour = muted),
        axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA))

ggsave(file.path(dirs$phase, sprintf("phase5_%spredicted_transitions.png", TAG)), p_heat,
       width = 8.0, height = 5.2, dpi = 200)

# The severity gradient is where the clinical content is, and a single-profile
# heatmap hides it entirely.
grad <- predicted_all[predicted_all$to %in%
                        c("no fentanyl", "continuous only", "extubated", "died"), ]
grad$profile <- factor(grad$profile, levels = PROFILES$label)
grad$from <- factor(grad$from, levels = LEVELS)
grad$to   <- factor(grad$to, levels = LEVELS)

p_grad <- ggplot(grad, aes(profile, probability, group = from, colour = from)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  facet_wrap(~ to, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c(
    "no fentanyl" = "#898781", "continuous only" = "#14427e",
    "bolus only" = "#eb6834", "continuous + bolus" = "#4a8bd8",
    "extubated" = "#7fb069")) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(title = "The same transitions across a severity gradient",
       subtitle = sprintf(
         "NEE, oxygenation and lactate at their p10, p50 and p90. Panels are the DESTINATION; colour the origin.\nNEE %.2f / %.2f / %.2f mcg/kg/min; P/F %.0f / %.0f / %.0f.",
         qq(tp$nee,.10), qq(tp$nee,.50), qq(tp$nee,.90),
         qq(tp$oxygenation,.90), qq(tp$oxygenation,.50), qq(tp$oxygenation,.10)),
       x = "Severity profile", y = sprintf("P(destination) in the next %dh", WINDOW_H),
       colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(colour = ink, face = "bold"),
        plot.subtitle = element_text(colour = muted, margin = margin(b = 10)),
        legend.position = "top",
        legend.text = element_text(colour = muted, size = 9),
        axis.title = element_text(colour = muted),
        axis.text = element_text(colour = muted),
        axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid.minor = element_blank(),
        strip.text = element_text(colour = ink, face = "bold"),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA))

ggsave(file.path(dirs$phase, sprintf("phase5_%sseverity_gradient.png", TAG)), p_grad,
       width = 9.5, height = 4.6, dpi = 200)

# ---- 10b. Transition hazards over the ventilation course ---------------------
# The discrete-time analogue of Lyons et al., Crit Care Explor 2022;4:e0784,
# Figure 3: one panel per ORIGIN state, every competing destination overlaid,
# read against a clock. Their estimator is Nelson-Aalen on 12h-gridded event
# times; ours is the discrete-time cause-specific hazard
#     h_ij(t) = P(X_{t+1} = j | X_t = i)
# which is the same quantity on a different clock (h ~ 1 - exp(-alpha*W)) and is
# exactly what multinom fits. Fitted probabilities are invariant to the choice
# of REF, so this figure does not depend on that decision.
#
# TIME AXIS: hours since intubation, NOT Lyons's time-since-entry-to-state.
# Their clock-reset is unusable for two of our origins -- MEASURED 2026-09-09,
# at-risk windows by hours in the current state:
#     bolus only          8,670 -> 667 at 12h -> 74 at 24h
#     continuous + bolus  6,441 -> 170 at 12h -> 17 at 24h
# because their states are durable conditions (KDIGO stage persists until the
# kidney changes) and two of ours are punctate events (a bolus has no duration).
# On hours-since-intubation every origin holds >= 298 at-risk windows across the
# whole 0-68h span. Time-in-state is retained as a covariate, so the clock-reset
# information is in the model even though it is not the axis.
#
# TWO CURVES, deliberately:
#   observed  -- empirical proportions, no model. Stands on its own if a reader
#                rejects the multinomial entirely. This is the Lyons estimator.
#   adjusted  -- the fitted model averaged over a case-mix held FIXED across
#                windows (g-computation / marginal standardisation). The two
#                diverge exactly where the surviving cohort's case-mix shifts,
#                which is the point of showing both: `observed` confounds "what
#                happens at hour 60" with "who is still ventilated at hour 60".

HAZ_BOOT <- as.integer(Sys.getenv("HAZARD_BOOT_REPS", "500"))
STD_N    <- 2000L   # standardisation sample per origin

# -- observed: empirical cause-specific hazard, cluster-bootstrapped over patients
# `grid` is passed in, not read off `d`: a bootstrap resample can miss a window
# entirely, and a table that silently loses a row breaks the alignment with the
# point estimate rather than erroring anywhere near the cause.
haz_counts <- function(d, grid) {
  tb <- table(factor(d$window_start_hr, levels = grid),
              factor(d$state_next, levels = LEVELS))
  sweep(tb, 1, pmax(rowSums(tb), 1), "/")
}

haz_observed <- do.call(rbind, lapply(ORIGINS, function(s) {
  d <- tp[tp$state == s, ]
  grid <- sort(unique(d$window_start_hr))
  p <- haz_counts(d, grid)
  den <- as.integer(table(factor(d$window_start_hr, levels = grid)))

  # resample PATIENTS, not windows -- the same clustering the coefficient
  # bootstrap uses. No refitting here, so 500 replicates costs seconds.
  ptn <- unique(d$patient_id)
  ix  <- split(seq_len(nrow(d)), d$patient_id)
  bs <- vapply(seq_len(HAZ_BOOT), function(b) {
    take <- sample(ptn, length(ptn), replace = TRUE)
    as.numeric(haz_counts(d[unlist(ix[as.character(take)], use.names = FALSE), ], grid))
  }, numeric(length(p)))

  q <- function(a) matrix(apply(bs, 1, stats::quantile, a, na.rm = TRUE),
                          nrow = nrow(p), dimnames = dimnames(p))
  lo <- q(0.025); hi <- q(0.975)

  data.frame(
    curve = "observed", origin = s,
    window_start_hr = as.integer(rep(rownames(p), times = ncol(p))),
    destination = rep(colnames(p), each = nrow(p)),
    n_at_risk = rep(den, times = ncol(p)),
    hazard = round(as.numeric(p), 5),
    lower  = round(as.numeric(lo), 5),
    upper  = round(as.numeric(hi), 5),
    row.names = NULL)
}))

# -- adjusted: model-averaged over a case-mix frozen at the origin's own
# distribution, so only the clock moves. Sampling the standardisation set once
# per origin (not per window) is what freezes it.
haz_adjusted <- do.call(rbind, lapply(ORIGINS, function(s) {
  d <- tp[tp$state == s, ]
  std <- d[sample(nrow(d), min(STD_N, nrow(d))), COVS, drop = FALSE]
  hrs <- sort(unique(d$window_start_hr))
  do.call(rbind, lapply(hrs, function(w) {
    nd <- std; nd$window_start_hr <- w
    pr <- predict(fits[[s]], newdata = nd, type = "probs")
    if (is.null(dim(pr))) pr <- matrix(pr, ncol = 1,
                                       dimnames = list(NULL, fits[[s]]$lev[2]))
    m <- setNames(rep(0, length(LEVELS)), LEVELS)
    m[colnames(pr)] <- colMeans(pr)
    data.frame(curve = "adjusted", origin = s, window_start_hr = w,
               destination = LEVELS, n_at_risk = sum(d$window_start_hr == w),
               hazard = round(as.numeric(m), 5),
               lower = NA_real_, upper = NA_real_, row.names = NULL)
  }))
}))

hazards <- rbind(haz_observed, haz_adjusted)

# Windows nobody occupies are dropped, not reported as zero: at hour 0 every
# episode is ventilated by construction, so `extubated` has no at-risk set and a
# row of zeros there would read as "extubated patients never move", which is the
# opposite of the truth.
hazards <- hazards[hazards$n_at_risk > 0, ]

# Check the distributions BEFORE suppression -- afterwards a fully suppressed
# window sums to 0 and the assertion would fire on its own tidying.
chk_h <- tapply(hazards$hazard,
                paste(hazards$curve, hazards$origin, hazards$window_start_hr),
                sum, na.rm = TRUE)
stopifnot(
  "each origin-window must be a distribution over destinations" =
    all(abs(chk_h - 1) < 1e-3),
  "adjusted curves must cover the same grid as observed" =
    nrow(haz_adjusted) > 0 && all(unique(haz_observed$origin) %in% haz_adjusted$origin))

# Structural zeros are real (you do not go from ventilated straight home), but a
# hazard estimated off a handful of windows is not reportable.
n_sup <- sum(hazards$n_at_risk < MIN_CELL)
hazards$hazard[hazards$n_at_risk < MIN_CELL] <- NA_real_
if (n_sup) cat(sprintf("  %d hazard cells suppressed (< %d at risk)\n", n_sup, MIN_CELL))

cat(sprintf("\nTransition hazards: %d origin x window x destination rows, %d bootstrap reps\n",
            nrow(hazards), HAZ_BOOT))
cat("  minimum at-risk windows per origin: ")
cat(paste(sprintf("%s %s", ORIGINS,
                  vapply(ORIGINS, function(s)
                    format(min(haz_observed$n_at_risk[haz_observed$origin == s]),
                           big.mark = ","), character(1))), collapse = "; "), "\n")

# -- figure: Lyons Figure 3, one panel per origin, self-transitions dropped
# (they dominate the axis at 0.80-0.97 and carry no information the others lack)
hz <- hazards[as.character(hazards$destination) != as.character(hazards$origin) &
                !is.na(hazards$hazard), ]
hz$origin      <- factor(hz$origin, levels = LEVELS)
hz$destination <- factor(hz$destination, levels = LEVELS)
hz$curve       <- factor(hz$curve, levels = c("observed", "adjusted"))

# Categorical for routes, an ordered ramp for intensity bands -- a reader should
# be able to see escalation as a gradient, and routes are not ordered.
STATE_COL <- if (identical(TAG, "")) {
  c("no fentanyl" = "#898781", "continuous only" = "#14427e",
    "bolus only" = "#eb6834", "continuous + bolus" = "#4a8bd8",
    "extubated" = "#7fb069", "discharged alive" = "#c9a227", "died" = "#a03030")
} else {
  setNames(c("#898781", "#a8c8ee", "#4a8bd8", "#14427e",
             "#7fb069", "#c9a227", "#a03030"), LEVELS)
}

p_haz <- ggplot(hz, aes(window_start_hr, hazard,
                        colour = destination, fill = destination)) +
  geom_ribbon(data = hz[hz$curve == "observed", ],
              aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(aes(linetype = curve), linewidth = 0.8) +
  facet_wrap(~ origin, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = STATE_COL) +
  scale_fill_manual(values = STATE_COL, guide = "none") +
  scale_linetype_manual(values = c(observed = "solid", adjusted = "22")) +
  scale_x_continuous(breaks = seq(0, EXTENT_H, 12)) +
  labs(title = sprintf("Probability of each transition in the next %dh, across the ventilation course",
                       WINDOW_H),
       subtitle = paste0(
         "Panels are the CURRENT state; colour the next state. Staying put is omitted.\n",
         "Solid = observed proportions with a 95% cluster-bootstrap band; dashed = model-adjusted to a fixed case-mix."),
       x = "Hours since intubation", y = sprintf("P(transition) per %dh window",
                                                 WINDOW_H),
       colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(colour = ink, face = "bold"),
        plot.subtitle = element_text(colour = muted, margin = margin(b = 10)),
        legend.position = "top",
        legend.text = element_text(colour = muted, size = 9),
        axis.title = element_text(colour = muted),
        axis.text = element_text(colour = muted),
        panel.grid.minor = element_blank(),
        strip.text = element_text(colour = ink, face = "bold"),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA))

ggsave(file.path(dirs$phase, sprintf("phase5_%stransition_hazards.png", TAG)), p_haz,
       width = 10.5, height = 7.0, dpi = 200)

  list(tag = TAG, label = LABEL, levels = LEVELS, reference = REF,
       tm = tm, counts = counts, coefs = coefs, predicted_all = predicted_all,
       model_fit = model_fit, cc = cc, missing_by_cov = missing_by_cov,
       hazards = hazards, fits = fits, formula = FORM, origins = ORIGINS,
       boot = boot)
}



# ---- 11. Run both definitions, then write ------------------------------------

RESULTS <- list(
  run_definition(derive_states(long), "", "Delivery route", "no fentanyl"),
  run_definition(derive_dose_states(long, DOSE_CUTS, DOSE_LAB), "dose_",
                 sprintf("Intensity band (0 / <=%s / <=%s / >%s mcg per %dh window)",
                         DOSE_CUTS[1], DOSE_CUTS[2], DOSE_CUTS[2], WINDOW_H),
                 DOSE_LAB[1]))

write_out <- function(x, name) {
  write.csv(x, file.path(dirs$phase, name), row.names = FALSE)
  cat(sprintf("written: %s\n", name))
}

cat("\n")
for (R in RESULTS) {
  g <- function(stem) sprintf("phase5_%s%s", R$tag, stem)
  write_out(R$tm,            g("transition_matrix.csv"))
  write_out(R$counts,        g("transition_counts.csv"))
  write_out(R$coefs,         g("model_coefficients.csv"))
  write_out(R$predicted_all, g("predicted_transitions.csv"))
  write_out(R$hazards,       g("transition_hazards.csv"))
  write_out(R$model_fit,     g("model_fit.csv"))
  write_out(R$cc,            g("measurement_by_origin.csv"))
  write_out(R$missing_by_cov, g("covariate_missingness.csv"))
}

# PHI: the fits carry fitted values per episode-window. Both definitions in one
# object, keyed by tag, so a reader cannot pick up one and think it is the other.
saveRDS(setNames(RESULTS, vapply(RESULTS, function(R)
                   if (nzchar(R$tag)) sub("_$", "", R$tag) else "route",
                   character(1))),
        file.path(dirs$out_phi, "transition_model.rds"))
cat("written: transition_model.rds (PHI, both definitions)\n")

write_json(prov, file.path(dirs$phase, "phase5_provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)


# ---- 12. Provenance ----------------------------------------------------------
# Which package versions produced these numbers?

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/06_transition_model.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "06_transition_model_sessioninfo.txt"))
