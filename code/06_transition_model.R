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
COV      <- fromJSON(here("config", "covariates.json"), simplifyVector = FALSE)
FENT_U   <- COV$exposure$units

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

OWNED <- list(
  out_phi = c("transition_model.rds"),
  phase   = c("phase5_transition_matrix.csv", "phase5_transition_counts.csv",
              "phase5_model_coefficients.csv", "phase5_predicted_transitions.csv",
              "phase5_predicted_transitions.png", "phase5_severity_gradient.png", "phase5_model_fit.csv",
              "phase5_provenance.json"))
n_cleared <- clear_owned_outputs(dirs, OWNED)
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
                 "total_dose", "died", "age", "sex", "cci", "sofa_total", "nee")))

stopifnot("window count on disk disagrees with the config grid" =
            length(unique(long$window_idx)) == EXTENT_H / WINDOW_H)

long <- derive_states(long)          # shared definition, code/utils/states.R
long <- add_history(long, window_h = WINDOW_H)

n_ep <- length(unique(long$encounter_block))
cat(sprintf("\nAnalytic cohort: %s ventilation episodes x %d windows\n",
            format(n_ep, big.mark = ","), length(unique(long$window_idx))))


# ---- 6. Transition pairs -----------------------------------------------------

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

COVS <- c("hours_in_state", "cumulative_dose", "window_start_hr",
          "sofa_total", "nee", "age", "cci", "sex")
FORM <- as.formula(paste("state_next ~", paste(COVS, collapse = " + ")))
REF  <- "no fentanyl"

tp$state <- droplevels(tp$state)
tp <- tp[stats::complete.cases(tp[, c("state", "state_next", COVS)]), ]
ORIGINS <- levels(tp$state)

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
# contribute more than one episode (section 11), so naive multinom standard
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

profile_at <- function(sofa, nee, label) {
  data.frame(
    label           = label,
    hours_in_state  = median(tp$hours_in_state,  na.rm = TRUE),
    cumulative_dose = median(tp$cumulative_dose, na.rm = TRUE),
    window_start_hr = median(tp$window_start_hr, na.rm = TRUE),
    sofa_total      = sofa,
    nee             = nee,
    age             = median(tp$age, na.rm = TRUE),
    cci             = median(tp$cci, na.rm = TRUE),
    sex             = names(sort(table(tp$sex), decreasing = TRUE))[1],
    stringsAsFactors = FALSE)
}

qq <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
PROFILES <- rbind(
  profile_at(qq(tp$sofa_total, .10), qq(tp$nee, .10), "least sick (p10)"),
  profile_at(qq(tp$sofa_total, .50), qq(tp$nee, .50), "median (p50)"),
  profile_at(qq(tp$sofa_total, .90), qq(tp$nee, .90), "sickest (p90)"))

predict_all <- function(prof) {
  do.call(rbind, lapply(ORIGINS, function(s) {
    pr <- predict(fits[[s]], newdata = prof, type = "probs")
    if (is.null(dim(pr))) pr <- setNames(as.numeric(pr), fits[[s]]$lev)
    p <- setNames(rep(0, length(STATE_LEVELS)), STATE_LEVELS)
    p[names(pr)] <- as.numeric(pr)
    data.frame(profile = prof$label, from = s, to = STATE_LEVELS,
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
predicted$to <- factor(predicted$to, levels = STATE_LEVELS)
predicted$from <- factor(as.character(predicted$from), levels = STATE_LEVELS)

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

ggsave(file.path(dirs$phase, "phase5_predicted_transitions.png"), p_heat,
       width = 8.0, height = 5.2, dpi = 200)

# The severity gradient is where the clinical content is, and a single-profile
# heatmap hides it entirely.
grad <- predicted_all[predicted_all$to %in%
                        c("no fentanyl", "continuous only", "extubated", "died"), ]
grad$profile <- factor(grad$profile, levels = PROFILES$label)
grad$from <- factor(grad$from, levels = STATE_LEVELS)
grad$to   <- factor(grad$to, levels = STATE_LEVELS)

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
         "SOFA and NEE at their p10, p50 and p90. Panels are the DESTINATION; colour the origin.\nSOFA p10/p50/p90 = %.0f / %.0f / %.0f.",
         qq(tp$sofa_total, .10), qq(tp$sofa_total, .50), qq(tp$sofa_total, .90)),
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

ggsave(file.path(dirs$phase, "phase5_severity_gradient.png"), p_grad,
       width = 9.5, height = 4.6, dpi = 200)


# ---- 11. Write ---------------------------------------------------------------

write_out <- function(x, name) {
  f <- file.path(dirs$phase, name)
  write.csv(x, f, row.names = FALSE)
  cat(sprintf("written: %s\n", name))
}

cat("\n")
write_out(tm,        "phase5_transition_matrix.csv")
write_out(counts,    "phase5_transition_counts.csv")
write_out(coefs,     "phase5_model_coefficients.csv")
write_out(predicted_all, "phase5_predicted_transitions.csv")
write_out(model_fit, "phase5_model_fit.csv")

# PHI: the fit carries fitted values per episode-window.
saveRDS(list(fits = fits, formula = FORM, origins = ORIGINS,
             reference = REF, boot = boot, coefficients = coefs),
        file.path(dirs$out_phi, "transition_model.rds"))
cat("written: transition_model.rds (PHI)\n")

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
