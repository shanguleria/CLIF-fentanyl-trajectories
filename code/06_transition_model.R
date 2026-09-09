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
              "phase5_predicted_transitions.png", "phase5_model_fit.csv",
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


# ---- 7. Model ----------------------------------------------------------------
# multinom estimates, for every state other than the reference, the log-odds of
# moving THERE rather than to the reference, given the covariates. History enters
# as covariates -- hours_in_state, cumulative_dose, prior_state -- which is what
# makes this a non-Markov model rather than a relaxation of one.
#
# Reference level for the OUTCOME is "no fentanyl": the modal destination and the
# clinically natural comparator ("relative to fentanyl being stopped").

tp$state_next <- relevel(droplevels(tp$state_next), ref = "no fentanyl")
tp$state      <- droplevels(tp$state)

FORM <- state_next ~ state + hours_in_state + cumulative_dose +
                     window_start_hr + sofa_total + nee + age + cci + sex

cat("\nFitting multinomial transition model ...\n")
t0 <- Sys.time()
fit <- multinom(FORM, data = tp, maxit = 500, trace = FALSE)
cat(sprintf("  fitted in %.1f min  (%d observations, %d coefficients)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            nrow(fit$fitted.values), length(coef(fit))))

model_fit <- data.frame(
  n_transitions = nrow(tp),
  n_episodes    = length(unique(tp$encounter_block)),
  n_patients    = length(unique(tp$patient_id)),
  deviance      = round(fit$deviance, 1),
  aic           = round(fit$AIC, 1),
  n_coefficients = length(coef(fit)),
  converged     = as.integer(fit$convergence == 0))
print(model_fit, row.names = FALSE)
if (fit$convergence != 0) {
  warning("multinom did not converge; raise maxit before reading the coefficients")
}


# ---- 8. Cluster bootstrap ----------------------------------------------------
# Resample PATIENTS with replacement, refit, and take percentile intervals. This
# is the honest interval here: 17 transitions from one episode are correlated,
# and 932 of 13,627 patients contribute more than one episode (section 11).
# Naive multinom SEs ignore both and are too narrow.

pt <- unique(tp$patient_id)
cat(sprintf("\nCluster bootstrap: %d replicates over %s patients ...\n",
            BOOT_REPS, format(length(pt), big.mark = ",")))
t0 <- Sys.time()
idx <- split(seq_len(nrow(tp)), tp$patient_id)

boot <- vapply(seq_len(BOOT_REPS), function(b) {
  take <- sample(pt, length(pt), replace = TRUE)
  s <- tp[unlist(idx[as.character(take)], use.names = FALSE), ]
  m <- try(multinom(FORM, data = s, maxit = 500, trace = FALSE), silent = TRUE)
  if (inherits(m, "try-error")) rep(NA_real_, length(coef(fit)))
  else as.numeric(coef(m))
}, numeric(length(coef(fit))))

cat(sprintf("  %.1f min; %d of %d replicates converged\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            sum(!is.na(boot[1, ])), BOOT_REPS))

est <- coef(fit)
coefs <- data.frame(
  destination = rep(rownames(est), times = ncol(est)),
  term        = rep(colnames(est), each = nrow(est)),
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
# More legible than a coefficient table: hold the covariates at the cohort median
# and read off where each state goes next.

ref <- data.frame(
  hours_in_state  = median(tp$hours_in_state, na.rm = TRUE),
  cumulative_dose = median(tp$cumulative_dose, na.rm = TRUE),
  window_start_hr = median(tp$window_start_hr, na.rm = TRUE),
  sofa_total      = median(tp$sofa_total, na.rm = TRUE),
  nee             = median(tp$nee, na.rm = TRUE),
  age             = median(tp$age, na.rm = TRUE),
  cci             = median(tp$cci, na.rm = TRUE),
  sex             = names(sort(table(tp$sex), decreasing = TRUE))[1])

grid <- do.call(rbind, lapply(levels(tp$state), function(st)
  cbind(data.frame(state = factor(st, levels = levels(tp$state))), ref)))
pred <- as.data.frame(predict(fit, newdata = grid, type = "probs"))
pred$state <- grid$state

predicted <- do.call(rbind, lapply(seq_len(nrow(pred)), function(i)
  data.frame(from = pred$state[i],
             to = setdiff(names(pred), "state"),
             probability = round(as.numeric(pred[i, setdiff(names(pred), "state")]), 4))))

cat("\nPredicted next-state probability at the cohort median covariate profile\n")
print(reshape(predicted, idvar = "from", timevar = "to", direction = "wide"),
      row.names = FALSE)


# ---- 10. Figure --------------------------------------------------------------

ink <- "#0b0b0b"; muted <- "#898781"; gridline <- "#e1e0d9"
predicted$to <- factor(predicted$to, levels = STATE_LEVELS)
predicted$from <- factor(as.character(predicted$from), levels = STATE_LEVELS)

p_heat <- ggplot(predicted, aes(to, from, fill = probability)) +
  geom_tile(colour = "#fcfcfb", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", probability),
                colour = probability > 0.5), size = 3.1, show.legend = FALSE) +
  scale_fill_gradient(low = "#f2f6fb", high = "#14427e", limits = c(0, 1)) +
  # a discrete y axis puts level 1 at the BOTTOM; reverse it so the matrix reads
  # top-to-bottom in state order, matching the printed table
  scale_y_discrete(limits = rev) +
  scale_colour_manual(values = c(`FALSE` = ink, `TRUE` = "white")) +
  labs(title = sprintf("Where an episode goes in the next %d hours", WINDOW_H),
       subtitle = paste0(
         "Predicted probability at the cohort median covariate profile.\n",
         "Rows are the current state, columns the next. Terminal states never originate."),
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
write_out(predicted, "phase5_predicted_transitions.csv")
write_out(model_fit, "phase5_model_fit.csv")

# PHI: the fit carries fitted values per episode-window.
saveRDS(list(fit = fit, formula = FORM, boot = boot),
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
