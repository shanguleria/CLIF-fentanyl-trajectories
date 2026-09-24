# Fentanyl delivery states, shared by Phase 1 (description) and Phase 5 (model).
#
# ONE definition, in one place. Phase 1 draws the alluvial and Phase 5 models the
# transitions; if each restated the rule they would drift, and the model would
# stop describing the figure.
#
# Seven mutually exclusive, exhaustive states per episode-window. Four describe
# HOW fentanyl was delivered while ventilated. `extubated` is transient (patients
# are reintubated, or leave); only `discharged alive` and `died` absorb.
# States are defined by the delivery ROUTE, not by a threshold on dose, so there
# are no cut points to defend -- which is the point, given that Phase 3/4 found
# dose level to be continuous and its latent classes a discretisation of it.

STATE_LEVELS <- c("no fentanyl", "continuous only", "bolus only",
                  "continuous + bolus", "extubated", "discharged alive", "died")
# ABSORBING states only. `extubated` is NOT one: 3.12% of extubated windows move
# on -- 1.28% back to a ventilated fentanyl state (reintubation) and 1.84% out of
# the hospital. Measured 2026-09-09 over 60,685 extubated windows. Treating it as
# absorbing dropped every extubated-origin row from the transition model, which
# removed reintubation and, worse, removed the extubated -> died / -> discharged
# transitions where most hospital exits actually happen.
STATE_ABSORBING <- c("discharged alive", "died")

# Needs: alive_admitted, imv_status, inf_dose, bolus_dose, died
derive_states <- function(d) {
  need <- c("alive_admitted", "imv_status", "inf_dose", "bolus_dose", "died")
  missing <- setdiff(need, names(d))
  if (length(missing)) {
    stop("derive_states needs columns the frame does not carry: ",
         paste(missing, collapse = ", "),
         ". Re-run code/01_build_cohort.py, or widen the col_select.",
         call. = FALSE)
  }
  vent <- !is.na(d$imv_status) & d$imv_status == 1
  s <- ifelse(
    !d$alive_admitted, ifelse(d$died, "died", "discharged alive"),
    ifelse(!vent, "extubated",
      ifelse(d$inf_dose > 0 & d$bolus_dose > 0, "continuous + bolus",
        ifelse(d$inf_dose > 0, "continuous only",
          ifelse(d$bolus_dose > 0, "bolus only", "no fentanyl")))))
  d$ventilated <- vent
  d$state <- factor(s, levels = STATE_LEVELS)
  stopifnot("every episode-window must land in exactly one state" =
              !any(is.na(d$state)))
  d
}

# Consecutive (state at t, state at t+1) pairs, one row per episode-window pair.
# Rows whose CURRENT state is ABSORBING are dropped: nothing transitions out of
# one, and leaving them in would train the model on rows whose outcome is
# deterministic. `extubated` is deliberately kept -- it is transient here, and
# the extubated -> discharged / -> died / -> reintubated moves are the ones the
# liberation story turns on.
transition_pairs <- function(d, id = "encounter_block", time = "window_idx") {
  d <- d[order(d[[id]], d[[time]]), ]
  # Levels come off the FACTOR, not a global: a data frame built by
  # derive_dose_states carries the dose levels, and hardcoding STATE_LEVELS here
  # would silently NA every dose band. The factor already knows its own levels.
  lv <- levels(d$state)
  nxt <- ave(as.character(d$state), d[[id]], FUN = function(v) c(v[-1], NA))
  d$state_next <- factor(nxt, levels = lv)
  out <- d[!is.na(d$state_next) & !(d$state %in% STATE_ABSORBING), ]
  droplevels(out, except = which(names(out) %in% c("state", "state_next")))
}

# State prevalence per window, as a percentage of the rows PASSED IN. The
# denominator is therefore the caller's choice of at-risk set, and the two that
# matter differ by more than half the cohort at 72h:
#   prevalence_by_window(long)                  -- all episodes; the bands shrink
#     as patients are extubated, discharged or die, which is the liberation view.
#   prevalence_by_window(long[long$ventilated,]) -- still ventilated only; the
#     bands answer "of the patients still on the vent, what are they getting?"
# One function rather than two so the two figures can never drift in their
# handling, which is why states.R exists at all.
prevalence_by_window <- function(d) {
  lv <- levels(d$state)
  do.call(rbind, lapply(sort(unique(d$window_idx)), function(w) {
    x <- d[d$window_idx == w, ]
    tb <- table(x$state)[lv]
    data.frame(window_idx = w, window_start_hr = x$window_start_hr[1],
               # A FACTOR, carrying the source level order. As a character
               # column ggplot re-sorts it alphabetically, which stacked the
               # area charts in a different order from the alluvial drawn off
               # the same states -- the two figures then disagree about which
               # band is which. Fixed 2026-09-24.
               state = factor(lv, levels = lv), n = as.integer(tb),
               pct = round(100 * as.numeric(tb) / nrow(x), 2))
  }))
}

# Row-percent transition matrix, with the at-risk count per origin state.
transition_matrix <- function(d) {
  lv <- levels(d$state)
  stopifnot("state and state_next must share a level set" =
              identical(lv, levels(d$state_next)))
  tm <- table(factor(d$state, levels = lv), factor(d$state_next, levels = lv))
  keep <- rowSums(tm) > 0
  out <- as.data.frame.matrix(round(100 * prop.table(tm[keep, , drop = FALSE], 1), 2))
  cbind(from = rownames(out), n_at_risk = as.integer(rowSums(tm[keep, , drop = FALSE])),
        out)
}

# History features, so the model is not forced into a Markov assumption. Each is
# computed from the past only -- nothing here may look forward.
add_history <- function(d, id = "encounter_block", time = "window_idx",
                        dose = "total_dose", window_h = 4) {
  d <- d[order(d[[id]], d[[time]]), ]
  # consecutive windows already spent in the CURRENT state, before this one
  run <- ave(as.character(d$state), d[[id]], FUN = function(v) {
    r <- rle(v); rep(sequence(r$lengths) - 1L, 1)
  })
  d$hours_in_state <- as.integer(run) * window_h
  d$cumulative_dose <- ave(d[[dose]], d[[id]],
                           FUN = function(v) cumsum(v) - v) * window_h
  d$prior_state <- factor(
    ave(as.character(d$state), d[[id]], FUN = function(v) c(NA, v[-length(v)])),
    levels = STATE_LEVELS)
  d
}


# ---- Second state definition: delivery INTENSITY ------------------------------
#
# The states above describe HOW fentanyl was delivered. These describe HOW MUCH.
# Both are run over the same cohort, the way Lyons et al. (Crit Care Explor
# 2022;4:e0784) fit two multistate models -- AKI stage, then AKI x IMV -- on one
# cohort to triangulate. The terminal three states are shared, so everything
# downstream (transition_pairs, transition_matrix, add_history, the per-origin
# multinomial, the hazard curves) works on either without modification.
#
# WHY A DECLARED BAND IS NOT THE LATENT CLASS WE REJECTED. Phases 3-4 found dose
# level continuously distributed with no gaps. That kills a LATENT class claim --
# an assertion that k kinds of patient exist, which gbmt made and could not
# support -- but says nothing about a DECLARED band, which asserts nothing and
# only labels. KDIGO stage, the state definition in Lyons, is itself a cut on a
# continuous creatinine ratio.
#
# WHY window_mcg AND NOT total_dose. total_dose is a RATE (mcg/hr); window_mcg is
# the AMOUNT actually delivered (mcg), summed over the hourly cells the patient
# was present for. They differ only in the 632 windows (0.34%) that straddle
# discharge or death, where rate * window_hours would credit a patient with drug
# they were not there to receive. See covariates.json exposure._two_scales.

DOSE_LABELS <- c("zero", "low", "medium", "high")

# Cuts are NOT defaulted. They are federation-critical (covariates.json
# exposure.dose_states._cuts_MUST_BE_ABSOLUTE) and a default here is exactly the
# hardcoded-duplicate-that-drifts failure: a site editing the config would keep
# silently cutting at the old values.
derive_dose_states <- function(d, cuts, labels = DOSE_LABELS) {
  need <- c("alive_admitted", "imv_status", "window_mcg", "died")
  missing <- setdiff(need, names(d))
  if (length(missing)) {
    stop("derive_dose_states needs columns the frame does not carry: ",
         paste(missing, collapse = ", "),
         ". Re-run code/01_build_cohort.py -- window_mcg was added 2026-09-09.",
         call. = FALSE)
  }
  stopifnot(
    "one cut fewer than the positive bands: labels = zero + one per interval" =
      length(cuts) == length(labels) - 2L,
    "cuts must be positive and strictly increasing" =
      all(cuts > 0) && !is.unsorted(cuts, strictly = TRUE),
    "window_mcg must never be negative" = all(d$window_mcg >= 0, na.rm = TRUE))

  lv <- dose_state_levels(labels)
  # breaks = (-Inf, 0], (0, c1], (c1, c2], (c2, Inf) -- zero is its own band, not
  # the bottom of `low`: giving no fentanyl is a different decision from giving a
  # little, and 50.8% of ventilated windows are in it.
  band <- as.character(cut(d$window_mcg, breaks = c(-Inf, 0, cuts, Inf),
                           labels = labels, right = TRUE))
  vent <- !is.na(d$imv_status) & d$imv_status == 1
  s <- ifelse(!d$alive_admitted, ifelse(d$died, "died", "discharged alive"),
              ifelse(!vent, "extubated", band))
  d$ventilated <- vent
  d$state <- factor(s, levels = lv)
  stopifnot("every episode-window must land in exactly one dose state" =
              !any(is.na(d$state)))
  d
}

dose_state_levels <- function(labels = DOSE_LABELS) {
  c(labels, "extubated", "discharged alive", "died")
}

# Both definitions share the terminal states, so STATE_ABSORBING governs both.
stopifnot("the two definitions must share their absorbing states" =
            all(STATE_ABSORBING %in% dose_state_levels()))
