# Fentanyl delivery states, shared by Phase 1 (description) and Phase 5 (model).
#
# ONE definition, in one place. Phase 1 draws the alluvial and Phase 5 models the
# transitions; if each restated the rule they would drift, and the model would
# stop describing the figure. See docs/design_notes.md section 10 Phase 1.
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
  nxt <- ave(as.character(d$state), d[[id]], FUN = function(v) c(v[-1], NA))
  d$state_next <- factor(nxt, levels = STATE_LEVELS)
  out <- d[!is.na(d$state_next) & !(d$state %in% STATE_ABSORBING), ]
  droplevels(out, except = which(names(out) %in% c("state", "state_next")))
}

# Row-percent transition matrix, with the at-risk count per origin state.
transition_matrix <- function(d) {
  tm <- table(factor(d$state, levels = STATE_LEVELS),
              factor(d$state_next, levels = STATE_LEVELS))
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
