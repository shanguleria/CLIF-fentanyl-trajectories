# Fentanyl delivery states, shared by Phase 1 (description) and Phase 5 (model).
#
# ONE definition, in one place. Phase 1 draws the alluvial and Phase 5 models the
# transitions; if each restated the rule they would drift, and the model would
# stop describing the figure. See docs/design_notes.md section 10 Phase 1.
#
# Seven mutually exclusive, exhaustive states per episode-window. Four describe
# HOW fentanyl was delivered while ventilated; three are terminal and absorb.
# States are defined by the delivery ROUTE, not by a threshold on dose, so there
# are no cut points to defend -- which is the point, given that Phase 3/4 found
# dose level to be continuous and its latent classes a discretisation of it.

STATE_LEVELS <- c("no fentanyl", "continuous only", "bolus only",
                  "continuous + bolus", "extubated", "discharged alive", "died")
STATE_TERMINAL <- c("extubated", "discharged alive", "died")

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
# Rows whose CURRENT state is terminal are dropped: nothing transitions out of
# an absorbing state, and leaving them in would train the model on rows whose
# outcome is deterministic.
transition_pairs <- function(d, id = "encounter_block", time = "window_idx") {
  d <- d[order(d[[id]], d[[time]]), ]
  nxt <- ave(as.character(d$state), d[[id]], FUN = function(v) c(v[-1], NA))
  d$state_next <- factor(nxt, levels = STATE_LEVELS)
  out <- d[!is.na(d$state_next) & !(d$state %in% STATE_TERMINAL), ]
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
