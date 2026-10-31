# Fentanyl delivery states, shared by the descriptive scripts (02, 04) and the
# tabled transition model.
# ONE definition in one place, or the model stops describing the figure.
#
# Seven mutually exclusive, exhaustive states per episode-window; four describe
# HOW fentanyl was delivered while ventilated, by ROUTE rather than by a dose
# threshold. `extubated` is transient; only `discharged alive` and `died` absorb.

STATE_LEVELS <- c("no fentanyl", "continuous only", "bolus only",
                  "continuous + bolus", "extubated", "discharged alive", "died")
# ABSORBING only. `extubated` is NOT: extubated windows move on, to reintubation
# or out of the hospital, and absorbing it drops those rows from the model.
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

# Consecutive (state at t, state at t+1) pairs. Rows whose CURRENT state is
# ABSORBING are dropped; `extubated` is kept, being transient here.
transition_pairs <- function(d, id = "encounter_block", time = "window_idx") {
  d <- d[order(d[[id]], d[[time]]), ]
  # Levels off the FACTOR, not a global: a frame from derive_dose_states carries
  # the dose levels, which a hardcoded STATE_LEVELS would NA.
  lv <- levels(d$state)
  nxt <- ave(as.character(d$state), d[[id]], FUN = function(v) c(v[-1], NA))
  d$state_next <- factor(nxt, levels = lv)
  out <- d[!is.na(d$state_next) & !(d$state %in% STATE_ABSORBING), ]
  droplevels(out, except = which(names(out) %in% c("state", "state_next")))
}

# State prevalence per window, as a percentage of the rows PASSED IN -- the
# denominator is the caller's at-risk set, and the two differ by half the cohort
# at 72h: all episodes (bands shrink as patients leave; the liberation view), or
# long[long$ventilated, ] (of those still on the vent, what are they getting).
prevalence_by_window <- function(d) {
  lv <- levels(d$state)
  do.call(rbind, lapply(sort(unique(d$window_idx)), function(w) {
    x <- d[d$window_idx == w, ]
    tb <- table(x$state)[lv]
    data.frame(window_idx = w, window_start_hr = x$window_start_hr[1],
               # A FACTOR, carrying the source level order: as character, ggplot
               # re-sorts it and the areas stack unlike the alluvial.
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

# The same transitions LONG, carrying COUNTS, which is what pools: sum n over
# (from, to), then recompute the percentage. The matrix above ships rounded.
transition_counts <- function(d) {
  lv <- levels(d$state)
  stopifnot("state and state_next must share a level set" =
              identical(lv, levels(d$state_next)))
  tm <- table(factor(d$state, levels = lv), factor(d$state_next, levels = lv))
  tm <- tm[rowSums(tm) > 0, , drop = FALSE]
  out <- as.data.frame(as.table(tm), stringsAsFactors = FALSE)
  names(out) <- c("from", "to", "n")
  out$n_at_risk <- as.integer(rowSums(tm)[out$from])
  out$pct <- round(100 * out$n / out$n_at_risk, 2)
  # Factors in the source level order, so a reader sorts as the figures stack.
  out$from <- factor(out$from, levels = lv[lv %in% rownames(tm)])
  out$to   <- factor(out$to, levels = lv)
  out <- out[order(out$from, out$to), c("from", "to", "n", "n_at_risk", "pct")]
  rownames(out) <- NULL
  out
}

# History features, so the model is not forced to be Markov. Past only.
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
# The states above describe HOW fentanyl was delivered; these, HOW MUCH. Both run
# over the same cohort (as Lyons et al., Crit Care Explor 2022;4:e0784, fit two
# multistate models on one), and the terminal three states are shared, so
# everything downstream works on either definition unmodified.
#
# A DECLARED band is not the latent class the tabled gbmt/lcmm work rejected:
# it asserts no
# subpopulations, it only labels, as KDIGO stage cuts a continuous ratio.
#
# window_mcg, NOT total_dose: total_dose is a RATE (mcg/hr), window_mcg the
# AMOUNT (mcg) over the hourly cells the patient was present for -- they differ
# in windows straddling discharge or death. See protocol.json exposure._two_scales.

DOSE_LABELS <- c("zero", "low", "medium", "high")

# Cuts are NOT defaulted: they are federation-critical (protocol.json
# exposure.dose_states._cuts_MUST_BE_ABSOLUTE) and would drift from an edited
# config.
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
  # breaks = (-Inf, 0], (0, c1], (c1, c2], (c2, Inf): zero is its own band, not
  # the bottom of `low`.
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

# ---- Predominant intensity, for stratification --------------------------------
# Over an episode's AT-RISK time -- every ventilated window -- one declared band
# holds the most time and labels the episode. A DECLARED rule, not a latent one.
# NOT A TRAJECTORY: modal time ignores ORDER, so "started high and weaned" and
# "stayed medium throughout" can share a label. One row per episode: the band
# (ORDERED, so a trend test is available), its share of at-risk time, whether
# that share tied, and the at-risk duration -- a coin-flip label is no stratum.
predominant_band <- function(d, cuts, labels = DOSE_LABELS, window_h = 4,
                             tie_break = c("higher", "lower")) {
  tie_break <- match.arg(tie_break)
  need <- c("encounter_block", "imv_status", "window_mcg")
  missing <- setdiff(need, names(d))
  if (length(missing)) {
    stop("predominant_band needs columns the frame does not carry: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  stopifnot("one cut fewer than the positive bands" =
              length(cuts) == length(labels) - 2L)

  vent <- d[!is.na(d$imv_status) & d$imv_status == 1, ]
  stopifnot("no ventilated windows to classify" = nrow(vent) > 0)
  band <- cut(vent$window_mcg, breaks = c(-Inf, 0, cuts, Inf), labels = labels,
              right = TRUE, ordered_result = TRUE)

  tab <- table(vent$encounter_block, band)
  counts <- as.matrix(tab)
  at_risk <- rowSums(counts)

  # Ties resolve toward the HIGHER band, never toward less exposure: max.col
  # "last" takes the rightmost maximum and the columns ascend.
  idx <- max.col(counts, ties.method = if (tie_break == "higher") "last" else "first")
  top <- counts[cbind(seq_len(nrow(counts)), idx)]
  # tied = at least two bands share the maximum
  n_at_max <- rowSums(counts == top)

  data.frame(
    encounter_block = rownames(counts),
    band = factor(labels[idx], levels = labels, ordered = TRUE),
    modal_share = round(top / at_risk, 3),
    tied = n_at_max > 1L,
    at_risk_windows = as.integer(at_risk),
    at_risk_hours = as.integer(at_risk) * window_h,
    row.names = NULL, stringsAsFactors = FALSE)
}
