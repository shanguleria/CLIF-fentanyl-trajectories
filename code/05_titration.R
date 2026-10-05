# ==============================================================================
# 05_titration.R  --  bolus co-administration at an infusion rate change
#
# Purpose : When the fentanyl infusion rate is increased, how often is a bolus given within +/-30 minutes? A rate change alone approaches the new steady state over roughly four to five half-lives; a bolus at the moment of uptitration gets there in minutes.
# Author  : Shan Guleria
# Created : 2026-09-24
# Inputs  : output/intermediate_phi/titration_{rate,bolus}_events.parquet, trajectory_long.parquet
# Outputs : output/final_no_phi/05_titration/ : the files listed in OWNED
#
# DESCRIPTIVE ONLY. Associations between adherence and average rate, cumulative
# dose or time to extubation are deliberately deferred (SG, 2026-09-24). The
# per-encounter adherence this script writes is shaped so that work needs no
# rebuild -- see section 8.
#
# RAW TIMESTAMPS, NEVER THE HOURLY GRID. The grid bins to whole hours, which
# would move a rate change up to an hour from the bolus that accompanied it and
# make a 30-minute window meaningless. Phase 0 exports the charted times for
# exactly this reason.
# ==============================================================================

# Run this in a FRESH R session (RStudio: Cmd+Shift+F10).

# ---- 1. Packages -------------------------------------------------------------

pkgs <- c("here", "jsonlite", "ggplot2", "arrow")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

source(here("code", "utils", "paths.R"))
source(here("code", "utils", "figures.R"))
# pool_hist() and boot_group(), shared with 02 and 06 so the three cannot drift
# in how they bin a distribution or bootstrap a clustered proportion.
source(here("code", "utils", "pooling.R"))


# ---- 2. Config ---------------------------------------------------------------

config <- fromJSON(here("config", "config.json"), simplifyVector = FALSE)
SEED <- config$model$seed
set.seed(SEED)

WINDOW_H <- config$cohort$window_hours
EXTENT_H <- config$cohort$granular_extent_hours
MIN_CELL <- config$reporting$small_cell_min_den

COV  <- fromJSON(here("config", "covariates.json"), simplifyVector = FALSE)
# Read from the config, never defaulted here: the threshold and the window are
# federation-critical, and a local default is how a site keeps pairing on the
# old numbers after the consortium moves them.
SPEC <- COV$titration
stopifnot("covariates.json must declare a `titration` block" = !is.null(SPEC))
MIN_DELTA <- as.numeric(SPEC$min_rate_change_mcg_hr)
WIN_MIN   <- as.numeric(SPEC$bolus_window_minutes)
WIN_SENS  <- sort(unique(c(as.numeric(unlist(SPEC$window_sensitivity_minutes)), WIN_MIN)))
VENT_ONLY <- isTRUE(SPEC$ventilated_only)
ADH_LM_H  <- as.numeric(SPEC$adherence_landmark_hours)
MIN_EV    <- as.numeric(SPEC$min_events_for_adherence)
BOOT_REPS <- as.integer(SPEC$bootstrap_resamples)
# Histogram bins. Absolute and config-declared, never local quantiles -- a
# site-specific binning cannot be summed with anyone else's.
HSPEC <- COV$pooling$histograms$variables
stopifnot("covariates.json must declare a `titration` block" = !is.null(SPEC),
          "covariates.json must declare pooling.histograms.variables" = !is.null(HSPEC),
          "the primary window must be one of the sensitivity widths" = WIN_MIN %in% WIN_SENS,
          "the threshold must be positive" = MIN_DELTA > 0,
          "bootstrap_resamples must be a positive whole number" = BOOT_REPS > 0)
# This script owns the histograms the config assigns to it, and only those.
HIST_MINE <- names(HSPEC)[vapply(HSPEC, function(s)
  identical(as.character(s$owner), "05"), logical(1))]

EVENT_LEVELS <- c("initiation", "uptitration", "downtitration")

# The indication question: was the dose change prompted by anything documented?
IND <- COV$indication
stopifnot("covariates.json must declare an `indication` block" = !is.null(IND))
PAIN_THRESH <- unlist(IND$pain_instruments)        # NVPS 4, CPOT 3
SED_THRESH  <- unlist(IND$sedation_instruments)    # RASS 1
IND_WIN_H   <- as.numeric(IND$indication_window_hours)
IND_SWEEP_H <- sort(unique(c(as.numeric(unlist(IND$indication_window_sensitivity_hours)),
                             IND_WIN_H)))
stopifnot(
  "at least one pain instrument must be declared" = length(PAIN_THRESH) >= 1,
  "at least one sedation instrument must be declared" = length(SED_THRESH) >= 1,
  "the primary indication window must be one of the sweep widths" =
    IND_WIN_H %in% IND_SWEEP_H)


# ---- 3. Paths and provenance -------------------------------------------------

dirs <- site_dirs()
dirs$phase <- phase_dir(dirs, "05_titration")
prov <- provenance(config)

message(sprintf("[05_titration] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))


# ---- 4. Guards ---------------------------------------------------------------

manifest <- require_manifest(dirs, here())
message(sprintf("  reading Phase 0 outputs from code %s, generated %s",
                manifest$code_version, manifest$generated))

OWNED <- list(
  out_phi = c("titration_adherence.parquet",
              "titration_events_classified.parquet",
              "indication_events.parquet"),
  phase   = c("coadministration.csv", "window_sensitivity.csv",
              "indication.csv", "indication_window_sensitivity.csv",
              "indication.png",
              "charting_precision.csv", "rate_change_magnitude.csv",
              "charting_agreement.csv", "adherence_distribution.csv",
              "pooling_histograms.csv",
              "coadministration.png", "window_sensitivity.png",
              "provenance.json", "captions.md"))
RETIRED <- character(0)          # new script on 2026-09-24
n_cleared <- clear_owned_outputs(dirs, OWNED, retired = RETIRED)
if (n_cleared) message(sprintf("  cleared %d output(s) from a previous run", n_cleared))


# ---- 5. Read -----------------------------------------------------------------

rate <- as.data.frame(read_parquet(file.path(dirs$out_phi, "titration_rate_events.parquet")))
bol  <- as.data.frame(read_parquet(file.path(dirs$out_phi, "titration_bolus_events.parquet")))

# If every timestamp were a whole hour the grid would have leaked in, and a
# 30-minute window would be meaningless. Asserted at export and again here.
frac_sub <- mean(rate$t_hr %% 1 != 0)
stopifnot(
  "rate events look hour-binned, not charted -- the grid has leaked in" = frac_sub > 0.5,
  "the clock must be relative hours" = is.numeric(rate$t_hr) && min(rate$t_hr) >= 0)

asm <- as.data.frame(read_parquet(file.path(dirs$out_phi, "assessment_events.parquet")))
# Scales are NOT comparable across instruments (NVPS 0-10, CPOT 0-8, RASS -5..+4)
# and must never be pooled as raw values. Each is cut at its OWN declared
# threshold; after the cut every instrument contributes the same yes/no.
declared <- c(names(PAIN_THRESH), names(SED_THRESH))
absent <- setdiff(declared, unique(asm$instrument))
if (length(absent))
  cat(sprintf("  NOTE declared but not charted at this site: %s\n",
              paste(absent, collapse = ", ")))
stopifnot("no declared indication instrument is present in the export" =
            length(setdiff(declared, absent)) >= 1)

long <- as.data.frame(read_parquet(
  file.path(dirs$out_phi, "trajectory_long.parquet"),
  col_select = c("encounter_block", "window_idx", "imv_status")))
long$ventilated <- !is.na(long$imv_status) & long$imv_status == 1

cat(sprintf("\n%s rate record(s), %s bolus record(s), %s episodes\n",
            format(nrow(rate), big.mark = ","), format(nrow(bol), big.mark = ","),
            format(length(unique(rate$encounter_block)), big.mark = ",")))


# ---- 6. Classify the events --------------------------------------------------
# THE RATE DELTA DEFINES THE EVENT; THE ACTION LABEL VALIDATES IT. 86.7% of
# fentanyl infusion records are `verify` or `going` -- re-chartings of an
# unchanged rate. They contribute a delta of zero and drop out on their own, so
# the label is not needed to exclude them. Defining events BY the label would
# miss a real change charted under `verify` and invent one from a `dose_change`
# charted at the same rate.
#
# `stop` and `not_administered` set the running rate to zero but are not events:
# stopping an infusion is not titrating it. A later restart from zero is then
# correctly an initiation.

rate <- rate[order(rate$encounter_block, rate$t_hr), ]
prev_rate  <- ave(rate$rate, rate$encounter_block, FUN = function(v) c(NA, v[-length(v)]))
rate$delta <- rate$rate - prev_rate
rate$first_of_block <- is.na(prev_rate)

is_stop <- rate$not_given | rate$action == "stop"
rate$kind <- NA_character_
# An initiation is a rise from a KNOWN zero. The first record of a block has no
# predecessor, so a block that opens mid-infusion is not called an initiation --
# we did not see it start.
rate$kind[!is_stop & !rate$first_of_block &
            prev_rate == 0 & rate$rate > 0]                    <- "initiation"
rate$kind[!is_stop & !rate$first_of_block &
            prev_rate > 0 & rate$delta >= MIN_DELTA]           <- "uptitration"
rate$kind[!is_stop & !rate$first_of_block &
            rate$delta <= -MIN_DELTA & rate$rate > 0]          <- "downtitration"

ev <- rate[!is.na(rate$kind), ]
ev$kind <- factor(ev$kind, levels = EVENT_LEVELS)
stopifnot(
  "an initiation must rise from zero" =
    all(ev$rate[ev$kind == "initiation"] > 0),
  "a titration must clear the declared threshold" =
    all(abs(ev$delta[ev$kind != "initiation"]) >= MIN_DELTA),
  "a stop is never an event" = !any(ev$not_given | ev$action == "stop"))

# Restrict to ventilated windows: the paper describes ventilated patients and
# every other figure uses that denominator. The unrestricted count is reported
# alongside so the restriction is visible rather than silent.
ev$window_idx <- pmin(floor(ev$t_hr / WINDOW_H), EXTENT_H / WINDOW_H - 1)
ev <- merge(ev, long[, c("encounter_block", "window_idx", "ventilated")],
            by = c("encounter_block", "window_idx"), all.x = TRUE)
ev$ventilated[is.na(ev$ventilated)] <- FALSE
n_all <- nrow(ev)
if (VENT_ONLY) ev <- ev[ev$ventilated, ]
cat(sprintf("  events: %s of %s in a ventilated window (%.1f%%)\n",
            format(nrow(ev), big.mark = ","), format(n_all, big.mark = ","),
            100 * nrow(ev) / max(n_all, 1)))


# ---- 7. Pair each event with a bolus -----------------------------------------
# Is there at least one fentanyl bolus within +/- w minutes of the event?

paired_within <- function(events, boluses, minutes) {
  w <- minutes / 60
  b <- split(boluses$t_hr, boluses$encounter_block)
  vapply(seq_len(nrow(events)), function(i) {
    tb <- b[[as.character(events$encounter_block[i])]]
    !is.null(tb) && any(abs(tb - events$t_hr[i]) <= w)
  }, logical(1))
}

ev$paired <- paired_within(ev, bol, WIN_MIN)
ev$on_hour <- (ev$t_hr %% 1) == 0      # the charting-precision stratum

# The unit of analysis is the EVENT, not the encounter: one episode contributes
# several rate changes, so n_events far exceeds the number of patients and the
# percentages are event-level. n_episodes travels beside it so a reader can see
# the clustering rather than having to infer it -- and so nobody reads a
# percentage of events as a percentage of patients.
summarise_pairing <- function(d, by) {
  do.call(rbind, lapply(split(d, d[[by]], drop = TRUE), function(x) {
    if (!nrow(x)) return(NULL)
    data.frame(group = as.character(x[[by]][1]),
               n_events = nrow(x),
               n_episodes = length(unique(x$encounter_block)),
               events_per_episode = round(nrow(x) / length(unique(x$encounter_block)), 2),
               n_paired = sum(x$paired),
               pct_paired = round(100 * mean(x$paired), 1))
  }))
}

# Interval and design effect for ONE clustered proportion. Events cluster
# within episodes, so a binomial interval is far too narrow here; the deff is
# what lets a coordinating centre widen the POOLED interval instead of assuming
# independence across sites. The resampling universe is stated at every call --
# the episodes contributing events to `d` and no others (lessons.md #26).
boot_one <- function(d, outcome_col) {
  eps <- sort(unique(d$encounter_block))
  x <- data.frame(encounter_block = d$encounter_block, g = "all",
                  y = as.integer(d[[outcome_col]]))
  r <- boot_group(x, "g", SEED, outcome = "y", episodes = eps,
                  n_resamples = BOOT_REPS)
  data.frame(ci_lo = round(r$ci_lo[1], 1), ci_hi = round(r$ci_hi[1], 1),
             deff = round(r$deff[1], 3))
}

coad <- summarise_pairing(ev, "kind")
inc <- ev[ev$kind != "downtitration", ]
coad <- rbind(
  coad,
  data.frame(group = "any increase", n_events = nrow(inc),
             n_episodes = length(unique(inc$encounter_block)),
             events_per_episode = round(nrow(inc) / length(unique(inc$encounter_block)), 2),
             n_paired = sum(inc$paired),
             pct_paired = round(100 * mean(inc$paired), 1)))
coad <- cbind(coad, do.call(rbind, lapply(coad$group, function(g) {
  boot_one(if (identical(g, "any increase")) inc else ev[ev$kind == g, ], "paired")
})))
coad$window_minutes <- WIN_MIN

cat(sprintf("\nBolus within +/-%g min of the event\n", WIN_MIN))
print(coad[, c("group", "n_events", "n_episodes", "events_per_episode",
               "n_paired", "pct_paired")], row.names = FALSE)
cat(sprintf("  percentages are EVENT-level; %s episodes contribute the %s increases\n",
            format(length(unique(inc$encounter_block)), big.mark = ","),
            format(nrow(inc), big.mark = ",")))

# The contrast that answers "is this deliberate?": pairing against window width,
# for uptitrations AND downtitrations. If the uptitration curve rises faster at
# short windows, pairing is intentional; if the two run parallel, the apparent
# rate is just bolus frequency.
sens <- do.call(rbind, lapply(WIN_SENS, function(w) {
  p <- paired_within(ev, bol, w)
  do.call(rbind, lapply(EVENT_LEVELS, function(k) {
    i <- ev$kind == k
    if (!any(i)) return(NULL)
    data.frame(window_minutes = w, group = k, n_events = sum(i),
               n_paired = sum(p[i]), pct_paired = round(100 * mean(p[i]), 1))
  }))
}))

cat("\nPairing by window width (%)\n")
sw <- reshape(sens[, c("window_minutes", "group", "pct_paired")],
              idvar = "window_minutes", timevar = "group", direction = "wide")
names(sw) <- sub("^pct_paired\\.", "", names(sw))
print(sw, row.names = FALSE)

# The charting-precision stratum. 15.7% of dose_change/start records land exactly
# on :00, and for those a true pairing could fall outside the window. If they
# pair less at the primary width but catch up as it widens while off-the-hour
# events do not, that gap is the artifact -- measured, not assumed.
prec <- do.call(rbind, lapply(WIN_SENS, function(w) {
  p <- paired_within(ev, bol, w)
  do.call(rbind, lapply(c(FALSE, TRUE), function(oh) {
    i <- ev$on_hour == oh & ev$kind != "downtitration"
    if (!any(i)) return(NULL)
    data.frame(window_minutes = w,
               charted = if (oh) "on the hour" else "off the hour",
               n_events = sum(i), n_paired = sum(p[i]),
               pct_paired = round(100 * mean(p[i]), 1))
  }))
}))

cat("\nPairing by charting precision, increases only (%)\n")
pw <- reshape(prec[, c("window_minutes", "charted", "pct_paired")],
              idvar = "window_minutes", timevar = "charted", direction = "wide")
names(pw) <- sub("^pct_paired\\.", "", names(pw))
print(pw, row.names = FALSE)

# Does a change at or above the threshold actually carry a dose_change/start
# label? Agreement is a finding about the site's charting, not an input.
agree <- as.data.frame(table(action = ev$action, kind = ev$kind))
agree <- agree[agree$Freq > 0, ]
names(agree)[3] <- "n_events"
cat("\nCharted action on events at or above the threshold\n")
print(reshape(agree, idvar = "action", timevar = "kind", direction = "wide"),
      row.names = FALSE)

# The magnitude distribution, so the 25 mcg/hr threshold is checkable against
# the site's own data rather than taken on trust.
d_all <- abs(rate$delta[!is.na(rate$delta) & rate$delta != 0 &
                          !is_stop & !rate$first_of_block])
# `value` carries two different units, so the unit is a column rather than
# being implied by the one it used to be named after (`mcg_per_hr`, which was a
# percentage on the share rows). `n` is the shared denominator -- every row is
# over the same set of non-zero rate changes -- and `n_at_or_above` is the
# numerator the share rows need in order to pool.
.thr <- c(5, 10, 25, 50, 100)
mag <- data.frame(
  statistic = c("p10", "p25", "p50", "p75", "p90", "max"),
  value = round(unname(quantile(d_all, c(.1, .25, .5, .75, .9, 1))), 1),
  unit = "mcg_per_hr",
  n = length(d_all),
  n_at_or_above = NA_integer_)
mag <- rbind(mag, data.frame(
  statistic = paste0("share >= ", .thr),
  value = round(100 * vapply(.thr, function(k) mean(d_all >= k), numeric(1)), 1),
  unit = "percent",
  n = length(d_all),
  n_at_or_above = vapply(.thr, function(k) sum(d_all >= k), integer(1))))
cat(sprintf("\nMagnitude of non-zero rate changes (n = %s), threshold %g mcg/hr\n",
            format(length(d_all), big.mark = ","), MIN_DELTA))
print(mag, row.names = FALSE)


# ---- 8. Per-encounter adherence ----------------------------------------------
# Patient-level, so the table stays PHI-side and only its distribution ships.
# Emitted TWICE -- overall and within the first ADH_LM_H hours. The second is an
# exposure measured strictly before outcome accrual, which is what a later
# association against time-to-extubation needs: an episode-long exposure summary
# against a time-to-event outcome is immortal-time bias.

adherence <- function(d, tag) {
  inc <- d[d$kind != "downtitration", ]
  if (!nrow(inc)) return(NULL)
  agg <- aggregate(cbind(n_events = rep(1, nrow(inc)), n_paired = as.integer(inc$paired))
                   ~ encounter_block, data = inc, FUN = sum)
  agg$adherence <- agg$n_paired / agg$n_events
  agg$scope <- tag
  agg
}
adh <- rbind(adherence(ev, "overall"),
             adherence(ev[ev$t_hr <= ADH_LM_H, ], sprintf("first_%gh", ADH_LM_H)))

adh_dist <- do.call(rbind, lapply(split(adh, adh$scope), function(x) {
  e <- x[x$n_events >= MIN_EV, ]
  data.frame(scope = x$scope[1],
             n_episodes = nrow(x),
             n_episodes_ge_min_events = nrow(e),
             median_adherence = round(median(e$adherence), 3),
             q1 = round(quantile(e$adherence, .25), 3),
             q3 = round(quantile(e$adherence, .75), 3),
             # The NUMERATORS behind the two percentages. Their denominator is
             # n_episodes_ge_min_events, NOT n_episodes -- without the counts a
             # site's rate cannot be pooled, only re-derived from a 1 dp value.
             n_never = sum(e$adherence == 0),
             n_always = sum(e$adherence == 1),
             pct_never = round(100 * mean(e$adherence == 0), 1),
             pct_always = round(100 * mean(e$adherence == 1), 1))
}))
cat(sprintf("\nPer-encounter adherence (episodes with >= %g increases)\n", MIN_EV))
print(adh_dist, row.names = FALSE)


# ---- 8b. Pooling histograms --------------------------------------------------
# A median has no closed form across sites; counts on FIXED absolute bins do.
# Config: covariates.json pooling.histograms. This script emits the two
# variables that block assigns to it, and asserts it emitted all of them -- a
# declared variable nobody writes is the config-integrity failure CLAUDE.md
# names first.
hist_rows <- list()
if ("rate_change_magnitude" %in% HIST_MINE) {
  hist_rows[["rate_change_magnitude"]] <-
    pool_hist("rate_change_magnitude", d_all, HSPEC$rate_change_magnitude)
}
if ("episode_adherence" %in% HIST_MINE) {
  e_all <- adh[adh$scope == "overall" & adh$n_events >= MIN_EV, ]
  hist_rows[["episode_adherence"]] <-
    pool_hist("episode_adherence", e_all$adherence, HSPEC$episode_adherence)
}
pool_histograms <- do.call(rbind, hist_rows)
stopifnot(
  "05 must emit every histogram the config assigns it" =
    setequal(unique(pool_histograms$variable), HIST_MINE))
cat(sprintf("\nPooling histograms: %s (%d bins total)\n",
            paste(HIST_MINE, collapse = ", "), nrow(pool_histograms)))


# ---- 8b. Was the dose change indicated? --------------------------------------
# The other half of the titration question. Coadministration above asks whether
# a bolus ACCOMPANIED an increase; this asks whether anything documented
# PROMPTED it -- a pain score at or above its instrument threshold, or agitation
# at or above the RASS threshold.
#
# A LEVEL, NOT A CHANGE. An NVPS of 4 held across three consecutive readings
# counts at every one of them; a rise from 0 to 3 counts at none. That is the
# right reading of "was there an indication" -- a patient in pain is in pain
# whether or not it is new -- but it is not what "change" means, and the caption
# says LEVEL.
#
# HOURS, NOT MINUTES. Assessments are charted hourly to q4h while titrations are
# charted to roughly the nearest 5 minutes, so the +/-30 min bolus window would
# find almost nothing here and the nothing would be an artifact of cadence.

# Times at which each class of instrument was ABOVE THRESHOLD, and at which it
# was recorded at all. Split once per episode; the same shape paired_within()
# uses for boluses.
qualifying_times <- function(a, thresholds) {
  hit <- rep(FALSE, nrow(a))
  for (nm in names(thresholds)) {
    i <- a$instrument == nm
    hit[i] <- a$value[i] >= thresholds[[nm]]
  }
  split(a$t_hr[hit], a$encounter_block[hit])
}
recorded_times <- function(a, thresholds) {
  i <- a$instrument %in% names(thresholds)
  split(a$t_hr[i], a$encounter_block[i])
}

# TRUE when the episode has at least one listed time within the window of the
# event. `back_only` restricts to times at or before the event, which is the
# literal causal reading -- an indication precedes its titration.
near_any <- function(events, times, hours, back_only = FALSE) {
  vapply(seq_len(nrow(events)), function(i) {
    tt <- times[[as.character(events$encounter_block[i])]]
    if (is.null(tt)) return(FALSE)
    d <- tt - events$t_hr[i]
    if (back_only) any(d <= 0 & d >= -hours) else any(abs(d) <= hours)
  }, logical(1))
}

PAIN_HIT <- qualifying_times(asm, PAIN_THRESH); PAIN_ANY <- recorded_times(asm, PAIN_THRESH)
SED_HIT  <- qualifying_times(asm, SED_THRESH);  SED_ANY  <- recorded_times(asm, SED_THRESH)
# "either" is not a merge of the two lists -- it is the same functions over the
# union of the threshold maps. Each instrument is still cut at its OWN value;
# merging the split lists by name would drop any episode present in one and not
# the other.
BOTH_THRESH <- c(PAIN_THRESH, SED_THRESH)
BOTH_HIT <- qualifying_times(asm, BOTH_THRESH); BOTH_ANY <- recorded_times(asm, BOTH_THRESH)

# Boluses are events here too: SG asked for increases AND boluses.
ev_sets <- list(increase = inc,
                bolus    = data.frame(encounter_block = bol$encounter_block,
                                      t_hr = bol$t_hr))

# TWO DENOMINATORS, the same distinction F3/F5 draw. Against ALL events the
# indicated fraction conflates "no score was taken" with "a score was taken and
# did not justify it"; against events WITH an assessment in the window it
# separates documentation from indication. Neither alone answers the question.
indication_row <- function(d, hours, event_type, class_name, hit, anyrec) {
  ind  <- near_any(d, hit, hours)
  asr  <- near_any(d, anyrec, hours)
  back <- near_any(d, hit, hours, back_only = TRUE)
  data.frame(
    event_type = event_type, class = class_name, window_hours = hours,
    n_events = nrow(d),
    n_episodes = length(unique(d$encounter_block)),
    n_assessed = sum(asr),
    pct_assessed = round(100 * mean(asr), 1),
    n_indicated = sum(ind),
    pct_indicated_all = round(100 * mean(ind), 1),
    pct_indicated_of_assessed = if (sum(asr)) round(100 * sum(ind & asr) / sum(asr), 1) else NA_real_,
    pct_indicated_backward_only = round(100 * mean(back), 1))
}

ind_rows <- list()
for (w in IND_SWEEP_H) {
  for (etype in names(ev_sets)) {
    d <- ev_sets[[etype]]
    ind_rows[[length(ind_rows) + 1]] <- rbind(
      indication_row(d, w, etype, "pain", PAIN_HIT, PAIN_ANY),
      indication_row(d, w, etype, "sedation", SED_HIT, SED_ANY),
      indication_row(d, w, etype, "either", BOTH_HIT, BOTH_ANY))
  }
}
ind_sens <- do.call(rbind, ind_rows)
ind_sens$class <- factor(ind_sens$class, levels = c("pain", "sedation", "either"))
ind <- ind_sens[ind_sens$window_hours == IND_WIN_H, ]

# Bootstrap the PRIMARY window only. The sweep's 24 rows would say nothing
# further about the clustering and would cost 24 x BOOT_REPS to learn it. The
# interval is on pct_indicated_all, which is why the columns say so -- this row
# carries four different percentages and a bare `ci_lo` would be ambiguous.
IND_HIT <- list(pain = PAIN_HIT, sedation = SED_HIT, either = BOTH_HIT)
ind <- cbind(ind, do.call(rbind, lapply(seq_len(nrow(ind)), function(i) {
  d <- ev_sets[[as.character(ind$event_type[i])]]
  b <- boot_one(data.frame(encounter_block = d$encounter_block,
                           y = near_any(d, IND_HIT[[as.character(ind$class[i])]],
                                        IND_WIN_H)), "y")
  names(b) <- paste0("indicated_", names(b))
  b
})))

# The accounting must close, and a wider window can only ADD qualifying scores.
# A fall as the window widens means the join is wrong, not that practice changed.
stopifnot(
  "indicated exceeds the event total" = all(ind_sens$n_indicated <= ind_sens$n_events),
  "assessed exceeds the event total"  = all(ind_sens$n_assessed <= ind_sens$n_events),
  "backward-only exceeds symmetric"   =
    all(ind_sens$pct_indicated_backward_only <= ind_sens$pct_indicated_all + 1e-9))
for (k in split(ind_sens, list(ind_sens$event_type, ind_sens$class), drop = TRUE)) {
  k <- k[order(k$window_hours), ]
  stopifnot("the indicated fraction fell as the window widened" =
              all(diff(k$n_indicated) >= 0))
}

cat("\n", sprintf("Documented indication within +/-%gh of the dose change\n", IND_WIN_H))
print(ind[, c("event_type", "class", "n_events", "pct_assessed",
              "pct_indicated_all", "pct_indicated_of_assessed",
              "pct_indicated_backward_only")], row.names = FALSE)


# ---- 9. Figures --------------------------------------------------------------

register_caption("coadministration.png",
  "Bolus co-administration at a fentanyl infusion rate change",
  paste(
  sprintf(paste0("Share of rate changes accompanied by at least one fentanyl ",
                 "bolus within +/-%g minutes, on raw charted timestamps rather ",
                 "than the hourly grid. A change must be at least %g mcg/hr to ",
                 "count. Initiations are a rise from a documented zero; ",
                 "uptitrations are a rise in a running infusion. ",
                 "DOWNTITRATIONS ARE A NEGATIVE CONTROL: a bolus at a rate ",
                 "decrease has no pharmacologic rationale, so its rate is the ",
                 "floor for coincidental pairing and the increase rates should ",
                 "be read against it, not against zero. %s"),
          WIN_MIN, MIN_DELTA,
          if (VENT_ONLY) "Restricted to ventilated windows." else ""),
    sprintf(paste0("THE DENOMINATOR IS EVENTS, NOT PATIENTS: %s episodes ",
                   "contribute the %s increases, a mean of %.2f each, so these ",
                   "percentages are not the share of patients and carry no valid ",
                   "confidence interval without accounting for clustering within ",
                   "episode."),
            format(length(unique(inc$encounter_block)), big.mark = ","),
            format(nrow(inc), big.mark = ","),
            nrow(inc) / length(unique(inc$encounter_block)))))

p_coad <- house(
  ggplot(coad[coad$group %in% EVENT_LEVELS, ],
         aes(factor(group, levels = EVENT_LEVELS), pct_paired,
             fill = factor(group, levels = EVENT_LEVELS))) +
    geom_col(width = 0.62) +
    geom_text(aes(label = sprintf("%.1f%%\n%s events", pct_paired,
                                  format(n_events, big.mark = ","))),
              vjust = -0.25, size = 3, colour = INK) +
    # No legend: the x axis already names each bar, so a legend would only
    # repeat it. Identity is never carried by colour alone here.
    scale_fill_manual(values = c(initiation = STATE_COLOURS[["continuous only"]],
                                 uptitration = STATE_COLOURS[["bolus only"]],
                                 downtitration = STATE_COLOURS[["discharged alive"]]),
                      guide = "none") +
    scale_x_discrete(labels = c(initiation = "Initiation", uptitration = "Uptitration",
                                downtitration = "Downtitration")) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.18))) +
    # The denominator is EVENTS, not patients -- one episode contributes several
    # rate changes. Naming it on the axis is where a reader is certain to look
    # before quoting the number.
    labs(x = NULL, y = sprintf("%% of Events Paired with Bolus (+/-%g min)", WIN_MIN)))
ggsave(file.path(dirs$phase, "coadministration.png"), p_coad,
       width = 7.5, height = 4.4, dpi = 200)

register_caption("window_sensitivity.png",
  "Does the pairing survive a narrower window?",
  sprintf(paste0("Share of rate changes with a bolus within the window, against ",
                 "window width. The comparison that matters is the SEPARATION ",
                 "between the increase curves and the downtitration control: a ",
                 "gap that is already present at the narrowest widths is ",
                 "deliberate pairing, while curves that rise together are ",
                 "measuring bolus frequency. The primary window is %g minutes. ",
                 "A change must be at least %g mcg/hr to count."),
          WIN_MIN, MIN_DELTA))

p_sens <- house(
  ggplot(sens, aes(window_minutes, pct_paired,
                   colour = factor(group, levels = EVENT_LEVELS))) +
    geom_vline(xintercept = WIN_MIN, colour = MUTED, linetype = "22", linewidth = 0.4) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
    scale_colour_manual(values = c(initiation = STATE_COLOURS[["continuous only"]],
                                   uptitration = STATE_COLOURS[["bolus only"]],
                                   downtitration = STATE_COLOURS[["discharged alive"]]),
                        labels = c(initiation = "Initiation", uptitration = "Uptitration",
                                   downtitration = "Downtitration"),
                        name = NULL) +
    scale_x_continuous(breaks = WIN_SENS) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Pairing window, +/- minutes",
         y = "% of rate-change EVENTS with a bolus"))
ggsave(file.path(dirs$phase, "window_sensitivity.png"), p_sens,
       width = 7.5, height = 4.4, dpi = 200)


register_caption("indication.png",
  "Was the dose change indicated? Documented pain or agitation at a fentanyl increase or bolus",
  paste(
    sprintf(paste0("Share of events with a qualifying score within +/-W hours, ",
                   "against window width. Pain is %s; sedation is %s. Each ",
                   "instrument is cut at ITS OWN threshold -- the scales are not ",
                   "comparable and raw values are never pooled -- after which ",
                   "every instrument contributes the same yes/no."),
            paste(sprintf("%s >= %g", names(PAIN_THRESH), PAIN_THRESH), collapse = ", "),
            paste(sprintf("%s >= %g", names(SED_THRESH), SED_THRESH), collapse = ", ")),
    paste0("THIS IS A LEVEL, NOT A CHANGE. A score held above threshold across ",
           "consecutive readings qualifies at every one of them; a rise that ",
           "stays below threshold qualifies at none."),
    sprintf(paste0("THE WINDOW IS HOURS, NOT MINUTES, because assessments are ",
                   "charted hourly to q4h while titrations are charted to ",
                   "roughly the nearest 5 minutes -- the +/-%g min bolus window ",
                   "would measure charting cadence rather than indication. The ",
                   "primary window is %gh."), WIN_MIN, IND_WIN_H),
    paste0("The denominator is ALL events. indication.csv also reports the ",
           "fraction among events that had any assessment in the window, which ",
           "separates documentation from indication, and a backward-only ",
           "variant restricted to scores at or before the event."),
    if (length(absent))
      sprintf("Declared but not charted at this site: %s.", paste(absent, collapse = ", "))
    else ""))

p_ind <- house(
  ggplot(ind_sens[ind_sens$class != "either", ],
         aes(window_hours, pct_indicated_all, colour = class,
             linetype = event_type)) +
    geom_vline(xintercept = IND_WIN_H, colour = MUTED, linetype = "22",
               linewidth = 0.4) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
    scale_colour_manual(
      values = c(pain = STATE_COLOURS[["bolus only"]],
                 sedation = STATE_COLOURS[["continuous only"]]),
      labels = c(pain = "Pain score", sedation = "Agitation (RASS)"), name = NULL) +
    scale_linetype_manual(values = c(increase = "solid", bolus = "42"),
                          labels = c(increase = "Rate increase", bolus = "Bolus"),
                          name = NULL) +
    scale_x_continuous(breaks = IND_SWEEP_H) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    # ORDER IS REQUIRED, not cosmetic. This is the only figure in the live
    # pipeline with two legends, and with both left at the default order = 0
    # ggplot2 placed them in an order that VARIED BETWEEN R SESSIONS: the same
    # data and the same code produced two different PNGs, the colour and
    # linetype boxes swapping sides. Measured 2026-10-05 -- 5 sessions gave one
    # byte-identical file with these two lines, two distinct files without them.
    # A figure a site cannot reproduce by checksum is a federation problem.
    guides(colour = guide_legend(order = 1),
           linetype = guide_legend(order = 2)) +
    labs(x = "Window, +/- hours",
         y = "% of EVENTS with a qualifying score"))
ggsave(file.path(dirs$phase, "indication.png"), p_ind,
       width = 7.5, height = 4.4, dpi = 200)


# ---- 10. Write ---------------------------------------------------------------

write_out <- function(x, name) {
  write.csv(x, file.path(dirs$phase, name), row.names = FALSE)
  cat(sprintf("written: %s\n", name))
}

report_small_cells <- function(x, name) {
  if (!"n_events" %in% names(x)) return(invisible(NULL))
  small <- x[x$n_events > 0 & x$n_events < MIN_CELL, ]
  if (nrow(small)) {
    cat(sprintf("  NOTE %s: %d cell(s) below reporting.small_cell_min_den = %d\n",
                name, nrow(small), MIN_CELL))
  }
}
cat("\n")
for (nm in c("coadministration", "window_sensitivity", "charting_precision")) {
  report_small_cells(get(switch(nm, coadministration = "coad",
                                window_sensitivity = "sens",
                                charting_precision = "prec")), paste0(nm, ".csv"))
}

write_out(coad, "coadministration.csv")
write_out(sens, "window_sensitivity.csv")
write_out(prec, "charting_precision.csv")
write_out(mag,  "rate_change_magnitude.csv")
write_out(agree, "charting_agreement.csv")
write_out(adh_dist, "adherence_distribution.csv")
write_out(pool_histograms, "pooling_histograms.csv")
write_out(ind, "indication.csv")
write_out(ind_sens, "indication_window_sensitivity.csv")

write_parquet(adh, file.path(dirs$out_phi, "titration_adherence.parquet"))
cat("written: titration_adherence.parquet (PHI -- per encounter)\n")

# The classified events themselves, for code/06_unit_variation.R. Handed off
# rather than re-derived: the event definition (section 6) and the pairing window
# are this script's, and a second implementation of "what counts as an
# uptitration" is a drift hazard with no upside. 06 splits these events by unit
# and by year; it must not redefine them.
#
# A relative clock, and nothing else -- the same contract Phase 0 holds for the
# raw events this was built from.
ev_out <- ev[, c("encounter_block", "t_hr", "kind", "paired", "on_hour")]
stopifnot(
  "classified events carry a datetime" =
    !any(vapply(ev_out, inherits, logical(1), what = c("POSIXct", "POSIXt", "Date"))),
  "classified events fall outside the grid" =
    min(ev_out$t_hr) >= 0 && max(ev_out$t_hr) <= EXTENT_H,
  "an event is unclassified" = all(ev_out$kind %in% EVENT_LEVELS)
)
write_parquet(ev_out, file.path(dirs$out_phi, "titration_events_classified.parquet"))
cat(sprintf("written: titration_events_classified.parquet (PHI -- %d events, %d episodes)\n",
            nrow(ev_out), length(unique(ev_out$encounter_block))))

# The per-event indication flags, for code/06_unit_variation.R. Handed off
# rather than re-derived, for the same reason the classified events are: the
# thresholds and the window are this script's, and a second implementation of
# "was this indicated" would drift from the first with nothing to catch it.
#
# The PRIMARY window only. The 1-4h sweep is this script's question; 06 asks a
# per-unit question at one window, and shipping all four would invite a reader
# to pick the flattering one.
#
# Boluses travel in the SAME file with event_type, and deliberately NOT inside
# titration_events_classified.parquet: 06 selects increases from that file with
# `kind != "downtitration"`, so a bolus row added there would be silently
# counted as a rate increase.
ind_ev <- rbind(
  data.frame(encounter_block = inc$encounter_block, t_hr = inc$t_hr,
             event_type = "increase", kind = as.character(inc$kind),
             stringsAsFactors = FALSE),
  data.frame(encounter_block = bol$encounter_block, t_hr = bol$t_hr,
             event_type = "bolus", kind = NA_character_,
             stringsAsFactors = FALSE))
ind_ev$pain_indicated     <- near_any(ind_ev, PAIN_HIT, IND_WIN_H)
ind_ev$sedation_indicated <- near_any(ind_ev, SED_HIT,  IND_WIN_H)
ind_ev$any_indicated      <- near_any(ind_ev, BOTH_HIT, IND_WIN_H)
ind_ev$any_assessed       <- near_any(ind_ev, BOTH_ANY, IND_WIN_H)

stopifnot(
  "indication events carry a datetime" =
    !any(vapply(ind_ev, inherits, logical(1), what = c("POSIXct", "POSIXt", "Date"))),
  "indication events fall outside the grid" =
    min(ind_ev$t_hr) >= 0 && max(ind_ev$t_hr) <= EXTENT_H,
  "an event is neither an increase nor a bolus" =
    all(ind_ev$event_type %in% c("increase", "bolus")),
  # any_indicated is the union of the two classes, so it can never be smaller.
  "any_indicated is not the union of pain and sedation" =
    all(ind_ev$any_indicated >= (ind_ev$pain_indicated | ind_ev$sedation_indicated)),
  "an indicated event was never assessed" =
    all(!ind_ev$any_indicated | ind_ev$any_assessed))

# The handoff must reproduce the table this script just printed, or 06 and 05
# would report different numbers for the same quantity.
for (et in c("increase", "bolus")) {
  here <- 100 * mean(ind_ev$any_indicated[ind_ev$event_type == et])
  there <- ind$pct_indicated_all[ind$event_type == et & ind$class == "either"]
  stopifnot("the handoff disagrees with indication.csv" = abs(here - there) < 0.05)
}

write_parquet(ind_ev, file.path(dirs$out_phi, "indication_events.parquet"))
cat(sprintf("written: indication_events.parquet (PHI -- %s events, %s increases, %s boluses)\n",
            format(nrow(ind_ev), big.mark = ","),
            format(sum(ind_ev$event_type == "increase"), big.mark = ","),
            format(sum(ind_ev$event_type == "bolus"), big.mark = ",")))

write_json(prov, file.path(dirs$phase, "provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)
cat("written: provenance.json\n")

write_captions(file.path(dirs$phase, "captions.md"), "05_titration.R",
               grep("\\.png$", OWNED$phase, value = TRUE), prov)
cat("written: captions.md\n")


# ---- 11. Provenance ----------------------------------------------------------

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/05_titration.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "05_titration_sessioninfo.txt")
)
