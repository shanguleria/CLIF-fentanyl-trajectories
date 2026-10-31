# ==============================================================================
# 05_titration.R  --  bolus co-administration at an infusion rate change
#
# Purpose : When the fentanyl infusion rate is increased, how often is a bolus
#           given within +/-30 minutes?
# Author  : Shan Guleria
# Created : 2026-09-24
# Inputs  : output/intermediate_phi/titration_{rate,bolus}_events.parquet, trajectory_long.parquet
# Outputs : output/final_no_phi/05_titration/            reported directly
#           output/final_no_phi/05_titration/diagnostics/ method checks
#           OWNED decides which folder each file is written to
#
# Descriptive only; associations are deferred -- see protocol.json.
# Raw charted timestamps, never the hourly grid: hour-binning would move a rate
# change up to an hour from its bolus.
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
# pool_hist() and boot_group(), shared with 02 and 06 so neither can drift.
source(here("code", "utils", "pooling.R"))


# ---- 2. Config ---------------------------------------------------------------

config <- fromJSON(here("config", "config.json"), simplifyVector = FALSE)
COV  <- fromJSON(here("config", "protocol.json"), simplifyVector = FALSE)
SEED <- COV$model$seed
set.seed(SEED)

WINDOW_H <- COV$cohort$window_hours
EXTENT_H <- COV$cohort$granular_extent_hours
MIN_CELL <- COV$reporting$small_cell_min_den

# Federation-critical: read from the config, never defaulted here.
SPEC <- COV$titration
stopifnot("protocol.json must declare a `titration` block" = !is.null(SPEC))
MIN_DELTA <- as.numeric(SPEC$min_rate_change_mcg_hr)
WIN_MIN   <- as.numeric(SPEC$bolus_window_minutes)
WIN_SENS  <- sort(unique(c(as.numeric(unlist(SPEC$window_sensitivity_minutes)), WIN_MIN)))
VENT_ONLY <- isTRUE(SPEC$ventilated_only)
ADH_LM_H  <- as.numeric(SPEC$adherence_landmark_hours)
# Every window a per-encounter adherence scope is emitted for. Declared rather
# than derived from ADH_LM_H alone, because 07 refits the primary model at a
# landmark of its own and needs the exposure measured over THAT window -- see
# titration._adherence_scope_hours_why. Asserted to include the reported one, so
# a config that drops it fails here rather than silently shipping fewer scopes.
ADH_SCOPES <- sort(unique(as.numeric(unlist(SPEC$adherence_scope_hours))))
stopifnot(
  "titration.adherence_scope_hours must include adherence_landmark_hours" =
    ADH_LM_H %in% ADH_SCOPES)
MIN_EV    <- as.numeric(SPEC$min_events_for_adherence)
BOOT_REPS <- as.integer(SPEC$bootstrap_resamples)
# Absolute config-declared bins; site-local quantiles cannot be pooled.
HSPEC <- COV$pooling$histograms$variables
stopifnot("protocol.json must declare a `titration` block" = !is.null(SPEC),
          "protocol.json must declare pooling.histograms.variables" = !is.null(HSPEC),
          "the primary window must be one of the sensitivity widths" = WIN_MIN %in% WIN_SENS,
          "the threshold must be positive" = MIN_DELTA > 0,
          "bootstrap_resamples must be a positive whole number" = BOOT_REPS > 0)
HIST_MINE <- names(HSPEC)[vapply(HSPEC, function(s)
  identical(as.character(s$owner), "05"), logical(1))]

EVENT_LEVELS <- c("initiation", "uptitration", "downtitration")

IND <- COV$indication
stopifnot("protocol.json must declare an `indication` block" = !is.null(IND))
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
dirs$phi <- phi_dir(dirs, "05_titration")   # patient-level, same subdivision
# Method checks and sensitivity analyses; everything here still ships.
dirs$diagnostics <- phase_dir(dirs, file.path("05_titration", "diagnostics"))
prov <- provenance(config)

message(sprintf("[05_titration] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))


# ---- 4. Guards ---------------------------------------------------------------

manifest <- require_manifest(dirs, here())
message(sprintf("  reading cohort-build outputs from code %s, generated %s",
                manifest$code_version, manifest$generated))

OWNED <- list(
  phi = c("titration_adherence.parquet",
              "titration_events_classified.parquet",
              "indication_events.parquet"),
  # Reported directly.
  phase   = c("coadministration.csv", "indication.csv",
              "adherence_distribution.csv", "assessment_frequency.csv",
              "rate_change_magnitude.csv", "pooling_histograms.csv",
              "coadministration.png", "indication.png",
              "provenance.json", "captions.md"),
  # Method checks, cited in Methods; none is a result.
  diagnostics = c("window_sensitivity.csv", "indication_window_sensitivity.csv",
                  "charting_precision.csv", "charting_agreement.csv",
                  "window_sensitivity.png"))
# Stale twins left in the phase directory from before the diagnostics/ split.
# Flat intermediate_phi locations, before it was subdivided by writing script
# on 2026-10-30. Declared so a stale twin is deleted, not read as current.
RETIRED_PHI <- file.path("output", "intermediate_phi",
                         c("titration_adherence.parquet",
                           "titration_events_classified.parquet",
                           "indication_events.parquet"))

RETIRED <- file.path("output", "final_no_phi", "05_titration",
                     c("window_sensitivity.csv", "indication_window_sensitivity.csv",
                       "charting_precision.csv", "charting_agreement.csv",
                       "window_sensitivity.png"))
n_cleared <- clear_owned_outputs(dirs, OWNED, retired = c(RETIRED, RETIRED_PHI))
if (n_cleared) message(sprintf("  cleared %d output(s) from a previous run", n_cleared))


# ---- 5. Read -----------------------------------------------------------------

rate <- as.data.frame(read_parquet(file.path(phi_dir(dirs, "01_cohort"), "titration_rate_events.parquet")))
bol  <- as.data.frame(read_parquet(file.path(phi_dir(dirs, "01_cohort"), "titration_bolus_events.parquet")))

# Whole-hour timestamps throughout would mean the hourly grid leaked in.
frac_sub <- mean(rate$t_hr %% 1 != 0)
stopifnot(
  "rate events look hour-binned, not charted -- the grid has leaked in" = frac_sub > 0.5,
  "the clock must be relative hours" = is.numeric(rate$t_hr) && min(rate$t_hr) >= 0)

asm <- as.data.frame(read_parquet(file.path(phi_dir(dirs, "01_cohort"), "assessment_events.parquet")))
# Scales are NOT comparable (NVPS 0-10, CPOT 0-8, RASS -5..+4); each is cut at its OWN threshold.
declared <- c(names(PAIN_THRESH), names(SED_THRESH))
absent <- setdiff(declared, unique(asm$instrument))
if (length(absent))
  cat(sprintf("  NOTE declared but not charted at this site: %s\n",
              paste(absent, collapse = ", ")))
stopifnot("no declared indication instrument is present in the export" =
            length(setdiff(declared, absent)) >= 1)

long <- as.data.frame(read_parquet(
  file.path(phi_dir(dirs, "01_cohort"), "trajectory_long.parquet"),
  col_select = c("encounter_block", "window_idx", "imv_status")))
long$ventilated <- !is.na(long$imv_status) & long$imv_status == 1

cat(sprintf("\n%s rate record(s), %s bolus record(s), %s episodes\n",
            format(nrow(rate), big.mark = ","), format(nrow(bol), big.mark = ","),
            format(length(unique(rate$encounter_block)), big.mark = ",")))


# ---- 5b. Cohort strata -------------------------------------------------------
# `cohort` is a KEY, not a label: landmark is a SUBSET of all and must never be summed into it.
tte <- as.data.frame(read_parquet(
  file.path(phi_dir(dirs, "01_cohort"), "time_to_event.parquet"),
  col_select = c("encounter_block", "landmark_eligible")))
ELIG <- tte$encounter_block[tte$landmark_eligible]
# NULL means no restriction, keeping the unrestricted pass byte-identical.
COHORTS <- list(all = NULL, landmark = ELIG)

keep_rows <- function(d, keep)
  if (is.null(keep)) d else d[d$encounter_block %in% keep, , drop = FALSE]

by_cohort <- function(build) {
  do.call(rbind, lapply(names(COHORTS), function(ck) {
    out <- build(COHORTS[[ck]])
    if (is.null(out) || !nrow(out)) return(NULL)
    rownames(out) <- NULL
    cbind(cohort = ck, out, stringsAsFactors = FALSE)
  }))
}
cat(sprintf("  landmark cohort: %s of %s episodes with a titration record\n",
            format(length(intersect(ELIG, unique(rate$encounter_block))), big.mark = ","),
            format(length(unique(rate$encounter_block)), big.mark = ",")))


# ---- 6. Classify the events --------------------------------------------------
# The rate delta defines the event, not the action label; a stop is not an event.

rate <- rate[order(rate$encounter_block, rate$t_hr), ]
prev_rate  <- ave(rate$rate, rate$encounter_block, FUN = function(v) c(NA, v[-length(v)]))
rate$delta <- rate$rate - prev_rate
rate$first_of_block <- is.na(prev_rate)

is_stop <- rate$not_given | rate$action == "stop"
rate$kind <- NA_character_
# An initiation rises from a KNOWN zero; a block opening mid-infusion is not one.
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

# Ventilated windows only; the unrestricted count prints alongside.
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

# Unit of analysis is the EVENT, not the encounter; n_episodes travels beside it.
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

# Interval and design effect for ONE clustered proportion; a binomial interval is far
# too narrow. The resampling universe is passed, never inferred (lessons.md #26).
boot_one <- function(d, outcome_col) {
  eps <- sort(unique(d$encounter_block))
  x <- data.frame(encounter_block = d$encounter_block, g = "all",
                  y = as.integer(d[[outcome_col]]))
  r <- boot_group(x, "g", SEED, outcome = "y", episodes = eps,
                  n_resamples = BOOT_REPS)
  data.frame(ci_lo = round(r$ci_lo[1], 1), ci_hi = round(r$ci_hi[1], 1),
             deff = round(r$deff[1], 3))
}

inc <- ev[ev$kind != "downtitration", ]

build_coad <- function(keep) {
  e <- keep_rows(ev, keep)
  i <- e[e$kind != "downtitration", ]
  if (!nrow(i)) return(NULL)
  out <- rbind(
    summarise_pairing(e, "kind"),
    data.frame(group = "any increase", n_events = nrow(i),
               n_episodes = length(unique(i$encounter_block)),
               events_per_episode = round(nrow(i) / length(unique(i$encounter_block)), 2),
               n_paired = sum(i$paired),
               pct_paired = round(100 * mean(i$paired), 1)))
  out <- cbind(out, do.call(rbind, lapply(out$group, function(g) {
    boot_one(if (identical(g, "any increase")) i else e[e$kind == g, ], "paired")
  })))
  out$window_minutes <- WIN_MIN
  out
}
coad <- by_cohort(build_coad)

cat(sprintf("\nBolus within +/-%g min of the event\n", WIN_MIN))
print(coad[, c("cohort", "group", "n_events", "n_episodes", "events_per_episode",
               "n_paired", "pct_paired")], row.names = FALSE)
cat(sprintf("  percentages are EVENT-level; %s episodes contribute the %s increases\n",
            format(length(unique(inc$encounter_block)), big.mark = ","),
            format(nrow(inc), big.mark = ",")))

# Pairing vs window width, increases AND the downtitration control: a gap already at
# narrow widths is deliberate pairing, not bolus frequency.
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

# Charting-precision stratum: for an event charted on :00 a true pairing can fall outside the window.
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

# Agreement between the threshold-defined event and the charted action label.
agree <- as.data.frame(table(action = ev$action, kind = ev$kind))
agree <- agree[agree$Freq > 0, ]
names(agree)[3] <- "n_events"
cat("\nCharted action on events at or above the threshold\n")
print(reshape(agree, idvar = "action", timevar = "kind", direction = "wide"),
      row.names = FALSE)

# The magnitude distribution, so the threshold is checkable against the site's data.
d_all <- abs(rate$delta[!is.na(rate$delta) & rate$delta != 0 &
                          !is_stop & !rate$first_of_block])
# `value` carries two units, so the unit is a column. `n` is the shared denominator
# over all non-zero rate changes; `n_at_or_above` is what the share rows pool on.
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
# Patient-level, so only the distribution ships. Emitted twice: overall, and within
# the first ADH_LM_H hours, measured before outcome accrual (immortal-time bias).

adherence <- function(d, tag) {
  inc <- d[d$kind != "downtitration", ]
  if (!nrow(inc)) return(NULL)
  agg <- aggregate(cbind(n_events = rep(1, nrow(inc)), n_paired = as.integer(inc$paired))
                   ~ encounter_block, data = inc, FUN = sum)
  agg$adherence <- agg$n_paired / agg$n_events
  agg$scope <- tag
  agg
}
# adh stays UNRESTRICTED: the pooling histogram below is a whole-cohort distribution.
scoped <- function(d) do.call(rbind, c(
  list(adherence(d, "overall")),
  lapply(ADH_SCOPES, function(h)
    adherence(d[d$t_hr <= h, ], sprintf("first_%gh", h)))))
adh <- scoped(ev)

summarise_adherence <- function(a) do.call(rbind, lapply(split(a, a$scope), function(x) {
  e <- x[x$n_events >= MIN_EV, ]
  data.frame(scope = x$scope[1],
             n_episodes = nrow(x),
             n_episodes_ge_min_events = nrow(e),
             median_adherence = round(median(e$adherence), 3),
             q1 = round(quantile(e$adherence, .25), 3),
             q3 = round(quantile(e$adherence, .75), 3),
             # Denominator is n_episodes_ge_min_events, NOT n_episodes.
             n_never = sum(e$adherence == 0),
             n_always = sum(e$adherence == 1),
             pct_never = round(100 * mean(e$adherence == 0), 1),
             pct_always = round(100 * mean(e$adherence == 1), 1))
}))

adh_dist <- by_cohort(function(keep) {
  e <- keep_rows(ev, keep)
  a <- scoped(e)
  if (is.null(a)) return(NULL)
  summarise_adherence(a)
})
cat(sprintf("\nPer-encounter adherence (episodes with >= %g increases)\n", MIN_EV))
print(adh_dist, row.names = FALSE)


# ---- 8b. Pooling histograms --------------------------------------------------
# A median has no closed form across sites; counts on FIXED absolute bins do.
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


# ---- 8a. How often is a score documented at all? -----------------------------
# Separates "no pain" from "no score charted". Denominator is VENTILATED windows, as in indication.csv.

asm_w <- asm
asm_w$window_idx <- floor(asm_w$t_hr / WINDOW_H)
vent <- long[long$ventilated, c("encounter_block", "window_idx")]

# The denominator restricts with the numerator: landmark windows for the landmark pass.
build_asm_freq <- function(keep) {
  a <- keep_rows(asm_w, keep)
  v <- keep_rows(vent, keep)
  vent_keys <- paste(v$encounter_block, v$window_idx, sep = "\r")
  n_vent_windows <- length(vent_keys)
  vent_days <- n_vent_windows * WINDOW_H / 24
  if (!nrow(a) || !n_vent_windows) return(NULL)
  out <- do.call(rbind, lapply(sort(unique(a$instrument)), function(ins) {
  d <- a[a$instrument == ins, ]
  # Gaps are WITHIN an episode; one spanning two episodes is a discharge.
  gaps <- unlist(lapply(split(d$t_hr, d$encounter_block), function(t) {
    t <- sort(t); if (length(t) < 2L) numeric(0) else diff(t)
  }), use.names = FALSE)
  in_vent <- paste(d$encounter_block, d$window_idx, sep = "\r") %in% vent_keys
  q <- function(p) if (!length(gaps)) NA_real_ else
    round(as.numeric(stats::quantile(gaps, p, names = FALSE)), 2)
  data.frame(
    instrument = ins,
    class = if (ins %in% names(PAIN_THRESH)) "pain" else "sedation",
    n_records = nrow(d),
    n_episodes_with_any = length(unique(d$encounter_block)),
    # The poolable pair: these two sum across sites, the quantiles below do not.
    n_ventilated_windows = n_vent_windows,
    n_windows_with_any = length(unique(
      paste(d$encounter_block, d$window_idx, sep = "\r")[in_vent])),
    records_per_ventilated_day = round(nrow(d) / vent_days, 2),
    n_gaps = length(gaps),
    gap_p50 = q(0.50), gap_p75 = q(0.75), gap_p90 = q(0.90), gap_p95 = q(0.95),
    stringsAsFactors = FALSE)
  }))
  out$pct_windows_with_any <-
    round(100 * out$n_windows_with_any / out$n_ventilated_windows, 1)
  out
}
asm_freq <- by_cohort(build_asm_freq)
stopifnot(
  "a documented window count exceeds the ventilated windows it is counted over" =
    all(asm_freq$n_windows_with_any <= asm_freq$n_ventilated_windows))

cat("\nAssessment documentation frequency (ventilated windows)\n")
print(asm_freq[, c("cohort", "instrument", "class", "n_records",
                   "records_per_ventilated_day", "pct_windows_with_any",
                   "gap_p50", "gap_p90")], row.names = FALSE)


# ---- 8b. Was the dose change indicated? --------------------------------------
# Was the increase PROMPTED by anything documented? A LEVEL, not a change, and in
# HOURS not minutes -- assessments are charted hourly to q4h.

# Times above threshold, and times recorded at all; split per episode.
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

# TRUE when the episode has a listed time in the window; `back_only` keeps times at or before it.
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
# "either" re-runs over the UNION of the threshold maps; merging the split lists by name would drop episodes.
BOTH_THRESH <- c(PAIN_THRESH, SED_THRESH)
BOTH_HIT <- qualifying_times(asm, BOTH_THRESH); BOTH_ANY <- recorded_times(asm, BOTH_THRESH)

ev_sets <- list(increase = inc,
                bolus    = data.frame(encounter_block = bol$encounter_block,
                                      t_hr = bol$t_hr))

# Two denominators: all events conflates "no score taken" with "score taken, not justifying"; assessed events separate them.
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

# The sweep stays whole-cohort: a method check about the measurement, not the population.
IND_HIT <- list(pain = PAIN_HIT, sedation = SED_HIT, either = BOTH_HIT)

# Bootstrap the PRIMARY window only; the interval is on pct_indicated_all, hence the prefix.
ind <- by_cohort(function(keep) {
  es <- lapply(ev_sets, keep_rows, keep = keep)
  rows <- do.call(rbind, lapply(names(es), function(etype) {
    d <- es[[etype]]
    if (!nrow(d)) return(NULL)
    rbind(indication_row(d, IND_WIN_H, etype, "pain", PAIN_HIT, PAIN_ANY),
          indication_row(d, IND_WIN_H, etype, "sedation", SED_HIT, SED_ANY),
          indication_row(d, IND_WIN_H, etype, "either", BOTH_HIT, BOTH_ANY))
  }))
  if (is.null(rows)) return(NULL)
  rows$class <- factor(rows$class, levels = c("pain", "sedation", "either"))
  cbind(rows, do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
    d <- es[[as.character(rows$event_type[i])]]
    b <- boot_one(data.frame(encounter_block = d$encounter_block,
                             y = near_any(d, IND_HIT[[as.character(rows$class[i])]],
                                          IND_WIN_H)), "y")
    names(b) <- paste0("indicated_", names(b))
    b
  })))
})
# Two code paths, one quantity: these rows must match the sweep at the same window.
local({
  a <- ind[ind$cohort == "all", ]
  b <- ind_sens[ind_sens$window_hours == IND_WIN_H, ]
  k <- c("event_type", "class", "n_events", "n_assessed", "n_indicated")
  a <- a[order(a$event_type, a$class), k]; b <- b[order(b$event_type, b$class), k]
  stopifnot("the per-cohort indication rows disagree with the window sweep" =
              isTRUE(all.equal(a, b, check.attributes = FALSE)))
})

# A wider window can only ADD qualifying scores; a fall means the join is wrong.
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
print(ind[, c("cohort", "event_type", "class", "n_events", "pct_assessed",
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
  # Whole cohort only: unfiltered this draws two indistinguishable bars per group.
  ggplot(coad[coad$cohort == "all" & coad$group %in% EVENT_LEVELS, ],
         aes(factor(group, levels = EVENT_LEVELS), pct_paired,
             fill = factor(group, levels = EVENT_LEVELS))) +
    geom_col(width = 0.62) +
    geom_text(aes(label = sprintf("%.1f%%\n%s events", pct_paired,
                                  format(n_events, big.mark = ","))),
              vjust = -0.25, size = 3, colour = INK) +
    # No legend: the x axis already names each bar, so colour carries no identity.
    scale_fill_manual(values = c(initiation = STATE_COLOURS[["continuous only"]],
                                 uptitration = STATE_COLOURS[["bolus only"]],
                                 downtitration = STATE_COLOURS[["discharged alive"]]),
                      guide = "none") +
    scale_x_discrete(labels = c(initiation = "Initiation", uptitration = "Uptitration",
                                downtitration = "Downtitration")) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.18))) +
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
ggsave(file.path(dirs$diagnostics, "window_sensitivity.png"), p_sens,
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
    # Order is REQUIRED: two legends both at the default order = 0 came out in a
    # different order between R sessions, so the PNG was not reproducible.
    guides(colour = guide_legend(order = 1),
           linetype = guide_legend(order = 2)) +
    labs(x = "Window, +/- hours",
         y = "% of EVENTS with a qualifying score"))
ggsave(file.path(dirs$phase, "indication.png"), p_ind,
       width = 7.5, height = 4.4, dpi = 200)


# ---- 10. Write ---------------------------------------------------------------

# Destination comes from OWNED, so a file cannot be written outside its folder.
write_out <- function(x, name) {
  key <- if (name %in% OWNED$diagnostics) "diagnostics" else "phase"
  stopifnot("write_out called for a file OWNED does not declare" =
              name %in% c(OWNED$phase, OWNED$diagnostics))
  write.csv(x, file.path(dirs[[key]], name), row.names = FALSE)
  cat(sprintf("written: %s%s\n", if (key == "diagnostics") "diagnostics/" else "",
              name))
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
write_out(asm_freq, "assessment_frequency.csv")
write_out(ind, "indication.csv")
write_out(ind_sens, "indication_window_sensitivity.csv")

write_parquet(adh, file.path(phi_dir(dirs, "05_titration"), "titration_adherence.parquet"))
cat("written: titration_adherence.parquet (PHI -- per encounter)\n")

# Classified events for 06, which splits them and must not redefine them. Relative clock only.
ev_out <- ev[, c("encounter_block", "t_hr", "kind", "paired", "on_hour")]
stopifnot(
  "classified events carry a datetime" =
    !any(vapply(ev_out, inherits, logical(1), what = c("POSIXct", "POSIXt", "Date"))),
  "classified events fall outside the grid" =
    min(ev_out$t_hr) >= 0 && max(ev_out$t_hr) <= EXTENT_H,
  "an event is unclassified" = all(ev_out$kind %in% EVENT_LEVELS)
)
write_parquet(ev_out, file.path(phi_dir(dirs, "05_titration"), "titration_events_classified.parquet"))
cat(sprintf("written: titration_events_classified.parquet (PHI -- %d events, %d episodes)\n",
            nrow(ev_out), length(unique(ev_out$encounter_block))))

# Per-event indication flags for 06, at the PRIMARY window. Boluses carry event_type
# and stay OUT of the classified-events file, which 06 filters by `kind`.
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
  "any_indicated is not the union of pain and sedation" =
    all(ind_ev$any_indicated >= (ind_ev$pain_indicated | ind_ev$sedation_indicated)),
  "an indicated event was never assessed" =
    all(!ind_ev$any_indicated | ind_ev$any_assessed))

# The handoff must reproduce the table printed above.
for (et in c("increase", "bolus")) {
  here <- 100 * mean(ind_ev$any_indicated[ind_ev$event_type == et])
  # The handoff is the WHOLE cohort, so compare it against the `all` rows.
  there <- ind$pct_indicated_all[ind$cohort == "all" & ind$event_type == et
                                 & ind$class == "either"]
  stopifnot("the handoff disagrees with indication.csv" = abs(here - there) < 0.05)
}

write_parquet(ind_ev, file.path(phi_dir(dirs, "05_titration"), "indication_events.parquet"))
cat(sprintf("written: indication_events.parquet (PHI -- %s events, %s increases, %s boluses)\n",
            format(nrow(ind_ev), big.mark = ","),
            format(sum(ind_ev$event_type == "increase"), big.mark = ","),
            format(sum(ind_ev$event_type == "bolus"), big.mark = ",")))

write_json(prov, file.path(dirs$phase, "provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)
cat("written: provenance.json\n")

write_captions(file.path(dirs$phase, "captions.md"), "05_titration.R",
               grep("\\.png$", c(OWNED$phase, OWNED$diagnostics), value = TRUE),
               prov)
cat("written: captions.md\n")


# ---- 11. Provenance ----------------------------------------------------------

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/05_titration.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "05_titration_sessioninfo.txt")
)
