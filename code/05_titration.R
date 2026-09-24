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


# ---- 2. Config ---------------------------------------------------------------

config <- fromJSON(here("config", "config.json"), simplifyVector = FALSE)
set.seed(config$model$seed)

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
stopifnot("the primary window must be one of the sensitivity widths" = WIN_MIN %in% WIN_SENS,
          "the threshold must be positive" = MIN_DELTA > 0)

EVENT_LEVELS <- c("initiation", "uptitration", "downtitration")


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
  out_phi = c("titration_adherence.parquet"),
  phase   = c("coadministration.csv", "window_sensitivity.csv",
              "charting_precision.csv", "rate_change_magnitude.csv",
              "charting_agreement.csv", "adherence_distribution.csv",
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

summarise_pairing <- function(d, by) {
  do.call(rbind, lapply(split(d, d[[by]], drop = TRUE), function(x) {
    if (!nrow(x)) return(NULL)
    data.frame(group = as.character(x[[by]][1]), n_events = nrow(x),
               n_paired = sum(x$paired),
               pct_paired = round(100 * mean(x$paired), 1))
  }))
}

coad <- summarise_pairing(ev, "kind")
coad <- rbind(
  coad,
  data.frame(group = "any increase", n_events = sum(ev$kind != "downtitration"),
             n_paired = sum(ev$paired & ev$kind != "downtitration"),
             pct_paired = round(100 * mean(ev$paired[ev$kind != "downtitration"]), 1)))
coad$window_minutes <- WIN_MIN

cat(sprintf("\nBolus within +/-%g min of the event\n", WIN_MIN))
print(coad[, c("group", "n_events", "n_paired", "pct_paired")], row.names = FALSE)

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
mag <- data.frame(
  quantile = c("p10", "p25", "p50", "p75", "p90", "max"),
  mcg_per_hr = round(unname(quantile(d_all, c(.1, .25, .5, .75, .9, 1))), 1))
mag <- rbind(mag, data.frame(
  quantile = paste0("share >= ", c(5, 10, 25, 50, 100)),
  mcg_per_hr = round(100 * vapply(c(5, 10, 25, 50, 100),
                                  function(k) mean(d_all >= k), numeric(1)), 1)))
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
             pct_never = round(100 * mean(e$adherence == 0), 1),
             pct_always = round(100 * mean(e$adherence == 1), 1))
}))
cat(sprintf("\nPer-encounter adherence (episodes with >= %g increases)\n", MIN_EV))
print(adh_dist, row.names = FALSE)


# ---- 9. Figures --------------------------------------------------------------

register_caption("coadministration.png",
  "Bolus co-administration at a fentanyl infusion rate change",
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
          if (VENT_ONLY) "Restricted to ventilated windows." else ""))

p_coad <- house(
  ggplot(coad[coad$group %in% EVENT_LEVELS, ],
         aes(factor(group, levels = EVENT_LEVELS), pct_paired,
             fill = factor(group, levels = EVENT_LEVELS))) +
    geom_col(width = 0.62) +
    geom_text(aes(label = sprintf("%.1f%%\n(%s)", pct_paired,
                                  format(n_events, big.mark = ","))),
              vjust = -0.25, size = 3, colour = INK) +
    scale_fill_manual(values = c(initiation = STATE_COLOURS[["continuous only"]],
                                 uptitration = STATE_COLOURS[["bolus only"]],
                                 downtitration = STATE_COLOURS[["discharged alive"]]),
                      labels = c(initiation = "Initiation (0 -> on)",
                                 uptitration = "Uptitration",
                                 downtitration = "Downtitration (negative control)"),
                      name = NULL) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = sprintf("%% with a bolus within +/-%g min", WIN_MIN)))
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
                        labels = c(initiation = "Initiation (0 -> on)",
                                   uptitration = "Uptitration",
                                   downtitration = "Downtitration (negative control)"),
                        name = NULL) +
    scale_x_continuous(breaks = WIN_SENS) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Pairing window, +/- minutes",
         y = "% of rate changes with a bolus"))
ggsave(file.path(dirs$phase, "window_sensitivity.png"), p_sens,
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

write_parquet(adh, file.path(dirs$out_phi, "titration_adherence.parquet"))
cat("written: titration_adherence.parquet (PHI -- per encounter)\n")

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
