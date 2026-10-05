# ==============================================================================
# 06_unit_variation.R  --  bolus co-administration by ICU unit and by calendar year
#
# Purpose : Split 05_titration.R's single pooled number by WHERE and WHEN. Which ICU was the patient in when the rate change was charted, and what calendar year was it? Caterpillar and funnel views per unit, and the metric across the years the extract spans.
# Author  : Shan Guleria
# Created : 2026-09-30
# Inputs  : output/intermediate_phi/titration_events_classified.parquet, icu_intervals.parquet, trajectory_long.parquet
# Outputs : output/final_no_phi/06_unit_variation/ : the files listed in OWNED
#
# NOT A QUALITY MEASURE. No evidence establishes that pairing a bolus with an
# up-titration is better care. This describes VARIATION IN PRACTICE. A unit with
# a lower rate is not thereby worse, and tests/test_covariates.py bans the words
# quality, performance and benchmark from this file for that reason.
#
# THE EVENTS ARE NOT REDEFINED HERE. 05_titration.R owns what counts as an
# initiation, an uptitration and a pairing; this script reads its classified
# output and only splits it. A second implementation of "what is an uptitration"
# would be a drift hazard with no upside.
#
# THE INTERVALS DESCRIBE, THEY DO NOT ADJUST. The bootstrap resamples EPISODES
# because events cluster within them. It says how much a unit's observed rate
# would wobble if the care process re-ran. It does NOT adjust for case mix -- a
# unit's rate reflects who it admits as much as how it practises -- and the
# association between unit and adherence is deliberately not modelled (SG,
# 2026-09-30). See covariates.json unit_variation._STATUS_model for what the
# deferred method is and how its estimand differs from this one.
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
# boot_group() lives in pooling.R so 05 and 06 cannot drift in how they
# bootstrap the same clustered proportion.
source(here("code", "utils", "pooling.R"))


# ---- 2. Config ---------------------------------------------------------------

config <- fromJSON(here("config", "config.json"), simplifyVector = FALSE)
SEED <- config$model$seed
set.seed(SEED)

WINDOW_H <- config$cohort$window_hours
EXTENT_H <- config$cohort$granular_extent_hours
MIN_CELL <- config$reporting$small_cell_min_den

COV <- fromJSON(here("config", "covariates.json"), simplifyVector = FALSE)
# Read from the config, never defaulted here: the attribution rule and the
# resample count are federation-critical, and a local default is how a site
# keeps computing on the old numbers after the consortium moves them.
SPEC <- COV$unit_variation
stopifnot("covariates.json must declare a `unit_variation` block" = !is.null(SPEC))
ATTR_CATS <- unlist(SPEC$attributed_location_categories)
UNIT_KEYS <- unlist(SPEC$unit_keys)
MIN_UNIT  <- as.numeric(SPEC$min_events_per_unit)
B_RESAMP  <- as.integer(SPEC$bootstrap_resamples)
YEAR_FROM <- as.character(SPEC$year_from)

TSPEC <- COV$titration
stopifnot("covariates.json must declare a `titration` block" = !is.null(TSPEC))
WIN_MIN   <- as.numeric(TSPEC$bolus_window_minutes)
MIN_DELTA <- as.numeric(TSPEC$min_rate_change_mcg_hr)

# Read from the protocol, not inherited from 05_titration.R: the flags in
# indication_events.parquet were computed at this window, and a caption naming
# a different one would be wrong in the direction nobody checks.
ISPEC <- COV$indication
stopifnot("covariates.json must declare an `indication` block" = !is.null(ISPEC))
IND_WIN_H <- as.numeric(ISPEC$indication_window_hours)

stopifnot(
  "at least one unit key must be declared" = length(UNIT_KEYS) >= 1,
  "min_events_per_unit must clear the disclosure floor" = MIN_UNIT >= MIN_CELL,
  "bootstrap_resamples must be a positive whole number" = B_RESAMP > 0
)

EVENT_LEVELS <- c("initiation", "uptitration", "downtitration")
CONTROL_KIND <- "downtitration"      # the negative control, excluded from "increases"

# Two-sided normal deviates for the funnel's control limits. 95% and 99.8% are
# Spiegelhalter's conventional pair for institutional comparison.
Z_INNER <- qnorm(0.975)
Z_OUTER <- qnorm(0.999)


# ---- 3. Paths and provenance -------------------------------------------------

dirs <- site_dirs()
dirs$phase <- phase_dir(dirs, "06_unit_variation")
prov <- provenance(config)

message(sprintf("[06_unit_variation] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))


# ---- 4. Guards ---------------------------------------------------------------

manifest <- require_manifest(dirs, here())
message(sprintf("  reading Phase 0 outputs from code %s, generated %s",
                manifest$code_version, manifest$generated))

# Two metrics, named so each file says which question it answers. Figure names
# are derived from METRICS x UNIT_KEYS rather than written out, so adding either
# cannot leave a figure unowned and therefore uncleared. CLAUDE.md's first
# config-integrity rule: iterate the keys rather than restating them.
METRICS <- c("coadmin", "indication")
FIGS <- c(as.vector(outer(METRICS, UNIT_KEYS,
                          function(m, k) paste0(m, "_caterpillar_", k, ".png"))),
          as.vector(outer(METRICS, UNIT_KEYS,
                          function(m, k) paste0(m, "_funnel_", k, ".png"))),
          paste0(METRICS, "_year.png"))
OWNED <- list(
  phase = c("unit_adherence.csv", "year_adherence.csv", "attribution_funnel.csv",
            FIGS, "provenance.json", "captions.md"))
# Renamed 2026-09-30 when the indication metric joined: "unit_caterpillar" did
# not say a caterpillar OF WHAT, which was merely terse with one metric and
# ambiguous with two. Built with paste0 so a name-level find-and-replace cannot
# reach into the list of old names -- lessons.md #13.
# Two rounds. The `unit_` prefix went on 2026-09-30 when the indication metric
# joined and "a caterpillar of what" became ambiguous. The `location_type`
# suffix went the same day, when the ED joined the ICU as an attributed setting
# and the key became care_setting -- location_type is null for every non-ICU
# row, so it could not name the new label.
OLD_KEYS <- c("location_type", "location_name")
RETIRED <- file.path("output", "final_no_phi", "06_unit_variation",
                     c(paste0(paste0("unit_", "caterpillar_"), OLD_KEYS, ".png"),
                       paste0(paste0("unit_", "funnel_"), OLD_KEYS, ".png"),
                       paste0("year_", "trend.png"),
                       as.vector(outer(METRICS, paste0("location_", "type"),
                                       function(m, k) paste0(m, "_caterpillar_", k, ".png"))),
                       as.vector(outer(METRICS, paste0("location_", "type"),
                                       function(m, k) paste0(m, "_funnel_", k, ".png")))))
n_cleared <- clear_owned_outputs(dirs, OWNED, retired = RETIRED)
if (n_cleared) message(sprintf("  cleared %d output(s) from a previous run", n_cleared))


# ---- 5. Read -----------------------------------------------------------------

ev <- as.data.frame(read_parquet(
  file.path(dirs$out_phi, "titration_events_classified.parquet")))
ind_ev <- as.data.frame(read_parquet(
  file.path(dirs$out_phi, "indication_events.parquet")))
icu <- as.data.frame(read_parquet(
  file.path(dirs$out_phi, "care_intervals.parquet")))
yrs <- as.data.frame(read_parquet(
  file.path(dirs$out_phi, "trajectory_long.parquet"),
  col_select = c("encounter_block", "anchor_year")))
yrs <- unique(yrs)

stopifnot(
  "classified events are missing a column" =
    all(c("encounter_block", "t_hr", "kind", "paired") %in% names(ev)),
  "an event kind is not one 05_titration.R produces" = all(ev$kind %in% EVENT_LEVELS),
  "indication events are missing a column" =
    all(c("encounter_block", "t_hr", "event_type", "any_indicated") %in% names(ind_ev)),
  "an indication event is neither an increase nor a bolus" =
    all(ind_ev$event_type %in% c("increase", "bolus")),
  "care_intervals is missing a declared unit key" = all(UNIT_KEYS %in% names(icu)),
  "care_intervals carries a location_category the config does not attribute" =
    all(tolower(icu$location_category) %in% tolower(ATTR_CATS)),
  "a care interval is zero-width or inverted" = all(icu$t_end_hr > icu$t_start_hr),
  "anchor_year is not one row per episode" =
    !anyDuplicated(yrs$encounter_block),
  "an episode has no anchor year" = all(!is.na(yrs$anchor_year))
)

# ---- 6. Attribution ----------------------------------------------------------
# Each event is credited to the unit the patient was PHYSICALLY IN when it was
# charted, not to a single label for the whole episode. 6.4% of attributed
# episodes touch more than one unit inside the 72h window (11.0% of whole ICU
# hospitalizations do, over their longer span -- two denominators, not one
# number). An episode-level label would credit every one of that episode's
# events to whichever unit held the most time.
# See covariates.json unit_variation._why_event_level.
#
# Written once and applied per metric: a bolus carries an encounter_block and a
# t_hr exactly as a rate change does, so the same interval join credits it.

attribute <- function(d, what) {
  d$row_id <- seq_len(nrow(d))
  n_events <- nrow(d)
  n_episodes <- length(unique(d$encounter_block))

  pairs <- merge(d[, c("row_id", "encounter_block", "t_hr")], icu,
                 by = "encounter_block", all.x = FALSE)
  pairs <- pairs[pairs$t_hr >= pairs$t_start_hr & pairs$t_hr < pairs$t_end_hr, ]

  # Phase 0 refuses to emit overlapping ICU intervals, so an event can match at
  # most one. Re-asserted rather than assumed: a duplicated event would be
  # counted twice and inflate whichever units overlapped.
  stopifnot("an event matched more than one care interval" = !anyDuplicated(pairs$row_id))

  att <- merge(d, pairs[, c("row_id", UNIT_KEYS)], by = "row_id", all.x = TRUE)
  stopifnot("the attribution join changed the event count" = nrow(att) == n_events)
  att <- merge(att, yrs, by = "encounter_block", all.x = TRUE)
  stopifnot("the year join changed the event count" = nrow(att) == n_events,
            "an event has no anchor year" = all(!is.na(att$anchor_year)))

  # UNATTRIBUTED IS A ROW, NOT A DROP. A ventilated patient can be off the unit
  # -- in a procedural or radiology location -- when the event is charted.
  # Dropping those would shrink the denominator in the direction that flatters
  # every unit.
  att$attributed <- !is.na(att[[UNIT_KEYS[1]]])
  n_unatt <- sum(!att$attributed)
  cat(sprintf("  %-11s %s event(s), %s episode(s); attributed %s (%.1f%%), unattributed %s (%.1f%%)\n",
              what, format(n_events, big.mark = ","),
              format(n_episodes, big.mark = ","),
              format(sum(att$attributed), big.mark = ","), 100 * mean(att$attributed),
              format(n_unatt, big.mark = ","), 100 * mean(!att$attributed)))

  funnel <- data.frame(
    metric = what,
    step = c("events from 05_titration.R", "inside an ICU interval",
             "unattributed -- ventilated but not in an ICU location"),
    n_events = c(n_events, sum(att$attributed), n_unatt),
    n_episodes = c(n_episodes,
                   length(unique(att$encounter_block[att$attributed])),
                   length(unique(att$encounter_block[!att$attributed]))))

  # The accounting must close. Watched failing before it was kept (lessons.md #6).
  stopifnot("attributed + unattributed does not equal the event total" =
              sum(att$attributed) + n_unatt == n_events)
  list(att = att, funnel = funnel, n_events = n_events, n_unatt = n_unatt)
}


# ---- 7. Cluster bootstrap ----------------------------------------------------
# Resample EPISODES, not events. 05_titration.R records that event-level
# percentages "carry no valid confidence interval without accounting for
# clustering within episode": episodes contribute several events each and the
# habit is correlated inside one. Resampling episodes keeps an episode's events
# together, which is the whole point.
#
# It measures the clustering rather than assuming a value for it, and it needs no
# package the pinned stack does not already have.

# THE EPISODE UNIVERSE IS PER METRIC, not shared. The coadmin metric runs on
# rate increases; the indication metric also covers boluses, and an episode can
# contribute a bolus without ever contributing an increase. Resampling one
# metric against the other's episode list would give those episodes no chance
# of being drawn -- match() would return NA and their events would silently
# carry NA weight. Passed in rather than inferred inside, so the caller has to
# say which universe it means.

# boot_group() is in code/utils/pooling.R -- shared with 05_titration.R.

# Two metrics, one machine. Each declares its events, the 0/1 column, its
# series, and which series is the COMBINED one -- the union of the others, used
# for the ordering, the reference rule and the funnel's single point per unit.
#
# Run as SUBSETS per series rather than one frame with duplicated rows: a
# duplicated row would double sum(n) inside boot_group() and halve the naive
# binomial variance, silently doubling the reported design effect.
# Increases only for the coadmin metric. The negative control belongs to
# 05_titration.R, where it establishes the coincidental-pairing floor; splitting
# it by unit would invite reading a control as a finding.
inc <- ev[ev$kind != CONTROL_KIND, ]

ANALYSES <- list(
  coadmin = list(
    events = inc, outcome = "paired", split_on = "kind",
    series = c("initiation", "uptitration", "any increase"),
    combined = "any increase",
    x_lab = sprintf("%% of EVENTS paired with a bolus (+/-%g min)", WIN_MIN),
    title = "Bolus co-administration at a fentanyl infusion increase"),
  indication = list(
    events = ind_ev, outcome = "any_indicated", split_on = "event_type",
    series = c("increase", "bolus", "any escalation"),
    combined = "any escalation",
    x_lab = sprintf("%% of EVENTS with a qualifying score (+/-%g h)", IND_WIN_H),
    title = "Documented pain or agitation at a fentanyl escalation")
)
stopifnot("METRICS and ANALYSES have drifted apart" =
            setequal(METRICS, names(ANALYSES)))

series_subset <- function(d, spec, s) {
  if (s == spec$combined) d else d[d[[spec$split_on]] == s, ]
}

cat("\n")
ATT <- list(); UNIT <- list(); YEAR <- list(); FUNNEL <- list()
for (m in METRICS) {
  spec <- ANALYSES[[m]]
  a <- attribute(spec$events, m)
  ATT[[m]] <- a$att; FUNNEL[[m]] <- a$funnel
  eps <- sort(unique(spec$events$encounter_block))

  u <- do.call(rbind, lapply(UNIT_KEYS, function(k) {
    do.call(rbind, lapply(spec$series, function(sname) {
      r <- boot_group(series_subset(a$att[a$att$attributed, ], spec, sname),
                      k, SEED, outcome = spec$outcome, episodes = eps,
                      n_resamples = B_RESAMP)
      r$metric <- m; r$unit_key <- k; r$series <- sname
      r$pooled_pct <- attr(r, "pooled_pct"); r$deff <- attr(r, "deff")
      r
    }))
  }))
  names(u)[names(u) == "group"] <- "unit"
  u$series <- factor(u$series, levels = spec$series)
  UNIT[[m]] <- u

  y <- boot_group(a$att, "anchor_year", SEED, outcome = spec$outcome,
                  episodes = eps, n_resamples = B_RESAMP)
  y$metric <- m
  y$pooled_pct <- attr(y, "pooled_pct"); y$deff <- attr(y, "deff")
  names(y)[names(y) == "group"] <- "year"
  y$year <- as.integer(y$year)
  YEAR[[m]] <- y[order(y$year), ]

  # A component series must reproduce its parent, or the subsetting has dropped
  # something. Watched failing before it was kept (lessons.md #6).
  comb <- u[u$series == spec$combined, ]
  chk <- aggregate(n_events ~ unit_key + unit,
                   data = u[u$series != spec$combined, ], FUN = sum)
  names(chk)[names(chk) == "n_events"] <- "parts"
  chk <- merge(chk, comb[, c("unit_key", "unit", "n_events")],
               by = c("unit_key", "unit"))
  stopifnot("the component series do not sum to the combined one" =
              all(chk$parts == chk$n_events))

  for (sname in spec$series) {
    r <- u[u$series == sname, ]
    cat(sprintf("  %-11s %-14s pooled %5.1f%%   design effect %.2f\n",
                m, sname, r$pooled_pct[1], r$deff[1]))
  }
}

unit_tbl <- do.call(rbind, UNIT)
year_tbl <- do.call(rbind, YEAR)
funnel_tbl <- do.call(rbind, FUNNEL)
rownames(unit_tbl) <- NULL; rownames(year_tbl) <- NULL; rownames(funnel_tbl) <- NULL

stopifnot(
  "a bootstrap interval does not contain its point estimate" =
    all(unit_tbl$ci_lo <= unit_tbl$pct_paired & unit_tbl$pct_paired <= unit_tbl$ci_hi) &&
    all(year_tbl$ci_lo <= year_tbl$pct_paired & year_tbl$pct_paired <= year_tbl$ci_hi),
  "the design effect is below 1, so clustering cannot be what widened it" =
    all(unit_tbl$deff >= 1) && all(year_tbl$deff >= 1)
)


# ---- 8. Disclosure -----------------------------------------------------------
# Suppression, not a printed note. 04 and 05 report small cells and publish them
# because their tables are prevalence and NA-ing a cell stops a column summing
# to 100. Here a suppressed unit simply does not appear, so the pooling.R
# convention applies instead: blank the RATE, keep the COUNT, flag the row.
#
# It keys on the unit's COMBINED event count and takes the whole unit with it,
# per metric. Suppressing series independently would leave one series published
# while another was withheld, and the combined row would let a reader recover
# the withheld one by subtraction.

unit_tbl$n_suppressed_small_cell <- 0L
for (m in METRICS) {
  cmb <- ANALYSES[[m]]$combined
  i <- unit_tbl$metric == m
  small <- unit_tbl[i & unit_tbl$series == cmb & unit_tbl$n_events < MIN_UNIT,
                    c("unit_key", "unit")]
  hit <- i & paste(unit_tbl$unit_key, unit_tbl$unit) %in%
             paste(small$unit_key, small$unit)
  if (any(hit)) {
    cat(sprintf("  NOTE %s: %d unit(s) below unit_variation.min_events_per_unit = %d, suppressed across all series\n",
                m, nrow(small), MIN_UNIT))
    unit_tbl$n_suppressed_small_cell[hit] <- 1L
    unit_tbl[hit, c("pct_paired", "ci_lo", "ci_hi", "n_paired")] <- NA
  }
}
plot_units <- unit_tbl[unit_tbl$n_suppressed_small_cell == 0, ]


# ---- 9. Figures --------------------------------------------------------------
# Journal style: no title, no subtitle, boxed legend at the foot, all type in
# INK. The n, the denominator, the seed and the method live in the caption.

# Deterministic label placement for the funnel, with a collision check.
#
# A label sits to the RIGHT of its point and flips LEFT exactly when drawing it
# right would run past the panel edge -- measured geometry, not a tuned
# fraction of the width, which happens to work for six units at one site and
# silently stops working for a longer unit name elsewhere.
#
# Nothing is nudged off its point: if two labels would collide the run FAILS
# rather than drawing them on top of each other, because a funnel with
# overlapping labels is worse than one with none -- it puts a unit somewhere it
# is not. A funnel's only legend is the control-limit linetype, so without
# these the reader cannot tell which dot is which unit; the figure shipped that
# way on 2026-09-30 and SG caught it.
#
# Extents are estimated in data units from the panel geometry: a 7.5 x 4.4in
# figure leaves roughly 6.3 x 3.2in of panel, and 2.5pt text is about 1.55mm
# per character wide and 3.0mm tall.
PANEL_W_MM <- 160; PANEL_H_MM <- 81
CHAR_W_MM  <- 1.55; LINE_H_MM <- 3.0

funnel_labels <- function(x, y, lab, xr, yr) {
  pad <- diff(xr) * 0.14
  xr <- c(xr[1] - pad, xr[2] + pad)
  x_per_mm <- diff(xr) / PANEL_W_MM
  y_per_mm <- diff(yr) / PANEL_H_MM
  gap <- 1.6 * x_per_mm
  w <- nchar(lab) * CHAR_W_MM * x_per_mm
  h <- LINE_H_MM * y_per_mm
  flip <- (x + gap + w) > xr[2]
  out <- data.frame(
    lab = lab, px = x, py = y,
    tx = ifelse(flip, x - gap, x + gap), ty = y,
    hjust = ifelse(flip, 1, 0), stringsAsFactors = FALSE)
  out$x0 <- ifelse(flip, out$tx - w, out$tx)
  out$x1 <- ifelse(flip, out$tx, out$tx + w)

  # Bounded vertical de-collision, per side. Two labels whose x extents overlap
  # cannot share a y band, so the upper one is pushed up until it clears. An
  # earlier version simply REFUSED on any collision, which was right when there
  # were six short ICU names and too brittle at seven settings with
  # "mixed_cardiothoracic_icu" among them -- it failed the whole run rather
  # than moving a label 2mm. Moving is honest as long as the reader can still
  # tell which point a label belongs to, which is what the leader line below is
  # for; refusing stays the behaviour when de-collision cannot converge.
  sep <- h * 1.15
  for (side in c(FALSE, TRUE)) {
    i <- which(flip == side)
    if (length(i) < 2) next
    i <- i[order(out$ty[i])]
    for (a in seq_along(i)[-1]) {
      cur <- i[a]
      for (b in seq_len(a - 1)) {
        prev <- i[b]
        x_overlap <- !(out$x1[cur] < out$x0[prev] || out$x1[prev] < out$x0[cur])
        if (x_overlap && (out$ty[cur] - out$ty[prev]) < sep)
          out$ty[cur] <- out$ty[prev] + sep
      }
    }
  }
  out$y0 <- out$ty - h / 2
  out$y1 <- out$ty + h / 2
  # Draw a connector only where the label has actually left its point.
  out$lead <- abs(out$ty - out$py) > h * 0.6

  olap <- function(a, b) !(a$x1 < b$x0 || b$x1 < a$x0 || a$y1 < b$y0 || b$y1 < a$y0)
  for (i in seq_len(nrow(out))) {
    for (j in seq_len(nrow(out))) {
      if (j <= i) next
      if (olap(out[i, ], out[j, ]))
        stop(sprintf("funnel labels '%s' and '%s' still overlap after de-collision",
                     out$lab[i], out$lab[j]))
    }
    for (j in seq_len(nrow(out))) {
      if (i == j) next
      if (out$px[j] >= out$x0[i] && out$px[j] <= out$x1[i] &&
          out$py[j] >= out$y0[i] && out$py[j] <= out$y1[i])
        stop(sprintf("funnel label '%s' covers the point for '%s'",
                     out$lab[i], out$lab[j]))
    }
    if (out$x0[i] < xr[1] || out$x1[i] > xr[2])
      stop(sprintf("funnel label '%s' runs outside the panel", out$lab[i]))
  }
  out
}

deff_note <- function(x) sprintf(
  paste0("Intervals are 2.5th-97.5th percentiles of %s bootstrap resamples of ",
         "EPISODES (seed %s), which keeps an episode's events together because ",
         "they cluster within it; the measured design effect is %.2f. The ",
         "interval says how much the observed rate would wobble if the care ",
         "process re-ran. IT DOES NOT ADJUST FOR CASE MIX: a unit's rate ",
         "reflects who it admits as much as how it practises, and no ",
         "association between unit and this metric is modelled here."),
  format(B_RESAMP, big.mark = ","), format(SEED), x)

METRIC_NOTE <- list(
  coadmin = sprintf(
    paste0("Share of fentanyl infusion INCREASES (initiation or uptitration of ",
           "at least %g mcg/hr) accompanied by a fentanyl bolus within +/-%g ",
           "minutes, on raw charted timestamps. THE DENOMINATOR IS EVENTS, NOT ",
           "PATIENTS. No evidence establishes that pairing a bolus with an ",
           "increase is better care: this is variation in practice, not a ranking."),
    MIN_DELTA, WIN_MIN),
  indication = sprintf(
    paste0("Share of fentanyl escalations -- a rate increase or a bolus -- with ",
           "a qualifying pain score or agitation within +/-%g h. Each instrument ",
           "is cut at ITS OWN threshold; the scales are not comparable and raw ",
           "values are never pooled. THIS IS A LEVEL, NOT A CHANGE. THE ",
           "DENOMINATOR IS EVENTS, NOT PATIENTS. A low share may mean scores are ",
           "not charted near the decision rather than that the decision was ",
           "unfounded -- 05_titration.R reports the fraction among events that ",
           "were assessed at all, which separates the two."),
    IND_WIN_H))

SERIES_NOTE <- list(
  coadmin = paste0("THREE SERIES PER UNIT: an **initiation** (a rise from a ",
                   "documented zero), an **uptitration** (a rise in a running ",
                   "infusion), and **any increase**, the two pooled rather than ",
                   "a third kind of event."),
  indication = paste0("TWO SERIES PER UNIT: a rate **increase** and a **bolus**, ",
                      "plus **any escalation**, the two pooled rather than a ",
                      "third kind of event."))

key_label <- function(k) if (k == "care_setting")
  paste0("care_setting -- the CLIF location_type where it exists, else ",
         "location_category, both mCIDE-controlled and so reproducible at ",
         "another site") else
  "site-local location_name (internal, not poolable)"

for (m in METRICS) {
  spec <- ANALYSES[[m]]
  pooled_m <- unit_tbl$pooled_pct[unit_tbl$metric == m &
                                  unit_tbl$series == spec$combined][1]
  deff_m <- unit_tbl$deff[unit_tbl$metric == m &
                          unit_tbl$series == spec$combined][1]
  n_unatt_m <- FUNNEL[[m]]$n_events[3]
  pal <- setNames(unname(SERIES_COLOURS)[seq_along(spec$series)], spec$series)

  for (k in UNIT_KEYS) {
    d <- plot_units[plot_units$metric == m & plot_units$unit_key == k, ]
    ord <- d[d$series == spec$combined, ]
    ord <- ord[order(ord$pct_paired), ]
    d$unit <- factor(d$unit, levels = ord$unit)
    d$series <- factor(as.character(d$series), levels = spec$series)
    n_units <- nrow(ord)
    dcomb <- d[d$series == spec$combined, ]

    # The event count rides on the AXIS, not inside the panel. Drawn beside the
    # interval it collided with the data at the busy end and forced a 30%
    # right-hand margin to be left empty for it; on the axis it reads as what it
    # is -- a property of the row, not a mark on the plot.
    # trim = TRUE: format() otherwise pads every number to the widest, which
    # puts a space inside the bracket ("( 4,445)") and reads as a typo. The
    # axis is right-aligned, so the brackets line up without the pad.
    axis_lab <- setNames(sprintf("%s  (%s)", ord$unit,
                                 format(ord$n_events, big.mark = ",", trim = TRUE)),
                         as.character(ord$unit))

    # The pooled series is the finding; the components are the decomposition.
    # Opaque and large for the first, lighter and smaller for the rest, so the
    # eye lands on the combined point and the components read as support.
    is_comb <- spec$series == spec$combined
    A_VAL <- setNames(ifelse(is_comb, 1.00, 0.55), spec$series)
    S_VAL <- setNames(ifelse(is_comb, 2.90, 1.70), spec$series)
    L_VAL <- setNames(ifelse(is_comb, 0.70, 0.40), spec$series)

    fn <- paste0(m, "_caterpillar_", k, ".png")
    register_caption(fn, sprintf("%s, by %s", spec$title, k),
      paste(METRIC_NOTE[[m]], SERIES_NOTE[[m]],
            paste0("The pooled series is drawn opaque and larger; the ",
                   "components are lighter and smaller, so the panel reads as ",
                   "one estimate per unit with its decomposition beside it."),
            sprintf(paste0("The pooled series necessarily lies between its ",
                           "components, at their event-weighted mean. Units are ",
                           "ordered by it and the dashed rule is its pooled rate ",
                           "(%.1f%%). %s units shown of %s; %s event(s) were ",
                           "charted while the patient was ventilated but in a ",
                           "care location this analysis does not attribute -- ",
                           "chiefly the operating and procedure suites, where ",
                           "these instruments are not charted under anaesthesia ",
                           "-- and are excluded from the unit views while ",
                           "remaining in attribution_funnel.csv. Labels are ",
                           "%s. The number in brackets on the axis is the ",
                           "unit's total event count; per-series counts are in ",
                           "unit_adherence.csv."),
                    pooled_m, n_units,
                    length(unique(unit_tbl$unit[unit_tbl$metric == m &
                                                unit_tbl$unit_key == k])),
                    format(n_unatt_m, big.mark = ","), key_label(k)),
            deff_note(deff_m)))

    # One layer per geom, with alpha/size/linewidth MAPPED rather than split into
    # separate layers per series: position_dodge divides the slot by the number
    # of groups present in a layer, so drawing the combined series on its own
    # would dodge it into a different position from its components.
    dodge <- position_dodge(width = 0.62)
    pl <- house(
      ggplot(d, aes(pct_paired, unit, colour = series, group = series)) +
        geom_vline(xintercept = pooled_m, colour = MUTED, linetype = "22",
                   linewidth = 0.4) +
        geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi, alpha = series,
                          linewidth = series),
                      orientation = "y", width = 0, position = dodge) +
        geom_point(aes(alpha = series, size = series), position = dodge) +
        scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
        # The legend carries identity by colour alone; alpha and size encode
        # emphasis, not a variable, so they must not appear in it.
        scale_alpha_manual(values = A_VAL, guide = "none") +
        scale_size_manual(values = S_VAL, guide = "none") +
        scale_linewidth_manual(values = L_VAL, guide = "none") +
        scale_y_discrete(labels = axis_lab) +
        scale_x_continuous(expand = expansion(mult = c(0.05, 0.06))) +
        labs(x = spec$x_lab, y = NULL))
    ggsave(file.path(dirs$phase, fn), pl,
           width = 7.5, height = max(3.0, 1.4 + 0.46 * n_units), dpi = 200)

    # The funnel stays on the COMBINED series: its point is one mark per unit
    # against that unit's volume, and overlapping funnels defeat the envelope.
    grid_n <- seq(max(1, min(dcomb$n_events) * 0.6), max(dcomb$n_events) * 1.15,
                  length.out = 200)
    se <- sqrt(pooled_m * (100 - pooled_m) / grid_n) * sqrt(deff_m)
    lim <- rbind(
      data.frame(n = grid_n, lo = pooled_m - Z_INNER * se,
                 hi = pooled_m + Z_INNER * se, band = "95%"),
      data.frame(n = grid_n, lo = pooled_m - Z_OUTER * se,
                 hi = pooled_m + Z_OUTER * se, band = "99.8%"))

    fn <- paste0(m, "_funnel_", k, ".png")
    register_caption(fn, sprintf("%s -- funnel against unit volume, by %s",
                                 spec$title, k),
      paste(METRIC_NOTE[[m]],
            sprintf(paste0("Each point is one unit, labelled beside it; the ",
                           "horizontal rule is the pooled rate (%.1f%%) and the ",
                           "envelopes are 95%% and 99.8%% control limits. The ",
                           "limits are binomial limits WIDENED BY ",
                           "sqrt(design effect) = %.2f, so they carry the same ",
                           "clustering correction as the caterpillar's intervals ",
                           "rather than assuming events are independent. Volume ",
                           "on the x axis is why this view sits beside the ",
                           "ordered one: a small unit sits far from the pooled ",
                           "rate more easily, and the envelope makes that ",
                           "visible. Only the pooled series is drawn -- ",
                           "overlapping funnels defeat the envelope reading. ",
                           "Labels are %s."),
                    pooled_m, sqrt(deff_m), key_label(k)),
            deff_note(deff_m)))

    lab <- funnel_labels(dcomb$n_events, dcomb$pct_paired, as.character(dcomb$unit),
                         xr = range(grid_n),
                         yr = range(c(lim$lo, lim$hi, dcomb$pct_paired)))
    pl <- house(
      ggplot(lim, aes(n)) +
        geom_line(aes(y = lo, linetype = band), colour = MUTED, linewidth = 0.4) +
        geom_line(aes(y = hi, linetype = band), colour = MUTED, linewidth = 0.4) +
        geom_hline(yintercept = pooled_m, colour = MUTED, linewidth = 0.4) +
        geom_point(data = dcomb, aes(x = n_events, y = pct_paired), size = 2.1,
                   colour = pal[[spec$combined]]) +
        geom_segment(data = lab[lab$lead, ],
                     aes(x = px, xend = tx, y = py, yend = ty),
                     inherit.aes = FALSE, colour = MUTED, linewidth = 0.25) +
        geom_text(data = lab, aes(x = tx, y = ty, label = lab, hjust = hjust),
                  inherit.aes = FALSE, size = 2.5, colour = INK) +
        scale_linetype_manual(values = c(`95%` = "22", `99.8%` = "11"), name = NULL) +
        scale_x_continuous(expand = expansion(mult = c(0.14, 0.14))) +
        labs(x = "Events contributed by the unit", y = spec$x_lab))
    ggsave(file.path(dirs$phase, fn), pl, width = 7.5, height = 4.4, dpi = 200)
  }

  # Year. Every event carries a year, attributed or not, so this view uses the
  # whole denominator rather than the ICU-attributed subset.
  yt <- year_tbl[year_tbl$metric == m, ]
  ypool <- yt$pooled_pct[1]
  fn <- paste0(m, "_year.png")
  register_caption(fn, sprintf("%s, by calendar year", spec$title),
    paste(METRIC_NOTE[[m]],
          sprintf(paste0("One point per calendar year of the ventilation ",
                         "anchor, %d to %d, with the pooled rate (%.1f%%) as a ",
                         "dashed rule. THE DENOMINATOR HERE IS ALL %s EVENTS, ",
                         "not only the ICU-attributed ones the unit figures use, ",
                         "because every event carries a year. Year counts follow ",
                         "the extract, not a study period: no date-based ",
                         "inclusion criterion is applied anywhere in this ",
                         "pipeline, so an end year may be partial."),
                  min(yt$year), max(yt$year), ypool,
                  format(sum(yt$n_events), big.mark = ",")),
          deff_note(yt$deff[1])))

  pl <- house(
    ggplot(yt, aes(year, pct_paired)) +
      geom_hline(yintercept = ypool, colour = MUTED, linetype = "22",
                 linewidth = 0.4) +
      geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0, linewidth = 0.5,
                    colour = INK) +
      geom_line(linewidth = 0.6, colour = pal[[spec$combined]]) +
      geom_point(size = 2.1, colour = pal[[spec$combined]]) +
      scale_x_continuous(breaks = yt$year) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
      labs(x = NULL, y = spec$x_lab))
  ggsave(file.path(dirs$phase, fn), pl, width = 7.5, height = 4.4, dpi = 200)
}


# ---- 10. Write ---------------------------------------------------------------

# The two files pool over DIFFERENT event sets, so the column cannot share a
# name. unit_adherence covers attributed events only -- an event in no ICU or ED
# interval cannot be credited to a unit -- while year_adherence covers all
# events, because every event has a year. At UCMC that is 13.17% against 13.53%:
# the same metric, two correct numbers, and a reader given one `pooled_pct`
# column has no way to tell which they are holding.
unit_tbl$pooled_pct_attributed <- unit_tbl$pooled_pct
year_tbl$pooled_pct_all_events <- year_tbl$pooled_pct

cols <- c("metric", "unit_key", "unit", "series", "n_events", "n_episodes",
          "n_paired", "pct_paired", "ci_lo", "ci_hi", "pooled_pct_attributed",
          "deff", "n_suppressed_small_cell")
write.csv(unit_tbl[order(unit_tbl$metric, unit_tbl$unit_key,
                         unit_tbl$unit, unit_tbl$series), cols],
          file.path(dirs$phase, "unit_adherence.csv"), row.names = FALSE)
write.csv(year_tbl[order(year_tbl$metric, year_tbl$year),
                   c("metric", "year", "n_events", "n_episodes", "n_paired",
                     "pct_paired", "ci_lo", "ci_hi", "pooled_pct_all_events",
                     "deff")],
          file.path(dirs$phase, "year_adherence.csv"), row.names = FALSE)
write.csv(funnel_tbl, file.path(dirs$phase, "attribution_funnel.csv"),
          row.names = FALSE)

write_json(prov, file.path(dirs$phase, "provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)
write_captions(file.path(dirs$phase, "captions.md"), "06_unit_variation.R",
               grep("\\.png$", OWNED$phase, value = TRUE), prov)

for (f in OWNED$phase) cat(sprintf("written: %s\n", f))

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/06_unit_variation.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "06_unit_variation_sessioninfo.txt")
)
