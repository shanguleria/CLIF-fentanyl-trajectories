# ==============================================================================
# 06_unit_variation.R  --  bolus co-administration by ICU unit and by calendar year
#
# Purpose : Split 05_titration.R's overall number by unit and by calendar year.
#           Caterpillar and funnel views per unit, plus a year series.
# Author  : Shan Guleria
# Created : 2026-09-30
# Inputs  : output/intermediate_phi/titration_events_classified.parquet, icu_intervals.parquet, trajectory_long.parquet
# Outputs : output/final_no_phi/06_unit_variation/ : the files listed in OWNED
#
# NOT A QUALITY MEASURE: this describes variation in practice, not a ranking.
# tests/test_protocol.py bans the words quality, performance and benchmark here.
# Events are NOT redefined here; 05_titration.R owns the definitions. The intervals
# describe sampling only -- see protocol.json unit_variation._STATUS_model.
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
# boot_group() lives in pooling.R so 05 and 06 cannot drift.
source(here("code", "utils", "pooling.R"))


# ---- 2. Config ---------------------------------------------------------------

config <- fromJSON(here("config", "config.json"), simplifyVector = FALSE)
COV <- fromJSON(here("config", "protocol.json"), simplifyVector = FALSE)
SEED <- COV$model$seed
set.seed(SEED)

WINDOW_H <- COV$cohort$window_hours
EXTENT_H <- COV$cohort$granular_extent_hours
MIN_CELL <- COV$reporting$small_cell_min_den

# Federation-critical: read from the config, never defaulted here.
SPEC <- COV$unit_variation
stopifnot("protocol.json must declare a `unit_variation` block" = !is.null(SPEC))
ATTR_CATS <- unlist(SPEC$attributed_location_categories)
UNIT_KEYS <- unlist(SPEC$unit_keys)
MIN_UNIT  <- as.numeric(SPEC$min_events_per_unit)
B_RESAMP  <- as.integer(SPEC$bootstrap_resamples)
YEAR_FROM <- as.character(SPEC$year_from)
OD_TRIM   <- as.numeric(SPEC$overdispersion_trim)

TSPEC <- COV$titration
stopifnot("protocol.json must declare a `titration` block" = !is.null(TSPEC))
WIN_MIN   <- as.numeric(TSPEC$bolus_window_minutes)
MIN_DELTA <- as.numeric(TSPEC$min_rate_change_mcg_hr)

# From the protocol, not from 05: the flags were computed at this window.
ISPEC <- COV$indication
stopifnot("protocol.json must declare an `indication` block" = !is.null(ISPEC))
IND_WIN_H <- as.numeric(ISPEC$indication_window_hours)

stopifnot(
  "at least one unit key must be declared" = length(UNIT_KEYS) >= 1,
  "min_events_per_unit must clear the disclosure floor" = MIN_UNIT >= MIN_CELL,
  "bootstrap_resamples must be a positive whole number" = B_RESAMP > 0,
  # Winsorising at 0.5 per end reports phi = 0 regardless of the data.
  "overdispersion_trim must be a fraction per end, in [0, 0.5)" =
    is.finite(OD_TRIM) && OD_TRIM >= 0 && OD_TRIM < 0.5
)

EVENT_LEVELS <- c("initiation", "uptitration", "downtitration")
CONTROL_KIND <- "downtitration"      # the negative control, excluded from "increases"

# Two-sided deviates for the funnel: Spiegelhalter's conventional 95% / 99.8% pair.
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

# Figure names derive from METRICS x UNIT_KEYS, so neither can leave a figure unowned.
METRICS <- c("coadmin", "indication")

# location_name is OPTIONAL in the CLIF adt schema, so availability is read from the
# SCHEMA before OWNED is derived from it, never asserted after the figures are written.
.ci_path <- file.path(dirs$out_phi, "care_intervals.parquet")
if (file.exists(.ci_path)) {
  .have <- names(arrow::open_dataset(.ci_path))
  .gone <- setdiff(UNIT_KEYS, .have)
  if (length(.gone)) {
    cat(sprintf(paste0("\n  NOTE unit key(s) not present in care_intervals and ",
                       "therefore SKIPPED: %s\n        present: %s\n",
                       "        This is expected where CLIF's optional ",
                       "location_name is not populated.\n"),
                paste(.gone, collapse = ", "), paste(.have, collapse = ", ")))
    UNIT_KEYS <- intersect(UNIT_KEYS, .have)
  }
  stopifnot("no declared unit key is present in care_intervals -- there is nothing to group by" =
              length(UNIT_KEYS) >= 1)
}

FIGS <- c(as.vector(outer(METRICS, UNIT_KEYS,
                          function(m, k) paste0(m, "_caterpillar_", k, ".png"))),
          as.vector(outer(METRICS, UNIT_KEYS,
                          function(m, k) paste0(m, "_funnel_", k, ".png"))),
          paste0(METRICS, "_year.png"))
OWNED <- list(
  phase = c("unit_adherence.csv", "year_adherence.csv", "attribution_funnel.csv",
            "overdispersion.csv",
            FIGS, "provenance.json", "captions.md"))
# Retired names, built with paste0 so find-and-replace cannot reach them (lessons.md #13).
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
  # UNIT_KEYS was narrowed to the schema above; this catches a different frame.
  "care_intervals is missing a unit key that its own schema advertised" =
    all(UNIT_KEYS %in% names(icu)),
  "care_intervals carries a location_category the config does not attribute" =
    all(tolower(icu$location_category) %in% tolower(ATTR_CATS)),
  "a care interval is zero-width or inverted" = all(icu$t_end_hr > icu$t_start_hr),
  "anchor_year is not one row per episode" =
    !anyDuplicated(yrs$encounter_block),
  "an episode has no anchor year" = all(!is.na(yrs$anchor_year))
)

# ---- 6. Attribution ----------------------------------------------------------
# Each event is credited to the unit the patient was PHYSICALLY IN when it was
# charted, not to one episode-level label -- see protocol.json
# unit_variation._why_event_level. A bolus joins the same way a rate change does.

attribute <- function(d, what) {
  d$row_id <- seq_len(nrow(d))
  n_events <- nrow(d)
  n_episodes <- length(unique(d$encounter_block))

  pairs <- merge(d[, c("row_id", "encounter_block", "t_hr")], icu,
                 by = "encounter_block", all.x = FALSE)
  pairs <- pairs[pairs$t_hr >= pairs$t_start_hr & pairs$t_hr < pairs$t_end_hr, ]

  # Phase 0 refuses overlapping intervals, so an event matches at most one.
  stopifnot("an event matched more than one care interval" = !anyDuplicated(pairs$row_id))

  att <- merge(d, pairs[, c("row_id", UNIT_KEYS)], by = "row_id", all.x = TRUE)
  stopifnot("the attribution join changed the event count" = nrow(att) == n_events)
  att <- merge(att, yrs, by = "encounter_block", all.x = TRUE)
  stopifnot("the year join changed the event count" = nrow(att) == n_events,
            "an event has no anchor year" = all(!is.na(att$anchor_year)))

  # Unattributed is a ROW, not a drop -- dropping off-unit events flatters every unit.
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

  # The accounting must close.
  stopifnot("attributed + unattributed does not equal the event total" =
              sum(att$attributed) + n_unatt == n_events)
  list(att = att, funnel = funnel, n_events = n_events, n_unatt = n_unatt)
}


# ---- 7. Cluster bootstrap ----------------------------------------------------
# Resample EPISODES, not events: the habit is correlated within an episode.

# The episode universe is PER METRIC: an episode can contribute a bolus without ever
# contributing an increase, and the other metric's list would give it NA weight.

# Each metric declares its events, the 0/1 column, its series, and which is COMBINED
# (the union of the others). Series run as SUBSETS -- duplicated rows would double
# sum(n) in boot_group() and double the reported deff. Increases only for coadmin.
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
      r$overall_pct <- attr(r, "overall_pct"); r$deff <- attr(r, "deff")
      r
    }))
  }))
  names(u)[names(u) == "group"] <- "unit"
  u$series <- factor(u$series, levels = spec$series)
  UNIT[[m]] <- u

  y <- boot_group(a$att, "anchor_year", SEED, outcome = spec$outcome,
                  episodes = eps, n_resamples = B_RESAMP)
  y$metric <- m
  y$overall_pct <- attr(y, "overall_pct"); y$deff <- attr(y, "deff")
  names(y)[names(y) == "group"] <- "year"
  y$year <- as.integer(y$year)
  YEAR[[m]] <- y[order(y$year), ]

  # A component series must reproduce its parent, or the subsetting dropped rows.
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
    cat(sprintf("  %-11s %-14s overall %5.1f%%   design effect %.2f\n",
                m, sname, r$overall_pct[1], r$deff[1]))
  }
}

unit_tbl <- do.call(rbind, UNIT)

# ---- the location_name -> care_setting crosswalk -----------------------------
# The mapping lives only in PHI-side care_intervals.parquet, so this column is what
# lets the pooled location_name panel be coloured by care setting.
.xw <- unique(icu[, c("location_name", "care_setting")])
.multi <- names(which(table(.xw$location_name) > 1L))
if (length(.multi))
  # A ward reclassified mid-period maps to two settings; reported, not resolved.
  cat(sprintf(paste0("\n  NOTE %d location_name(s) map to more than one ",
                     "care_setting and ship as 'multiple': %s\n"),
              length(.multi), paste(.multi, collapse = ", ")))
XWALK <- stats::setNames(
  vapply(split(.xw$care_setting, .xw$location_name),
         function(v) if (length(unique(v)) > 1L) "multiple" else unique(v)[1],
         character(1)),
  sort(unique(.xw$location_name)))
# care_setting rows are their own setting; location_name rows look it up.
unit_tbl$care_setting <- ifelse(
  unit_tbl$unit_key == "care_setting", unit_tbl$unit,
  unname(XWALK[unit_tbl$unit]))
stopifnot(
  "a location_name row has no care_setting -- the crosswalk missed a unit" =
    !any(is.na(unit_tbl$care_setting)))
year_tbl <- do.call(rbind, YEAR)
funnel_tbl <- do.call(rbind, FUNNEL)
rownames(unit_tbl) <- NULL; rownames(year_tbl) <- NULL; rownames(funnel_tbl) <- NULL

stopifnot(
  "a bootstrap interval does not contain its point estimate" =
    all(unit_tbl$ci_lo <= unit_tbl$pct_paired & unit_tbl$pct_paired <= unit_tbl$ci_hi) &&
    all(year_tbl$ci_lo <= year_tbl$pct_paired & year_tbl$pct_paired <= year_tbl$ci_hi),
  # deff is reported, never applied; below 1 means no clustering is being measured.
  # na.rm is required -- naive is zero at exactly 0% or 100% and NA aborts stopifnot().
  "the design effect is below 1, so the bootstrap is not measuring clustering" =
    all(unit_tbl$deff >= 1, na.rm = TRUE) && all(year_tbl$deff >= 1, na.rm = TRUE)
)
.nd <- sum(is.na(unit_tbl$deff)) + sum(is.na(year_tbl$deff))
if (.nd)
  cat(sprintf(paste0("\n  NOTE %d design effect(s) are undefined and ship as NA: ",
                     "the proportion is exactly 0%% or 100%%, so the binomial ",
                     "variance it divides by is zero. Not a fault; the ",
                     "coordinating centre must skip these when pooling.\n"), .nd))


# ---- 8. Disclosure -----------------------------------------------------------
# Per the pooling.R convention: blank the RATE, keep the COUNT, flag the row. Keys on
# the COMBINED count and takes the whole unit, so subtraction cannot recover a series.

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


# ---- 8b. Over-dispersion, measured and never applied -------------------------
# Spiegelhalter, Qual Saf Health Care 2005;14:347-51. z is a unit's deviation from the
# overall rate in null binomial SEs; Q = sum(z^2) is chi-square on I-1 df; phi = Q/I
# would widen the section 9 limits by sqrt(phi) if applied. BETWEEN-unit, unlike the
# WITHIN-episode deff. The ma_* columns overdispersion() also returns are pooled-only.
FUNNEL_COLS <- c("k", "q_stat", "q_df", "q_p", "phi", "phi_winsorised",
                 "phi_crit_95", "tau2", "mean_binomial_var",
                 "share_not_sampling")

OD <- list()
for (m in METRICS) {
  for (k in UNIT_KEYS) {
    for (sname in ANALYSES[[m]]$series) {
      d <- plot_units[plot_units$metric == m & plot_units$unit_key == k &
                        plot_units$series == sname, ]
      p0 <- unit_tbl$overall_pct[unit_tbl$metric == m &
                                  unit_tbl$series == sname][1]
      # Below three units Q has 0 or 1 df, so the row would not be a measurement.
      if (nrow(d) < 3) {
        cat(sprintf("  NOTE %s / %s / %s: %d unit(s), too few to test\n",
                    m, k, sname, nrow(d)))
        next
      }
      od <- overdispersion(d$pct_paired, d$n_events, p0, OD_TRIM)[, FUNNEL_COLS]
      # `k` is the shared function's grouping-agnostic name; here they are units.
      names(od)[names(od) == "k"] <- "n_units"
      OD[[length(OD) + 1L]] <- cbind(
        data.frame(metric = m, unit_key = k, series = sname,
                   stringsAsFactors = FALSE), od)
    }
  }
}
# Fewer than three attributed units anywhere produces NO rows, which is a true
# statement about the site. The empty frame's columns come from the function's OWN shape.
.od_shape <- overdispersion(c(10, 20), c(100, 100), 15, OD_TRIM)[, FUNNEL_COLS]
names(.od_shape)[names(.od_shape) == "k"] <- "n_units"
od_tbl <- if (length(OD)) do.call(rbind, OD) else
  cbind(data.frame(metric = character(), unit_key = character(),
                   series = character(), stringsAsFactors = FALSE),
        .od_shape[0, ])
rownames(od_tbl) <- NULL
if (!nrow(od_tbl))
  cat(paste0("\n  NOTE no metric x unit key has the >= 3 units the over-dispersion\n",
             "       diagnostic needs, so overdispersion.csv ships with headers and\n",
             "       no rows. Expected at a single-unit site; the funnel and\n",
             "       caterpillar panels are unaffected.\n"))

stopifnot(
  "phi cannot be negative" = all(od_tbl$phi >= 0),
  "tau2 cannot be negative" = all(od_tbl$tau2 >= 0),
  "share_not_sampling must be a proportion" =
    all(od_tbl$share_not_sampling >= 0 & od_tbl$share_not_sampling <= 1)
)

cat("\n  Over-dispersion (measured, applied nowhere; trim =",
    sprintf("%.2f per end)\n", OD_TRIM))
for (i in seq_len(nrow(od_tbl))) {
  r <- od_tbl[i, ]
  cat(sprintf(paste0("  %-11s %-14s %-15s I=%2d  phi=%6.2f  winsorised=%6.2f",
                     "  (95%% crit %.2f)  tau=%5.2f pp  %s\n"),
              r$metric, r$unit_key, r$series, r$n_units, r$phi,
              r$phi_winsorised, r$phi_crit_95, sqrt(r$tau2),
              if (r$phi_winsorised > r$phi_crit_95)
                "OVER-DISPERSED" else "within binomial"))
}
cat("\n")


# ---- 9. Figures --------------------------------------------------------------
# Journal style: no title, no subtitle, boxed legend at the foot, type in INK.

# Deterministic funnel label placement, with a collision check. A label sits RIGHT of
# its point and flips LEFT only when it would run past the panel edge. Extents in data
# units: ~6.3 x 3.2in of panel, 2.5pt text ~1.55mm x 3.0mm.
PANEL_W_MM <- 160; PANEL_H_MM <- 81
CHAR_W_MM  <- 1.55; LINE_H_MM <- 3.0

# Leading words EVERY label shares can cost half the panel. Cut to a word boundary.
strip_common_prefix <- function(lab) {
  if (length(lab) < 2) return(lab)
  w <- strsplit(lab, " ", fixed = TRUE)
  n <- 0
  while (n < min(lengths(w)) - 1 && length(unique(vapply(w, `[`, "", n + 1))) == 1) n <- n + 1
  if (n == 0) return(lab)
  vapply(w, function(x) paste(x[-seq_len(n)], collapse = " "), "")
}

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

  # Bounded vertical de-collision: labels whose x extents overlap cannot share a y
  # band, so the upper one is pushed clear and the leader line keeps it attached. ALL
  # labels are compared -- a left-flipped and a right-drawn label can still overlap.
  sep <- h * 1.15
  i <- order(out$ty)
  for (a in seq_along(i)[-1]) {
    cur <- i[a]
    for (b in seq_len(a - 1)) {
      prev <- i[b]
      x_overlap <- !(out$x1[cur] < out$x0[prev] || out$x1[prev] < out$x0[cur])
      if (x_overlap && (out$ty[cur] - out$ty[prev]) < sep)
        out$ty[cur] <- out$ty[prev] + sep
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

# `n_overall` is a REQUIRED argument, never taken from scope -- two nested loops (lessons.md #26).
deff_note <- function(x, n_overall) sprintf(
  paste0("Intervals are 2.5th-97.5th percentiles of %s bootstrap resamples of ",
         "EPISODES (seed %s), which keeps an episode's events together because ",
         "they cluster within it. The interval says how much the observed rate ",
         "would wobble if the care process re-ran. The design effect of %.2f is ",
         "a property of the OVERALL rate and of no interval drawn here -- each ",
         "unit's interval is its own bootstrap percentile carrying its own ",
         "clustering, and no bar's width comes from it. It is reported ",
         "because a coordinating centre combining this site's overall rate with ",
         "another site's needs it: on the overall rate, %s events carry the ",
         "precision of roughly %s independent ones. IT DOES NOT ADJUST FOR ",
         "CASE MIX: a unit's rate reflects who it admits as much as how it ",
         "practises, and no association between unit and this metric is ",
         "modelled here."),
  format(B_RESAMP, big.mark = ","), format(SEED), x,
  format(n_overall, big.mark = ","),
  format(round(n_overall / x), big.mark = ","))

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
                   "infusion), and **any increase**, the two combined rather than ",
                   "a third kind of event."),
  indication = paste0("TWO SERIES PER UNIT: a rate **increase** and a **bolus**, ",
                      "plus **any escalation**, the two combined rather than a ",
                      "third kind of event."))

key_label <- function(k) if (k == "care_setting")
  paste0("care_setting -- the CLIF location_type where it exists, else ",
         "location_category, both mCIDE-controlled and so reproducible at ",
         "another site") else
  "site-local location_name (internal, not poolable)"

for (m in METRICS) {
  spec <- ANALYSES[[m]]
  overall_m <- unit_tbl$overall_pct[unit_tbl$metric == m &
                                  unit_tbl$series == spec$combined][1]
  deff_m <- unit_tbl$deff[unit_tbl$metric == m &
                          unit_tbl$series == spec$combined][1]
  n_unatt_m <- FUNNEL[[m]]$n_events[3]
  # Events per episode: the funnel's limits assume independence within a unit.
  comb_ev <- series_subset(ATT[[m]][ATT[[m]]$attributed, ], spec, spec$combined)
  ev_per_ep <- nrow(comb_ev) / length(unique(comb_ev$encounter_block))
  pal <- setNames(unname(SERIES_COLOURS)[seq_along(spec$series)], spec$series)

  for (k in UNIT_KEYS) {
    d <- plot_units[plot_units$metric == m & plot_units$unit_key == k, ]
    ord <- d[d$series == spec$combined, ]
    ord <- ord[order(ord$pct_paired), ]
    d$unit <- factor(d$unit, levels = ord$unit)
    d$series <- factor(as.character(d$series), levels = spec$series)
    n_units <- nrow(ord)
    dcomb <- d[d$series == spec$combined, ]

    # trim = TRUE, or format() pads to the widest and puts a space inside the bracket.
    axis_lab <- setNames(sprintf("%s  (%s)", ord$unit,
                                 format(ord$n_events, big.mark = ",", trim = TRUE)),
                         as.character(ord$unit))

    # Combined series opaque and large, components lighter and smaller.
    is_comb <- spec$series == spec$combined
    A_VAL <- setNames(ifelse(is_comb, 1.00, 0.55), spec$series)
    S_VAL <- setNames(ifelse(is_comb, 2.90, 1.70), spec$series)
    L_VAL <- setNames(ifelse(is_comb, 0.70, 0.40), spec$series)

    # Ordered-adjacent pairs with overlapping intervals, from the run that drew the panel.
    n_adj_overlap <- if (n_units > 1L)
      sum(ord$ci_lo[-1] <= ord$ci_hi[-n_units]) else 0L

    fn <- paste0(m, "_caterpillar_", k, ".png")
    register_caption(fn, sprintf("%s, by %s", spec$title, k),
      paste(METRIC_NOTE[[m]], SERIES_NOTE[[m]],
            paste0("The combined series is drawn opaque and larger; the ",
                   "components are lighter and smaller, so the panel reads as ",
                   "one estimate per unit with its decomposition beside it."),
            sprintf(paste0("The combined series necessarily lies between its ",
                           "components, at their event-weighted mean. Units are ",
                           "ordered by it and the dashed rule is its overall rate ",
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
                    overall_m, n_units,
                    length(unique(unit_tbl$unit[unit_tbl$metric == m &
                                                unit_tbl$unit_key == k])),
                    format(n_unatt_m, big.mark = ","), key_label(k)),
            deff_note(deff_m, sum(dcomb$n_events)),
            sprintf(paste0("WHAT THE INTERVALS ARE FOR. Ordering by rate ",
                           "invites reading the rank itself as the result, and ",
                           "the intervals are what say how much of that ",
                           "ordering is resolvable: %d of the %d ordered-",
                           "adjacent pairs here OVERLAP. Read the panel as the ",
                           "extremes being far apart with the middle a cluster, ",
                           "not as a step per unit. Overlap is an informal ",
                           "read and NOT a pairwise test -- there is no ",
                           "multiplicity control across %d units and none is ",
                           "intended. The interval covers SAMPLING ONLY: it is ",
                           "silent on whether the charting behind the numerator ",
                           "is accurate, which is a larger source of doubt than ",
                           "its width suggests. A narrow interval means the ",
                           "rate is precisely measured, not that it reflects ",
                           "how the unit practises rather than whom it admits."),
                    n_adj_overlap, max(n_units - 1L, 0L), n_units),
            paste0("This panel and the funnel beside it ask different ",
                   "questions. The funnel asks whether the spread is explained ",
                   "by how much each unit contributes; this one asks how much ",
                   "of the resulting ordering can be resolved.")))

    # alpha/size/linewidth are MAPPED, not a layer per series: position_dodge divides
    # the slot by the groups in a layer, so a lone series would dodge to another spot.
    dodge <- position_dodge(width = 0.62)
    pl <- house(
      ggplot(d, aes(pct_paired, unit, colour = series, group = series)) +
        geom_vline(xintercept = overall_m, colour = MUTED, linetype = "22",
                   linewidth = 0.4) +
        geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi, alpha = series,
                          linewidth = series),
                      orientation = "y", width = 0, position = dodge) +
        geom_point(aes(alpha = series, size = series), position = dodge) +
        scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
        # alpha and size encode emphasis, not a variable: keep them out of the legend.
        scale_alpha_manual(values = A_VAL, guide = "none") +
        scale_size_manual(values = S_VAL, guide = "none") +
        scale_linewidth_manual(values = L_VAL, guide = "none") +
        scale_y_discrete(labels = axis_lab) +
        scale_x_continuous(expand = expansion(mult = c(0.05, 0.06))) +
        labs(x = spec$x_lab, y = NULL))
    ggsave(file.path(dirs$phase, fn), pl,
           width = 7.5, height = max(3.0, 1.4 + 0.46 * n_units), dpi = 200)

    # COMBINED series only; limits are UNADJUSTED binomial, over-dispersion never applied.
    grid_n <- seq(max(1, min(dcomb$n_events) * 0.6), max(dcomb$n_events) * 1.15,
                  length.out = 200)
    se <- sqrt(overall_m * (100 - overall_m) / grid_n)
    lim <- rbind(
      data.frame(n = grid_n, lo = overall_m - Z_INNER * se,
                 hi = overall_m + Z_INNER * se, band = "95%"),
      data.frame(n = grid_n, lo = overall_m - Z_OUTER * se,
                 hi = overall_m + Z_OUTER * se, band = "99.8%"))

    # Quoted in this panel's caption; omitted where section 8b skipped the row.
    od_mk <- od_tbl[od_tbl$metric == m & od_tbl$unit_key == k &
                      od_tbl$series == spec$combined, ]
    od_note <- if (nrow(od_mk) == 1L)
      sprintf(paste0("WHAT THIS PANEL ANSWERS, as a number: the volume axis ",
                     "exists to ask whether the spread is just low-volume ",
                     "units being noisy, and it is not -- %.1f%% of the ",
                     "variation across units is NOT sampling error ",
                     "(between-unit SD %.2f pp against a mean within-unit ",
                     "binomial SD of %.2f pp). The spread is real, and that is ",
                     "this figure's contribution; the caterpillar beside it ",
                     "then shows how much of the ordering can be resolved. ",
                     "Over-dispersion is correspondingly large -- Winsorised ",
                     "phi %.1f against a 95%% reference of %.2f -- and is ",
                     "reported in overdispersion.csv, NOT corrected for. ",
                     "Inflating the envelope by sqrt(phi) would widen it past ",
                     "the whole observed range and past 0%%, which removes the ",
                     "figure's ability to say anything rather than making it ",
                     "cautious."),
              100 * od_mk$share_not_sampling, sqrt(od_mk$tau2),
              sqrt(od_mk$mean_binomial_var), od_mk$phi_winsorised,
              od_mk$phi_crit_95) else ""

    fn <- paste0(m, "_funnel_", k, ".png")
    register_caption(fn, sprintf("%s -- funnel against unit volume, by %s",
                                 spec$title, k),
      paste(METRIC_NOTE[[m]],
            sprintf(paste0("Each point is one unit, labelled beside it; the ",
                           "horizontal rule is the overall rate (%.1f%%) and the ",
                           "envelopes are 95%% and 99.8%% control limits, the ",
                           "two-sided p < 0.05 and p < 0.002 pair of ",
                           "Spiegelhalter (Qual Saf Health Care ",
                           "2005;14:347-51). THE LIMITS ARE UNADJUSTED ",
                           "BINOMIAL -- p0 +/- z*sqrt(p0(100-p0)/n) over a grid ",
                           "of n, centred on the overall rate. They therefore ",
                           "assume events are independent within a unit, which ",
                           "is NOT true here: this series averages %.2f events ",
                           "per episode and the habit is correlated inside one. ",
                           "No over-dispersion factor and no clustering ",
                           "correction is applied, and the proportions are not ",
                           "risk-adjusted, so a unit outside the envelope ",
                           "differs from the overall rate by more than sampling ",
                           "alone -- it does not follow that it differs in ",
                           "practice. Both quantities are MEASURED and shipped ",
                           "rather than applied: the design effect in ",
                           "unit_adherence.csv, over-dispersion in ",
                           "overdispersion.csv. The caterpillar beside this ",
                           "figure uses cluster-bootstrap intervals, so the two ",
                           "panels imply DIFFERENT widths for the same data by ",
                           "design -- an envelope is a null drawn as a shape, a ",
                           "caterpillar interval is uncertainty about one unit. ",
                           "Volume on the x axis is why this view sits beside ",
                           "the ordered one: a small unit sits far from the ",
                           "overall rate more easily, and the envelope makes ",
                           "that visible. Only the combined series is drawn -- ",
                           "overlapping funnels defeat the envelope reading. ",
                           "Labels are %s."),
                    overall_m, ev_per_ep, key_label(k)),
            od_note))

    lab <- funnel_labels(dcomb$n_events, dcomb$pct_paired,
                         strip_common_prefix(as.character(dcomb$unit)),
                         xr = range(grid_n),
                         yr = range(c(lim$lo, lim$hi, dcomb$pct_paired)))
    pl <- house(
      ggplot(lim, aes(n)) +
        geom_line(aes(y = lo, linetype = band), colour = MUTED, linewidth = 0.4) +
        geom_line(aes(y = hi, linetype = band), colour = MUTED, linewidth = 0.4) +
        geom_hline(yintercept = overall_m, colour = MUTED, linewidth = 0.4) +
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

  # Every event carries a year, so this view uses the whole denominator.
  yt <- year_tbl[year_tbl$metric == m, ]
  ypool <- yt$overall_pct[1]
  fn <- paste0(m, "_year.png")
  register_caption(fn, sprintf("%s, by calendar year", spec$title),
    paste(METRIC_NOTE[[m]],
          sprintf(paste0("One point per calendar year of the ventilation ",
                         "anchor, %d to %d, with the overall rate (%.1f%%) as a ",
                         "dashed rule. THE DENOMINATOR HERE IS ALL %s EVENTS, ",
                         "not only the ICU-attributed ones the unit figures use, ",
                         "because every event carries a year. Year counts follow ",
                         "the extract, not a study period: no date-based ",
                         "inclusion criterion is applied anywhere in this ",
                         "pipeline, so an end year may be partial."),
                  min(yt$year), max(yt$year), ypool,
                  format(sum(yt$n_events), big.mark = ",")),
          deff_note(yt$deff[1], sum(yt$n_events))))

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

# Different event sets, so the column cannot share a name: attributed events vs all.
unit_tbl$overall_pct_attributed <- unit_tbl$overall_pct
year_tbl$overall_pct_all_events <- year_tbl$overall_pct

cols <- c("metric", "unit_key", "unit", "care_setting", "series",
          "n_events", "n_episodes", "n_paired", "pct_paired", "ci_lo", "ci_hi",
          "overall_pct_attributed", "deff", "n_suppressed_small_cell")
write.csv(unit_tbl[order(unit_tbl$metric, unit_tbl$unit_key,
                         unit_tbl$unit, unit_tbl$series), cols],
          file.path(dirs$phase, "unit_adherence.csv"), row.names = FALSE)
write.csv(year_tbl[order(year_tbl$metric, year_tbl$year),
                   c("metric", "year", "n_events", "n_episodes", "n_paired",
                     "pct_paired", "ci_lo", "ci_hi", "overall_pct_all_events",
                     "deff")],
          file.path(dirs$phase, "year_adherence.csv"), row.names = FALSE)
write.csv(funnel_tbl, file.path(dirs$phase, "attribution_funnel.csv"),
          row.names = FALSE)
# Rounded on write, unlike deff: phi and tau2 are read against a threshold.
od_out <- od_tbl
num <- vapply(od_out, is.numeric, logical(1)) & names(od_out) != "n_units" &
       names(od_out) != "q_df"
od_out[num] <- lapply(od_out[num], round, 6)
write.csv(od_out, file.path(dirs$phase, "overdispersion.csv"), row.names = FALSE)

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
