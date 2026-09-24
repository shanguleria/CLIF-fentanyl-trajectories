# The single-patient exemplar (F1) -- shared shape, gap handling, and a synthetic
# patient to develop against.
#
# WHY A SYNTHETIC PATIENT EXISTS. F1 is the only figure in this project that
# draws one identifiable clinical course in detail. Designing it normally means
# rendering it and looking, which would mean looking at a real patient. So the
# figure is built and iterated against a GENERATED patient that exercises every
# visual feature the real one has -- step risers up and down, a true off-drug
# gap, a staircase of near-simultaneous boluses, ordinal scores with a gap past
# the LOCF cap, and extubation before the window ends. Only when the drawing is
# settled does the real series go through the same code.
#
# The synthetic patient is not a test fixture in the usual sense -- this repo has
# no R test harness (tests/ is eight Python files, and renv.lock pins no
# testthat) -- it is the development surface, and it keeps the figure code
# runnable at a site with no data at all.

# One long frame, whatever the source. `series` is a FACTOR with declared levels
# so the panel order cannot be re-sorted alphabetically underneath us -- the same
# failure that reordered every stacked area before 2026-09-24 (lessons.md #15).
EXEMPLAR_SERIES <- c("infusion", "bolus", "rass", "nvps")

EXEMPLAR_LABELS <- c(infusion = "Fentanyl infusion\n(mcg/hr)",
                     bolus    = "Fentanyl bolus\n(mcg)",
                     rass     = "RASS",
                     nvps     = "NVPS")

# A step asserts the last value persisted. Across a gap longer than the LOCF cap
# that assertion is false, so the step must BREAK rather than carry a score over
# hours nobody documented. Baker draws missingness as nothing at all; this is the
# one place F1 deliberately departs (docs/references.md, F1 entry).
#
# Implemented by inserting an NA immediately after the last observation before
# each over-long gap: ggplot breaks a line or step on NA.
break_on_gap <- function(d, cap_hours) {
  if (nrow(d) < 2) return(d)
  d <- d[order(d$t_hr), ]
  gap <- c(diff(d$t_hr), 0)
  brk <- which(gap > cap_hours)
  if (!length(brk)) return(d)
  filler <- d[brk, ]
  filler$t_hr  <- filler$t_hr + 1e-6
  filler$value <- NA_real_
  out <- rbind(d, filler)
  out[order(out$t_hr), ]
}

# A step function must be CLOSED at both ends, so that "off drug" reads as the
# trace sitting at baseline rather than as the trace stopping. The infusion grid
# is hourly and complete -- absence_means_zero, so a zero is a measurement, not a
# gap -- and only the right-hand end needs carrying to the edge.
close_step <- function(d, extent_h) {
  if (!nrow(d)) return(d)
  d <- d[order(d$t_hr), ]
  last <- d[nrow(d), ]
  if (last$t_hr < extent_h) {
    last$t_hr <- extent_h
    d <- rbind(d, last)
  }
  d
}

# A generated patient. Deterministic given the seed: the same fake person every
# run, so a change in the figure is a change in the CODE and never in the data.
#
# Every feature below is deliberate, and each is something the real figure must
# survive. Change them only to make the development surface harder, never easier.
synthetic_exemplar <- function(seed = 1L, extent_h = 72, window_h = 4) {
  set.seed(seed)

  # Infusion: started late, escalated, weaned, stopped dead for 8h, restarted,
  # weaned to zero before extubation. Hourly, complete, zeros are real.
  rate <- rep(0, extent_h)
  rate[5:16]  <- 50     # started at hour 4
  rate[17:24] <- 100    # escalated
  rate[25:30] <- 75
  rate[31:38] <- 0      # a GENUINE off-drug gap, 8h at baseline
  rate[39:48] <- 50     # restarted
  rate[49:56] <- 25     # weaned
  extubation_hr <- 56
  rate[57:extent_h] <- 0

  infusion <- data.frame(t_hr = seq(0, extent_h - 1), series = "infusion",
                         value = rate)

  # Boluses, including a STAIRCASE of three inside forty minutes at the point of
  # escalation. That overplotting is real signal about how hard a patient was
  # being chased and must not be jittered away.
  bolus <- data.frame(
    t_hr   = c(4.2, 16.1, 16.4, 16.7, 28.5, 39.0, 50.2),
    series = "bolus",
    value  = c(50, 100, 100, 50, 75, 50, 25))

  # RASS: deep early, lightening as fentanyl is weaned, with a 12h documentation
  # gap that exceeds the 4h cap and must break the step.
  rass <- data.frame(
    t_hr   = c(1, 4, 8, 12, 16, 20, 24, 28, 42, 46, 50, 54),
    series = "rass",
    value  = c(-3, -4, -4, -5, -5, -4, -4, -3, -2, -2, -1, 0))

  nvps <- data.frame(
    t_hr   = c(2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 44, 48, 52, 55),
    series = "nvps",
    value  = c(2, 1, 0, 0, 3, 5, 4, 2, 6, 7, 3, 2, 4, 3))

  series <- rbind(infusion, bolus, rass, nvps)
  series$series <- factor(series$series, levels = EXEMPLAR_SERIES)

  list(series = series,
       meta = list(extent_h = extent_h, window_h = window_h,
                   extubation_hr = extubation_hr, synthetic = TRUE,
                   n_eligible = NA_integer_, criteria = NULL))
}

# Everything I am allowed to know about the REAL exemplar. Counts and spans, no
# values, no identifiers -- the figure itself is reviewed by a human, not here.
exemplar_diagnostics <- function(series, meta, cap_hours) {
  stopifnot("exemplar series must be a factor with the declared levels" =
              is.factor(series$series),
            "exemplar series carries an unexpected column" =
              setequal(names(series), c("t_hr", "series", "value")),
            "the exemplar clock must be relative and start at or after 0" =
              is.numeric(series$t_hr) && min(series$t_hr) >= 0,
            "the exemplar clock must not exceed the window" =
              max(series$t_hr) <= meta$extent_h)

  cat(sprintf("\nExemplar series (%s)\n",
              if (isTRUE(meta$synthetic)) "SYNTHETIC -- generated, not a patient"
              else "real, de-identified"))
  cat(sprintf("  span                 %.1f - %.1f h (window %d h)\n",
              min(series$t_hr), max(series$t_hr), meta$extent_h))
  if (!is.na(meta$extubation_hr)) {
    cat(sprintf("  extubated at         %.1f h\n", meta$extubation_hr))
  }
  for (s in EXEMPLAR_SERIES) {
    d <- series[series$series == s & !is.na(series$value), ]
    if (!nrow(d)) { cat(sprintf("  %-20s none\n", s)); next }
    gaps <- if (nrow(d) > 1) diff(sort(d$t_hr)) else 0
    cat(sprintf("  %-20s %3d points, longest gap %.1fh%s\n",
                s, nrow(d), max(gaps),
                if (s == "infusion")
                  sprintf(", %d step changes, %.0fh at zero",
                          sum(diff(d$value) != 0), sum(d$value == 0))
                else if (s == "bolus")
                  sprintf(", closest pair %.2fh apart", min(gaps))
                else
                  sprintf(", %d break(s) past the %gh cap", sum(gaps > cap_hours),
                          cap_hours)))
  }
  invisible(NULL)
}
