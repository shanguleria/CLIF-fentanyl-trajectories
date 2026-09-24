# ==============================================================================
# 05_exemplar.R  --  F1, one ventilation course in detail
#
# Purpose : A single patient's first 72h of ventilation at sub-hourly resolution -- fentanyl infusion rate, boluses, RASS and NVPS -- so that "continuous, then bolus, then off" reads as a person rather than as a marginal.
# Author  : Shan Guleria
# Created : 2026-09-24
# Inputs  : output/intermediate_phi/exemplar_series.parquet + exemplar_meta.json
#           (or a generated patient, see EXEMPLAR_SOURCE below)
# Outputs : output/final_no_phi/05_exemplar/ : the files listed in OWNED
#
# Modelled on Baker et al. Sci Rep 2020;10:10718 Figure 1, with four deliberate
# departures -- see docs/references.md, the F1 entry, which is the binding spec:
#   * stacked facets on ONE shared x axis, never Baker's three overlaid y-axes
#   * a rate and an amount never share a numeric axis, so boluses get their own
#     panel rather than riding on the infusion scale
#   * ordinal scores are points plus a step, never straight-line interpolation,
#     which would assert the patient passed through intermediate values
#   * the step BREAKS across a gap longer than the LOCF cap, rather than
#     asserting a score persisted through hours nobody documented
# And the one thing Baker's figure is weakest on, which this must improve:
# "Representative ICU admission" states no selection rule. Ours is declared in
# covariates.json, applied by code, and written into the caption.
#
# DEVELOPMENT. Set EXEMPLAR_SOURCE=synthetic to draw a generated patient into
# output/dev_synthetic/ instead. The figure was designed that way on purpose --
# see code/utils/exemplar.R.
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
source(here("code", "utils", "exemplar.R"))


# ---- 2. Config ---------------------------------------------------------------

config <- fromJSON(here("config", "config.json"), simplifyVector = FALSE)
set.seed(config$model$seed)

WINDOW_H <- config$cohort$window_hours
EXTENT_H <- config$cohort$granular_extent_hours

COV <- fromJSON(here("config", "covariates.json"), simplifyVector = FALSE)

# The selection rule is PROTOCOL, not a local preference: two sites picking
# exemplars by different rules produce F1s that cannot be compared. It therefore
# lives in covariates.json, which is tracked, and never in config.json, which is
# gitignored and site-local. Read, never defaulted.
EX_SPEC <- COV$exemplar
stopifnot("covariates.json must declare an `exemplar` block" = !is.null(EX_SPEC),
          "the exemplar draw must be the declared one" =
            identical(EX_SPEC$draw, "uniform_random_from_eligible"))

# The step-breaking threshold is the SAME number as the LOCF cap applied to
# rass/nvps in Phase 0 -- read from the same key, never restated here, or the
# panel would assert a persistence the analytic table does not.
CAP_H <- as.numeric(COV$time_varying$rass$locf$cap_hours)
stopifnot("rass and nvps must share one LOCF cap" =
            identical(CAP_H, as.numeric(COV$time_varying$nvps$locf$cap_hours)))

SOURCE <- Sys.getenv("EXEMPLAR_SOURCE", "real")
stopifnot("EXEMPLAR_SOURCE must be 'real' or 'synthetic'" =
            SOURCE %in% c("real", "synthetic"))


# ---- 3. Paths ----------------------------------------------------------------

dirs <- site_dirs()
prov <- provenance(config)

OWNED <- list(phase = c("exemplar.png", "provenance.json", "captions.md"))

if (SOURCE == "synthetic") {
  # Deliberately NOT under final_no_phi/: a generated figure must never be able
  # to reach the coordinating-centre upload set. Covered by .gitignore's bare
  # `output/`.
  dirs$phase <- file.path(dirname(dirs$out_final), "dev_synthetic")
  dir.create(dirs$phase, recursive = TRUE, showWarnings = FALSE)
  message("[05_exemplar] SYNTHETIC source -- writing to output/dev_synthetic/")
} else {
  dirs$phase <- phase_dir(dirs, "05_exemplar")
  manifest <- require_manifest(dirs, here())
  message(sprintf("[05_exemplar] site=%s  reading Phase 0 outputs from code %s",
                  config$site_name, manifest$code_version))
  n_cleared <- clear_owned_outputs(dirs, OWNED, retired = character(0))
  if (n_cleared) message(sprintf("  cleared %d output(s) from a previous run", n_cleared))
}


# ---- 4. Read -----------------------------------------------------------------

if (SOURCE == "synthetic") {
  ex <- synthetic_exemplar(seed = config$model$seed, extent_h = EXTENT_H,
                           window_h = WINDOW_H)
} else {
  f_series <- file.path(dirs$out_phi, "exemplar_series.parquet")
  f_meta   <- file.path(dirs$out_phi, "exemplar_meta.json")
  if (!file.exists(f_series)) {
    stop("no exemplar export found. Re-run code/01_build_cohort.py, which ",
         "selects the exemplar and writes ", basename(f_series),
         ". To design the figure without data, set EXEMPLAR_SOURCE=synthetic.",
         call. = FALSE)
  }
  series <- as.data.frame(read_parquet(f_series))
  series$series <- factor(series$series, levels = EXEMPLAR_SERIES)
  ex <- list(series = series, meta = fromJSON(f_meta, simplifyVector = TRUE))
  ex$meta$synthetic <- FALSE
}

series <- ex$series
meta   <- ex$meta

# The de-identification is asserted here as well as at export, because this is
# the script that turns the series into something a person looks at. Same
# pattern as the 100-patient raster in 03_delivery_states.R.
stopifnot(
  "the exemplar frame must carry no identifier" =
    !any(c("encounter_block", "patient_id", "hospitalization_id") %in% names(series)),
  "the exemplar clock must be relative hours, never a date" =
    is.numeric(series$t_hr),
  "every row must belong to a declared series" = !any(is.na(series$series)))

exemplar_diagnostics(series, meta, CAP_H)


# ---- 5. Figure ---------------------------------------------------------------
# ONE panel, three y scales -- fentanyl on the left, RASS and NVPS on two right
# spines -- mirroring Baker et al. Figure 1 (glucose left, insulin and D50
# right). SG, 2026-09-24, choosing compaction over the four stacked panels this
# script drew first. See docs/references.md, the F1 entry, which records what
# that trades away.
#
# The assessments are given a RESERVED BAND in the upper part of the panel
# rather than being stretched across its whole height. Baker's insulin and D50
# separated from glucose by accident of magnitude; ours would not -- an ordinal
# score mapped over 0-250 mcg would run straight through the infusion trace.
# The band is a layout choice; the right-hand ticks still read true values.

FMAX <- max(250, ceiling(max(series$value[series$series %in% c("infusion", "bolus")],
                             na.rm = TRUE) / 50) * 50)
if (FMAX > 250) {
  message(sprintf("  fentanyl axis extended to %g mcg: this episode exceeds the "
                  , FMAX),
          "nominal 250 ceiling, and clipping a dose would be a lie")
}
BAND <- c(0.58, 1.00) * (FMAX / 0.50)   # where the ordinal series live
TOP  <- FMAX / 0.50

# Ordinal value -> panel height, and back for the axis ticks. One pair of
# functions per series so a scale can never drift from the axis that reads it.
to_band   <- function(v, lo, hi) BAND[1] + (v - lo) / (hi - lo) * diff(BAND)
from_band <- function(y, lo, hi) lo + (y - BAND[1]) / diff(BAND) * (hi - lo)
RASS_LIM <- c(-5, 4)
NVPS_LIM <- c(0, 10)

SERIES_COL <- c(infusion = STATE_COLOURS[["continuous only"]],
                bolus    = STATE_COLOURS[["bolus only"]],
                rass     = STATE_COLOURS[["extubated"]],
                nvps     = STATE_COLOURS[["discharged alive"]])
# Shape carries identity too, so the four series are never colour-alone.
SERIES_SHP <- c(infusion = NA, bolus = 18, rass = 16, nvps = 17)
SERIES_LAB <- c(infusion = "Fentanyl infusion (mcg/hr)",
                bolus    = "Fentanyl bolus (mcg)",
                rass     = "RASS", nvps = "NVPS")

# The traces STOP at extubation rather than running on at zero. A zero infusion
# line past extubation would assert the patient was ventilated and off fentanyl,
# when they were not ventilated at all.
END_HR <- if (!is.null(meta$extubation_hr) && !is.na(meta$extubation_hr)) {
  min(meta$extubation_hr, EXTENT_H)
} else EXTENT_H
series <- series[series$t_hr <= END_HR, ]

# Extubation. After it the patient is no longer ventilated, which is why the
# traces stop rather than continue at zero. Labelled once -- an unlabelled rule
# is a puzzle, and a reference-line label is data, not narration.
vent_rule <- if (!is.null(meta$extubation_hr) && !is.na(meta$extubation_hr)) {
  list(geom_vline(xintercept = meta$extubation_hr, colour = MUTED,
                  linetype = "22", linewidth = 0.4))
} else list()
vent_label <- if (length(vent_rule)) {
  list(annotate("text", x = meta$extubation_hr, y = TOP, label = " extubated",
                hjust = 0, vjust = 1.2, colour = MUTED, size = 2.6))
} else list()

inf  <- close_step(series[series$series == "infusion", ], END_HR)
bol  <- series[series$series == "bolus", ]
rass <- break_on_gap(series[series$series == "rass", ], CAP_H)
nvps <- break_on_gap(series[series$series == "nvps", ], CAP_H)
rass$y <- to_band(rass$value, RASS_LIM[1], RASS_LIM[2])
nvps$y <- to_band(nvps$value, NVPS_LIM[1], NVPS_LIM[2])

p <- ggplot(mapping = aes(t_hr)) +
  # A rule under the assessment band, so the reader sees that the upper region
  # is a different quantity rather than more fentanyl.
  geom_hline(yintercept = BAND[1] - 0.02 * TOP, colour = GRIDLINE, linewidth = 0.4) +
  # Square risers, closed to zero at both ends: time off the drug is the trace
  # sitting at baseline, not the trace disappearing.
  geom_step(data = inf, aes(y = value, colour = "infusion"),
            direction = "hv", linewidth = 0.55) +
  # Instantaneous, never carried forward. Dose is read off POSITION. Overlapping
  # administrations are left overplotted -- that staircase is real signal about
  # how hard the patient was being chased.
  geom_point(data = bol, aes(y = value, colour = "bolus", shape = "bolus"),
             size = 1.9) +
  # Points plus a step, never a straight line: interpolating an ordinal score
  # asserts the patient passed through intermediate values nobody recorded.
  geom_step(data = rass, aes(y = y, colour = "rass"), direction = "hv",
            linewidth = 0.4, na.rm = TRUE) +
  geom_point(data = rass, aes(y = y, colour = "rass", shape = "rass"),
             size = 1.1, na.rm = TRUE) +
  geom_step(data = nvps, aes(y = y, colour = "nvps"), direction = "hv",
            linewidth = 0.4, linetype = "22", na.rm = TRUE) +
  geom_point(data = nvps, aes(y = y, colour = "nvps", shape = "nvps"),
             size = 1.1, na.rm = TRUE) +
  scale_colour_manual(values = SERIES_COL, labels = SERIES_LAB,
                      breaks = names(SERIES_LAB), name = NULL) +
  scale_shape_manual(values = SERIES_SHP[!is.na(SERIES_SHP)],
                     labels = SERIES_LAB[!is.na(SERIES_SHP)],
                     breaks = names(SERIES_SHP)[!is.na(SERIES_SHP)], name = NULL) +
  scale_x_continuous(limits = c(0, EXTENT_H), breaks = seq(0, EXTENT_H, by = 12),
                     expand = expansion(mult = 0.01)) +
  scale_y_continuous(
    limits = c(0, TOP), expand = expansion(mult = c(0.01, 0.02)),
    breaks = seq(0, FMAX, by = 50),
    name = "Fentanyl: mcg/hr infused, mcg per bolus",
    sec.axis = sec_axis(~ from_band(., RASS_LIM[1], RASS_LIM[2]), name = "RASS",
                        breaks = seq(RASS_LIM[1], RASS_LIM[2], by = 1))) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, order = 1),
         shape = "none") +
  labs(x = "Hours since first IMV episode")

if (length(vent_rule)) p <- p + vent_rule + vent_label

p <- house(p) + theme(panel.grid.major.x = element_blank())

# The third scale. Borrowed from a bare plot over the identical panel range, so
# its ticks land where NVPS actually is.
donor <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) + geom_blank() +
  scale_y_continuous(limits = c(0, TOP), expand = expansion(mult = c(0.01, 0.02)),
                     breaks = to_band(seq(NVPS_LIM[1], NVPS_LIM[2], by = 2),
                                      NVPS_LIM[1], NVPS_LIM[2]),
                     labels = seq(NVPS_LIM[1], NVPS_LIM[2], by = 2),
                     position = "right") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y.right = element_text(colour = MUTED),
        axis.title = element_blank())

if (SOURCE == "synthetic") {
  # The one sanctioned exception to journal style. A generated figure mistaken
  # for a result is the real hazard here, so it says so on its face.
  p <- p + labs(title = "SYNTHETIC -- GENERATED DATA, NOT A PATIENT") +
    theme(plot.title = element_text(colour = STATE_COLOURS[["died"]],
                                    face = "bold", size = 10))
}

g <- add_third_axis(p, donor, title = "NVPS")
ggsave(file.path(dirs$phase, "exemplar.png"), g,
       width = 7.5, height = 4.6, dpi = 200)
cat(sprintf("\nwritten: exemplar.png (%s)\n", SOURCE))


# ---- 6. Write ----------------------------------------------------------------

crit <- meta$criteria
register_caption("exemplar.png",
  "One ventilation course, from intubation to extubation",
  paste0(
    if (isTRUE(meta$synthetic))
      "GENERATED DATA, NOT A PATIENT -- a development render. " else "",
    sprintf(paste0("Fentanyl infusion rate, bolus administrations, RASS and ",
                   "NVPS over the first %dh of invasive ventilation, on a clock ",
                   "relative to the first IMV episode. "), EXTENT_H),
    if (!is.null(meta$extubation_hr) && !is.na(meta$extubation_hr))
      sprintf("The patient was extubated at hour %.0f; the traces end there. ",
              meta$extubation_hr) else "",
    sprintf(paste0("The infusion is a step function closed to zero at both ",
                   "ends, so time off the drug is the trace at baseline rather ",
                   "than absent. Boluses are instantaneous and are not carried ",
                   "forward; near-simultaneous administrations are left ",
                   "overplotted. Infusion and bolus share the left axis, which ",
                   "names both quantities because they are not the same one: ",
                   "the infusion is a RATE in mcg/hr and each bolus an AMOUNT ",
                   "in mcg, so a bolus plotted at a given height is not an ",
                   "infusion of equal size. RASS and NVPS are ordinal, share no ",
                   "scale with the fentanyl or with each other, and are drawn ",
                   "in a reserved band in the upper panel against their own ",
                   "right-hand axes; each is points with a ",
                   "last-value-carried-forward step, broken across any gap ",
                   "longer than the %gh cap rather than asserting a score ",
                   "persisted. "), CAP_H),
    if (!is.null(crit) && !is.na(meta$n_eligible))
      sprintf(paste0("SELECTION: this episode was drawn at random, seeded at ",
                     "%s, from the %s episodes meeting pre-specified criteria ",
                     "(%s). It was not chosen by inspection. "),
              format(config$model$seed, scientific = FALSE),
              format(meta$n_eligible, big.mark = ","),
              paste(sprintf("%s = %s", names(crit), unlist(crit)),
                    collapse = "; "))
    else "",
    "DE-IDENTIFICATION: relative hours only, no dates and no identifiers; ",
    "asserted in code at export and again at draw."))

write_json(prov, file.path(dirs$phase, "provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)
cat("written: provenance.json\n")

write_captions(file.path(dirs$phase, "captions.md"), "05_exemplar.R",
               grep("\\.png$", OWNED$phase, value = TRUE), prov)
cat("written: captions.md\n")


# ---- 7. Provenance -----------------------------------------------------------

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/05_exemplar.R"),
    paste("Source :", SOURCE),
    "",
    capture.output(sessionInfo())),
  here("logs", "05_exemplar_sessioninfo.txt")
)
