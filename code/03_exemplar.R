# ==============================================================================
# 03_exemplar.R  --  F1, one ventilation course in detail
#
# Purpose : A single patient's first 72h of ventilation at sub-hourly resolution -- fentanyl infusion rate, boluses, RASS and NVPS -- so that "continuous, then bolus, then off" reads as a person rather than as a marginal.
# Author  : Shan Guleria
# Created : 2026-09-24
# Inputs  : output/intermediate_phi/exemplar_series.parquet + exemplar_meta.json
#           (or a generated patient, see EXEMPLAR_SOURCE below)
# Outputs : output/final_no_phi/03_exemplar/ : the files listed in OWNED
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

# Renumbered 05 -> 03 on 2026-09-24, so the exemplar sits with the other
# granularity exhibit (the raster, in 04) and ahead of the aggregates -- Baker's
# Figure 1 is a proof-of-granularity exhibit BEFORE any aggregate. The old
# folder is retired explicitly: a site that pulls this update would otherwise
# keep a 05_exemplar/ nobody owns, full of files that look current. Built with
# paste0 so a name-level find-and-replace cannot reach into it (lessons.md #13).
RETIRED <- file.path("output", "final_no_phi", paste0("05_", "exemplar"),
                     c("exemplar.png", "provenance.json", "captions.md"))

if (SOURCE == "synthetic") {
  # Deliberately NOT under final_no_phi/: a generated figure must never be able
  # to reach the coordinating-centre upload set. Covered by .gitignore's bare
  # `output/`.
  dirs$phase <- file.path(dirname(dirs$out_final), "dev_synthetic")
  dir.create(dirs$phase, recursive = TRUE, showWarnings = FALSE)
  message("[03_exemplar] SYNTHETIC source -- writing to output/dev_synthetic/")
} else {
  dirs$phase <- phase_dir(dirs, "03_exemplar")
  manifest <- require_manifest(dirs, here())
  message(sprintf("[03_exemplar] site=%s  reading Phase 0 outputs from code %s",
                  config$site_name, manifest$code_version))
  n_cleared <- clear_owned_outputs(dirs, OWNED, retired = RETIRED)
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
# pattern as the 100-patient raster in 04_delivery_states.R.
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
# The assessments OVERLAY the fentanyl across the whole panel height, as Baker's
# insulin and dextrose overlay his glucose (SG, 2026-09-24). An earlier version
# reserved them a band in the upper panel to keep the traces apart; that is not
# what Baker does and not what was asked for. Crossing traces are the cost, and
# are why the series carry separate colours, separate marks and -- for the two
# ordinal ones -- their own axis spines.

FMAX <- max(250, ceiling(max(series$value[series$series %in% c("infusion", "bolus")],
                             na.rm = TRUE) / 50) * 50)
if (FMAX > 250) {
  message(sprintf("  fentanyl axis extended to %g mcg: this episode exceeds the ",
                  FMAX),
          "nominal 250 ceiling, and clipping a dose would be a lie")
}
TOP <- FMAX

# Ordinal value -> panel height, and back for the axis ticks. One pair of
# functions so a scale can never drift from the axis that reads it. Each ordinal
# series is stretched over the FULL height, which is what makes the right-hand
# axes readable at the same tick spacing as the left one.
to_band   <- function(v, lo, hi) (v - lo) / (hi - lo) * TOP
from_band <- function(y, lo, hi) lo + y / TOP * (hi - lo)
RASS_LIM <- c(-5, 4)
NVPS_LIM <- c(0, 10)

SERIES_COL <- c(infusion = STATE_COLOURS[["continuous only"]],
                bolus    = STATE_COLOURS[["bolus only"]],
                rass     = STATE_COLOURS[["extubated"]],
                nvps     = STATE_COLOURS[["discharged alive"]])
# Shape AND line style carry identity too, so the four series are never
# colour-alone -- and the fentanyl step stays the only solid line, so it reads as
# the primary quantity even where the ordinal traces cross it.
SERIES_SHP <- c(infusion = NA, bolus = 18, rass = 16, nvps = 17)
# Boluses are drawn markedly larger than the assessment points: they are the
# thing the figure is least able to show any other way -- an instantaneous event
# with no line to trace it -- so they have to carry on their own.
SERIES_SZ  <- c(infusion = 1.2, bolus = 3.0, rass = 1.2, nvps = 1.2)
SERIES_LTY <- c(infusion = "solid", bolus = "blank", rass = "13", nvps = "42")
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
                hjust = 0, vjust = 1.2, colour = INK, size = 2.6))
} else list()

inf  <- close_step(series[series$series == "infusion", ], END_HR)
bol  <- series[series$series == "bolus", ]
rass <- break_on_gap(series[series$series == "rass", ], CAP_H)
nvps <- break_on_gap(series[series$series == "nvps", ], CAP_H)
rass$y <- to_band(rass$value, RASS_LIM[1], RASS_LIM[2])
nvps$y <- to_band(nvps$value, NVPS_LIM[1], NVPS_LIM[2])

p <- ggplot(mapping = aes(t_hr)) +
  # Square risers, closed to zero at both ends: time off the drug is the trace
  # sitting at baseline, not the trace disappearing.
  geom_step(data = inf, aes(y = value, colour = "infusion"),
            direction = "hv", linewidth = 0.55, linetype = SERIES_LTY[["infusion"]]) +
  # Instantaneous, never carried forward. Dose is read off POSITION. Overlapping
  # administrations are left overplotted -- that staircase is real signal about
  # how hard the patient was being chased.
  # A PAPER halo behind each bolus, so that at this size a run of
  # near-simultaneous administrations still reads as separate marks rather than
  # one blob. The staircase is real signal and must stay countable; the halo is
  # the surface gap that keeps overlapping marks legible. Drawn as its own layer
  # with no colour mapping, so it stays out of the legend.
  geom_point(data = bol, aes(y = value), colour = PAPER,
             shape = SERIES_SHP[["bolus"]], size = SERIES_SZ[["bolus"]] + 1.2,
             show.legend = FALSE) +
  geom_point(data = bol, aes(y = value, colour = "bolus"),
             shape = SERIES_SHP[["bolus"]], size = SERIES_SZ[["bolus"]]) +
  # Points plus a step, never a straight line: interpolating an ordinal score
  # asserts the patient passed through intermediate values nobody recorded.
  geom_step(data = rass, aes(y = y, colour = "rass"), direction = "hv",
            linewidth = 0.45, linetype = SERIES_LTY[["rass"]], na.rm = TRUE) +
  geom_point(data = rass, aes(y = y, colour = "rass"),
             shape = SERIES_SHP[["rass"]], size = SERIES_SZ[["rass"]], na.rm = TRUE) +
  geom_step(data = nvps, aes(y = y, colour = "nvps"), direction = "hv",
            linewidth = 0.45, linetype = SERIES_LTY[["nvps"]], na.rm = TRUE) +
  geom_point(data = nvps, aes(y = y, colour = "nvps"),
             shape = SERIES_SHP[["nvps"]], size = SERIES_SZ[["nvps"]], na.rm = TRUE) +
  scale_colour_manual(values = SERIES_COL, labels = SERIES_LAB,
                      breaks = names(SERIES_LAB), name = NULL) +
  scale_x_continuous(limits = c(0, EXTENT_H), breaks = seq(0, EXTENT_H, by = 12),
                     expand = expansion(mult = 0.01)) +
  scale_y_continuous(
    limits = c(0, TOP), expand = expansion(mult = c(0.01, 0.02)),
    breaks = seq(0, FMAX, by = 50),
    name = "Fentanyl: mcg/hr infused, mcg per bolus",
    sec.axis = sec_axis(~ from_band(., RASS_LIM[1], RASS_LIM[2]), name = "RASS",
                        breaks = seq(RASS_LIM[1], RASS_LIM[2], by = 1))) +
  # ONE legend, its keys forced to the right glyphs. Mapping shape and linetype
  # as aesthetics instead produced three legends side by side -- ggplot merges
  # guides only when the scales agree key for key, and a shape scale carrying an
  # NA for the line-only series never will.
  guides(colour = guide_legend(
    nrow = 2, byrow = TRUE,
    override.aes = list(linetype = unname(SERIES_LTY[names(SERIES_LAB)]),
                        shape    = unname(SERIES_SHP[names(SERIES_LAB)]),
                        size     = unname(SERIES_SZ[names(SERIES_LAB)])))) +
  labs(x = "Hours since first IMV episode")

if (length(vent_rule)) p <- p + vent_rule + vent_label

# Three scales need three visible spines with tick marks, or a reader cannot
# tell which axis a trace belongs to. The SPINE and TICKS take the series
# colour -- they are marks, and a mark is what may carry identity -- while every
# label reads black at the shared axis-text size.
AXIS_PT <- 8.8   # theme_minimal(base_size = 12)'s axis.text size
p <- house(p) +
  theme(panel.grid.major.x = element_blank(),
        axis.text            = element_text(size = AXIS_PT),
        axis.line.y.left     = element_line(colour = INK, linewidth = 0.4),
        axis.ticks.y.left    = element_line(colour = INK, linewidth = 0.4),
        axis.line.y.right    = element_line(colour = SERIES_COL[["rass"]],
                                            linewidth = 0.4),
        axis.ticks.y.right   = element_line(colour = SERIES_COL[["rass"]],
                                            linewidth = 0.4),
        # theme_minimal() blanks ticks, so these are re-enabled by name; and they
        # must be long enough to read AS hatches against their own spine rather
        # than disappearing into it.
        axis.ticks.length    = grid::unit(4.5, "pt"))

# The third scale. Borrowed from a bare plot over the identical panel range, so
# its ticks land where NVPS actually is.
donor <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) + geom_blank() +
  scale_y_continuous(limits = c(0, TOP), expand = expansion(mult = c(0.01, 0.02)),
                     breaks = to_band(seq(NVPS_LIM[1], NVPS_LIM[2], by = 1),
                                      NVPS_LIM[1], NVPS_LIM[2]),
                     labels = seq(NVPS_LIM[1], NVPS_LIM[2], by = 1),
                     position = "right") +
  theme_minimal(base_size = 12) +
  theme(axis.title        = element_blank(),
        axis.text.y.right = element_text(colour = INK, size = AXIS_PT),
        axis.line.y.right = element_line(colour = SERIES_COL[["nvps"]],
                                         linewidth = 0.4),
        axis.ticks.y.right = element_line(colour = SERIES_COL[["nvps"]],
                                          linewidth = 0.4),
        axis.ticks.length = grid::unit(4.5, "pt"))

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
                   "infusion of equal size. RASS and NVPS are ordinal and ",
                   "share no scale with the fentanyl or with each other; each ",
                   "is drawn across the full panel height against its own ",
                   "right-hand axis, so the three traces overlay as in the ",
                   "source figure. Each is points with a ",
                   "last-value-carried-forward step, broken across any gap ",
                   "longer than the %gh cap rather than asserting a score ",
                   "persisted. "), CAP_H),
    if (!is.null(crit) && !is.na(meta$n_eligible))
      sprintf(paste0("SELECTION: this episode was drawn at random, seeded at ",
                     "%s, from the %s episodes meeting pre-specified criteria ",
                     "(%s). It was not chosen by inspection. "),
              # From the config, not the run metadata: draw_seed is not a
              # criterion and was removed from that list, which left the
              # caption reading "seeded at NULL".
              format(EX_SPEC$draw_seed, scientific = FALSE),
              format(meta$n_eligible, big.mark = ","),
              # vapply per element, never unlist(): unlist() on a list mixing
              # numbers and logicals coerces TRUE to 1, so a declared rule read
              # "require_survived_hospitalization = 1" in the caption.
              paste(vapply(names(crit), function(k) {
                v <- crit[[k]]
                sprintf("%s = %s", k,
                        if (is.logical(v)) (if (isTRUE(v)) "yes" else "no")
                        else format(v, scientific = FALSE))
              }, character(1)), collapse = "; "))
    else "",
    "DE-IDENTIFICATION: relative hours only, no dates and no identifiers; ",
    "asserted in code at export and again at draw."))

write_json(prov, file.path(dirs$phase, "provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)
cat("written: provenance.json\n")

write_captions(file.path(dirs$phase, "captions.md"), "03_exemplar.R",
               grep("\\.png$", OWNED$phase, value = TRUE), prov)
cat("written: captions.md\n")


# ---- 7. Provenance -----------------------------------------------------------

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/03_exemplar.R"),
    paste("Source :", SOURCE),
    "",
    capture.output(sessionInfo())),
  here("logs", "03_exemplar_sessioninfo.txt")
)
