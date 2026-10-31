# ==============================================================================
# 03_exemplar.R  --  F1, one ventilation course in detail
#
# Purpose : A single patient's first 72h of ventilation at sub-hourly resolution
#           -- fentanyl infusion rate, boluses, RASS and the site's pain score.
# Author  : Shan Guleria
# Created : 2026-09-24
# Inputs  : output/intermediate_phi/exemplar_series.parquet + exemplar_meta.json
# Outputs : output/final_no_phi/03_exemplar/ : the files listed in OWNED
#
# Modelled on Baker et al. Sci Rep 2020;10:10718 Figure 1. The selection rule is
# declared in protocol.json exemplar. EXEMPLAR_SOURCE=synthetic draws a generated
# patient into output/dev_synthetic/ instead; see code/utils/exemplar.R.
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
COV <- fromJSON(here("config", "protocol.json"), simplifyVector = FALSE)
set.seed(COV$model$seed)

WINDOW_H <- COV$cohort$window_hours
EXTENT_H <- COV$cohort$granular_extent_hours


# The selection rule is protocol: tracked protocol.json, never site-local config.json.
EX_SPEC <- COV$exemplar
stopifnot("protocol.json must declare an `exemplar` block" = !is.null(EX_SPEC),
          "the exemplar draw must be the declared one" =
            identical(EX_SPEC$draw, "uniform_random_from_eligible"))

# The pain instrument varies by site (protocol.json indication.pain_instruments), so
# instrument, display name and scale all come from config -- a hardcoded 0-10 axis
# would stretch a 0-8 trace over ticks no value can occupy, with no error anywhere.
PAIN_CATS <- names(COV$indication$pain_instruments)
PAIN_VARS <- Filter(function(v)
  identical(COV$time_varying[[v]]$source$table, "assessments") &&
    COV$time_varying[[v]]$source$category %in% PAIN_CATS,
  grep("^_", names(COV$time_varying), value = TRUE, invert = TRUE))
stopifnot("protocol.json declares no pain instrument as a covariate" =
            length(PAIN_VARS) >= 1)

# The step-breaking threshold is the cohort build's LOCF cap itself, read from the same key.
CAP_H <- as.numeric(COV$time_varying$rass$locf$cap_hours)
for (.v in PAIN_VARS)
  stopifnot("rass and every pain instrument must share one LOCF cap" =
              identical(CAP_H, as.numeric(COV$time_varying[[.v]]$locf$cap_hours)))

# Axis scales come from outlier_config analysis_unit_bounds, never restated here.
OUTLIERS <- fromJSON(here("config", "outlier_config.json"), simplifyVector = FALSE)

SOURCE <- Sys.getenv("EXEMPLAR_SOURCE", "real")
stopifnot("EXEMPLAR_SOURCE must be 'real' or 'synthetic'" =
            SOURCE %in% c("real", "synthetic"))


# ---- 3. Paths ----------------------------------------------------------------

dirs <- site_dirs()
prov <- provenance(config)

OWNED <- list(phase = c("exemplar.png", "provenance.json", "captions.md"))

# Retire the old 05_exemplar/ explicitly, or a site that pulls this keeps stale files
# that look current. Built with paste0 so a find-and-replace cannot reach the name.
RETIRED <- file.path("output", "final_no_phi", paste0("05_", "exemplar"),
                     c("exemplar.png", "provenance.json", "captions.md"))

if (SOURCE == "synthetic") {
  # Never under final_no_phi/: a generated figure must not reach the upload set.
  dirs$phase <- file.path(dirname(dirs$out_final), "dev_synthetic")
  dir.create(dirs$phase, recursive = TRUE, showWarnings = FALSE)
  message("[03_exemplar] SYNTHETIC source -- writing to output/dev_synthetic/")
} else {
  dirs$phase <- phase_dir(dirs, "03_exemplar")
  manifest <- require_manifest(dirs, here())
  message(sprintf("[03_exemplar] site=%s  reading cohort-build outputs from code %s",
                  config$site_name, manifest$code_version))
  n_cleared <- clear_owned_outputs(dirs, OWNED, retired = RETIRED)
  if (n_cleared) message(sprintf("  cleared %d output(s) from a previous run", n_cleared))
}


# ---- 4. Read -----------------------------------------------------------------

if (SOURCE == "synthetic") {
  ex <- synthetic_exemplar(seed = COV$model$seed, extent_h = EXTENT_H,
                           window_h = WINDOW_H)
  PAIN_VAR <- PAIN_VARS[1]
} else {
  f_series <- file.path(phi_dir(dirs, "01_cohort"), "exemplar_series.parquet")
  f_meta   <- file.path(phi_dir(dirs, "01_cohort"), "exemplar_meta.json")
  if (!file.exists(f_series)) {
    stop("no exemplar export found. Re-run code/01_build_cohort.py, which ",
         "selects the exemplar and writes ", basename(f_series),
         ". To design the figure without data, set EXEMPLAR_SOURCE=synthetic.",
         call. = FALSE)
  }
  series <- as.data.frame(read_parquet(f_series))
  # The cohort build exports the trace under its variable name; relabel to `pain` before the
  # factor, or the levels below turn the whole trace into NA at a non-NVPS site.
  .found <- intersect(PAIN_VARS, unique(as.character(series$series)))
  stopifnot("the exemplar export carries no declared pain instrument" =
              length(.found) >= 1)
  PAIN_VAR <- .found[1]
  series$series[as.character(series$series) %in% PAIN_VARS] <- "pain"
  series$series <- factor(series$series, levels = EXEMPLAR_SERIES)
  ex <- list(series = series, meta = fromJSON(f_meta, simplifyVector = TRUE))
  ex$meta$synthetic <- FALSE
}

series <- ex$series
meta   <- ex$meta

stopifnot(
  "the exemplar frame must carry no identifier" =
    !any(c("encounter_block", "patient_id", "hospitalization_id") %in% names(series)),
  "the exemplar clock must be relative hours, never a date" =
    is.numeric(series$t_hr),
  "every row must belong to a declared series" = !any(is.na(series$series)))

exemplar_diagnostics(series, meta, CAP_H)


# ---- 5. Figure ---------------------------------------------------------------
# One panel, three y scales: fentanyl left, RASS and the pain score on two right
# spines, with the ordinal traces overlaying the full panel height. The left axis
# carries both a rate and an amount, so its title names both quantities.

FMAX <- max(250, ceiling(max(series$value[series$series %in% c("infusion", "bolus")],
                             na.rm = TRUE) / 50) * 50)
if (FMAX > 250) {
  message(sprintf("  fentanyl axis extended to %g mcg: this episode exceeds the ",
                  FMAX),
          "nominal 250 ceiling, and clipping a dose would be a lie")
}
TOP <- FMAX

# Ordinal value -> panel height and back, one pair so scale and axis cannot drift.
to_band   <- function(v, lo, hi) (v - lo) / (hi - lo) * TOP
from_band <- function(y, lo, hi) lo + y / TOP * (hi - lo)
PAIN_CAT <- COV$time_varying[[PAIN_VAR]]$source$category
PAIN_LAB <- COV$indication$instrument_display_names[[PAIN_CAT]]
if (is.null(PAIN_LAB)) PAIN_LAB <- PAIN_CAT
RASS_LIM <- unlist(OUTLIERS$analysis_unit_bounds$patient_assessments$RASS)
PAIN_LIM <- unlist(OUTLIERS$analysis_unit_bounds$patient_assessments[[PAIN_CAT]])
stopifnot("the pain instrument has no analysis_unit_bounds entry to scale its axis" =
            length(PAIN_LIM) == 2 && PAIN_LIM[2] > PAIN_LIM[1])

SERIES_COL <- c(infusion = STATE_COLOURS[["continuous only"]],
                bolus    = STATE_COLOURS[["bolus only"]],
                rass     = STATE_COLOURS[["extubated"]],
                pain     = STATE_COLOURS[["discharged alive"]])
# Shape and line style carry identity too, so no series is colour-alone.
SERIES_SHP <- c(infusion = NA, bolus = 18, rass = 16, pain = 17)
# Boluses are drawn larger: an instantaneous event with no line to trace it.
SERIES_SZ  <- c(infusion = 1.2, bolus = 3.0, rass = 1.2, pain = 1.2)
SERIES_LTY <- c(infusion = "solid", bolus = "blank", rass = "13", pain = "42")
SERIES_LAB <- c(infusion = "Fentanyl infusion (mcg/hr)",
                bolus    = "Fentanyl bolus (mcg)",
                rass     = "RASS", pain = PAIN_LAB)

# Traces stop at extubation; a zero line past it would assert continued ventilation.
END_HR <- if (!is.null(meta$extubation_hr) && !is.na(meta$extubation_hr)) {
  min(meta$extubation_hr, EXTENT_H)
} else EXTENT_H
series <- series[series$t_hr <= END_HR, ]

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
pain <- break_on_gap(series[series$series == "pain", ], CAP_H)
rass$y <- to_band(rass$value, RASS_LIM[1], RASS_LIM[2])
pain$y <- to_band(pain$value, PAIN_LIM[1], PAIN_LIM[2])

p <- ggplot(mapping = aes(t_hr)) +
  # Closed to zero at both ends: time off the drug is the trace at baseline.
  geom_step(data = inf, aes(y = value, colour = "infusion"),
            direction = "hv", linewidth = 0.55, linetype = SERIES_LTY[["infusion"]]) +
  # Instantaneous, never carried forward; overlaps are left overplotted. The PAPER halo
  # is an unmapped layer, out of the legend, keeping near-simultaneous boluses
  # countable.
  geom_point(data = bol, aes(y = value), colour = PAPER,
             shape = SERIES_SHP[["bolus"]], size = SERIES_SZ[["bolus"]] + 1.2,
             show.legend = FALSE) +
  geom_point(data = bol, aes(y = value, colour = "bolus"),
             shape = SERIES_SHP[["bolus"]], size = SERIES_SZ[["bolus"]]) +
  # Points plus a step: interpolating an ordinal score asserts values nobody recorded.
  geom_step(data = rass, aes(y = y, colour = "rass"), direction = "hv",
            linewidth = 0.45, linetype = SERIES_LTY[["rass"]], na.rm = TRUE) +
  geom_point(data = rass, aes(y = y, colour = "rass"),
             shape = SERIES_SHP[["rass"]], size = SERIES_SZ[["rass"]], na.rm = TRUE) +
  geom_step(data = pain, aes(y = y, colour = "pain"), direction = "hv",
            linewidth = 0.45, linetype = SERIES_LTY[["pain"]], na.rm = TRUE) +
  geom_point(data = pain, aes(y = y, colour = "pain"),
             shape = SERIES_SHP[["pain"]], size = SERIES_SZ[["pain"]], na.rm = TRUE) +
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
  # One legend, keys forced by override.aes: mapping shape and linetype as aesthetics
  # gives three legends, because the scales cannot merge key for key.
  guides(colour = guide_legend(
    nrow = 2, byrow = TRUE,
    override.aes = list(linetype = unname(SERIES_LTY[names(SERIES_LAB)]),
                        shape    = unname(SERIES_SHP[names(SERIES_LAB)]),
                        size     = unname(SERIES_SZ[names(SERIES_LAB)])))) +
  labs(x = "Hours since first IMV episode")

if (length(vent_rule)) p <- p + vent_rule + vent_label

# Three scales need three visible spines; spine and ticks take the series colour.
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
        # theme_minimal() blanks ticks, so these are re-enabled by name.
        axis.ticks.length    = grid::unit(4.5, "pt"))

# The third scale, off a bare plot over the identical range so its ticks land right.
donor <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) + geom_blank() +
  scale_y_continuous(limits = c(0, TOP), expand = expansion(mult = c(0.01, 0.02)),
                     breaks = to_band(seq(PAIN_LIM[1], PAIN_LIM[2], by = 1),
                                      PAIN_LIM[1], PAIN_LIM[2]),
                     labels = seq(PAIN_LIM[1], PAIN_LIM[2], by = 1),
                     position = "right") +
  theme_minimal(base_size = 12) +
  theme(axis.title        = element_blank(),
        axis.text.y.right = element_text(colour = INK, size = AXIS_PT),
        axis.line.y.right = element_line(colour = SERIES_COL[["pain"]],
                                         linewidth = 0.4),
        axis.ticks.y.right = element_line(colour = SERIES_COL[["pain"]],
                                          linewidth = 0.4),
        axis.ticks.length = grid::unit(4.5, "pt"))

if (SOURCE == "synthetic") {
  # The one sanctioned exception to journal style: a generated figure says so.
  p <- p + labs(title = "SYNTHETIC -- GENERATED DATA, NOT A PATIENT") +
    theme(plot.title = element_text(colour = STATE_COLOURS[["died"]],
                                    face = "bold", size = 10))
}

g <- add_third_axis(p, donor, title = PAIN_LAB)
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
                   "%s over the first %dh of invasive ventilation, on a clock ",
                   "relative to the first IMV episode. "), PAIN_LAB, EXTENT_H),
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
                   "infusion of equal size. RASS and ", PAIN_LAB, " are ordinal and ",
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
              # From the config, not the run metadata: draw_seed is not a criterion.
              format(EX_SPEC$draw_seed, scientific = FALSE),
              format(meta$n_eligible, big.mark = ","),
              # vapply, never unlist(): unlist() coerces TRUE to 1 in a mixed list.
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
