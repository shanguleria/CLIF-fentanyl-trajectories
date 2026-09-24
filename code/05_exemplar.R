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
# Four panels, one shared x axis, stacked with the same machinery the risk
# tables use. Separate plots rather than facets because ggplot gives every
# `scales = "free_y"` continuous facet the same height, and these four do not
# want it: the infusion trace carries the most information and the bolus spikes
# the least.

X <- list(scale_x_continuous(limits = c(0, EXTENT_H),
                             breaks = seq(0, EXTENT_H, by = 12),
                             expand = expansion(mult = 0.01)))

# The traces STOP at extubation rather than running on at zero. A zero infusion
# line past extubation would assert the patient was ventilated and off fentanyl,
# when in fact they were not ventilated at all -- the same distinction the state
# definitions make between `no fentanyl` and `extubated`. The empty span to the
# right of the rule is the patient no longer being on the ventilator, and the x
# axis still runs to 72h so this panel is comparable with F2-F5.
END_HR <- if (!is.null(meta$extubation_hr) && !is.na(meta$extubation_hr)) {
  min(meta$extubation_hr, EXTENT_H)
} else EXTENT_H
series <- series[series$t_hr <= END_HR, ]

# Extubation, drawn on every panel so the eye carries it down the stack. After
# it the patient is no longer ventilated, which is why the traces stop rather
# than continue at zero.
vent_rule <- if (!is.null(meta$extubation_hr) && !is.na(meta$extubation_hr)) {
  list(geom_vline(xintercept = meta$extubation_hr, colour = MUTED,
                  linetype = "22", linewidth = 0.4))
} else list()

# An unlabelled rule is a puzzle. Labelled once, on the top panel only, so the
# eye carries it down the stack without the word repeating four times.
vent_label <- if (length(vent_rule)) {
  list(annotate("text", x = meta$extubation_hr, y = Inf, label = " extubated",
                hjust = 0, vjust = 1.4, colour = MUTED, size = 2.6))
} else list()

panel <- function(p, ylab, bottom = FALSE) {
  p <- house(p + X + vent_rule + labs(x = NULL, y = ylab)) +
    theme(panel.grid.major.x = element_blank(),
          plot.margin = margin(t = 2, r = 6, b = 2, l = 6))
  if (bottom) {
    p + labs(x = "Hours since first IMV episode")
  } else {
    p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
}

inf <- close_step(series[series$series == "infusion", ], END_HR)
p_inf <- panel(
  ggplot(inf, aes(t_hr, value)) +
    # Square risers, closed to zero at both ends, so "off drug" is the trace
    # sitting at baseline rather than the trace disappearing.
    geom_step(direction = "hv", colour = STATE_COLOURS[["continuous only"]],
              linewidth = 0.55) +
    vent_label +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))),
  EXEMPLAR_LABELS[["infusion"]])

bol <- series[series$series == "bolus", ]
p_bol <- panel(
  ggplot(bol, aes(t_hr, value)) +
    # Spikes from zero plus a reserved marker at the tip. Dose is read off
    # POSITION, never off marker size -- and near-simultaneous boluses are left
    # overplotted, because that staircase is real signal about how hard the
    # patient was being chased.
    geom_segment(aes(xend = t_hr, yend = 0),
                 colour = STATE_COLOURS[["bolus only"]], linewidth = 0.4) +
    geom_point(shape = 18, size = 1.9,
               colour = STATE_COLOURS[["bolus only"]]) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))),
  EXEMPLAR_LABELS[["bolus"]])

ordinal_panel <- function(name, breaks, limits, rule_at = NULL) {
  d <- break_on_gap(series[series$series == name, ], CAP_H)
  layers <- list()
  if (!is.null(rule_at)) {
    # MUTED and dashed, not GRIDLINE: drawn in the gridline colour it was
    # indistinguishable from the gridlines and carried no meaning at all.
    layers <- c(layers, list(geom_hline(yintercept = rule_at, colour = MUTED,
                                        linetype = "22", linewidth = 0.4)))
  }
  ggplot(d, aes(t_hr, value)) + layers +
    # Points plus a step, never a straight line between them: interpolating
    # ordinal scores asserts the patient passed through intermediate values at
    # times nobody recorded.
    geom_step(direction = "hv", colour = INK, linewidth = 0.4, na.rm = TRUE) +
    geom_point(colour = INK, size = 1.1, na.rm = TRUE) +
    # The FULL ordinal range, always. Letting the axis shrink to the
    # observed values would make "deeply sedated throughout" look like the
    # normal spread of the scale, and would stop two exemplars being
    # comparable at a glance.
    scale_y_continuous(breaks = breaks, limits = limits)
}

p_rass <- panel(ordinal_panel("rass", seq(-5, 4, by = 1), c(-5, 4), rule_at = 0),
                EXEMPLAR_LABELS[["rass"]])
p_nvps <- panel(ordinal_panel("nvps", seq(0, 10, by = 2), c(0, 10)),
                EXEMPLAR_LABELS[["nvps"]], bottom = TRUE)

if (SOURCE == "synthetic") {
  # The one sanctioned exception to journal style. A generated figure mistaken
  # for a result is the real hazard here, so it says so on its face.
  p_inf <- p_inf + labs(title = "SYNTHETIC -- GENERATED DATA, NOT A PATIENT") +
    theme(plot.title = element_text(colour = STATE_COLOURS[["died"]],
                                    face = "bold", size = 10))
}

g <- stack_panels(list(p_inf, p_bol, p_rass, p_nvps),
                  c(NA, 0.70, 1.30, 1.10))
ggsave(file.path(dirs$phase, "exemplar.png"), g,
       width = 7.5, height = 6.0, dpi = 200)
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
                   "overplotted. RASS and NVPS are ordinal and are drawn as ",
                   "points with a last-value-carried-forward step, broken ",
                   "across any gap longer than the %gh cap rather than ",
                   "asserting a score persisted. "), CAP_H),
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
