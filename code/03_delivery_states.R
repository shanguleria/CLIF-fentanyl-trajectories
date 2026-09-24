# ==============================================================================
# 03_delivery_states.R  --  fentanyl delivery states: prevalence and transitions
#
# Purpose : The seven delivery states over the ventilation course -- prevalence per window, transitions window to window, and the figures that draw them. Run twice: delivery ROUTE and dose INTENSITY band.
# Author  : Shan Guleria
# Created : 2026-09-05
# Split   : 2026-09-23, out of 02_descriptive_trajectory.R -- one subject per
#           script, so the figure work has somewhere to go.
# Inputs  : output/intermediate_phi/trajectory_long.parquet, time_to_event.parquet
# Outputs : output/final_no_phi/03_states/ : the CSVs and figures listed in OWNED
# ==============================================================================

# Run this in a FRESH R session (RStudio: Cmd+Shift+F10).
# Do not use rm(list = ls()) -- it does not unload packages or reset options,
# so it only gives the appearance of a clean slate.

# ---- 1. Packages -------------------------------------------------------------

pkgs <- c("here", "jsonlite", "ggplot2", "arrow", "ggalluvial")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

source(here("code", "utils", "paths.R"))
source(here("code", "utils", "states.R"))
source(here("code", "utils", "figures.R"))


# ---- 2. Config (never setwd(); here() anchors to the .Rproj) -----------------

config <- fromJSON(here("config", "config.json"), simplifyVector = FALSE)
set.seed(config$model$seed)

WINDOW_H  <- config$cohort$window_hours
EXTENT_H  <- config$cohort$granular_extent_hours
PANELS    <- unlist(config$cohort$balanced_panel_hours)
LANDMARK  <- config$cohort$landmark_hours

COV       <- fromJSON(here("config", "covariates.json"), simplifyVector = FALSE)
# Dose-band definition. Read from the config, never defaulted here: the cuts are
# federation-critical and a local default is how a site keeps cutting at the old
# values after the consortium moves them.
DS_SPEC   <- COV$exposure$dose_states
DOSE_CUTS <- as.numeric(unlist(DS_SPEC$cuts))
DOSE_LAB  <- unlist(DS_SPEC$labels)
stopifnot("covariates.json exposure.dose_states must declare cuts and labels" =
            length(DOSE_CUTS) > 0 && length(DOSE_LAB) == length(DOSE_CUTS) + 2,
          "dose states must be defined on window_mcg" =
            identical(DS_SPEC$variable, "window_mcg"))
SED_UNITS <- COV$exposure$sedatives$units

# Read the drug list FROM the config rather than restating it. A hardcoded list
# here is how a declared sedative becomes a no-op: adding dexmedetomidine to
# covariates.json would have changed nothing.
SED_COLS <- unlist(COV$exposure$sedatives$columns)
DRUGS <- c(fentanyl = "total_dose",
           setNames(SED_COLS, sub("_dose$", "", SED_COLS)))
UNITS <- c(fentanyl = COV$exposure$units,
           setNames(vapply(SED_COLS, function(c) SED_UNITS[[c]], character(1)),
                    sub("_dose$", "", SED_COLS)))
SEDATIVES <- setdiff(names(DRUGS), "fentanyl")

MIN_CELL <- config$reporting$small_cell_min_den

# Race is collapsed for DISPLAY only; the full CLIF granularity still reaches
# pooling_categorical.csv. Map lives in covariates.json, not here.
RACE_COLLAPSE <- COV$time_invariant$race$reporting_collapse




# ---- 3. Paths and provenance -------------------------------------------------
# One site, one output tree. site_dirs() creates them and labels the PHI ones.

dirs <- site_dirs()
dirs$phase <- phase_dir(dirs, "03_states")   # shareable outputs, subdivided by script
prov <- provenance(config)

message(sprintf("[03_delivery_states] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))


# ---- 4. Guards ---------------------------------------------------------------
# Phase 0 (Python) wrote what this script reads. Refuse to analyse tables that a
# different code version or a different config produced.

manifest <- require_manifest(dirs, here())
message(sprintf("  reading Phase 0 outputs from code %s, generated %s",
                manifest$code_version, manifest$generated))

OWNED <- list(phase = c(
  "state_prevalence.csv", "state_transitions.csv",
  "dose_state_prevalence.csv", "dose_state_transitions.csv",
  "state_prevalence.png", "state_alluvial.png",
  "dose_state_prevalence.png", "dose_state_alluvial.png",
  "provenance.json"))

# This script is new on 2026-09-23; it has no retired stems of its own. The
# copies of these outputs that 02 used to write into 02_descriptive/ are retired
# by that script, which still owns that folder.
RETIRED <- character(0)
n_cleared <- clear_owned_outputs(dirs, OWNED, retired = RETIRED)
if (n_cleared) message(sprintf("  cleared %d output(s) from a previous run", n_cleared))


# ---- 5. Read -----------------------------------------------------------------

long <- as.data.frame(read_parquet(
  file.path(dirs$out_phi, "trajectory_long.parquet"),
  col_select = c("encounter_block", "patient_id", "window_idx", "window_start_hr",
                 "alive_admitted", "imv_status", "inf_dose", "bolus_dose",
                 "total_dose", "window_mcg", "propofol_dose", "midazolam_dose",
                 "dexmedetomidine_dose",
                 "age", "sex", "race", "cci", "bmi_admission", "weight_kg",
                 "sofa_total", "nee", "oxygenation", "lactate",
                 "first_imv_episode_hours", "n_imv_episodes",
                 "died")))          # terminal states, section 12b below

tte <- as.data.frame(read_parquet(
  file.path(dirs$out_phi, "time_to_event.parquet"),
  col_select = c("encounter_block", "landmark_eligible")))

n_windows_expected <- EXTENT_H / WINDOW_H
stopifnot(
  "window count on disk disagrees with the config grid" =
    length(unique(long$window_idx)) == n_windows_expected,
  "the grid does not span granular_extent_hours" =
    max(long$window_start_hr) + WINDOW_H == EXTENT_H
)

# The Phase 1 denominator is VENTILATED, not alive_admitted. The two differ by
# more than half the cohort by 72h.
long$ventilated <- !is.na(long$imv_status) & long$imv_status == 1
stopifnot(
  "a window cannot be ventilated without being alive and admitted" =
    all(long$ventilated <= long$alive_admitted)
)

missing_cols <- setdiff(unname(DRUGS), names(long))
if (length(missing_cols)) {
  stop("covariates.json declares dose columns that trajectory_long does not carry: ",
       paste(missing_cols, collapse = ", "),
       ". Re-run code/01_build_cohort.py rather than dropping them silently.",
       call. = FALSE)
}

anchor_n <- length(unique(long$encounter_block))
elig     <- tte$encounter_block[tte$landmark_eligible]

cat(sprintf("\nCohort: %s ventilation episodes; %s landmark-eligible at T=%dh\n",
            format(anchor_n, big.mark = ","), format(length(elig), big.mark = ","),
            LANDMARK))


# ---- 12b. Delivery states and transitions ------------------------------------
# A different description of the same 72h: not "what shape is the dose curve"
# but "what state is the episode in, and where does it go next".
#
# Seven states. Four describe HOW fentanyl was delivered while the patient was
# ventilated; three are terminal. Because they are defined by which route was
# used rather than by a threshold on a continuum, there are no cut points to
# defend -- which is the whole point after the Phase 3/4 finding that dose level
# is continuous and its classes were a discretisation of it.
#
# WHOLE ANALYTIC COHORT, not the landmark set. The landmark conditions on being
# ventilated at T, so inside [0, T] nobody dies, is discharged or is
# permanently extubated -- the three terminal states would be empty by
# construction, and the liberation pathway is exactly what makes this figure
# worth drawing. 8,169 of 14,897 episodes leave before 72h, median at 24h.

# TWO definitions, same machinery. The first describes HOW fentanyl was
# delivered (route); the second HOW MUCH (intensity band on window_mcg). Lyons
# et al. fit two multistate models on one cohort for the same reason -- AKI
# stage, then AKI x IMV -- and called it triangulation. Everything below is
# written once and run twice, so the two can never drift apart in their handling.

state_summary <- function(d, tag) {
  lv <- levels(d$state)
  prevalence <- do.call(rbind, lapply(sort(unique(d$window_idx)), function(w) {
    x <- d[d$window_idx == w, ]
    tb <- table(x$state)[lv]
    data.frame(window_idx = w, window_start_hr = x$window_start_hr[1],
               state = lv, n = as.integer(tb),
               pct = round(100 * as.numeric(tb) / nrow(x), 2))
  }))
  tp <- transition_pairs(d)
  # Terminal states must absorb; transition_pairs() drops rows starting in one,
  # so a terminal state appearing as an ORIGIN means the definition is wrong.
  stopifnot("an absorbing state must not originate a transition" =
              !any(tp$state %in% STATE_ABSORBING))

  cat(sprintf("\n%s states, %% of the analytic cohort\n", tag))
  pv <- reshape(prevalence[, c("window_start_hr", "state", "pct")],
                idvar = "window_start_hr", timevar = "state", direction = "wide")
  names(pv) <- sub("^pct\\.", "", names(pv))
  print(pv[pv$window_start_hr %in% c(0, 24, 48, 68), ], row.names = FALSE)

  tmx <- transition_matrix(tp)
  cat(sprintf("\n%s transition matrix (row = t, col = t+1, row %%)\n", tag))
  print(tmx, row.names = FALSE)
  cat("  terminal states verified absorbing\n")

  list(prevalence = prevalence, matrix = tmx, data = d)
}

long  <- derive_states(long)                              # route definition
route <- state_summary(long, "Delivery")

# The intensity definition. derive_dose_states() overwrites `state`, so it runs
# on a COPY -- `long` keeps the route states for the rest of the script.
dose <- state_summary(
  derive_dose_states(long, DOSE_CUTS, DOSE_LAB),
  sprintf("Intensity (0 / <=%s / <=%s / >%s mcg per %dh window)",
          DOSE_CUTS[1], DOSE_CUTS[2], DOSE_CUTS[2], WINDOW_H))

state_prevalence      <- route$prevalence
transition_matrix_tbl <- route$matrix


# ---- 13. Figures -------------------------------------------------------------

# --- Delivery-state figures ---------------------------------------------------
# Alluvial: every episode is a ribbon, its width the number of episodes moving
# between states. Ribbons entering a terminal state never leave it, so
# liberation and death are visible as the flow drains out of the fentanyl states.

STATE_PAL <- c(
  "no fentanyl"        = "#d8d6cf",
  "continuous only"    = "#14427e",
  "bolus only"         = "#eb6834",
  "continuous + bolus" = "#4a8bd8",
  "extubated"          = "#7fb069",
  "discharged alive"   = "#a8c8ee",
  "died"               = "#8c2f18")

# A sequential ramp for intensity, because the bands are ORDERED: a reader
# should be able to see escalation as a colour gradient. The route palette above
# is categorical, because routes are not ordered.
DOSE_PAL <- setNames(
  c("#e8e6df", "#a8c8ee", "#4a8bd8", "#14427e", "#7fb069", "#a8c8ee", "#8c2f18"),
  dose_state_levels(DOSE_LAB))
DOSE_PAL[["discharged alive"]] <- "#c9a227"     # must not collide with `low`

# Alluvial at every second window keeps the ribbons legible; 18 axes is a smear.
# The x axis MUST be discrete: with a continuous x, ggalluvial draws the strata
# but no flows at all, because it cannot tell which axes are adjacent.
draw_states <- function(res, pal, tag, title_flow, title_area, sub_area) {
  d <- res$data
  alv <- d[d$window_idx %% 2 == 0,
           c("encounter_block", "window_start_hr", "state")]
  alv$hr <- factor(alv$window_start_hr)
  stopifnot("alluvial data must be in lodes form" =
              is_lodes_form(alv, key = hr, value = state,
                            id = encounter_block, silent = TRUE))

  p_alluvial <- house(
    ggplot(alv, aes(x = hr, stratum = state, alluvium = encounter_block,
                    fill = state)) +
      geom_flow(alpha = 0.55, width = 0.42) +
      geom_stratum(width = 0.42, colour = "#fcfcfb", linewidth = 0.25) +
      scale_fill_manual(values = pal, drop = FALSE) +
      guides(fill = guide_legend(nrow = 2)) +
      labs(title = title_flow,
           subtitle = sprintf(
             "All %s ventilation episodes from the first IMV episode. Terminal states absorb.\nShown every %dh for legibility; the underlying grid is %dh.",
             format(anchor_n, big.mark = ","), 2 * WINDOW_H, WINDOW_H),
           x = "Hours since first IMV episode", y = "Ventilation episodes",
           fill = NULL))
  ggsave(file.path(dirs$phase, sprintf("%salluvial.png", tag)),
         p_alluvial, width = 9.0, height = 5.8, dpi = 200)

  p_states <- house(
    ggplot(res$prevalence, aes(window_start_hr, pct, fill = state)) +
      geom_area(colour = "#fcfcfb", linewidth = 0.2) +
      scale_fill_manual(values = pal, drop = FALSE) +
      scale_x_continuous(breaks = seq(0, EXTENT_H, by = 12)) +
      guides(fill = guide_legend(nrow = 2)) +
      labs(title = title_area, subtitle = sub_area,
           x = "Hours since first IMV episode", y = "% of episodes", fill = NULL))
  ggsave(file.path(dirs$phase, sprintf("%sprevalence.png", tag)),
         p_states, width = 7.5, height = 4.8, dpi = 200)
}

draw_states(route, STATE_PAL, "state_",
  "How fentanyl is delivered, and how episodes leave",
  "Delivery state over the first 72h of ventilation",
  sprintf("%s ventilation episodes; states are mutually exclusive and exhaustive.",
          format(anchor_n, big.mark = ",")))

draw_states(dose, DOSE_PAL, "dose_state_",
  "How much fentanyl, and how episodes leave",
  "Fentanyl intensity over the first 72h of ventilation",
  sprintf("%s episodes. Bands are mcg delivered per %dh window: 0 / <=%s / <=%s / >%s.",
          format(anchor_n, big.mark = ","), WINDOW_H,
          DOSE_CUTS[1], DOSE_CUTS[2], DOSE_CUTS[2]))


# ---- Write -------------------------------------------------------------------

write_out <- function(x, name) {
  f <- file.path(dirs$phase, name)
  write.csv(x, f, row.names = FALSE)
  cat(sprintf("written: %s\n", name))
}

cat("\n")
write_out(state_prevalence, "state_prevalence.csv")
write_out(transition_matrix_tbl, "state_transitions.csv")
write_out(dose$prevalence, "dose_state_prevalence.csv")
write_out(dose$matrix,     "dose_state_transitions.csv")

write_json(prov, file.path(dirs$phase, "provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)
cat("written: provenance.json\n")


# ---- Provenance --------------------------------------------------------------
# Which package versions produced these numbers?

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/03_delivery_states.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "03_delivery_states_sessioninfo.txt")
)
