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
  "state_prevalence.csv", "state_transitions.csv", "state_prevalence_at_risk.csv",
  "dose_state_prevalence.csv", "dose_state_transitions.csv",
  "state_prevalence.png", "state_alluvial.png", "state_prevalence_at_risk.png",
  "state_raster.png",
  "dose_state_prevalence.png", "dose_state_alluvial.png",
  "provenance.json", "captions.md"))

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
  # Denominator = the whole analytic cohort, terminal states included. The
  # at-risk twin (ventilated only) is built further down from the same helper.
  prevalence <- prevalence_by_window(d)
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


# ---- 12c. At-risk prevalence -------------------------------------------------
# The same states on the denominator the question actually asks about.
#
# state_prevalence divides by the whole analytic cohort, so every band shrinks
# as patients are extubated, discharged or die -- by hour 68 more than half the
# cohort has left (14,897 -> 6,728 still ventilated). That is the right picture
# for liberation and the wrong one for delivery: a band falling because patients
# went home reads identically to one falling because fentanyl was weaned.
#
# Conditioning on still being ventilated separates the two. Among ventilated
# windows only the four ROUTE states can occur -- the three terminal states each
# require !alive_admitted or !ventilated -- so this is a denominator change, not
# a second definition, and droplevels() is what makes that explicit.
#
# BOTH ship. Neither is a corrected version of the other (SG, 2026-09-24).

at_risk <- long[long$ventilated, ]

# The classification and the dose column are derived independently --
# derive_states() reads inf_dose and bolus_dose, never total_dose -- so this is
# a real cross-check on the route definition, not a tautology.
# 01_build_cohort.py:712 sets total_dose = inf_dose + bolus_dose, both >= 0.
stopifnot(
  "a ventilated window cannot be in a terminal state" =
    !any(at_risk$state %in% c("extubated", STATE_ABSORBING)),
  "route state disagrees with total_dose about who received fentanyl" =
    all((at_risk$state != "no fentanyl") == (at_risk$total_dose > 0)))

at_risk$state <- droplevels(at_risk$state)
state_prevalence_at_risk <- prevalence_by_window(at_risk)

cat("\nAt-risk delivery states, % of episodes STILL VENTILATED in that window\n")
pv_ar <- reshape(state_prevalence_at_risk[, c("window_start_hr", "state", "pct")],
                 idvar = "window_start_hr", timevar = "state", direction = "wide")
names(pv_ar) <- sub("^pct\\.", "", names(pv_ar))
print(pv_ar[pv_ar$window_start_hr %in% c(0, 24, 48, 68), ], row.names = FALSE)

# Cross-script check. 02_descriptive_cohort.R:221 builds dose_summary on the
# identical subset (`window_idx == w & ventilated`), so the three exposed states
# must sum to its pct_receiving_any in every window. Precedent for a later
# script checking an earlier one's shareable CSV: 04_landmark_cohort.R:228.
ds_file <- file.path(dirs$out_final, "02_descriptive", "dose_summary.csv")
if (file.exists(ds_file)) {
  ds <- read.csv(ds_file)
  ds <- ds[ds$drug == "fentanyl" & ds$denominator == "all_ventilated",
           c("window_idx", "n", "pct_receiving_any")]
  exposed <- aggregate(
    pct ~ window_idx,
    data = state_prevalence_at_risk[state_prevalence_at_risk$state != "no fentanyl", ],
    FUN = sum)
  chk <- merge(ds, exposed, by = "window_idx")
  n_here <- aggregate(n ~ window_idx, data = state_prevalence_at_risk, FUN = sum)
  chk <- merge(chk, n_here, by = "window_idx", suffixes = c("_02", "_here"))
  stopifnot(
    "every window must appear in both tables" = nrow(chk) == nrow(ds),
    "the at-risk denominator disagrees with dose_summary's all_ventilated n" =
      all(chk$n_02 == chk$n_here),
    # 0.1pp: dose_summary rounds pct_receiving_any to 1dp, and this sums three
    # values each rounded to 2dp.
    "exposed states do not sum to dose_summary's pct_receiving_any" =
      max(abs(chk$pct - chk$pct_receiving_any)) <= 0.1)
  cat(sprintf("  cross-check vs 02_descriptive/dose_summary.csv: %d windows agree, max gap %.2fpp\n",
              nrow(chk), max(abs(chk$pct - chk$pct_receiving_any))))
} else {
  cat("  cross-check SKIPPED: 02_descriptive/dose_summary.csv not found -- run 02 first\n")
}


# ---- 13. Figures -------------------------------------------------------------
#
# JOURNAL STYLE: no figure carries a title or a subtitle. A journal figure is a
# panel -- axes, legend, data -- and the n, the denominator, the seed and the
# sampling rule belong in the caption, which is manuscript text. Nothing is
# discarded to achieve that: every draw call registers its caption here and the
# set is written to captions.md beside the figures.


# --- Delivery-state figures ---------------------------------------------------
# Alluvial: every episode is a ribbon, its width the number of episodes moving
# between states. Ribbons entering a terminal state never leave it, so
# liberation and death are visible as the flow drains out of the fentanyl states.

# Both palettes live in code/utils/figures.R. They were declared here and
# again there, and the two copies had drifted on six of seven states before
# 2026-09-24 -- the alluvial and the raster must agree on colour or a reader
# cannot move between them.
STATE_PAL <- STATE_COLOURS
DOSE_PAL  <- dose_palette(DOSE_LAB)
stopifnot("the dose palette must cover exactly the declared band levels" =
            identical(names(DOSE_PAL), dose_state_levels(DOSE_LAB)))

# Alluvial at every second window keeps the ribbons legible; 18 axes is a smear.
# The x axis MUST be discrete: with a continuous x, ggalluvial draws the strata
# but no flows at all, because it cannot tell which axes are adjacent.
draw_states <- function(res, pal, tag, title_flow, title_area, note_area) {
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
      geom_stratum(width = 0.42, colour = PAPER, linewidth = 0.25) +
      scale_fill_manual(values = pal, drop = FALSE) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
      labs(x = "Hours since first IMV episode", y = "Ventilation episodes",
           fill = NULL))
  register_caption(sprintf("%salluvial.png", tag), title_flow, sprintf(
    "All %s ventilation episodes from the first IMV episode. Terminal states absorb. Shown every %dh for legibility; the underlying grid is %dh.",
    format(anchor_n, big.mark = ","), 2 * WINDOW_H, WINDOW_H))
  ggsave(file.path(dirs$phase, sprintf("%salluvial.png", tag)),
         p_alluvial, width = 9.0, height = 5.8, dpi = 200)

  draw_prevalence_area(res$prevalence, pal, sprintf("%sprevalence.png", tag),
                       title_area, note_area, "% of episodes",
                       denom_label = "All episodes")
}

# The stacked-area panel, factored out so the whole-cohort figure and its
# at-risk twin below are the same drawing on two denominators. A thin PAPER
# stroke between bands is the surface gap that keeps adjacent fills separable
# without adding a second encoding.
draw_prevalence_area <- function(prevalence, pal, file, title, note, ylab,
                                 denom_label) {
  wins <- risk_table_windows(prevalence$window_start_hr, WINDOW_H)

  p_main <- house(
    ggplot(prevalence, aes(window_start_hr, pct, fill = state)) +
      geom_area(colour = PAPER, linewidth = 0.2) +
      scale_fill_manual(values = pal, drop = FALSE) +
      scale_x_continuous(breaks = seq(0, EXTENT_H, by = 12)) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
      labs(x = NULL, y = ylab, fill = NULL))

  p_tab <- risk_table_panel(prevalence, wins, denom_label,
                            x_title = "Hours since first IMV episode")
  g <- stack_with_risk_table(p_main, p_tab)

  register_caption(file, title, paste0(
    note, " Table: n (% of the denominator) for window 0 and for the windows ",
    "ending at ", paste(wins$label[-1], collapse = ", "), " hours."))
  ggsave(file.path(dirs$phase, file), g,
         width = 7.5, height = 4.4 + 0.20 * attr(p_tab, "n_rows") + 0.3, dpi = 200)
}

# Windows the risk table reports: window 0, then the windows ENDING at 24, 48
# and 72 hours, so the columns are evenly spaced in elapsed time. Everywhere
# else in this script a window is named by its START, so the one ending at 24h
# starts at 20h -- the table carries its own header naming the END hour, which
# is the quantity a reader wants. Window 0 has no 24h-spaced end of its own and
# is reported as the baseline, labelled 0.
#
# Kept only where the grid actually carries them: at a site running a different
# window_hours the start may not be a window start, and a column that silently
# vanishes is worse than a shorter table.
risk_table_windows <- function(hrs, window_h, ends = c(24, 48, 72)) {
  starts <- ends - window_h
  keep   <- starts %in% hrs
  if (any(!keep)) {
    message(sprintf("  risk table: window(s) ending at %s are not on the %dh grid, dropped",
                    paste(ends[!keep], collapse = ", "), window_h))
  }
  stopifnot("window 0 must exist" = 0 %in% hrs)
  data.frame(start_hr = c(0, starts[keep]),
             label    = c("0", as.character(ends[keep])))
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

# --- At-risk prevalence -------------------------------------------------------
# The twin of the figure above on the denominator the question asks about. The
# numbers in the subtitle are read off the data, never typed: a hardcoded
# denominator is one re-run away from describing a cohort that no longer exists.
n_at_risk <- aggregate(n ~ window_start_hr, data = state_prevalence_at_risk, FUN = sum)
last_hr   <- max(n_at_risk$window_start_hr)

draw_prevalence_area(
  state_prevalence_at_risk, STATE_PAL, "state_prevalence_at_risk.png",
  "What the patients still on the ventilator are getting",
  # Caption text, not panel text -- so it can run to full sentences. The
  # denominators are read off the data, never typed: a hardcoded count is one
  # re-run away from describing a cohort that no longer exists.
  sprintf(paste0("Denominator is ventilation episodes STILL VENTILATED in the ",
                 "window, not the whole cohort: %s at the start of ventilation ",
                 "falling to %s in the window ending at %dh. Bands are stacked ",
                 "in state order, with `no fentanyl` uppermost. The x axis is ",
                 "labelled by window START; the risk table is labelled by ",
                 "window END, so its %dh column is the %dh window covering ",
                 "hours %d-%d."),
          format(n_at_risk$n[n_at_risk$window_start_hr == 0], big.mark = ","),
          format(n_at_risk$n[n_at_risk$window_start_hr == last_hr], big.mark = ","),
          last_hr + WINDOW_H, last_hr + WINDOW_H, WINDOW_H, last_hr,
          last_hr + WINDOW_H),
  "% of episodes still ventilated",
  denom_label = "Still ventilated")

# --- Per-patient raster -------------------------------------------------------
# The one figure here that draws individual courses. Everything else is a
# marginal; this is what makes "continuous, then bolus, then off" legible as a
# sequence rather than as four curves that happen to move. Iyer et al. Figure 1B
# is the precedent for the form -- not for the row ordering, which that paper
# does not state, or for its ragged right edge, which we do not copy.
#
# Rows are a SIMPLE random sample. Not stratified on the exit: stratifying would
# over-represent death and discharge relative to the cohort, and the panel is
# meant to read as "what 100 episodes look like", not as a balanced design.

N_RASTER    <- 100
EXIT_STATES <- c("extubated", STATE_ABSORBING)

# Re-seeded here rather than relied on from the top of the script: any RNG in
# between would silently change which patients a reader sees.
set.seed(config$model$seed)
blocks <- unique(long$encounter_block)
stopifnot("not enough episodes to fill the raster" = length(blocks) >= N_RASTER)
picked <- sample(blocks, N_RASTER)

r <- long[long$encounter_block %in% picked,
          c("encounter_block", "window_idx", "window_start_hr", "state")]

# A ragged right edge would make the raster and the alluvial disagree about who
# is still in the picture. The grid is already rectangular, because the terminal
# states persist in the source data rather than being carried forward here --
# assert that rather than trust it.
stopifnot("every sampled episode must span the full window grid" =
            all(table(r$encounter_block) == n_windows_expected))

# Hierarchical row order: state at the last window, then time to the first exit
# descending, then count of fentanyl-exposed windows. Every tie is broken, so
# the panel is reproducible from the seed alone.
key <- do.call(rbind, lapply(split(r, r$encounter_block), function(e) {
  exit <- which(e$state %in% EXIT_STATES)
  data.frame(encounter_block = e$encounter_block[1],
             final_state  = e$state[which.max(e$window_idx)],
             t_first_exit = if (length(exit)) min(e$window_idx[exit]) else Inf,
             n_exposed    = sum(!(e$state %in% c("no fentanyl", EXIT_STATES))))
}))
key <- key[order(key$final_state, -key$t_first_exit, -key$n_exposed), ]
key$row_index <- seq_len(nrow(key))

r <- merge(r, key[, c("encounter_block", "row_index")], by = "encounter_block")

# De-identification asserted in code, not left to review. The plotted frame
# carries an anonymous row index, hours RELATIVE to the first IMV episode, and a
# state -- nothing that identifies or links a patient.
raster <- r[, c("row_index", "window_start_hr", "state")]
stopifnot(
  "the raster frame must carry no identifier" =
    !any(c("encounter_block", "patient_id", "hospitalization_id") %in% names(raster)),
  "rows must be an anonymous 1..N index" =
    identical(sort(unique(raster$row_index)), seq_len(N_RASTER)),
  "the x axis must be relative hours, never a date" =
    is.numeric(raster$window_start_hr) && min(raster$window_start_hr) == 0)

p_raster <- house(
  ggplot(raster, aes(window_start_hr, row_index, fill = state)) +
    geom_raster() +
    scale_fill_manual(values = STATE_PAL, drop = FALSE) +
    scale_x_continuous(breaks = seq(0, EXTENT_H, by = 12), expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(x = "Hours since first IMV episode", y = NULL, fill = NULL)) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

register_caption("state_raster.png", "One hundred ventilation episodes, end to end",
  sprintf(paste0("A random %d of %s episodes, seeded at %s; one row per episode, ",
                 "one cell per %dh window. Ordered by state at hour %d, then by ",
                 "time to first extubation, death or discharge, then by number of ",
                 "fentanyl-exposed windows. Relative hours only; no identifiers."),
          N_RASTER, format(anchor_n, big.mark = ","),
          format(config$model$seed, scientific = FALSE), WINDOW_H, last_hr))
ggsave(file.path(dirs$phase, "state_raster.png"),
       p_raster, width = 7.5, height = 5.5, dpi = 200)


# ---- Write -------------------------------------------------------------------

write_out <- function(x, name) {
  f <- file.path(dirs$phase, name)
  write.csv(x, f, row.names = FALSE)
  cat(sprintf("written: %s\n", name))
}

# Small cells in the shareable set. config.json declares
# reporting.small_cell_min_den, and until 2026-09-24 this script read the value
# and applied it to nothing -- a threshold that is declared but never consumed
# reads as a policy and is not one. Reported per file rather than suppressed:
# these tables are prevalence, so NA-ing a cell stops the column summing to 100,
# and that trade is a disclosure decision, not a code decision. Whoever runs
# this sees the count; at another site the counts will differ.
report_small_cells <- function(x, name) {
  if (!all(c("n", "state", "window_start_hr") %in% names(x))) return(invisible(0L))
  small <- x[x$n > 0 & x$n < MIN_CELL, ]
  if (!nrow(small)) return(invisible(0L))
  cat(sprintf("  NOTE %s: %d cell(s) below reporting.small_cell_min_den = %d --\n",
              name, nrow(small), MIN_CELL))
  for (i in seq_len(nrow(small))) {
    cat(sprintf("       hour %-3s %-18s n = %d\n",
                small$window_start_hr[i], small$state[i], small$n[i]))
  }
  invisible(nrow(small))
}

cat("\nSmall-cell check against reporting.small_cell_min_den\n")
n_small <- sum(
  report_small_cells(state_prevalence,         "state_prevalence.csv"),
  report_small_cells(state_prevalence_at_risk, "state_prevalence_at_risk.csv"),
  report_small_cells(dose$prevalence,          "dose_state_prevalence.csv"))
if (!n_small) cat(sprintf("  all cells at or above %d\n", MIN_CELL))

cat("\n")
write_out(state_prevalence, "state_prevalence.csv")
write_out(transition_matrix_tbl, "state_transitions.csv")
write_out(dose$prevalence, "dose_state_prevalence.csv")
write_out(dose$matrix,     "dose_state_transitions.csv")
write_out(state_prevalence_at_risk, "state_prevalence_at_risk.csv")

write_json(prov, file.path(dirs$phase, "provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)
cat("written: provenance.json\n")

# The figures carry no title or subtitle -- they are journal panels. What used
# to sit on them is written here instead, so the n, the denominator, the seed
# and the ordering rules survive in a form a manuscript can lift directly.
# Every figure this script writes must appear; a panel whose provenance exists
# only in the code is a panel nobody can caption later.
write_captions(file.path(dirs$phase, "captions.md"), "03_delivery_states.R",
               grep("\\.png$", OWNED$phase, value = TRUE), prov)
cat("written: captions.md\n")


# ---- Provenance --------------------------------------------------------------
# Which package versions produced these numbers?

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/03_delivery_states.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "03_delivery_states_sessioninfo.txt")
)
