# ==============================================================================
# 04_delivery_states.R  --  fentanyl delivery states: prevalence and transitions
#
# Purpose : The seven delivery states over the ventilation course -- prevalence
#           per window, transitions window to window, and the figures that draw
#           them. Run twice: delivery ROUTE and dose INTENSITY band.
# Author  : Shan Guleria
# Created : 2026-09-05
# Inputs  : output/intermediate_phi/trajectory_long.parquet, time_to_event.parquet
# Outputs : output/final_no_phi/04_states/ : the CSVs and figures listed in OWNED
# ==============================================================================

# Run this in a FRESH R session; rm(list = ls()) does not unload packages.

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
COV       <- fromJSON(here("config", "protocol.json"), simplifyVector = FALSE)
set.seed(COV$model$seed)

WINDOW_H  <- COV$cohort$window_hours
EXTENT_H  <- COV$cohort$granular_extent_hours
PANELS    <- unlist(COV$cohort$balanced_panel_hours)
LANDMARK  <- COV$cohort$landmark_hours

# Dose-band cuts are federation-critical: read from the config, never defaulted here.
DS_SPEC   <- COV$exposure$dose_states
DOSE_CUTS <- as.numeric(unlist(DS_SPEC$cuts))
DOSE_LAB  <- unlist(DS_SPEC$labels)
stopifnot("protocol.json exposure.dose_states must declare cuts and labels" =
            length(DOSE_CUTS) > 0 && length(DOSE_LAB) == length(DOSE_CUTS) + 2,
          "dose states must be defined on window_mcg" =
            identical(DS_SPEC$variable, "window_mcg"))
SED_UNITS <- COV$exposure$sedatives$units

# Drug list comes from the config; a hardcoded one makes a declared sedative a no-op.
SED_COLS <- unlist(COV$exposure$sedatives$columns)
DRUGS <- c(fentanyl = "total_dose",
           setNames(SED_COLS, sub("_dose$", "", SED_COLS)))
UNITS <- c(fentanyl = COV$exposure$units,
           setNames(vapply(SED_COLS, function(c) SED_UNITS[[c]], character(1)),
                    sub("_dose$", "", SED_COLS)))
SEDATIVES <- setdiff(names(DRUGS), "fentanyl")

MIN_CELL <- COV$reporting$small_cell_min_den

# Collapsed for display only; full CLIF granularity still reaches
# pooling_categorical.csv.
RACE_COLLAPSE <- COV$time_invariant$race$reporting_collapse




# ---- 3. Paths and provenance -------------------------------------------------

dirs <- site_dirs()
dirs$phase <- phase_dir(dirs, "04_states")   # shareable outputs, subdivided by script
prov <- provenance(config)

message(sprintf("[04_delivery_states] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))


# ---- 4. Guards ---------------------------------------------------------------

manifest <- require_manifest(dirs, here())
message(sprintf("  reading Phase 0 outputs from code %s, generated %s",
                manifest$code_version, manifest$generated))

OWNED <- list(phase = c(
  "state_prevalence.csv", "state_transitions.csv", "state_prevalence_at_risk.csv",
  "state_transition_counts.csv", "dose_state_transition_counts.csv",
  "dose_state_prevalence.csv", "dose_state_transitions.csv",
  "state_prevalence.png", "state_alluvial.png", "state_prevalence_at_risk.png",
  "state_raster.png",
  "dose_state_prevalence.png", "dose_state_alluvial.png",
  "provenance.json", "captions.md"))

# Retire the old 03_states/ explicitly; paste0 so a find-and-replace cannot reach it.
RETIRED <- file.path("output", "final_no_phi", paste0("03_", "states"),
                     OWNED$phase)
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

# The Phase 1 denominator is ventilated, not alive_admitted.
long$ventilated <- !is.na(long$imv_status) & long$imv_status == 1
stopifnot(
  "a window cannot be ventilated without being alive and admitted" =
    all(long$ventilated <= long$alive_admitted)
)

missing_cols <- setdiff(unname(DRUGS), names(long))
if (length(missing_cols)) {
  stop("protocol.json declares dose columns that trajectory_long does not carry: ",
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
# Seven states: four routes of delivery while ventilated, three terminal. The
# denominator is the WHOLE analytic cohort, not the landmark set -- the landmark
# conditions on being ventilated at T, so the terminal states would be empty.

state_summary <- function(d, tag) {
  prevalence <- prevalence_by_window(d)
  tp <- transition_pairs(d)
  # Terminal states must absorb: a terminal state as an ORIGIN means the
  # definition is wrong.
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

  list(prevalence = prevalence, matrix = tmx, counts = transition_counts(tp),
       data = d)
}

long  <- derive_states(long)                              # route definition
route <- state_summary(long, "Delivery")

# derive_dose_states() overwrites `state`, so it runs on a copy; `long` keeps
# the routes.
dose <- state_summary(
  derive_dose_states(long, DOSE_CUTS, DOSE_LAB),
  sprintf("Intensity (0 / <=%s / <=%s / >%s mcg per %dh window)",
          DOSE_CUTS[1], DOSE_CUTS[2], DOSE_CUTS[2], WINDOW_H))

state_prevalence      <- route$prevalence
transition_matrix_tbl <- route$matrix


# ---- 12c. At-risk prevalence -------------------------------------------------
# The same states conditioned on still being ventilated: a denominator change, not a
# second definition, which is what droplevels() makes explicit. Both views ship.

at_risk <- long[long$ventilated, ]

# derive_states() reads inf_dose and bolus_dose, never total_dose: not a tautology.
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

# Cross-script check: 02 builds dose_summary on the identical subset, so the exposed
# states must sum to its n_receiving in every window.
ds_file <- file.path(dirs$out_final, "02_descriptive", "dose_summary.csv")
if (file.exists(ds_file)) {
  ds <- read.csv(ds_file)
  ds <- ds[ds$drug == "fentanyl" & ds$denominator == "all_ventilated",
           c("window_idx", "n", "n_receiving")]
  exposed <- aggregate(
    n ~ window_idx,
    data = state_prevalence_at_risk[state_prevalence_at_risk$state != "no fentanyl", ],
    FUN = sum)
  names(exposed)[names(exposed) == "n"] <- "n_exposed"
  chk <- merge(ds, exposed, by = "window_idx")
  n_here <- aggregate(n ~ window_idx, data = state_prevalence_at_risk, FUN = sum)
  names(n_here)[names(n_here) == "n"] <- "n_here"
  chk <- merge(chk, n_here, by = "window_idx")
  stopifnot(
    "every window must appear in both tables" = nrow(chk) == nrow(ds),
    "the at-risk denominator disagrees with dose_summary's all_ventilated n" =
      all(chk$n == chk$n_here),
    # Exact: dose_summary ships the numerator, so the counts must agree outright.
    "exposed states do not sum to dose_summary's receiving count" =
      all(chk$n_exposed == chk$n_receiving))
  cat(sprintf("  cross-check vs 02_descriptive/dose_summary.csv: %d windows agree exactly\n",
              nrow(chk)))
} else {
  cat("  cross-check SKIPPED: 02_descriptive/dose_summary.csv not found -- run 02 first\n")
}


# ---- 13. Figures -------------------------------------------------------------
# Journal style: no titles or subtitles on the panels. Every draw call registers its
# caption, and the set is written to captions.md beside the figures.


# --- Delivery-state figures ---------------------------------------------------
# Alluvial: one ribbon per episode, its width the episodes moving between states.

# Palettes live in code/utils/figures.R only; the alluvial and the raster must match.
STATE_PAL <- STATE_COLOURS
DOSE_PAL  <- dose_palette(DOSE_LAB)
stopifnot("the dose palette must cover exactly the declared band levels" =
            identical(names(DOSE_PAL), dose_state_levels(DOSE_LAB)))

# Every second window keeps the ribbons legible. The x axis must be discrete: with a
# continuous x, ggalluvial draws the strata but no flows at all.
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

# One drawing on two denominators: the whole-cohort figure and its at-risk twin.
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

# Window 0, then the windows ENDING at 24, 48 and 72h. Everywhere else a window is
# named by its START, so the one ending at 24h starts at 20h; the table header names
# the end hour. Dropped where a site's window_hours puts them off the grid.
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
n_at_risk <- aggregate(n ~ window_start_hr, data = state_prevalence_at_risk, FUN = sum)
last_hr   <- max(n_at_risk$window_start_hr)

draw_prevalence_area(
  state_prevalence_at_risk, STATE_PAL, "state_prevalence_at_risk.png",
  "What the patients still on the ventilator are getting",
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
# Individual courses, after Iyer et al. Figure 1B. Rows are a simple random sample,
# not stratified on the exit, so the panel reads as the cohort does.

N_RASTER    <- 100
EXIT_STATES <- c("extubated", STATE_ABSORBING)

# Re-seeded here: any RNG in between would change which patients a reader sees.
set.seed(COV$model$seed)
blocks <- unique(long$encounter_block)
stopifnot("not enough episodes to fill the raster" = length(blocks) >= N_RASTER)
picked <- sample(blocks, N_RASTER)

r <- long[long$encounter_block %in% picked,
          c("encounter_block", "window_idx", "window_start_hr", "state")]

# The grid is already rectangular -- terminal states persist in the source data.
stopifnot("every sampled episode must span the full window grid" =
            all(table(r$encounter_block) == n_windows_expected))

# Row order: final state, time to first exit descending, then exposed windows;
# ties all broken.
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
          format(COV$model$seed, scientific = FALSE), WINDOW_H, last_hr))
ggsave(file.path(dirs$phase, "state_raster.png"),
       p_raster, width = 7.5, height = 5.5, dpi = 200)


# ---- Write -------------------------------------------------------------------

write_out <- function(x, name) {
  f <- file.path(dirs$phase, name)
  write.csv(x, f, row.names = FALSE)
  cat(sprintf("written: %s\n", name))
}

# Reported, not suppressed: these tables are prevalence, so NA-ing a cell stops the
# column summing to 100 -- a disclosure decision, not a code decision.
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
write_out(route$counts, "state_transition_counts.csv")
write_out(dose$counts,  "dose_state_transition_counts.csv")
write_out(dose$prevalence, "dose_state_prevalence.csv")
write_out(dose$matrix,     "dose_state_transitions.csv")
write_out(state_prevalence_at_risk, "state_prevalence_at_risk.csv")

write_json(prov, file.path(dirs$phase, "provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)
cat("written: provenance.json\n")

write_captions(file.path(dirs$phase, "captions.md"), "04_delivery_states.R",
               grep("\\.png$", OWNED$phase, value = TRUE), prov)
cat("written: captions.md\n")


# ---- Provenance --------------------------------------------------------------

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/04_delivery_states.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "04_delivery_states_sessioninfo.txt")
)
