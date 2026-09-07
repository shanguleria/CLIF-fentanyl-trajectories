# ==============================================================================
# 02_descriptive_trajectory.R  --  Phase 1 -- cohort dose curves + balanced panels
#
# Purpose : Whole-cohort dose trajectory over the granular grid, with balanced-panel overlays to separate real dose change from cohort composition change, plus the retention table that chooses the landmark T and the zero fraction that chooses the model family.
# Author  : Shan Guleria
# Created : 2026-09-05
# Inputs  : output/intermediate_phi/trajectory_long.parquet, time_to_event.parquet
# Outputs : output/final_no_phi/ : phase1_*.csv, phase1_*.png, phase1_provenance.json
#
# Spec: docs/design_notes.md section 10.
# ==============================================================================

# Run this in a FRESH R session (RStudio: Cmd+Shift+F10).
# Do not use rm(list = ls()) -- it does not unload packages or reset options,
# so it only gives the appearance of a clean slate.


# ---- 1. Packages -------------------------------------------------------------

pkgs <- c("here", "jsonlite", "ggplot2", "arrow")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

source(here("code", "utils", "paths.R"))


# ---- 2. Config (never setwd(); here() anchors to the .Rproj) -----------------

config <- fromJSON(here("config", "config.json"), simplifyVector = FALSE)
set.seed(config$model$seed)

WINDOW_H  <- config$cohort$window_hours
EXTENT_H  <- config$cohort$granular_extent_hours
PANELS    <- unlist(config$cohort$balanced_panel_hours)
LANDMARK  <- config$cohort$landmark_hours

COV       <- fromJSON(here("config", "covariates.json"), simplifyVector = FALSE)
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


# ---- 3. Paths and provenance -------------------------------------------------
# One site, one output tree. site_dirs() creates them and labels the PHI ones.

dirs <- site_dirs()
prov <- provenance(config)

message(sprintf("[02_descriptive_trajectory] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))


# ---- 4. Guards ---------------------------------------------------------------
# Phase 0 (Python) wrote what this script reads. Refuse to analyse tables that a
# different code version or a different config produced.

manifest <- require_manifest(dirs, here())
message(sprintf("  reading Phase 0 outputs from code %s, generated %s",
                manifest$code_version, manifest$generated))

OWNED <- list(out_final = c(
  "phase1_baseline_characteristics.csv", "phase1_retention.csv",
  "phase1_dose_summary.csv", "phase1_dose_distribution.csv",
  "phase1_balanced_panels.csv", "phase1_zero_fraction.csv",
  "phase1_imv_episodes.csv", "phase1_choosing_T.csv", "phase1_provenance.json",
  "phase1_fentanyl_curves.png", "phase1_fentanyl_balanced_panels.png",
  "phase1_fentanyl_distribution.png", "phase1_sedative_curves.png"))

# Figures renamed 2026-09-07 when fentanyl became the primary view. A rename
# leaves a stale twin the owned list no longer names.
RETIRED <- file.path("output", "final_no_phi",
                     c("phase1_dose_curves.png", "phase1_pct_receiving.png",
                       "phase1_dose_distribution.png", "phase1_balanced_panels.png"))
n_cleared <- clear_owned_outputs(dirs, OWNED, retired = RETIRED)
if (n_cleared) message(sprintf("  cleared %d output(s) from a previous run", n_cleared))


# ---- 5. Read -----------------------------------------------------------------

long <- as.data.frame(read_parquet(
  file.path(dirs$out_phi, "trajectory_long.parquet"),
  col_select = c("encounter_block", "patient_id", "window_idx", "window_start_hr",
                 "alive_admitted", "imv_status", "inf_dose", "bolus_dose",
                 "total_dose", "propofol_dose", "midazolam_dose",
                 "dexmedetomidine_dose",
                 "age", "sex", "race", "cci", "bmi_admission", "weight_kg",
                 "sofa_total", "nee", "oxygenation", "lactate",
                 "first_imv_episode_hours", "n_imv_episodes")))

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
# more than half the cohort by 72h; design_notes.md section 10 works the example.
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


# ---- 6. Retention ------------------------------------------------------------
# at-risk = STILL VENTILATED. alive_admitted is carried beside it so the gap is
# visible rather than implicit.

retention <- do.call(rbind, lapply(sort(unique(long$window_idx)), function(w) {
  x <- long[long$window_idx == w, ]
  data.frame(window_idx = w,
             window_start_hr = x$window_start_hr[1],
             n_ventilated = sum(x$ventilated),
             n_alive_admitted = sum(x$alive_admitted),
             pct_of_anchor_cohort = round(100 * sum(x$ventilated) / anchor_n, 1))
}))

# Section 8's "Choosing T" tabulation: state at each candidate landmark, which is
# the window whose span ENDS at that hour.
choosing_T <- do.call(rbind, lapply(c(12, 24, 48, 72), function(h) {
  w <- h / WINDOW_H - 1
  if (!(w %in% retention$window_idx)) return(NULL)
  r <- retention[retention$window_idx == w, ]
  data.frame(T_hours = h, n_ventilated = r$n_ventilated,
             pct_of_intubated_cohort = r$pct_of_anchor_cohort,
             windows_available = h / WINDOW_H)
}))

cat("\nRetention (at-risk = still ventilated)\n")
print(retention, row.names = FALSE)
cat("\nChoosing T (design_notes.md section 8)\n")
print(choosing_T, row.names = FALSE)


# ---- 7. Dose summaries -------------------------------------------------------
# Three curves, not one: weaning-to-zero among the still-ventilated and dropout
# of low-dose patients move the overall median in opposite directions.


summarise_dose <- function(v, drug, denom, w, hr) {
  v <- v[!is.na(v)]
  if (!length(v)) return(NULL)
  q <- unname(quantile(v, c(0.25, 0.5, 0.75)))
  data.frame(drug = drug, unit = UNITS[[drug]], denominator = denom,
             window_idx = w, window_start_hr = hr, n = length(v),
             median = round(q[2], 3), q1 = round(q[1], 3), q3 = round(q[3], 3),
             mean = round(mean(v), 3), sd = round(stats::sd(v), 3),
             pct_receiving_any = round(100 * mean(v > 0), 1))
}

dose_summary <- do.call(rbind, lapply(names(DRUGS), function(drug) {
  col <- DRUGS[[drug]]
  do.call(rbind, lapply(sort(unique(long$window_idx)), function(w) {
    x <- long[long$window_idx == w & long$ventilated, ]
    rbind(summarise_dose(x[[col]], drug, "all_ventilated", w, x$window_start_hr[1]),
          summarise_dose(x[[col]][x[[col]] > 0], drug, "receivers_only", w,
                         x$window_start_hr[1]))
  }))
}))

cat("\nFentanyl dose by window (all ventilated, zeros included)\n")
print(dose_summary[dose_summary$drug == "fentanyl" &
                     dose_summary$denominator == "all_ventilated",
                   c("window_start_hr", "n", "median", "q1", "q3",
                     "pct_receiving_any")], row.names = FALSE)


# ---- 8. Distribution shape ---------------------------------------------------
# The zero fraction says how much of the mass sits at zero; it says nothing about
# whether the rest is skewed or multimodal, which is the other half of the
# gbmt-vs-crimCV question.

dose_distribution <- do.call(rbind, lapply(names(DRUGS), function(drug) {
  v <- long[[DRUGS[[drug]]]][long$ventilated]
  v <- v[!is.na(v)]
  nz <- v[v > 0]
  probs <- seq(0.1, 0.9, by = 0.1)
  rbind(
    data.frame(drug = drug, unit = UNITS[[drug]], population = "all_ventilated",
               n = length(v), decile = probs * 10,
               value = round(unname(quantile(v, probs)), 3)),
    if (length(nz)) data.frame(
      drug = drug, unit = UNITS[[drug]], population = "non_zero_only",
      n = length(nz), decile = probs * 10,
      value = round(unname(quantile(nz, probs)), 3))
  )
}))


# ---- 9. Balanced panels ------------------------------------------------------
# The all-ventilated curve answers a different question at every timepoint,
# because its denominator keeps changing. A frozen denominator makes movement
# mean within-patient change. Each panel spans only the window it is balanced for.

still_vent_at <- function(h) {
  w <- h / WINDOW_H - 1
  unique(long$encounter_block[long$window_idx == w & long$ventilated])
}

balanced_panels <- do.call(rbind, lapply(PANELS, function(h) {
  ids <- still_vent_at(h)
  if (!length(ids)) return(NULL)
  x <- long[long$encounter_block %in% ids & long$window_start_hr < h, ]
  do.call(rbind, lapply(sort(unique(x$window_idx)), function(w) {
    y <- x[x$window_idx == w, ]
    s <- summarise_dose(y$total_dose, "fentanyl", sprintf("balanced_ge_%dh", h),
                        w, y$window_start_hr[1])
    if (is.null(s)) return(NULL)
    s$panel_hours <- h
    s$panel_n <- length(ids)
    s
  }))
}))

cat("\nBalanced panels (fentanyl, frozen denominators)\n")
for (h in PANELS) {
  p <- balanced_panels[balanced_panels$panel_hours == h, ]
  if (!nrow(p)) next
  cat(sprintf("  >=%2dh  n=%5s  median %.2f -> %.2f %s (within-panel change %+.1f%%)\n",
              h, format(p$panel_n[1], big.mark = ","), p$median[1],
              p$median[nrow(p)], UNITS[["fentanyl"]],
              100 * (p$median[nrow(p)] / p$median[1] - 1)))
}


# ---- 10. Zero fraction -------------------------------------------------------
# Governs gbmt (continuous) versus crimCV (zero-inflated Poisson), section 9.
# The number that decides it is the LANDMARK-cohort one: that is the population
# Phases 2-6 actually model.

zero_fraction <- do.call(rbind, lapply(names(DRUGS), function(drug) {
  col <- DRUGS[[drug]]
  do.call(rbind, lapply(list(
    list("whole_cohort", rep(TRUE, nrow(long))),
    list("landmark_cohort", long$encounter_block %in% elig)
  ), function(p) {
    v <- long[[col]][long$ventilated & p[[2]]]
    v <- v[!is.na(v)]
    data.frame(drug = drug, unit = UNITS[[drug]], population = p[[1]],
               n_ventilated_windows = length(v),
               n_zero = sum(v == 0),
               pct_zero = round(100 * mean(v == 0), 1))
  }))
}))

cat("\nZero fraction of dose across ventilated windows\n")
print(zero_fraction, row.names = FALSE)


# ---- 11. IMV episodes --------------------------------------------------------

blocks <- long[!duplicated(long$encounter_block),
               c("encounter_block", "patient_id", "first_imv_episode_hours",
                 "n_imv_episodes")]
blocks$landmark_eligible <- blocks$encounter_block %in% elig

q_row <- function(v, name, unit) {
  v <- v[!is.na(v)]
  q <- unname(quantile(v, c(0.25, 0.5, 0.75, 0.9, 0.95)))
  data.frame(measure = name, unit = unit, n = length(v),
             median = round(q[2], 2), q1 = round(q[1], 2), q3 = round(q[3], 2),
             p90 = round(q[4], 2), p95 = round(q[5], 2),
             min = round(min(v), 2), max = round(max(v), 2))
}

per_patient <- table(blocks$patient_id)
imv_episodes <- rbind(
  q_row(blocks$first_imv_episode_hours, "first_imv_episode_hours", "hours"),
  q_row(blocks$n_imv_episodes, "n_imv_episodes_per_block", "count"),
  q_row(as.numeric(per_patient), "blocks_per_patient", "count")
)

cat("\nVentilation episodes\n")
print(imv_episodes, row.names = FALSE)
cat(sprintf("  patients contributing more than one block: %s of %s (%.1f%%)\n",
            format(sum(per_patient > 1), big.mark = ","),
            format(length(per_patient), big.mark = ","),
            100 * mean(per_patient > 1)))


# ---- 12. Baseline characteristics --------------------------------------------
# Strata are MUTUALLY EXCLUSIVE (eligible vs not) so the p-value is meaningful;
# a whole-cohort column against its own subset would not be. This doubles as the
# included-versus-excluded disclosure section 8 requires.
# Shape is gtsummary-like: __bold__ parent labels, indented sub-levels, N in the
# strata headers. Rendering to HTML is a separate step.

base <- long[long$window_idx == 0, ]
base$landmark_eligible <- base$encounter_block %in% elig
grp <- ifelse(base$landmark_eligible, "eligible", "not_eligible")

hdr <- function(label, n) sprintf("%s (N = %s)", label, format(n, big.mark = ","))
COLS <- c(hdr("Overall", nrow(base)),
          hdr("Landmark-eligible", sum(grp == "eligible")),
          hdr("Not eligible", sum(grp == "not_eligible")))

fmt_p <- function(p) {
  if (is.na(p)) return("--")
  if (p < 0.001) "<0.001" else sprintf("%.3f", p)
}

fmt_iqr <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return("--")
  q <- unname(quantile(v, c(0.25, 0.5, 0.75)))
  sprintf("%.1f (%.1f, %.1f)", q[2], q[1], q[3])
}
fmt_pct <- function(k, n) sprintf("%s (%.1f%%)", format(k, big.mark = ","),
                                  100 * k / max(n, 1))

row_continuous <- function(label, v) {
  p <- tryCatch(stats::wilcox.test(v[grp == "eligible"], v[grp == "not_eligible"])$p.value,
                error = function(e) NA_real_)
  data.frame(characteristic = sprintf("__%s__, median (IQR)", label),
             overall = fmt_iqr(v),
             eligible = fmt_iqr(v[grp == "eligible"]),
             not_eligible = fmt_iqr(v[grp == "not_eligible"]),
             p_value = fmt_p(p), stringsAsFactors = FALSE)
}

row_categorical <- function(label, v) {
  v <- ifelse(is.na(v), "Missing", as.character(v))
  lv <- sort(unique(v))
  tab <- table(v, grp)
  keep <- rownames(tab)[rownames(tab) != "Missing"]
  p <- if (length(keep) > 1 && all(dim(tab[keep, , drop = FALSE]) > 1))
    tryCatch(stats::chisq.test(tab[keep, , drop = FALSE])$p.value,
             error = function(e) NA_real_) else NA_real_
  head <- data.frame(characteristic = sprintf("__%s__, n (%%)", label),
                     overall = "", eligible = "", not_eligible = "",
                     p_value = fmt_p(p), stringsAsFactors = FALSE)
  body <- do.call(rbind, lapply(lv, function(l) data.frame(
    characteristic = paste0("    ", l),
    overall = fmt_pct(sum(v == l), length(v)),
    eligible = fmt_pct(sum(v == l & grp == "eligible"), sum(grp == "eligible")),
    not_eligible = fmt_pct(sum(v == l & grp == "not_eligible"), sum(grp == "not_eligible")),
    p_value = "", stringsAsFactors = FALSE)))
  rbind(head, body)
}

baseline <- rbind(
  row_continuous("Age, years", base$age),
  row_categorical("Sex", base$sex),
  row_categorical("Race", base$race),
  row_continuous("Charlson Comorbidity Index", base$cci),
  row_continuous("BMI at admission, kg/m2", base$bmi_admission),
  row_continuous("Weight, kg", base$weight_kg),
  row_continuous("SOFA, first window", base$sofa_total),
  row_continuous("Norepinephrine equivalent, mcg/kg/min", base$nee),
  row_continuous("P/F ratio, first window", base$oxygenation),
  row_continuous("Lactate, mmol/L", base$lactate),
  row_continuous(sprintf("Fentanyl dose, first window (%s)", UNITS[["fentanyl"]]),
                 base$total_dose),
  do.call(rbind, lapply(SEDATIVES, function(d) row_continuous(
    sprintf("%s dose, first window (%s)",
            paste0(toupper(substring(d, 1, 1)), substring(d, 2)), UNITS[[d]]),
    base[[DRUGS[[d]]]]))),
  row_continuous("First IMV episode, hours", base$first_imv_episode_hours),
  row_continuous("IMV episodes per block", base$n_imv_episodes)
)
names(baseline) <- c("Characteristic", COLS, "p-value")

stopifnot(
  "the two strata must partition the cohort" =
    sum(grp == "eligible") + sum(grp == "not_eligible") == nrow(base),
  "window 0 must hold every analytic block exactly once" =
    nrow(base) == anchor_n
)

cat("\nBaseline characteristics\n")
print(baseline, row.names = FALSE)


# ---- 13. Figures -------------------------------------------------------------

ink <- "#0b0b0b"; muted <- "#898781"; gridline <- "#e1e0d9"

house <- function(p) {
  p + theme_minimal(base_size = 12) +
    theme(plot.title = element_text(colour = ink, face = "bold"),
          plot.subtitle = element_text(colour = muted, margin = margin(b = 10)),
          legend.position = "top",
          legend.text = element_text(colour = muted, size = 9),
          axis.title = element_text(colour = muted),
          axis.text = element_text(colour = muted),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = gridline, linewidth = 0.4),
          plot.background = element_rect(fill = "#fcfcfb", colour = NA),
          strip.text = element_text(colour = ink, face = "bold"))
}

# FENTANYL IS THE STUDY. It gets the primary figures; propofol and midazolam are
# companions and share one secondary figure.
FENT_U <- UNITS[["fentanyl"]]
DOSE_COLS <- c("All ventilated, zeros included" = "#14427e",
               "Receivers only" = "#4a8bd8")

curves <- dose_summary
curves$series <- ifelse(curves$denominator == "all_ventilated",
                        "All ventilated, zeros included", "Receivers only")

# --- Primary figure: fentanyl, all three curves of section 10 -----------------
# Two stacked FACETS rather than a secondary axis: a dose and a proportion are
# different quantities, and sec_axis applies ONE transform to every facet, which
# draws the proportion against the wrong scale wherever the scales are free.
f <- curves[curves$drug == "fentanyl", ]
f_pct <- f[f$denominator == "all_ventilated", ]

SERIES <- c("Median dose, all ventilated (zeros included)" = "#14427e",
            "Median dose, receivers only"                  = "#4a8bd8",
            "% receiving any fentanyl"                     = "#eb6834")

fent_df <- rbind(
  data.frame(window_start_hr = f$window_start_hr,
             value = f$median, lo = f$q1, hi = f$q3,
             series = ifelse(f$denominator == "all_ventilated",
                             names(SERIES)[1], names(SERIES)[2]),
             quantity = sprintf("Dose (%s), median and IQR", FENT_U)),
  data.frame(window_start_hr = f_pct$window_start_hr,
             value = f_pct$pct_receiving_any, lo = NA_real_, hi = NA_real_,
             series = names(SERIES)[3],
             quantity = "% of ventilated episodes receiving any")
)
fent_df$series <- factor(fent_df$series, levels = names(SERIES))
fent_df$quantity <- factor(fent_df$quantity, levels = unique(fent_df$quantity))

p_fent <- house(
  ggplot(fent_df, aes(window_start_hr, value, colour = series, fill = series)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.13, colour = NA,
                na.rm = TRUE, show.legend = FALSE) +
    geom_line(linewidth = 1.0) + geom_point(size = 1.6) +
    facet_wrap(~ quantity, scales = "free_y", ncol = 1) +
    scale_colour_manual(values = SERIES) +
    scale_fill_manual(values = SERIES) +
    guides(fill = "none",
           colour = guide_legend(nrow = 2, override.aes = list(fill = NA))) +
    labs(title = "Fentanyl dose over the first 72h of ventilation",
         subtitle = paste0(
           "Denominator is episodes STILL VENTILATED in each window.\n",
           "The two medians diverge: exposure narrows to fewer episodes rather than falling within them."),
         x = "Hours since first IMV episode", y = NULL, colour = NULL))

ggsave(file.path(dirs$out_final, "phase1_fentanyl_curves.png"), p_fent,
       width = 7.5, height = 7.6, dpi = 200)

# --- Primary figure: fentanyl balanced panels --------------------------------
# On the MEAN, not the median. Over half of ventilated windows are exactly zero,
# so the median sits on the floor from h24 and every panel collapses onto the
# same line -- it stops discriminating exactly where the cohort starts shrinking
# fastest. The mean uses the zeros without being pinned by them, and a
# proportion has no floor at all.
bp <- balanced_panels
bp$curve <- sprintf("Ventilated >=%dh (n=%s)", bp$panel_hours,
                    trimws(format(bp$panel_n, big.mark = ",")))
allc <- dose_summary[dose_summary$drug == "fentanyl" &
                       dose_summary$denominator == "all_ventilated", ]
allc$curve <- "All ventilated (changing denominator)"

keep <- c("window_start_hr", "mean", "pct_receiving_any", "curve")
pdf_ <- rbind(bp[, keep], allc[, keep])
lv <- c(sort(unique(bp$curve)), "All ventilated (changing denominator)")
pdf_$curve <- factor(pdf_$curve, levels = lv)
pal <- setNames(c("#a8c8ee", "#4a8bd8", "#14427e", "#eb6834")[seq_along(lv)], lv)

panel_df <- rbind(
  data.frame(pdf_[, c("window_start_hr", "curve")], value = pdf_$mean,
             quantity = sprintf("Mean fentanyl dose (%s)", FENT_U)),
  data.frame(pdf_[, c("window_start_hr", "curve")], value = pdf_$pct_receiving_any,
             quantity = "% receiving any fentanyl")
)
panel_df$quantity <- factor(panel_df$quantity, levels = unique(panel_df$quantity))

p_panels <- house(
  ggplot(panel_df, aes(window_start_hr, value, colour = curve)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.3) +
    facet_wrap(~ quantity, scales = "free_y", ncol = 1) +
    scale_colour_manual(values = pal) +
    guides(colour = guide_legend(nrow = 2)) +
    labs(title = "Real dose change, or a changing mix of patients?",
         subtitle = paste0(
           "Panels freeze the denominator at a ventilation duration. Shown on the mean,\n",
           "because the median is pinned at zero from h24 and stops discriminating."),
         x = "Hours since first IMV episode", y = NULL, colour = NULL))

ggsave(file.path(dirs$out_final, "phase1_fentanyl_balanced_panels.png"), p_panels,
       width = 7.5, height = 7.0, dpi = 200)

# --- Primary figure: fentanyl distribution -----------------------------------
# The zero spike is excluded because it is a different kind of observation from
# the continuous part, and it is what decides crimCV. Tail clipped at p99.5 --
# otherwise a handful of extreme windows stretch the axis and the shape that
# matters occupies a tenth of the panel.
fent <- long$total_dose[long$ventilated & !is.na(long$total_dose)]
nz <- fent[fent > 0]
zero_pct <- 100 * mean(fent == 0)
clip <- unname(quantile(nz, 0.995))
n_clipped <- sum(nz > clip)

p_dist <- house(
  ggplot(data.frame(dose = nz[nz <= clip]), aes(dose)) +
    geom_histogram(bins = 60, fill = "#14427e", colour = NA) +
    labs(title = "Fentanyl dose distribution across ventilated windows",
         subtitle = sprintf(
           "Non-zero windows only; %.1f%% of ventilated windows are exactly zero.\nTail clipped at p99.5 = %.0f %s (%s windows above it).",
           zero_pct, clip, FENT_U, trimws(format(n_clipped, big.mark = ","))),
         x = sprintf("Dose (%s)", FENT_U), y = "Windows"))

ggsave(file.path(dirs$out_final, "phase1_fentanyl_distribution.png"), p_dist,
       width = 7.5, height = 4.8, dpi = 200)

# --- Secondary figure: the companion sedatives -------------------------------
# Separate units per drug, so free_y and a label carrying the unit.
sed <- curves[curves$drug %in% SEDATIVES, ]
# Facet order follows the config, not the alphabet, so the panels stay put when
# a drug is added.
sed$facet <- factor(sprintf("%s (%s)", sed$drug, sed$unit),
                    levels = sprintf("%s (%s)", SEDATIVES, UNITS[SEDATIVES]))

# Prevalence per drug, so the caption states which of these are actually used
# here rather than naming one by hand.
zf_w <- zero_fraction[zero_fraction$population == "whole_cohort", ]
prev <- vapply(SEDATIVES, function(d)
  100 - zf_w$pct_zero[zf_w$drug == d], numeric(1))
prev_txt <- paste(sprintf("%s %.1f%%", SEDATIVES, prev), collapse = ", ")
# Wrap by hand: ggplot does not wrap a subtitle, it clips it at the canvas edge.
prev_txt <- paste(strwrap(prev_txt, width = 66), collapse = "\n")

p_sed <- house(
  ggplot(sed, aes(window_start_hr, median, colour = series, fill = series)) +
    geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.13, colour = NA,
                show.legend = FALSE) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.4) +
    facet_wrap(~ facet, scales = "free_y", ncol = 1) +
    scale_colour_manual(values = DOSE_COLS) +
    scale_fill_manual(values = DOSE_COLS) +
    guides(fill = "none", colour = guide_legend(override.aes = list(fill = NA))) +
    labs(title = "Companion sedatives",
         subtitle = sprintf(
           "Secondary to the fentanyl exposure. Infusions only.\nShare of ventilated windows with any drug:\n%s",
           prev_txt),
         x = "Hours since first IMV episode", y = "Median dose (band = IQR)",
         colour = NULL))

ggsave(file.path(dirs$out_final, "phase1_sedative_curves.png"), p_sed,
       width = 7.5, height = 2.2 * length(SEDATIVES) + 2.2, dpi = 200)


# ---- 14. Write ---------------------------------------------------------------

write_out <- function(x, name) {
  f <- file.path(dirs$out_final, name)
  write.csv(x, f, row.names = FALSE)
  cat(sprintf("written: %s\n", name))
}

cat("\n")
write_out(baseline, "phase1_baseline_characteristics.csv")
write_out(retention, "phase1_retention.csv")
write_out(choosing_T, "phase1_choosing_T.csv")
write_out(dose_summary, "phase1_dose_summary.csv")
write_out(dose_distribution, "phase1_dose_distribution.csv")
write_out(balanced_panels, "phase1_balanced_panels.csv")
write_out(zero_fraction, "phase1_zero_fraction.csv")
write_out(imv_episodes, "phase1_imv_episodes.csv")

write_json(prov, file.path(dirs$out_final, "phase1_provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)
cat("written: phase1_provenance.json\n")


# ---- 15. Provenance ----------------------------------------------------------
# Which package versions produced these numbers?

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/02_descriptive_trajectory.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "02_descriptive_trajectory_sessioninfo.txt")
)
