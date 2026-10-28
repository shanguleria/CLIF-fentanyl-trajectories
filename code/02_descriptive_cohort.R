# ==============================================================================
# 02_descriptive_cohort.R  --  cohort description: who, how long, how much fentanyl
#
# Purpose : Dose trajectory over the granular grid with balanced-panel overlays,
#           the retention table that supplies the at-risk denominator, Table 1,
#           and the pooling exports.
# Author  : Shan Guleria
# Created : 2026-09-05
# Inputs  : output/intermediate_phi/trajectory_long.parquet, time_to_event.parquet
# Outputs : output/final_no_phi/02_descriptive/ : the CSVs and figures listed in OWNED
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
source(here("code", "utils", "pooling.R"))
source(here("code", "utils", "figures.R"))


# ---- 2. Config (never setwd(); here() anchors to the .Rproj) -----------------

config <- fromJSON(here("config", "config.json"), simplifyVector = FALSE)
COV       <- fromJSON(here("config", "protocol.json"), simplifyVector = FALSE)
set.seed(COV$model$seed)

WINDOW_H  <- COV$cohort$window_hours
EXTENT_H  <- COV$cohort$granular_extent_hours
PANELS    <- unlist(COV$cohort$balanced_panel_hours)
LANDMARK  <- COV$cohort$landmark_hours

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

# Where a categorical ships in both a full and a collapsed form, the full form is
# published at `overall` only: both views over the same strata subtract.
FULL_SCOPE <- COV$disclosure$full_levels_strata_scope
stopifnot("protocol.json must declare disclosure.full_levels_strata_scope" =
            !is.null(FULL_SCOPE))

# Absolute config-declared bins, never local quantiles -- those cannot be pooled.
HSPEC <- COV$pooling$histograms$variables
stopifnot("protocol.json must declare pooling.histograms.variables" =
            !is.null(HSPEC))
HIST_MINE <- names(HSPEC)[vapply(HSPEC, function(s)
  identical(as.character(s$owner), "02"), logical(1))]

POOL <- list()
POOL_CAT <- list()



# ---- 3. Paths and provenance -------------------------------------------------

dirs <- site_dirs()
dirs$phase <- phase_dir(dirs, "02_descriptive")   # shareable outputs, subdivided by script
prov <- provenance(config)

message(sprintf("[02_descriptive_cohort] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))


# ---- 4. Guards ---------------------------------------------------------------

manifest <- require_manifest(dirs, here())
message(sprintf("  reading Phase 0 outputs from code %s, generated %s",
                manifest$code_version, manifest$generated))

OWNED <- list(phase = c(
  "baseline_characteristics.csv", "retention.csv",
  "dose_summary.csv", "balanced_panels.csv",
  "imv_encounter_duration.csv",
  "pooling_continuous.csv", "pooling_categorical.csv",
  "pooling_histograms.csv",
  "provenance.json",
  "fentanyl_curves.png", "fentanyl_balanced_panels.png",
  "sedative_curves.png", "captions.md"))

# Renames leave stale twins that clear_owned_outputs no longer names; retire them here.
RETIRED <- c(
  file.path("output", "final_no_phi",
            paste0("phase1_", c("dose_curves.png", "pct_receiving.png",
                                "dose_distribution.png", "balanced_panels.png"))),
  file.path("output", "final_no_phi", "02_descriptive",
            paste0("phase1_", c(
              "baseline_characteristics.csv", "retention.csv",
              "dose_summary.csv", "dose_distribution.csv",
              "balanced_panels.csv", "zero_fraction.csv",
              "imv_episodes.csv", "choosing_T.csv",
              "pooling_continuous.csv", "pooling_categorical.csv",
              "provenance.json",
              "fentanyl_curves.png", "fentanyl_balanced_panels.png",
              "fentanyl_distribution.png", "sedative_curves.png",
              "state_prevalence.csv", "state_transitions.csv",
              "state_alluvial.png", "state_prevalence.png",
              "dose_state_prevalence.csv", "dose_state_transitions.csv",
              "dose_state_alluvial.png", "dose_state_prevalence.png"))),
  file.path("output", "final_no_phi", "02_descriptive",
            c("state_prevalence.csv", "state_transitions.csv",
              "dose_state_prevalence.csv", "dose_state_transitions.csv",
              "state_prevalence.png", "state_alluvial.png",
              "dose_state_prevalence.png", "dose_state_alluvial.png")))
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


# ---- 6. Retention ------------------------------------------------------------

retention <- do.call(rbind, lapply(sort(unique(long$window_idx)), function(w) {
  x <- long[long$window_idx == w, ]
  data.frame(window_idx = w,
             window_start_hr = x$window_start_hr[1],
             n_ventilated = sum(x$ventilated),
             n_alive_admitted = sum(x$alive_admitted),
             pct_of_anchor_cohort = round(100 * sum(x$ventilated) / anchor_n, 1))
}))

cat("\nRetention (at-risk = still ventilated)\n")
print(retention, row.names = FALSE)


# ---- 7. Dose summaries -------------------------------------------------------

summarise_dose <- function(v, drug, denom, w, hr) {
  v <- v[!is.na(v)]
  if (!length(v)) return(NULL)
  q <- unname(quantile(v, c(0.25, 0.5, 0.75)))
  # n_receiving is the numerator behind pct_receiving_any, so the proportion
  # pools exactly.
  k <- sum(v > 0)
  # receivers_only is already filtered to v > 0, so the percentage is structurally 100.
  receivers <- identical(denom, "receivers_only")
  data.frame(drug = drug, unit = UNITS[[drug]], denominator = denom,
             window_idx = w, window_start_hr = hr, n = length(v),
             n_receiving = k,
             median = round(q[2], 3), q1 = round(q[1], 3), q3 = round(q[3], 3),
             mean = round(mean(v), 3), sd = round(stats::sd(v), 3),
             pct_receiving_any = if (receivers) NA_real_
                                 else round(100 * k / length(v), 1))
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

TV_POOL <- c(sofa_total = "SOFA", nee = "Norepinephrine equivalent (mcg/kg/min)",
             oxygenation = "P/F ratio", lactate = "Lactate (mmol/L)")

for (drug in names(DRUGS)) {
  col <- DRUGS[[drug]]
  for (w in sort(unique(long$window_idx))) {
    x <- long[long$window_idx == w & long$ventilated, ]
    hr <- x$window_start_hr[1]
    POOL[[length(POOL) + 1]] <- pool_row(
      "by_window", sprintf("%s dose", drug), UNITS[[drug]], "all_ventilated",
      w, hr, x[[col]], MIN_CELL)
    POOL[[length(POOL) + 1]] <- pool_row(
      "by_window", sprintf("%s dose", drug), UNITS[[drug]], "receivers_only",
      w, hr, x[[col]][!is.na(x[[col]]) & x[[col]] > 0], MIN_CELL)
  }
}
for (v in names(TV_POOL)) {
  for (w in sort(unique(long$window_idx))) {
    x <- long[long$window_idx == w & long$ventilated, ]
    POOL[[length(POOL) + 1]] <- pool_row(
      "by_window", TV_POOL[[v]], NA_character_, "ventilated",
      w, x$window_start_hr[1], x[[v]], MIN_CELL)
  }
}

cat("\nFentanyl dose by window (all ventilated, zeros included)\n")
print(dose_summary[dose_summary$drug == "fentanyl" &
                     dose_summary$denominator == "all_ventilated",
                   c("window_start_hr", "n", "median", "q1", "q3",
                     "pct_receiving_any")], row.names = FALSE)


# ---- 9. Balanced panels ------------------------------------------------------
# Frozen denominator: movement means within-patient change. Each panel spans
# only its own window.

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
imv_encounter_duration <- rbind(
  q_row(blocks$first_imv_episode_hours, "first_imv_episode_hours", "hours"),
  q_row(blocks$n_imv_episodes, "n_imv_episodes_per_block", "count"),
  q_row(as.numeric(per_patient), "blocks_per_patient", "count")
)

cat("\nVentilation episode duration and encounter structure\n")
print(imv_encounter_duration, row.names = FALSE)
cat(sprintf("  patients contributing more than one block: %s of %s (%.1f%%)\n",
            format(sum(per_patient > 1), big.mark = ","),
            format(length(per_patient), big.mark = ","),
            100 * mean(per_patient > 1)))


# ---- 12. Baseline characteristics --------------------------------------------
# gtsummary-like shape: __bold__ parent labels, indented sub-levels, N in the headers.

base <- long[long$window_idx == 0, ]

# Overall, Landmark (the outcome subset) and its complement. The complement is what
# makes the test legal. See protocol.json table1._why_the_complement_column.
T1 <- COV$table1
stopifnot("protocol.json must declare table1" = !is.null(T1),
          "Table 1 must be stratified on cohort membership" =
            identical(T1$stratify_on, "cohort_membership"),
          "table1.strata must name exactly two mutually exclusive groups" =
            length(unlist(T1$strata)) == 2)
STRATA <- unlist(T1$strata)
# The declared test names the operation and labels the column; never hardcode either.
T1_TESTS <- list(wilcoxon = list(label = "difference (Wilcoxon)",
                                 fn = function(v, g) stats::wilcox.test(v ~ g, exact = FALSE)$p.value),
                 `chi-square` = list(label = "heterogeneity (chi-square)",
                                     fn = function(tab) stats::chisq.test(tab)$p.value))
stopifnot("table1 declares a test this script cannot run" =
            all(c(T1$test_continuous, T1$test_categorical) %in% names(T1_TESTS)))
T1_CONT <- T1_TESTS[[T1$test_continuous]]
T1_CAT  <- T1_TESTS[[T1$test_categorical]]
base$landmark_eligible <- base$encounter_block %in% elig
grp <- factor(ifelse(base$landmark_eligible, STRATA[1], STRATA[2]),
              levels = STRATA)
cat(sprintf("\nTable 1 stratified on cohort membership (%s of %s blocks reach T=%dh)\n",
            format(sum(base$landmark_eligible), big.mark = ","),
            format(nrow(base), big.mark = ","), LANDMARK))
print(table(grp))

hdr <- function(label, n) sprintf("%s (N = %s)", label, format(n, big.mark = ","))
COLS <- c(hdr("Overall", nrow(base)),
          vapply(STRATA, function(st) hdr(st, sum(grp == st)), character(1)))

fmt_p <- function(p) {
  if (is.na(p)) return("--")
  if (p < 0.001) "<0.001" else sprintf("%.3f", p)
}

# `digits` is per-measure, not per-table; default 1 leaves an unspecified row as it was.
fmt_iqr <- function(v, digits = 1) {
  v <- v[!is.na(v)]
  if (!length(v)) return("--")
  q <- unname(quantile(v, c(0.25, 0.5, 0.75)))
  f <- sprintf("%%.%df (%%.%df, %%.%df)", digits, digits, digits)
  sprintf(f, q[2], q[1], q[3])
}
fmt_pct <- function(k, n) sprintf("%s (%.1f%%)", format(k, big.mark = ","),
                                  100 * k / max(n, 1))

row_continuous <- function(label, v, unit = NA_character_, digits = 1) {
  for (st in c("overall", STRATA)) {
    x <- if (st == "overall") v else v[grp == st]
    POOL[[length(POOL) + 1]] <<- pool_row("baseline", label, unit, st, NA, NA, x, MIN_CELL)
  }
  p <- tryCatch(T1_CONT$fn(v, grp), error = function(e) NA_real_)
  out <- data.frame(characteristic = sprintf("__%s__, median (IQR)", label),
                    overall = fmt_iqr(v, digits), stringsAsFactors = FALSE)
  for (st in STRATA) out[[st]] <- fmt_iqr(v[grp == st], digits)
  out$p_value <- fmt_p(p)
  out$test <- T1_CONT$label
  out
}

row_categorical <- function(label, v, collapse = NULL, pool_raw = TRUE) {
  raw <- ifelse(is.na(v), "Missing", as.character(v))
  # With no collapsed twin this is the only view of itself, so it ships per stratum.
  if (pool_raw) POOL_CAT[[length(POOL_CAT) + 1]] <<-
    pool_cat(label, raw, grp, MIN_CELL,
             strata_scope = if (is.null(collapse)) "all" else FULL_SCOPE)
  v <- if (is.null(collapse)) raw else collapse_levels(raw, collapse)
  if (!identical(v, raw)) POOL_CAT[[length(POOL_CAT) + 1]] <<-
    pool_cat(paste(label, "(collapsed)"), v, grp, MIN_CELL)
  lv <- sort(unique(v))
  tab <- table(v, grp)
  keep <- rownames(tab)[rownames(tab) != "Missing"]
  p <- if (length(keep) > 1 && all(dim(tab[keep, , drop = FALSE]) > 1))
    tryCatch(T1_CAT$fn(tab[keep, , drop = FALSE]), error = function(e) NA_real_)
    else NA_real_
  head <- data.frame(characteristic = sprintf("__%s__, n (%%)", label),
                     overall = "", stringsAsFactors = FALSE)
  for (st in STRATA) head[[st]] <- ""
  head$p_value <- fmt_p(p)
  head$test <- T1_CAT$label
  body <- do.call(rbind, lapply(lv, function(l) {
    r <- data.frame(characteristic = paste0("    ", l),
                    overall = fmt_pct(sum(v == l), length(v)),
                    stringsAsFactors = FALSE)
    for (st in STRATA) r[[st]] <- fmt_pct(sum(v == l & grp == st), sum(grp == st))
    r$p_value <- ""
    r$test <- ""
    r
  }))
  rbind(head, body)
}

baseline <- rbind(
  row_continuous("Age, years", base$age, "years"),
  row_categorical("Sex", base$sex),
  row_categorical("Race", base$race, collapse = RACE_COLLAPSE),
  row_continuous("Charlson Comorbidity Index", base$cci, "index", digits = 0),
  row_continuous("BMI at admission, kg/m2", base$bmi_admission, "kg/m2"),
  row_continuous("Weight, kg", base$weight_kg, "kg"),
  row_continuous("SOFA, first window", base$sofa_total, "points"),
  row_continuous("Norepinephrine equivalent, mcg/kg/min", base$nee, "mcg/kg/min",
                 digits = 2),
  row_continuous("P/F ratio, first window", base$oxygenation, "mmHg", digits = 0),
  row_continuous("Lactate, mmol/L", base$lactate, "mmol/L"),
  row_continuous(sprintf("Fentanyl dose, first window (%s)", UNITS[["fentanyl"]]),
                 base$total_dose, UNITS[["fentanyl"]]),
  do.call(rbind, lapply(SEDATIVES, function(d) row_continuous(
    sprintf("%s dose, first window (%s)",
            paste0(toupper(substring(d, 1, 1)), substring(d, 2)), UNITS[[d]]),
    base[[DRUGS[[d]]]], UNITS[[d]]))),
  row_continuous("First IMV episode, hours", base$first_imv_episode_hours, "hours")
)
names(baseline) <- c("Characteristic", COLS, "p-value", "test")

stopifnot(
  "the strata must partition the cohort" = sum(table(grp)) == nrow(base),
  "every stratum must be non-empty -- check that time_to_event.parquet loaded" =
    all(table(grp) > 0),
  "the Landmark column must match the landmark cohort exactly" =
    sum(grp == STRATA[1]) == length(intersect(elig, base$encounter_block)),
  "window 0 must hold every analytic block exactly once" =
    nrow(base) == anchor_n
)

cat("\nBaseline characteristics\n")
print(baseline, row.names = FALSE)


# ---- 12b. Pooling histograms -------------------------------------------------
# Counts on fixed absolute bins pool where a median does not. See protocol.json
# pooling.histograms.
pool_histograms <- do.call(rbind, lapply(HIST_MINE, function(v) {
  spec <- HSPEC[[v]]
  col <- sub("^trajectory_long\\.", "", as.character(spec$source))
  stopifnot("a declared histogram names a column trajectory_long does not carry" =
              col %in% names(long))
  x <- if (identical(as.character(spec$restrict), "one row per episode")) {
    # One value per episode: off the window grid unfiltered it would weight long stays.
    long[[col]][!duplicated(long$encounter_block)]
  } else {
    long[[col]][long$ventilated]
  }
  pool_hist(v, x, spec)
}))
stopifnot(
  "02 must emit every histogram the config assigns it" =
    setequal(unique(pool_histograms$variable), HIST_MINE))
cat(sprintf("\nPooling histograms: %d variables, %d bins total\n",
            length(HIST_MINE), nrow(pool_histograms)))


# ---- 13. Figures -------------------------------------------------------------

FENT_U <- UNITS[["fentanyl"]]
DOSE_COLS <- c("All ventilated, zeros included" = "#14427e",
               "Receivers only" = "#4a8bd8")

curves <- dose_summary
curves$series <- ifelse(curves$denominator == "all_ventilated",
                        "All ventilated, zeros included", "Receivers only")

# --- Primary figure: fentanyl, all three dose curves -------------------------
# Stacked facets, not sec_axis: sec_axis applies one transform to every facet and
# would draw the proportion against the wrong scale under free scales.
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
    labs(x = "Hours since first IMV episode", y = NULL, colour = NULL))

register_caption("fentanyl_curves.png",
  "Fentanyl dose over the first 72h of ventilation",
  paste0("Denominator is episodes STILL VENTILATED in each window. The two ",
         "medians diverge: exposure narrows to fewer episodes rather than ",
         "falling within them."))
ggsave(file.path(dirs$phase, "fentanyl_curves.png"), p_fent,
       width = 7.5, height = 7.6, dpi = 200)

# --- Primary figure: fentanyl balanced panels --------------------------------
# Mean as well as median: most ventilated windows are exactly zero, so the
# median floors.
bp <- balanced_panels
bp$curve <- sprintf("Ventilated >=%dh (n=%s)", bp$panel_hours,
                    trimws(format(bp$panel_n, big.mark = ",")))
allc <- dose_summary[dose_summary$drug == "fentanyl" &
                       dose_summary$denominator == "all_ventilated", ]
allc$curve <- "All ventilated (changing denominator)"

keep <- c("window_start_hr", "mean", "median", "pct_receiving_any", "curve")
pdf_ <- rbind(bp[, keep], allc[, keep])
lv <- c(sort(unique(bp$curve)), "All ventilated (changing denominator)")
pdf_$curve <- factor(pdf_$curve, levels = lv)
pal <- setNames(c("#a8c8ee", "#4a8bd8", "#14427e", "#eb6834")[seq_along(lv)], lv)

panel_df <- rbind(
  data.frame(pdf_[, c("window_start_hr", "curve")], value = pdf_$mean,
             quantity = sprintf("MEAN dose (%s)", FENT_U)),
  data.frame(pdf_[, c("window_start_hr", "curve")], value = pdf_$median,
             quantity = sprintf("MEDIAN dose (%s)", FENT_U)),
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
    labs(x = "Hours since first IMV episode", y = NULL, colour = NULL))

register_caption("fentanyl_balanced_panels.png",
  "Real dose change, or a changing mix of patients?",
  paste0("Panels freeze the denominator at a ventilation duration. Mean and ",
         "median are shown together: over half of ventilated windows are ",
         "exactly zero, so the median falls onto the floor at h24 and every ",
         "panel collapses onto one line."))
ggsave(file.path(dirs$phase, "fentanyl_balanced_panels.png"), p_panels,
       width = 7.5, height = 9.2, dpi = 200)

# --- Secondary figure: the companion sedatives -------------------------------
sed <- curves[curves$drug %in% SEDATIVES, ]
# Facet order follows the config, not the alphabet.
sed$facet <- factor(sprintf("%s (%s)", sed$drug, sed$unit),
                    levels = sprintf("%s (%s)", SEDATIVES, UNITS[SEDATIVES]))

prev <- vapply(SEDATIVES, function(d) {
  v <- long[[DRUGS[[d]]]][long$ventilated]
  100 * mean(v > 0, na.rm = TRUE)
}, numeric(1))
prev_txt <- paste(sprintf("%s %.1f%%", SEDATIVES, prev), collapse = ", ")
# Wrap by hand: ggplot does not wrap a subtitle, it clips it at the canvas edge.
prev_txt <- paste(strwrap(prev_txt, width = 66), collapse = "\n")

sed_long <- rbind(
  data.frame(sed[, c("window_start_hr", "series", "facet")],
             value = sed$median, statistic = "median"),
  data.frame(sed[, c("window_start_hr", "series", "facet")],
             value = sed$mean, statistic = "mean")
)
sed_long$statistic <- factor(sed_long$statistic, levels = c("mean", "median"))

p_sed <- house(
  ggplot(sed_long, aes(window_start_hr, value, colour = series)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.2) +
    facet_grid(facet ~ statistic, scales = "free_y", switch = "y") +
    scale_colour_manual(values = DOSE_COLS) +
    labs(x = "Hours since first IMV episode", y = NULL, colour = NULL) +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90)))

register_caption("sedative_curves.png", "Companion sedatives",
  sprintf(paste0("Secondary to the fentanyl exposure. Infusions only. Share of ",
                 "ventilated windows with any drug: %s"),
          gsub("\n", "; ", trimws(prev_txt))))
ggsave(file.path(dirs$phase, "sedative_curves.png"), p_sed,
       width = 7.5, height = 2.2 * length(SEDATIVES) + 2.2, dpi = 200)


# ---- Write -------------------------------------------------------------------

write_out <- function(x, name) {
  f <- file.path(dirs$phase, name)
  write.csv(x, f, row.names = FALSE)
  cat(sprintf("written: %s\n", name))
}

cat("\n")
write_out(baseline, "baseline_characteristics.csv")
write_out(retention, "retention.csv")
write_out(dose_summary, "dose_summary.csv")
write_out(balanced_panels, "balanced_panels.csv")
write_out(imv_encounter_duration, "imv_encounter_duration.csv")

pooling_continuous <- do.call(rbind, POOL)
pooling_categorical <- do.call(rbind, POOL_CAT)
write_out(pooling_continuous, "pooling_continuous.csv")
write_out(pooling_categorical, "pooling_categorical.csv")
write_out(pool_histograms, "pooling_histograms.csv")

write_captions(file.path(dirs$phase, "captions.md"), "02_descriptive_cohort.R",
               grep("\\.png$", OWNED$phase, value = TRUE), prov)
cat("written: captions.md\n")

write_json(prov, file.path(dirs$phase, "provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)
cat("written: provenance.json\n")


# ---- Provenance --------------------------------------------------------------

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/02_descriptive_cohort.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "02_descriptive_cohort_sessioninfo.txt")
)
