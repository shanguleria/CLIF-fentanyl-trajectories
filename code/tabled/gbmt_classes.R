# ==============================================================================
# gbmt_classes.R  --  Phase 3 -- gbmt classes on combined dose
#
# Purpose : Group-based trajectory model on total fentanyl dose among episodes that received ANY fentanyl in [0, T]; sweep ng and select by the full criterion conjunction.
# Author  : Shan Guleria
# Created : 2026-09-07
# Inputs  : output/intermediate_phi/landmark_cohort.parquet
# Outputs : output/intermediate_phi/gbmt_fits.rds, gbmt_group_assignments.csv; IC table, scree plot, trajectory figure
#
# ESTIMAND: trajectory classes among ventilation episodes that received any
# fentanyl while intubated in [0, T]. Episodes with zero dose in every window
# are a named stratum, not a class -- see config model.min_exposed_windows.
#
# ==============================================================================

# Run this in a FRESH R session (RStudio: Cmd+Shift+F10).
# Do not use rm(list = ls()) -- it does not unload packages or reset options,
# so it only gives the appearance of a clean slate.


# ---- 1. Packages -------------------------------------------------------------

pkgs <- c("here", 
          "jsonlite", 
          "gbmt", 
          "arrow",
          "tidyverse",
          "ggplot2")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

source(here("code", "utils", "paths.R"))
source(here("code", "utils", "pooling.R"))


# ---- 2. Config (never setwd(); here() anchors to the .Rproj) -----------------

config <- fromJSON(here("config", "config.json"), simplifyVector = FALSE)
set.seed(config$model$seed)


# ---- 3. Paths and provenance -------------------------------------------------
# One site, one output tree. site_dirs() creates them and labels the PHI ones.

dirs <- site_dirs()
dirs$phase <- phase_dir(dirs, "gbmt")   # shareable outputs, subdivided by script
prov <- provenance(config)

message(sprintf("[gbmt_classes] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))


# ---- 4. Guards ---------------------------------------------------------------
# Confirm that the Phase 0 data on file is consistent with most recent config

manifest <- require_manifest(dirs, here())
message(sprintf("  reading Phase 0 outputs from code %s, generated %s",
                manifest$code_version, manifest$generated))

# gbmt_fits.rds is deliberately NOT cleared since it's so computationally expensive
# gbmt_fits_provenance.json IS cleared, which makes the pairing checkable:
#   sidecar missing        -> the fits are from an interrupted or older run
#   digests disagree       -> the fits predate the current config
OWNED <- list(
  out_phi   = c("gbmt_group_assignments.csv", "gbmt_ic_comparison.csv",
                "gbmt_fits_provenance.json"),
  phase = c("gbmt_exposure_stratum.csv", "gbmt_model_selection.csv",
                "bic_ng_plot.png",
                "gbmt_trajectories.png"))
# Taken off the number line on 2026-09-24: the live pipeline reclaimed 04 and
# 05 when the exemplar moved ahead of the states script, so two scripts were
# claiming each number. These are parked, not part of the sequence, and their
# folder now says so. The numbered folder is retired -- built with paste0 so a
# name-level find-and-replace cannot reach the old names (lessons.md #13).
# Existing results were MOVED, not discarded; this only clears the twin a
# re-run would otherwise leave behind.
RETIRED <- file.path("output", "final_no_phi", paste0("04_", "gbmt"),
                     OWNED$phase)
n_cleared <- clear_owned_outputs(dirs, OWNED, retired = RETIRED)
if (n_cleared) message(sprintf("  cleared %d output(s) from a previous run", n_cleared))


# ---- 5. gbmt ANALYSIS part 1 ----
# read from  : dirs$out_phi
# PHI out    : dirs$out_phi
# aggregate  : dirs$phase   (output/final_no_phi/gbmt/)

## ---- Data ----
# Read landmark_cohort parquet as df foy gbmt
df_all <- as.data.frame(read_parquet(file.path(dirs$out_phi, "landmark_cohort.parquet")))
df_all <- df_all[, c("id_num", "window_idx", "total_dose")]

## ---- Exposure stratum ----
# Episodes with zero dose in EVERY window have no trajectory to model -- their
# within-group variance is exactly zero, which makes the Gaussian mixture
# likelihood unbounded and the IC comparison invalid across ng. Hold them out as
# a named stratum and fit classes among the exposed. Threshold from config.

MIN_EXPOSED <- config$model$min_exposed_windows
DEG <- config$model$polynomial_degree
stopifnot("config model.min_exposed_windows is missing" =
            "min_exposed_windows" %in% names(config$model))

n_exposed_windows <- tapply(df_all$total_dose > 0, df_all$id_num, sum)
exposed_ids  <- as.integer(names(n_exposed_windows)[n_exposed_windows >= MIN_EXPOSED])
never_ids    <- as.integer(names(n_exposed_windows)[n_exposed_windows <  MIN_EXPOSED])

stratum <- data.frame(
  stratum = c("exposed (modelled)", "never exposed (stratum, not modelled)"),
  definition = c(sprintf(">= %d of %d windows with dose > 0", MIN_EXPOSED,
                         length(unique(df_all$window_idx))),
                 sprintf("< %d windows with dose > 0", MIN_EXPOSED)),
  n = c(length(exposed_ids), length(never_ids)),
  pct = round(100 * c(length(exposed_ids), length(never_ids)) /
                length(n_exposed_windows), 1))
print(stratum, row.names = FALSE)
write.csv(stratum, file.path(dirs$phase, "gbmt_exposure_stratum.csv"),
          row.names = FALSE)

# Filter out encounters with no fentanyl exposure over entire landmark window
df <- df_all[df_all$id_num %in% exposed_ids, ]
stopifnot("the exposed set is empty" = nrow(df) > 0,
          "the split must partition the cohort" =
            length(exposed_ids) + length(never_ids) == length(n_exposed_windows))
message(sprintf("fitting on %s exposed episodes; %s never-exposed held out",
                format(length(exposed_ids), big.mark = ","),
                format(length(never_ids), big.mark = ",")))

## ---- Fit Models ----
# Expected structure: 
# gbmt(x.names, unit, time, ng=1, d=2, data,  
#      scaling=0, pruning=TRUE, delete.empty=FALSE,
#      nstart=NULL, tol=1e-4, maxit=1000, quiet=FALSE)

# Read the sweep range FROM the config rather than restating it here. A
# hardcoded range is how config and run drift apart: model.ng_range said [1, 6]
# while this line said 2:5, so the config described a sweep that never happened.
ng_range <- unlist(config$model$ng_range)
stopifnot("config model.ng_range is missing" =
            "ng_range" %in% names(config$model),
          "model.ng_range must be [min, max] with min >= 2" =
            length(ng_range) == 2 && ng_range[1] >= 2 && ng_range[2] >= ng_range[1])
n_groups <- seq(ng_range[1], ng_range[2])
message(sprintf("sweeping ng = %s (config model.ng_range)",
                paste(range(n_groups), collapse = ":")))

# Storing the fits in a named list rather than individual objects
fits <- lapply(n_groups, function(g) {
  message("Fitting ng = ", g, " ...")
  gbmt(x.names = "total_dose", 
       unit = "id_num", 
       time = "window_idx",
       d = DEG,     # config$model$polynomial_degree
       ng = g, 
       data = df, 
       scaling = 0, # Scaling must be 0 for this project
       nstart = NULL,
       quiet = TRUE)
})
names(fits) <- paste0("ng", n_groups)

# Save fits
saveRDS(fits, file.path(dirs$out_phi, "gbmt_fits.rds"))

# Written only after the save succeeds, so its presence certifies the .rds. It
# records everything that would make those fits stale.
write_json(list(
  saved            = format(Sys.time(), tz = config$timezone, usetz = TRUE),
  code_version     = prov$code_version,
  phase0_code      = manifest$code_version,
  phase0_generated = manifest$generated,
  config_digests   = manifest$config_digests,
  ng_range            = ng_range,
  min_exposed_windows = MIN_EXPOSED,
  scaling             = config$model$scaling,
  polynomial_degree   = DEG,
  seed                = config$model$seed,
  n_units          = length(exposed_ids),
  n_windows        = length(unique(df$window_idx)),
  models           = names(fits)),
  file.path(dirs$out_phi, "gbmt_fits_provenance.json"),
  auto_unbox = TRUE, pretty = TRUE)

message("saved gbmt_fits.rds (", length(fits), " fits) + provenance sidecar")

## ---- Compare and choose models ----
# $ic holds the information criteria; 
# $appa is the average posterior probability of assignment 
#  (a common rule of thumb wants APPA > 0.7)

ic_table <- do.call(rbind, lapply(names(fits), function(nm) {
  f <- fits[[nm]]
  data.frame(model    = nm,
             ng       = length(f$prior),
             t(f$ic),
             appa_min = min(f$appa),   # weakest group is what you report
             smallest_group = min(f$prior))
}))
stopifnot(nrow(ic_table) == length(fits))

print(ic_table)

write.csv(ic_table,
          file.path(dirs$out_phi, "gbmt_ic_comparison.csv"),
          row.names = FALSE)

# Chen et al. layout: statistics as rows, one column per candidate ng, and a
# %classK row per class so every solution's class sizes are visible at once
# rather than only the smallest. Follows BMC Infect Dis 2026;26:950 Table 1.
# Aggregate, so this one is shareable.
sel <- selection_table(
  stats = ic_table[, c("ng", "aic", "bic", "caic", "ssbic", "hqic", "appa_min")],
  sizes = lapply(fits, function(f) as.integer(table(f$assign))))
print(sel, row.names = FALSE)
write.csv(sel, file.path(dirs$phase, "gbmt_model_selection.csv"), row.names = FALSE)

# Plot BIC to visualize elbow

best_ic <- ic_table[which.min(ic_table$bic), ]

blue <- "#2a78d6"; accent <- "darkblue"
ink  <- "#0b0b0b"; muted  <- "#898781"; gridline <- "#e1e0d9"

bic_ng_plot <- ggplot (ic_table, aes(ng, bic)) +
  geom_line(linewidth = 0.7, color = blue) +
  geom_point(size = 2.6, color = blue) +
  geom_point(data = best_ic, size = 4.4, shape = 21,
             fill = accent, color = "#fcfcfb", stroke = 1.1) +
  geom_text(data = best_ic, aes(label = sprintf("ng = %d\nBIC = %.0f", ng, bic)),
            vjust = 1.9, size = 3.4, color = ink, lineheight = 1.1) +
  scale_x_continuous(breaks = ic_table$ng) +
  scale_y_continuous(expand = expansion(mult = c(0.22, 0.08))) +
  labs(title    = "BIC by number of trajectory groups",
       x = "Number of groups (ng)", y = "BIC") +
  theme_bw(base_size=12) + theme(plot.title         = element_text(face = "bold", color = ink),
                                 plot.subtitle      = element_text(color = muted, margin = margin(b = 12)),
                                 axis.title         = element_text(color = muted),
                                 axis.text          = element_text(color = muted),
                                 panel.grid.minor   = element_blank(),
                                 panel.grid.major.x = element_blank(),
                                 panel.grid.major.y = element_line(color = gridline, linewidth = 0.4),
                                 plot.background    = element_rect(fill = "#fcfcfb", color = NA))

ggsave(file.path(dirs$phase, "bic_ng_plot.png"), bic_ng_plot,
       width = 6.5, height = 4.2, dpi = 200)

# ---- 6. gbmt ANALYSIS part 2 ----

## ---- Inspect chosen model ----
chosen_model <- fits[["ng4"]]   # insert chosen model here

# Summarize chosen model's assignments
print(data.frame(group = names(table(chosen_model$assign)),
                 n     = as.integer(table(chosen_model$assign)),
                 pct   = round(100 * as.numeric(chosen_model$prior), 1),
                 appa  = round(chosen_model$appa, 3), row.names = NULL))

# Every landmark episode gets a level: a trajectory class, or the named stratum.
# Phase 6 needs all of them, with "never exposed" the natural reference.
assignments <- rbind(
  data.frame(unit  = as.integer(names(chosen_model$assign)),
             group = as.character(chosen_model$assign),
             stratum = "exposed", row.names = NULL),
  data.frame(unit = never_ids, group = "never_exposed",
             stratum = "never exposed", row.names = NULL))
stopifnot("assignments must cover every landmark episode" =
            nrow(assignments) == length(n_exposed_windows))

write.csv(assignments,
          file.path(dirs$out_phi, "gbmt_group_assignments.csv"),
          row.names = FALSE)

## ---- Plot chosen model ----
# Fitted class curve with the observed group mean overlaid

fitted_curves <- function(fit, window_h) {
  do.call(rbind, lapply(names(fit$fitted), function(g) {
    y <- fit$fitted[[g]][[1]]
    data.frame(group = g, window_idx = seq_along(y) - 1L,
               hr = (seq_along(y) - 1L) * window_h, fitted = y)
  }))
}

# observed mean per group per window pulled from the analytic table
observed_curves <- function(dat, assign, window_h, col = "total_dose") {
  a <- data.frame(id_num = as.integer(names(assign)),
                  group  = as.character(assign))
  m <- merge(dat, a, by = "id_num")
  agg <- aggregate(m[[col]], by = list(group = m$group, window_idx = m$window_idx),
                   FUN = mean, na.rm = TRUE)
  names(agg)[3] <- "observed"
  agg$hr <- agg$window_idx * window_h
  agg
}

plot_trajectories <- function(fit, dat, window_h, unit_label, title, subtitle) {
  fc <- fitted_curves(fit, window_h)
  oc <- observed_curves(dat, fit$assign, window_h)
  n  <- table(fit$assign)

  # order classes by level so the ramp reads low -> high
  lvl  <- names(sort(tapply(fc$fitted, fc$group, mean)))
  labs <- sprintf("Class %s  (n=%s, %.1f%%)", lvl,
                  format(as.integer(n[lvl]), big.mark = ","),
                  100 * as.integer(n[lvl]) / sum(n))
  fc$group <- factor(fc$group, levels = lvl, labels = labs)
  oc$group <- factor(oc$group, levels = lvl, labels = labs)
  pal <- setNames(colorRampPalette(c("#a8c8ee", "#14427e"))(length(lvl)), labs)

  ggplot(fc, aes(hr, fitted, colour = group)) +
    geom_line(data = oc, aes(y = observed), linetype = "22", linewidth = 0.6,
              alpha = 0.8, show.legend = FALSE) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 1.5) +
    scale_colour_manual(values = pal) +
    guides(colour = guide_legend(nrow = 2)) +
    labs(title = title, subtitle = subtitle,
         x = "Hours since first IMV episode",
         y = sprintf("Fentanyl dose (%s)", unit_label), colour = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title         = element_text(colour = ink, face = "bold"),
          plot.subtitle      = element_text(colour = muted, margin = margin(b = 10)),
          legend.position    = "top",
          legend.text        = element_text(colour = muted, size = 9),
          axis.title         = element_text(colour = muted),
          axis.text          = element_text(colour = muted),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = gridline, linewidth = 0.4),
          plot.background    = element_rect(fill = "#fcfcfb", colour = NA))
}

COV <- fromJSON(here("config", "covariates.json"), simplifyVector = FALSE)

traj_plot <- plot_trajectories(
  chosen_model, df,
  window_h   = config$cohort$window_hours,
  unit_label = COV$exposure$units,
  title      = "Fentanyl trajectory classes over the landmark window",
  subtitle   = sprintf(
    "%s episodes that received any fentanyl while intubated in [0, %dh]; %s never-exposed held out.\nSolid = fitted class curve, dashed = observed group mean.",
    format(length(chosen_model$assign), big.mark = ","),
    config$cohort$landmark_hours,
    format(length(never_ids), big.mark = ",")))

ggsave(file.path(dirs$phase, "gbmt_trajectories.png"), traj_plot,
       width = 7.5, height = 5.4, dpi = 200)

# ---- 7. Provenance -----------------------------------------------------------
# Which package versions produced these numbers?

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/tabled/gbmt_classes.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "gbmt_classes_sessioninfo.txt")
)
