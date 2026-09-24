# ==============================================================================
# 05_lcmm_classes.R  --  Phase 4 -- lcmm classes on the same data, ARI vs Phase 3
#
# Purpose : Latent-class mixed model on total fentanyl dose over the identical landmark panel Phase 3 used; sweep ng, evaluate by the same criterion conjunction, and compare the partition to gbmt by adjusted Rand index.
# Author  : Shan Guleria
# Created : 2026-09-08
# Inputs  : output/intermediate_phi/landmark_cohort.parquet, gbmt_group_assignments.csv
# Outputs : output/intermediate_phi/lcmm_fits.rds, lcmm_group_assignments.csv; IC table, ARI table, trajectory figure
#
# ESTIMAND: trajectory classes among ventilation episodes that received any 
# fentanyl while intubated in [0, T].
#
# METHOD: heterogeneous linear mixed effects (hlme) from lcmm package to identify latent classes
# Contains fixed effects for shared within-class features and random effects for individual variation
# ==============================================================================

# Run this in a FRESH R session (RStudio: Cmd+Shift+F10).
# Do not use rm(list = ls()) -- it does not unload packages or reset options,
# so it only gives the appearance of a clean slate.

# ---- 1. Packages ----

pkgs <- c("here", 
          "jsonlite", 
          "arrow", 
          "lcmm",
          "tidyverse",
          "ggplot2")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

source(here("code", "utils", "paths.R"))

# ---- 2. Config (never setwd(); here() anchors to the .Rproj) ----

config <- fromJSON(here("config", "config.json"), simplifyVector = FALSE)
set.seed(config$model$seed)

COV      <- fromJSON(here("config", "covariates.json"), simplifyVector = FALSE)
WINDOW_H <- config$cohort$window_hours
LANDMARK <- config$cohort$landmark_hours
DEG      <- config$model$polynomial_degree
MIN_EXPOSED <- config$model$min_exposed_windows

stopifnot("config model.ng_range is missing" =
            "ng_range" %in% names(config$model),
          "config model.min_exposed_windows is missing" =
            "min_exposed_windows" %in% names(config$model))
ng_range <- unlist(config$model$ng_range)
n_groups <- seq(ng_range[1], ng_range[2])

# ---- 3. Paths and provenance ----

dirs <- site_dirs()
dirs$phase <- phase_dir(dirs, "05_lcmm")   # shareable outputs, subdivided by script
prov <- provenance(config)

message(sprintf("[05_lcmm_classes] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))

# ---- 4. Guards ----

# Confirm that the Phase 0 data on file is consistent with most recent config
manifest <- require_manifest(dirs, here())
message(sprintf("  reading Phase 0 outputs from code %s, generated %s",
                manifest$code_version, manifest$generated))

OWNED <- list(
  out_phi   = c("lcmm_fits.rds", "lcmm_group_assignments.csv"),
  phase = c("lcmm_ic_comparison.csv", "lcmm_vs_gbmt_ari.csv",
                "lcmm_trajectories.png", "lcmm_bic_ng_plot.png",
                "lcmm_provenance.json"))
n_cleared <- clear_owned_outputs(dirs, OWNED)
if (n_cleared) message(sprintf("  cleared %d output(s) from a previous run", n_cleared))

# ---- 5. Data ----

df_all <- as.data.frame(read_parquet(file.path(dirs$out_phi, "landmark_cohort.parquet")))
df_all <- df_all[, c("id_num", "window_idx", "total_dose")]


# Filter out encounters with no fentanyl exposure over entire landmark window
n_exposed_windows <- tapply(df_all$total_dose > 0, df_all$id_num, sum)
exposed_ids <- as.integer(names(n_exposed_windows)[n_exposed_windows >= MIN_EXPOSED])
never_ids   <- as.integer(names(n_exposed_windows)[n_exposed_windows <  MIN_EXPOSED])

df <- df_all[df_all$id_num %in% exposed_ids, ]
df$t <- df$window_idx          # hlme wants the time variable by name in a formula

stopifnot("the exposed set is empty" = nrow(df) > 0,
          "panel must stay balanced" = length(unique(table(df$id_num))) == 1)
message(sprintf("  fitting on %s exposed episodes x %d windows",
                format(length(exposed_ids), big.mark = ","),
                length(unique(df$window_idx))))

# Cross-check against gbmt cohort so a silent cohort drift cannot pass.
gb_file <- file.path(dirs$out_phi, "gbmt_group_assignments.csv")
gb <- if (file.exists(gb_file)) read.csv(gb_file) else NULL
if (!is.null(gb)) {
  gb_exposed <- gb$unit[gb$group != "never_exposed"]
  stopifnot("lcmm and gbmt must fit the SAME episodes, or the ARI is meaningless" =
              setequal(as.integer(gb_exposed), exposed_ids))
  message("  cohort matches gbmt run exactly")
}

# ---- 6. hlme ANALYSIS part 1 ----

## ---- Fit ----
# hlme default pattern fits ng=1 model and derive every ng > 1 start from it (B = m1)

# fixed   = the mean trajectory shape
# mixture = which of those terms are allowed to differ BY CLASS (all of them)
# random  = the per-episode deviation. THIS is what gbmt does not have.

f_fixed  <- as.formula(sprintf("total_dose ~ poly(t, %d, raw = TRUE)", DEG))
  #DEG from config$model$polynomial_degree
f_mix    <- as.formula(sprintf("~ poly(t, %d, raw = TRUE)", DEG))
f_random <- ~ 1                     # random intercept; see the note below

message("Fitting ng = 1 (base model for starting values) ...")
t0 <- Sys.time()
m1 <- do.call(hlme, list(fixed = f_fixed, random = f_random, subject = "id_num",
                         ng = 1, data = quote(df), verbose = FALSE))
message(sprintf("  ng = 1 done in %.1f min", as.numeric(
  difftime(Sys.time(), t0, units = "mins"))))

fits <- list()
for (g in n_groups) {
  message("Fitting ng = ", g, " ...")
  t0 <- Sys.time()
  fits[[paste0("ng", g)]] <- do.call(hlme, list(
    fixed = f_fixed, mixture = f_mix, random = f_random,
    subject = "id_num", ng = g, data = quote(df), B = m1, verbose = FALSE))
  message(sprintf("  ng = %d done in %.1f min", g,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# Save as soon as they exist, PHI
saveRDS(list(m1 = m1, fits = fits), file.path(dirs$out_phi, "lcmm_fits.rds"))
message("saved lcmm_fits.rds")


## ---- Compare models ----
# hlme reports the posterior matrix in $pprob

class_of  <- function(m) setNames(m$pprob[[2]], as.character(m$pprob[[1]]))
postmat   <- function(m) as.matrix(m$pprob[, -(1:2), drop = FALSE])

criteria <- function(m, label) {
  P <- postmat(m); K <- ncol(P); cl <- class_of(m)
  prior <- as.numeric(table(factor(cl, levels = seq_len(K)))) / length(cl)
  appa  <- vapply(seq_len(K), function(k) mean(P[cl == k, k]), numeric(1))
  occ   <- (appa / (1 - appa)) / (prior / (1 - prior))
  pl    <- P * log(pmax(P, .Machine$double.eps))
  data.frame(model = label, ng = K,
             loglik = round(m$loglik, 1), npar = length(m$best),
             aic = round(m$AIC, 1), bic = round(m$BIC, 1),
             entropy = round(1 - (-sum(pl)) / (nrow(P) * log(K)), 4),
             appa_min = round(min(appa), 4),
             occ_min  = round(min(occ), 1),
             smallest = round(min(prior), 4),
             converged = m$conv)
}

ic_table <- do.call(rbind, lapply(names(fits), function(nm) criteria(fits[[nm]], nm)))
ic_table$dBIC <- c(NA, round(diff(ic_table$bic)))
ic_table$passes_section7 <- with(ic_table,
  entropy > 0.80 & appa_min >= 0.70 & occ_min > 5.0 & smallest >= 0.05)

print(ic_table, row.names = FALSE)

# conv == 1 is convergence; anything else means the fit did not settle and its
# information criteria are not interpretable.
if (any(ic_table$converged != 1)) {
  warning("these fits did NOT converge and must not be selected: ",
          paste(ic_table$model[ic_table$converged != 1], collapse = ", "))
}

write.csv(ic_table, file.path(dirs$phase, "lcmm_ic_comparison.csv"),
          row.names = FALSE)

# Plot BIC to visualize elbow
ink <- "#0b0b0b"; muted <- "#898781"; gridline <- "#e1e0d9"

bic_plot <- ggplot(ic_table, aes(ng, bic)) +
  geom_line(linewidth = 0.7, colour = "#2a78d6") +
  geom_point(size = 2.6, colour = "#2a78d6") +
  scale_x_continuous(breaks = ic_table$ng) +
  labs(title = "BIC by number of trajectory groups -- lcmm",
       x = "Number of groups (ng)", y = "BIC") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(colour = ink, face = "bold"),
        axis.title = element_text(colour = muted),
        axis.text = element_text(colour = muted),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = gridline, linewidth = 0.4),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA))

ggsave(file.path(dirs$phase, "lcmm_bic_ng_plot.png"), bic_plot,
       width = 6.5, height = 4.2, dpi = 200)

# ---- 7. hlme ANALYSIS part 2 ----

## ---- Choosing a model/ng ----
CHOSEN <- "ng2"                       # <-- set after deciding on ng
chosen_model <- fits[[CHOSEN]]
stopifnot("CHOSEN is not one of the fitted models" = !is.null(chosen_model))

assignments <- rbind(
  data.frame(unit = as.integer(names(class_of(chosen_model))),
             group = as.character(class_of(chosen_model)),
             stratum = "exposed", row.names = NULL),
  data.frame(unit = never_ids, group = "never_exposed",
             stratum = "never exposed", row.names = NULL))
stopifnot("assignments must cover every landmark episode" =
            nrow(assignments) == length(n_exposed_windows))
write.csv(assignments, file.path(dirs$out_phi, "lcmm_group_assignments.csv"),
          row.names = FALSE)

# Adjusted Rand index (ARI): agreement between two partitions of the same units,
# 1 = identical, 0 = chance.
ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  s  <- sum(choose(tab, 2))
  a1 <- sum(choose(rowSums(tab), 2)); b1 <- sum(choose(colSums(tab), 2))
  e  <- a1 * b1 / choose(n, 2)
  (s - e) / ((a1 + b1) / 2 - e)
}

if (!is.null(gb)) {
  gbx <- gb[gb$group != "never_exposed", ]
  ari_table <- do.call(rbind, lapply(names(fits), function(nm) {
    cl <- class_of(fits[[nm]])
    k  <- intersect(names(cl), as.character(gbx$unit))
    g  <- setNames(as.character(gbx$group), as.character(gbx$unit))
    data.frame(lcmm_model = nm, lcmm_ng = fits[[nm]]$ng,
               gbmt_reference = "chosen gbmt partition",
               n_units = length(k),
               ari = round(ari(cl[k], g[k]), 4))
  }))
  print(ari_table, row.names = FALSE)
  write.csv(ari_table, file.path(dirs$phase, "lcmm_vs_gbmt_ari.csv"),
            row.names = FALSE)
}


# ---- 8. Figure ----

pred <- predictY(chosen_model,
                 newdata = data.frame(t = sort(unique(df$t))),
                 var.time = "t")

stopifnot("predictY returned one column per class" =
            ncol(pred$pred) == chosen_model$ng)

fc <- do.call(rbind, lapply(seq_len(chosen_model$ng), function(k) {
  data.frame(group = as.character(k),
             hr = sort(unique(df$t)) * WINDOW_H,
             fitted = pred$pred[, k])
}))

cl <- class_of(chosen_model)
oc <- merge(df, data.frame(id_num = as.integer(names(cl)),
                           group = as.character(cl)), by = "id_num")
oc <- aggregate(oc$total_dose, by = list(group = oc$group, t = oc$t), FUN = mean)
names(oc)[3] <- "observed"; oc$hr <- oc$t * WINDOW_H

n    <- table(cl)
lvl  <- names(sort(tapply(fc$fitted, fc$group, mean)))
labs <- sprintf("Class %s  (n=%s, %.1f%%)", lvl,
                format(as.integer(n[lvl]), big.mark = ","),
                100 * as.integer(n[lvl]) / sum(n))
fc$group <- factor(fc$group, levels = lvl, labels = labs)
oc$group <- factor(oc$group, levels = lvl, labels = labs)
pal <- setNames(colorRampPalette(c("#a8c8ee", "#14427e"))(length(lvl)), labs)

traj_plot <- ggplot(fc, aes(hr, fitted, colour = group)) +
  geom_line(data = oc, aes(y = observed), linetype = "22", linewidth = 0.6,
            alpha = 0.8, show.legend = FALSE) +
  geom_line(linewidth = 1.1) + geom_point(size = 1.5) +
  scale_colour_manual(values = pal) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(title = "Fentanyl trajectory classes -- lcmm (random intercept)",
       subtitle = sprintf(
         "%s episodes that received any fentanyl while intubated in [0, %dh].\nSolid = fitted class curve, dashed = observed group mean.",
         format(length(cl), big.mark = ","), LANDMARK),
       x = "Hours since first IMV episode",
       y = sprintf("Fentanyl dose (%s)", COV$exposure$units), colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(colour = ink, face = "bold"),
        plot.subtitle = element_text(colour = muted, margin = margin(b = 10)),
        legend.position = "top",
        legend.text = element_text(colour = muted, size = 9),
        axis.title = element_text(colour = muted),
        axis.text = element_text(colour = muted),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = gridline, linewidth = 0.4),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA))

ggsave(file.path(dirs$phase, "lcmm_trajectories.png"), traj_plot,
       width = 7.5, height = 5.4, dpi = 200)

write_json(prov, file.path(dirs$phase, "lcmm_provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)


# ---- 9. Provenance ----
# Which package versions produced these numbers?

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/05_lcmm_classes.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "05_lcmm_classes_sessioninfo.txt")
)
