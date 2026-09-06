# ==============================================================================
# 03_scaling_experiments.R
#
# Purpose : Reproducible evidence for the scaling decisions recorded in
#           references/sedation_gbtm_design_notes.md. Three experiments on
#           synthetic, fentanyl-shaped data:
#             E1  what each `scaling` value erases
#             E2  what a near-all-zero second indicator does under scaling >= 1
#             E3  whether scaling = 0 is sensitive to the units of each indicator
#             E4  what an unbalanced panel (early extubation) does to the model
# Author  : Shan Guleria
# Created : 2026-09-04
# Inputs  : none (all data simulated)
# Outputs : output/final_no_phi/validation/scaling_experiments.csv
#           logs/03_scaling_experiments_sessioninfo.txt
# ==============================================================================

pkgs <- c("here", "gbmt")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

out_interim <- here("output", "final_no_phi", "validation")
dir.create(out_interim, recursive = TRUE, showWarnings = FALSE)


# ---- Adjusted Rand Index (see 02_indicator_sensitivity.R for the explanation)
adjusted_rand <- function(x, y) {
  tb <- table(x, y); ch2 <- function(k) k * (k - 1) / 2
  n <- sum(tb); idx <- sum(ch2(tb))
  ai <- sum(ch2(rowSums(tb))); bj <- sum(ch2(colSums(tb)))
  e <- ai * bj / ch2(n)
  (idx - e) / (0.5 * (ai + bj) - e)
}

# One simulated "patient": id, window index, and one or two dose columns.
sim_patients <- function(ids, f_inf, f_push = NULL, sd_inf = 8, sd_push = 1.5,
                         n_win = 10) {
  do.call(rbind, lapply(ids, function(i) {
    w <- seq_len(n_win)
    out <- data.frame(id = i, w = w,
                      inf = pmax(0, f_inf(w) + rnorm(n_win, 0, sd_inf)))
    if (!is.null(f_push)) out$push <- pmax(0, f_push(w) + rnorm(n_win, 0, sd_push))
    out
  }))
}

results <- list()


# ---- E1: what each scaling value erases --------------------------------------
# Three true phenotypes that differ in LEVEL and in SHAPE:
#   HI = flat high dose, LO = flat low dose, UP = rising.
# scaling >= 1 normalises each patient against their OWN mean/SD, so the two
# flat groups should become indistinguishable.

set.seed(42)
d1 <- rbind(
  sim_patients(sprintf("HI%02d", 1:4), function(t) rep(200, length(t))),
  sim_patients(sprintf("LO%02d", 1:4), function(t) rep( 50, length(t))),
  sim_patients(sprintf("UP%02d", 1:4), function(t) 40 + 17 * t)
)
truth1 <- substr(unique(d1$id), 1, 2)

for (s in 0:2) {
  f <- gbmt(x.names = "inf", unit = "id", time = "w", d = 1, ng = 3,
            data = d1, scaling = s, quiet = TRUE)
  results[[length(results) + 1]] <- data.frame(
    experiment = "E1_what_scaling_erases", setting = paste0("scaling=", s),
    metric = "ARI_vs_truth", value = adjusted_rand(truth1, f$assign[unique(d1$id)])
  )
}


# ---- E2: a near-all-zero second indicator ------------------------------------
# Push doses are zero for most patients in most windows. Under scaling >= 1 the
# within-patient SD of that column is exactly 0 for a patient who never received
# a bolus, and the normalisation divides by it.

set.seed(3)
d2 <- rbind(
  sim_patients(sprintf("P%02d", 1:6),  function(t) rep(150, length(t))),
  sim_patients(sprintf("P%02d", 7:12), function(t) 40 + 12 * t)
)
d2$push <- 0
d2$push[d2$id %in% c("P01", "P08") & d2$w %in% c(3, 7)] <- 50
truth2 <- c(rep("flat", 6), rep("rise", 6))

for (s in c(0, 2)) {
  f <- gbmt(x.names = c("inf", "push"), unit = "id", time = "w", d = 1, ng = 2,
            data = d2, scaling = s, quiet = TRUE)
  results[[length(results) + 1]] <- data.frame(
    experiment = "E2_zero_inflated_indicator", setting = paste0("scaling=", s),
    metric = "ARI_vs_truth", value = adjusted_rand(truth2, f$assign[unique(d2$id)])
  )
}


# ---- E3: is scaling = 0 sensitive to indicator units? ------------------------
# A Gaussian mixture with a freely estimated full covariance per group is
# theoretically equivariant to rescaling a variable. But gbmt initialises EM
# from a Ward hierarchical clustering, which is NOT scale-invariant -- so on
# poorly separated data the two runs can land in different local optima.
# Repeated over seeds because a single seed would not show the instability.

for (seed in 1:6) {
  set.seed(seed)
  d3 <- rbind(
    sim_patients(sprintf("A%02d", 1:8), function(t) rep(110, length(t)),
                 function(t) rep(2, length(t)), sd_inf = 45, sd_push = 4),
    sim_patients(sprintf("B%02d", 1:8), function(t) 55 + 8 * t,
                 function(t) rep(2, length(t)), sd_inf = 45, sd_push = 4),
    sim_patients(sprintf("C%02d", 1:8), function(t) rep(60, length(t)),
                 function(t) rep(7, length(t)), sd_inf = 45, sd_push = 4)
  )
  units <- unique(d3$id)
  fit_assign <- function(dd) {
    gbmt(x.names = c("inf", "push"), unit = "id", time = "w", d = 1, ng = 3,
         data = dd, scaling = 0, quiet = TRUE)$assign[units]
  }
  d3_rescaled <- d3; d3_rescaled$push <- d3_rescaled$push * 10

  base <- fit_assign(d3)
  results[[length(results) + 1]] <- data.frame(
    experiment = "E3_unit_sensitivity_scaling0", setting = paste0("seed=", seed),
    metric = "ARI_base_vs_push_x10", value = adjusted_rand(base, fit_assign(d3_rescaled))
  )
}


# ---- E4: unbalanced panels and the polynomial degree cap ---------------------
# Motivation: if patients extubated early are kept in the model with short
# trajectories rather than excluded at a landmark, the panel becomes unbalanced.
# gbmt permits this, but `d` cannot exceed (shortest unit's time points - 1),
# and the package enforces that SILENTLY -- it caps d and issues a warning
# rather than failing. One briefly-ventilated patient therefore constrains the
# trajectory shape available to the entire cohort.
# Evidence for references/sedation_gbtm_design_notes.md section 8.

set.seed(5)
d4 <- do.call(rbind, c(
  lapply(sprintf("L%02d", 1:8), function(i)
    sim_patients(i, function(t) rep(120, length(t)), n_win = 18)),
  lapply(sprintf("M%02d", 1:8), function(i)
    sim_patients(i, function(t) 150 - 6 * t, n_win = 18)),
  list(sim_patients("E01", function(t) rep(100, length(t)), n_win = 4),   # extubated early
       sim_patients("E02", function(t) rep( 90, length(t)), n_win = 9))
))
shortest <- min(table(d4$id))

for (dd in c(2, 4)) {
  w <- NULL
  f <- withCallingHandlers(
    gbmt(x.names = "inf", unit = "id", time = "w", d = dd, ng = 2,
         data = d4, scaling = 0, quiet = TRUE),
    warning = function(cond) { w <<- conditionMessage(cond); invokeRestart("muffleWarning") }
  )
  effective_d <- f$call$d          # gbmt records the d it ACTUALLY used
  results[[length(results) + 1]] <- data.frame(
    experiment = "E4_unbalanced_panel_degree_cap",
    setting = paste0("requested_d=", dd, " shortest_unit_windows=", shortest),
    metric = "effective_d_used", value = effective_d
  )
  cat(sprintf("E4: requested d=%d -> effective d=%d %s\n", dd, effective_d,
              if (is.null(w)) "(no warning)" else paste0("[warning: ", w, "]")))
}


# ---- Report ------------------------------------------------------------------

res <- do.call(rbind, results)
res$value <- round(res$value, 3)
print(res, row.names = FALSE)

e3 <- res$value[res$experiment == "E3_unit_sensitivity_scaling0"]
cat("\nE3 summary: mean ARI =", round(mean(e3), 3),
    "| runs where the partition changed:", sum(e3 < 0.999), "of", length(e3), "\n")

write.csv(res, file.path(out_interim, "scaling_experiments.csv"), row.names = FALSE)

writeLines(
  c(paste("Run at:", format(Sys.time(), usetz = TRUE)), capture.output(sessionInfo())),
  here("logs", "03_scaling_experiments_sessioninfo.txt")
)
