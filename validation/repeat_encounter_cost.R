# ==============================================================================
# repeat_encounter_cost.R
#
# Purpose : Measure the cost of using the ENCOUNTER BLOCK as the unit of
#           analysis when a patient can contribute more than one block
#           (design_notes.md section 11, "DECIDED: all qualifying blocks are kept").
#
#           The question this answers is NOT "does gbmt run" -- it does, and E0
#           below demonstrates that explicitly. It is: what does the repeat
#           structure do to CLASS ENUMERATION, which is the single most
#           consequential decision in the study (section 7).
#
#           E0  gbmt runs to convergence on a cohort containing repeat episodes
#           E1  BIC-selected ng, first-episodes-only vs all-episodes, over a
#               sweep of repeat fractions
#           E2  agreement (ARI) between the two solutions on the shared episodes
#           E3  design effect implied by the repeat structure
#
# Author  : Shan Guleria
# Created : 2026-09-05
# STATUS  : WRITTEN, NOT YET RUN TO COMPLETION. E0 was run separately and is
#           reported in design_notes.md section 11: gbmt converged in 5.9s on
#           168 episodes from 120 patients (48 contributing two). The E1-E3
#           sweep is slow (20 gbmt fits) and was deprioritised 2026-09-05 --
#           the decision it informs is reversible with a one-line filter, so it
#           is not blocking. Run it before the class solution is finalised.
# Inputs  : none (all data simulated)
# Outputs : output/final_no_phi/validation/repeat_encounter_cost.csv
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

set.seed(20260905)

N_PAT    <- 300L    # patients
N_WIN    <- 18L     # 4h windows to 72h
NG_RANGE <- 1:5
NSTART   <- 10L     # lower than the protocol's 50; this is a simulation, not a fit
DEG      <- 2L

# ---- Adjusted Rand Index ------------------------------------------------------
ari <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  s_ij <- sum(choose(tab, 2))
  s_i  <- sum(choose(rowSums(tab), 2))
  s_j  <- sum(choose(colSums(tab), 2))
  exp_ <- s_i * s_j / choose(n, 2)
  max_ <- (s_i + s_j) / 2
  if (isTRUE(all.equal(max_, exp_))) return(NA_real_)
  (s_ij - exp_) / (max_ - exp_)
}

# ---- Three fentanyl-shaped trajectory classes ---------------------------------
# escalating, stable-moderate, weaning. mcg/kg/hr, on the scale of real dosing.
shape <- function(cls, t) {
  switch(cls,
    "escalate" = 0.8 + 0.9 * (t / N_WIN) + 0.6 * (t / N_WIN)^2,
    "stable"   = 1.4 + 0.1 * (t / N_WIN),
    "wean"     = 2.2 - 1.7 * (t / N_WIN) - 0.3 * (t / N_WIN)^2
  )
}
CLASSES <- c("escalate", "stable", "wean")

# ---- Simulate one cohort ------------------------------------------------------
# repeat_frac of patients contribute a SECOND episode. The second episode shares
# the patient's random intercept and, with prob same_class_p, the same class --
# which is what makes the two episodes non-independent.
simulate <- function(repeat_frac, same_class_p = 0.7) {
  cls  <- sample(CLASSES, N_PAT, replace = TRUE)
  ri   <- rnorm(N_PAT, 0, 0.30)                 # patient-level random intercept
  n_rep <- round(N_PAT * repeat_frac)
  repeaters <- if (n_rep > 0) sample.int(N_PAT, n_rep) else integer(0)

  rows <- list()
  ep <- 0L
  for (i in seq_len(N_PAT)) {
    eps <- 1L + as.integer(i %in% repeaters)
    for (k in seq_len(eps)) {
      ep <- ep + 1L
      c_k <- if (k == 1L) cls[i] else
             if (runif(1) < same_class_p) cls[i] else sample(CLASSES, 1)
      t <- seq_len(N_WIN)
      y <- shape(c_k, t) + ri[i] + rnorm(N_WIN, 0, 0.25)
      rows[[length(rows) + 1L]] <- data.frame(
        block   = sprintf("b%04d", ep),
        patient = sprintf("p%04d", i),
        episode = k,
        time    = t,
        truth   = c_k,
        dose    = pmax(y, 0),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# ---- Fit gbmt across ng, return BIC and the selected solution -----------------
fit_sweep <- function(df, label) {
  d <- df[, c("block", "time", "dose")]
  res <- list(); bics <- rep(NA_real_, length(NG_RANGE))
  for (j in seq_along(NG_RANGE)) {
    g <- NG_RANGE[j]
    m <- try(gbmt(x.names = "dose", unit = "block", time = "time",
                  ng = g, d = DEG, data = d, scaling = 0,
                  nstart = NSTART, quiet = TRUE), silent = TRUE)
    if (inherits(m, "try-error")) {
      cat(sprintf("  %-22s ng=%d  FAILED\n", label, g)); next
    }
    bics[j] <- m$ic["bic"]
    res[[as.character(g)]] <- m
    cat(sprintf("  %-22s ng=%d  BIC=%12.1f\n", label, g, bics[j]))
  }
  list(bic = bics, fits = res,
       best = NG_RANGE[which.min(bics)])
}

assign_of <- function(m, df) {
  a <- m$assign.list
  lab <- setNames(rep(NA_integer_, length(unique(df$block))), unique(df$block))
  for (g in seq_along(a)) lab[a[[g]]] <- g
  lab
}

cat("\n================ E0: does gbmt run with repeat episodes? ================\n")
d40 <- simulate(0.40)
t0 <- Sys.time()
m0 <- try(gbmt(x.names = "dose", unit = "block", time = "time", ng = 3, d = DEG,
               data = d40[, c("block", "time", "dose")], scaling = 0,
               nstart = NSTART, quiet = TRUE), silent = TRUE)
cat(sprintf("  cohort: %d episodes from %d patients (%d repeaters)\n",
            length(unique(d40$block)), length(unique(d40$patient)),
            sum(table(d40$patient) > N_WIN)))
if (inherits(m0, "try-error")) {
  cat("  RESULT: gbmt ERRORED\n"); print(m0)
} else {
  cat(sprintf("  RESULT: gbmt converged in %.1fs; ng=3 fitted, BIC=%.1f\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")), m0$ic["bic"]))
  cat("  => Repeat episodes do NOT prevent gbmt from running. The cost is inferential.\n")
}

cat("\n================ E1/E2/E3: the cost, by repeat fraction ================\n")
out <- list()
for (rf in c(0.00, 0.10, 0.25, 0.40)) {
  cat(sprintf("\n-- repeat fraction %.2f --\n", rf))
  df   <- simulate(rf)
  firsts <- df[df$episode == 1L, ]

  s_all   <- fit_sweep(df,     sprintf("all-episodes rf=%.2f", rf))
  s_first <- fit_sweep(firsts, sprintf("first-only   rf=%.2f", rf))

  # ARI between the two solutions, on the first episodes they share
  a_all <- assign_of(s_all$fits[[as.character(s_all$best)]], df)
  a_fst <- assign_of(s_first$fits[[as.character(s_first$best)]], firsts)
  shared <- intersect(names(a_all), names(a_fst))
  ari_solutions <- if (length(shared) > 10) ari(a_all[shared], a_fst[shared]) else NA_real_
  # ARI of each against the truth, on first episodes
  truth_f <- tapply(firsts$truth, firsts$block, function(z) z[1])
  ari_all_truth <- ari(a_all[shared], truth_f[shared])
  ari_fst_truth <- ari(a_fst[shared], truth_f[shared])

  m_bar <- mean(table(df$patient) / N_WIN)          # mean episodes per patient
  deff  <- 1 + (m_bar - 1) * 0.5                    # ICC = 0.5, a pessimistic guess

  out[[length(out) + 1L]] <- data.frame(
    repeat_frac      = rf,
    n_patients       = length(unique(df$patient)),
    n_episodes       = length(unique(df$block)),
    ng_first_only    = s_first$best,
    ng_all_episodes  = s_all$best,
    ari_solutions    = ari_solutions,
    ari_first_truth  = ari_fst_truth,
    ari_all_truth    = ari_all_truth,
    mean_eps_per_pat = m_bar,
    deff_icc0.5      = deff,
    se_inflation_pct = 100 * (sqrt(deff) - 1)
  )
}

res <- do.call(rbind, out)
print(res, row.names = FALSE, digits = 3)
write.csv(res, file.path(out_interim, "repeat_encounter_cost.csv"), row.names = FALSE)
cat(sprintf("\nwritten: %s\n", file.path(out_interim, "repeat_encounter_cost.csv")))
