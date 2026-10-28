# Known-answer checks on the inverse-variance pool in code/utils/pooling.R.
#
# metafor is deliberately not a dependency, so the reference is hand computation:
# every expectation below is derived on paper and written as a literal, not
# produced by the code under test.
#
# Run standalone:  Rscript tests/test_coefficient_meta.R
suppressPackageStartupMessages(library(here))
source(here("code", "utils", "pooling.R"))

pass <- 0L; fail <- 0L
check <- function(name, ok, detail = "") {
  if (isTRUE(ok)) { pass <<- pass + 1L; cat(sprintf("  PASS  %s\n", name)) }
  else { fail <<- fail + 1L; cat(sprintf("  FAIL  %s\n        %s\n", name, detail)) }
}
near <- function(a, b, tol = 1e-6) isTRUE(abs(a - b) < tol)

SPEC <- list(method = "coefficient_meta", coefficient = "estimate",
             standard_error = "std_error")
frame <- function(site, est, se, term = "Pindicated_per10pp",
                  model = "primary", poolable = TRUE, n = 1000) {
  data.frame(site = site, model = model, term = term, estimate = est,
             std_error = se, poolable = poolable, n_admissions = n,
             stringsAsFactors = FALSE)
}
pool <- function(d, min_k = 5) suppressWarnings(
  pool_coefficient_meta(d, c("model", "term"), SPEC, "t.csv", min_k = min_k))

# ---- 1. identical sites pool to themselves, with SE / sqrt(k) ---------------
d <- rbind(frame("A", 0.07, 0.013), frame("B", 0.07, 0.013), frame("C", 0.07, 0.013))
r <- pool(d)
check("identical sites keep the estimate", near(r$estimate, 0.07),
      sprintf("got %.8f", r$estimate))
check("identical sites give SE/sqrt(k)", near(r$std_error, 0.013 / sqrt(3)),
      sprintf("got %.8f, want %.8f", r$std_error, 0.013 / sqrt(3)))
check("identical sites have Q = 0", near(r$Q, 0), sprintf("got %.8f", r$Q))
check("identical sites have I2 = 0", near(r$I2_pct, 0))
check("k_sites counts the sites", r$k_sites == 3L)

# ---- 2. hand-computed two-site pool -----------------------------------------
# beta = (0.1, 0.2), se = (0.1, 0.2) -> w = (100, 25), sum w = 125
#   est = (100*0.1 + 25*0.2)/125 = 15/125 = 0.12
#   se  = sqrt(1/125) = 0.0894427191
#   Q   = 100*(0.1-0.12)^2 + 25*(0.2-0.12)^2 = 0.04 + 0.16 = 0.20
d <- rbind(frame("A", 0.1, 0.1), frame("B", 0.2, 0.2))
r <- pool(d)
check("two-site estimate matches hand computation", near(r$estimate, 0.12),
      sprintf("got %.8f want 0.12", r$estimate))
check("two-site SE matches hand computation", near(r$std_error, sqrt(1/125)),
      sprintf("got %.8f want %.8f", r$std_error, sqrt(1/125)))
check("Cochran's Q matches hand computation", near(r$Q, 0.20),
      sprintf("got %.8f want 0.20", r$Q))
check("Q below its df gives I2 = 0", near(r$I2_pct, 0) && near(r$tau2, 0))
check("the odds ratio is exp of the pooled log-odds",
      near(r$odds_ratio, round(exp(0.12), 4)))
check("the interval is on the log-odds and exponentiated for reading",
      near(r$or_ci_lo, round(exp(0.12 - 1.959964 * sqrt(1/125)), 4)))

# ---- 3. a near-zero-weight site must not move the estimate ------------------
d <- rbind(frame("A", 0.10, 0.01), frame("B", 5.00, 100))
r <- pool(d)
check("a huge-SE site barely moves the pooled estimate",
      abs(r$estimate - 0.10) < 1e-3, sprintf("got %.8f", r$estimate))

# ---- 4. ORs are never averaged ----------------------------------------------
# exp() does not commute with a weighted mean: the mean of exp(0.1), exp(0.2)
# weighted 100/25 is 1.1280, while exp of the pooled log-odds is 1.1275.
d <- rbind(frame("A", 0.1, 0.1), frame("B", 0.2, 0.2))
r <- pool(d)
naive <- (100 * exp(0.1) + 25 * exp(0.2)) / 125
check("the pooled OR is not the weighted mean of the site ORs",
      !near(r$odds_ratio, round(naive, 4)),
      sprintf("pooled %.4f vs naive %.4f -- these must differ", r$odds_ratio, naive))

# ---- 5. non-poolable rows are excluded --------------------------------------
d <- rbind(frame("A", 0.07, 0.013), frame("B", 0.07, 0.013),
           frame("B", 9.99, 0.001, term = "icuN08S", poolable = FALSE))
r <- pool(d)
check("a poolable = FALSE row never reaches the pool",
      nrow(r) == 1 && near(r$estimate, 0.07),
      sprintf("%d row(s), estimate %.4f", nrow(r), r$estimate[1]))

# ---- 6. heterogeneity role follows heterogeneity_min_k -----------------------
d <- rbind(frame("A", 0.1, 0.1), frame("B", 0.2, 0.2), frame("C", 0.15, 0.15))
check("below min_k the heterogeneity is flagged descriptive_only",
      pool(d, min_k = 5)$heterogeneity_role == "descriptive_only")
check("at or above min_k it is reportable",
      pool(d, min_k = 3)$heterogeneity_role == "reportable")

# ---- 7. real heterogeneity is detected --------------------------------------
# beta = (0.0, 1.0) with se = 0.1 each -> w = 100, est = 0.5,
#   Q = 100*0.25 + 100*0.25 = 50 on 1 df, so I2 = 100*(50-1)/50 = 98%
d <- rbind(frame("A", 0.0, 0.1), frame("B", 1.0, 0.1))
r <- pool(d)
check("Q detects a genuine split", near(r$Q, 50), sprintf("got %.4f", r$Q))
check("I2 matches hand computation", near(r$I2_pct, round(100 * 49 / 50, 1)),
      sprintf("got %.1f want %.1f", r$I2_pct, 100 * 49 / 50))
check("tau2 is positive when Q exceeds its df", r$tau2 > 0)

# ---- 8. the pooled covariance diagonal reproduces the pooled variances ------
dc <- rbind(frame("A", 0.1, 0.1), frame("B", 0.2, 0.2),
            frame("A", -0.02, 0.05, term = "Pbolus_per10pp"),
            frame("B", -0.03, 0.08, term = "Pbolus_per10pp"))
rc <- pool(dc)
cov_in <- do.call(rbind, lapply(c("A", "B"), function(s) {
  se1 <- dc$std_error[dc$site == s & dc$term == "Pindicated_per10pp"]
  se2 <- dc$std_error[dc$site == s & dc$term == "Pbolus_per10pp"]
  data.frame(site = s, model = "primary",
             term_row = c("Pindicated_per10pp", "Pindicated_per10pp",
                          "Pbolus_per10pp", "Pbolus_per10pp"),
             term_col = c("Pindicated_per10pp", "Pbolus_per10pp",
                          "Pindicated_per10pp", "Pbolus_per10pp"),
             covariance = c(se1^2, 0.3 * se1 * se2, 0.3 * se1 * se2, se2^2),
             stringsAsFactors = FALSE)
}))
rv <- suppressWarnings(pool_covariance_meta(cov_in, SPEC, "cov.csv", rc))
check("the covariance pool returns a cell per term pair", !is.null(rv) && nrow(rv) == 4)
for (tm in c("Pindicated_per10pp", "Pbolus_per10pp")) {
  dg <- rv$covariance[rv$term_row == tm & rv$term_col == tm]
  se <- rc$std_error[rc$term == tm]
  check(sprintf("pooled var of %s equals its pooled SE squared", tm),
        near(sqrt(dg), se, 1e-5), sprintf("sqrt(var) %.8f vs SE %.8f", sqrt(dg), se))
}
check("the pooled covariance is symmetric",
      near(rv$covariance[rv$term_row == "Pindicated_per10pp" & rv$term_col == "Pbolus_per10pp"],
           rv$covariance[rv$term_row == "Pbolus_per10pp" & rv$term_col == "Pindicated_per10pp"]))

cat(sprintf("\n%d/%d passed\n", pass, pass + fail))
if (fail) quit(status = 1)
