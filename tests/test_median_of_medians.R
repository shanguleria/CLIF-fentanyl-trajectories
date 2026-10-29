# Known-answer checks on the weighted median of medians in code/utils/pooling.R.
#
# Reference is McGrath 2019 section 2.2 (references/ ... meta-analysis of
# medians.pdf) plus hand computation. metamedian and Hmisc are NOT dependencies,
# so every expectation below is derived on paper and written as a literal.
#
# Run standalone:  Rscript tests/test_median_of_medians.R
suppressPackageStartupMessages(library(here))
source(here("code", "utils", "pooling.R"))

pass <- 0L; fail <- 0L
check <- function(name, ok, detail = "") {
  if (isTRUE(ok)) { pass <<- pass + 1L; cat(sprintf("  PASS  %s\n", name)) }
  else { fail <<- fail + 1L; cat(sprintf("  FAIL  %s\n        %s\n", name, detail)) }
}
near <- function(a, b, tol = 1e-9) isTRUE(abs(a - b) < tol)

# ---- weighted_quantile ------------------------------------------------------
check("equal weights, odd k, q=0.5 is the middle value",
      near(weighted_quantile(c(1, 2, 3), c(1, 1, 1), 0.5), 2))
check("equal weights, even k, q=0.5 is the midpoint",
      near(weighted_quantile(c(1, 3), c(1, 1), 0.5), 2))
check("q=0 is the minimum", near(weighted_quantile(c(5, 1, 3), c(2, 9, 4), 0), 1))
check("q=1 is the maximum", near(weighted_quantile(c(5, 1, 3), c(2, 9, 4), 1), 5))
check("input order does not matter",
      near(weighted_quantile(c(3, 1, 2), c(1, 1, 1), 0.5),
           weighted_quantile(c(1, 2, 3), c(1, 1, 1), 0.5)))
# A heavy study pulls the weighted median TOWARD its value without pinning it
# there -- interpolation, not a step. Hand: w = (0.8, 0.1, 0.1), midpoint
# positions (0.40, 0.85, 0.95); q = 0.5 sits between the first two, so
# t = (0.5-0.40)/(0.85-0.40) = 2/9 and the value is 1 + (2/9)(2-1) = 1.2222.
check("a heavy study pulls the weighted median toward its own value",
      near(weighted_quantile(c(1, 2, 3), c(8, 1, 1), 0.5), 1 + 2/9, 1e-9),
      sprintf("got %.6f want %.6f", weighted_quantile(c(1, 2, 3), c(8, 1, 1), 0.5), 1 + 2/9))
check("and it lands below the unweighted median",
      weighted_quantile(c(1, 2, 3), c(8, 1, 1), 0.5) <
        weighted_quantile(c(1, 2, 3), c(1, 1, 1), 0.5))
# Weights scale-invariant: only relative size matters.
check("weights are scale-invariant",
      near(weighted_quantile(c(1, 2, 3), c(2, 2, 2), 0.5),
           weighted_quantile(c(1, 2, 3), c(50, 50, 50), 0.5)))
# A zero or negative weight is dropped, not silently treated as 1.
check("a zero-weight study is dropped",
      near(weighted_quantile(c(1, 2, 99), c(1, 1, 0), 0.5), 1.5),
      sprintf("got %.6f", weighted_quantile(c(1, 2, 99), c(1, 1, 0), 0.5)))
check("an NA value is dropped",
      near(weighted_quantile(c(1, 2, NA), c(1, 1, 5), 0.5), 1.5))

# Hand-computed, on midpoint positions: x = (10, 20), w = (0.25, 0.75);
# positions (0.125, 0.625); q = 0.5 -> t = (0.5-0.125)/0.5 = 0.75 -> 17.5
check("interpolation matches hand computation",
      near(weighted_quantile(c(10, 20), c(1, 3), 0.5), 17.5, 1e-9),
      sprintf("got %.9f want 17.5", weighted_quantile(c(10, 20), c(1, 3), 0.5)))

# ---- wm_of_medians: the interval degenerates below k = 5 --------------------
# z/(2 sqrt k) > 1/2 whenever sqrt(k) < 1.959964, i.e. k <= 3; at k = 4 it is
# 0.49, so the limits are the 1st and 99th weighted percentiles -- still
# essentially the range. This is the fact that makes "WM 95% CI" and "range of
# site medians" the SAME COLUMN at three sites.
# The interval is pinned to the range at every k we will plausibly have. The
# flag is MEASURED from the result, so it reports the truth rather than the
# half < 1/2 prediction, which understates how long the pinning lasts.
for (k in 2:5) {
  r <- wm_of_medians(seq_len(k) * 10, rep(100, k))
  check(sprintf("k=%d: the interval is the range, and is flagged as such", k),
        isTRUE(r$degenerate) && near(r$lo, r$min) && near(r$hi, r$max),
        sprintf("lo %.4f min %.4f hi %.4f max %.4f degenerate=%s",
                r$lo, r$min, r$hi, r$max, r$degenerate))
}
# It has to stop being degenerate eventually, or the flag means nothing.
kk <- Filter(function(k) isFALSE(wm_of_medians(seq_len(k), rep(10, k))$degenerate), 2:60)
check("the interval does narrow at some k, so the flag can be FALSE",
      length(kk) > 0, "degenerate at every k from 2 to 60 -- the flag is inert")
check("degeneracy ends at k = 6 on equal weights", length(kk) && kk[1] == 6L,
      sprintf("first non-degenerate k = %s", if (length(kk)) kk[1] else "none"))

# ---- wm_of_medians: the point estimate --------------------------------------
r <- wm_of_medians(c(60, 64), c(2471, 1599))
# w = (0.60712, 0.39288); midpoint positions (0.30356, 0.80356).
# q = 0.5 -> t = (0.5-0.30356)/0.5 = 0.39288 -> 60 + 0.39288*4 = 61.5715
check("two sites: the weighted median sits between, nearer the larger study",
      near(r$estimate, 60 + (0.5 - 2471/4070/2) / 0.5 * 4, 1e-6) &&
        r$estimate > 60 && r$estimate < 62,
      sprintf("got %.6f", r$estimate))
check("two sites: k and range are reported", r$k == 2L && r$min == 60 && r$max == 64)
check("argument order does not change the result",
      near(wm_of_medians(c(64, 60), c(1599, 2471))$estimate,
           wm_of_medians(c(60, 64), c(2471, 1599))$estimate))
check("all-missing input returns NA rather than erroring",
      is.na(wm_of_medians(c(NA, NA), c(10, 10))$estimate))

cat(sprintf("\n%d/%d passed\n", pass, pass + fail))
if (fail) quit(status = 1)
