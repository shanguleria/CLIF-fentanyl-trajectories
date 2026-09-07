# Federated-pooling exports, shared by Phases 1 and 2.
#
# A median cannot be pooled across sites; a mean can, exactly, from n and the two
# sums. Carrying sum and sum_sq rather than only mean and sd is what makes the
# pooled figures exact rather than an approximation that assumes equal variances:
#     mean_pooled = sum(sum) / sum(n)
#     var_pooled  = (sum(sum_sq) - sum(sum)^2 / sum(n)) / (sum(n) - 1)
# Medians travel alongside because the doses here are right-skewed with a large
# mass at zero: pool the means, report the medians.
#
# Cells below reporting.small_cell_min_den are suppressed -- a mean over n = 1 is
# that patient's value.

# One pooling row. `v` is the raw vector; everything else is labelling.
pool_row <- function(scope, variable, unit, stratum, w, hr, v, min_cell) {
  v <- v[!is.na(v)]
  n <- length(v)
  small <- n > 0 && n < min_cell
  blank <- function(x) if (small || n == 0) NA_real_ else round(x, 6)
  q <- if (n) unname(quantile(v, c(0.25, 0.5, 0.75))) else rep(NA_real_, 3)
  data.frame(
    scope = scope, variable = variable, unit = unit, stratum = stratum,
    window_idx = w, window_start_hr = hr,
    n = n, n_suppressed_small_cell = as.integer(small),
    mean = blank(if (n) mean(v) else NA_real_),
    sd   = blank(if (n > 1) stats::sd(v) else NA_real_),
    sum  = blank(if (n) sum(v) else NA_real_),
    sum_sq = blank(if (n) sum(v^2) else NA_real_),
    min = blank(if (n) min(v) else NA_real_),
    max = blank(if (n) max(v) else NA_real_),
    median = blank(q[2]), q1 = blank(q[1]), q3 = blank(q[3]),
    stringsAsFactors = FALSE)
}

# Level counts per stratum. `strata` is a vector as long as `v`; an "overall"
# row is always emitted alongside the named strata.
pool_cat <- function(variable, v, strata, min_cell) {
  levels_of <- c("overall", sort(unique(strata)))
  do.call(rbind, lapply(levels_of, function(st) {
    x <- if (st == "overall") v else v[strata == st]
    tb <- table(x)
    do.call(rbind, lapply(names(tb), function(l) {
      k <- as.integer(tb[[l]])
      small <- k > 0 && k < min_cell
      data.frame(variable = variable, level = l, stratum = st,
                 n = if (small) NA_integer_ else k,
                 denominator = length(x),
                 pct = if (small) NA_real_ else round(100 * k / length(x), 2),
                 n_suppressed_small_cell = as.integer(small),
                 stringsAsFactors = FALSE)
    }))
  }))
}

# Collapse category levels for display through a map declared in covariates.json.
collapse_levels <- function(v, map) {
  # Look up through a NAMED VECTOR, not list subsetting. `unlist(map[v])` drops
  # the NULLs for unmatched values, so the result is shorter than v and ifelse
  # recycles it -- which silently misaligns every row rather than erroring.
  named <- names(map)[!startsWith(names(map), "_")]
  lut <- unlist(map[named])
  out <- unname(lut[v])
  out[is.na(out)] <- map[["_default"]]
  out[v == "Missing"] <- "Missing"
  stopifnot("collapse_levels changed the vector length" = length(out) == length(v))
  as.character(out)
}
