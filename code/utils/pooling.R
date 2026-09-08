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


# ---- Model-selection table, Chen et al. layout -------------------------------
# Statistics as ROWS, one column per candidate class count, with a %classK row
# per class so every solution's class sizes are visible at once -- not just the
# smallest. Format follows Chen et al., BMC Infect Dis 2026;26:950
# (doi:10.1186/s12879-026-12981-9) Table 1.
#
# stats: a data.frame with one row per model, an `ng` column, and any numeric
#        columns to show as rows.
# sizes: a named list, one element per model, each a vector of class counts.
selection_table <- function(stats, sizes, rows = NULL) {
  stopifnot("stats needs an ng column" = "ng" %in% names(stats),
            "one size vector per model" = nrow(stats) == length(sizes))
  if (is.null(rows)) {
    rows <- setdiff(names(stats), c("model", "ng"))
    rows <- rows[vapply(stats[rows], is.numeric, logical(1))]
  }
  cols <- paste0("ng", stats$ng)

  body <- lapply(rows, function(r) {
    v <- stats[[r]]
    data.frame(statistic = r, t(setNames(as.character(round(v, 4)), cols)),
               check.names = FALSE, stringsAsFactors = FALSE)
  })

  kmax <- max(vapply(sizes, length, integer(1)))
  pct <- lapply(seq_len(kmax), function(k) {
    v <- vapply(sizes, function(s)
      if (length(s) >= k) sprintf("%.2f", 100 * s[k] / sum(s)) else NA_character_,
      character(1))
    data.frame(statistic = sprintf("%%class%d", k),
               t(setNames(v, cols)), check.names = FALSE,
               stringsAsFactors = FALSE)
  })

  n <- lapply(seq_len(kmax), function(k) {
    v <- vapply(sizes, function(s)
      if (length(s) >= k) as.character(s[k]) else NA_character_, character(1))
    data.frame(statistic = sprintf("n_class%d", k),
               t(setNames(v, cols)), check.names = FALSE,
               stringsAsFactors = FALSE)
  })

  do.call(rbind, c(body, pct, n))
}
