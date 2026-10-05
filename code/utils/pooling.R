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
# `strata_scope`: "all" emits overall plus every named stratum; "overall" emits
# the overall row only. Config: covariates.json disclosure.full_levels_strata_scope
# -- a categorical published in BOTH a full and a collapsed form ships the full
# form at overall only, because publishing the two views over the same strata
# makes each suppressed cell recoverable by subtraction.
pool_cat <- function(variable, v, strata, min_cell, strata_scope = "all") {
  # as.character FIRST: c("overall", <factor>) dispatches on the character and
  # coerces the factor to its INTEGER CODES, so every named stratum then matches
  # nothing and only the overall rows survive. Table 1's `grp` is a factor.
  strata <- as.character(strata)
  stopifnot("strata_scope must be 'all' or 'overall'" =
              strata_scope %in% c("all", "overall"))
  levels_of <- if (identical(strata_scope, "overall")) "overall"
               else c("overall", sort(unique(strata)))
  do.call(rbind, lapply(levels_of, function(st) {
    x <- if (st == "overall") v else v[strata == st]
    tb <- table(x)
    k <- as.integer(tb)
    names(k) <- names(tb)
    sup <- k > 0 & k < min_cell
    # A row of counts sums to its PUBLISHED denominator, so ONE suppressed cell
    # is just denominator minus the rest. Suppress the next-smallest as well.
    if (sum(sup) == 1L) {
      rest <- which(!sup & k > 0)
      if (length(rest)) sup[rest[which.min(k[rest])]] <- TRUE
    }
    do.call(rbind, lapply(names(tb), function(l) {
      small <- unname(sup[[l]])
      data.frame(variable = variable, level = l, stratum = st,
                 n = if (small) NA_integer_ else unname(k[[l]]),
                 denominator = length(x),
                 pct = if (small) NA_real_ else round(100 * k[[l]] / length(x), 2),
                 n_suppressed_small_cell = as.integer(small),
                 stringsAsFactors = FALSE)
    }))
  }))
}

# ---- Histograms: the only way a median pools --------------------------------
# A median has no closed form across sites, and for a skewed or zero-inflated
# variable the mean is not an acceptable substitute. Counts on FIXED bins do
# pool: sum them, then read any quantile off the pooled distribution, accurate
# to within one bin width. Config: covariates.json pooling.histograms.
#
# `spec` is one entry from that block. Returns one row per bin, lowest first,
# carrying numeric bin_lo/bin_hi so a pooler never has to parse a label.
pool_hist <- function(variable, v, spec) {
  v <- v[!is.na(v)]
  from <- as.numeric(spec$from)
  to   <- as.numeric(spec$to)
  by   <- as.numeric(spec$by)
  zero_bin <- isTRUE(spec$zero_bin)
  open_top <- isTRUE(spec$open_top)
  # round(): seq() accumulates float error at 0.05 steps, which would put an
  # edge at 0.30000000000000004 and make two sites' nominally identical bins
  # compare unequal. The bins are the federation contract, so they must be exact.
  edges <- round(seq(from, to, by = by), 10)
  stopifnot("a histogram needs at least two edges" = length(edges) >= 2)

  lo <- numeric(0); hi <- numeric(0); cnt <- integer(0); lab <- character(0)
  rest <- v
  if (zero_bin) {
    lo <- c(lo, from); hi <- c(hi, from)
    cnt <- c(cnt, sum(rest == from)); lab <- c(lab, sprintf("%g", from))
    rest <- rest[rest != from]
  }
  n_int <- length(edges) - 1L
  for (i in seq_len(n_int)) {
    a <- edges[i]; b <- edges[i + 1L]
    # Half-open [a, b), matching windows.interval_convention -- except the final
    # interval of a BOUNDED variable, which must close or its maximum falls out
    # of every bin and the completeness assert below fires.
    last_closed <- (i == n_int) && !open_top
    sel <- if (last_closed) rest >= a & rest <= b else rest >= a & rest < b
    lo <- c(lo, a); hi <- c(hi, b); cnt <- c(cnt, sum(sel))
    lab <- c(lab, sprintf(if (last_closed) "[%g, %g]" else "[%g, %g)",
                          if (zero_bin && i == 1L) from else a, b))
  }
  if (open_top) {
    top <- edges[length(edges)]
    lo <- c(lo, top); hi <- c(hi, Inf)
    cnt <- c(cnt, sum(rest >= top)); lab <- c(lab, sprintf("[%g, Inf)", top))
  }

  # A value that lands in no bin is silently lost exposure, which is the whole
  # failure mode this file exists to prevent.
  stopifnot("every value must land in exactly one bin" = sum(cnt) == length(v))
  data.frame(
    variable = variable,
    unit = if (is.null(spec$unit)) NA_character_ else as.character(spec$unit),
    bin = seq_along(cnt), bin_label = lab, bin_lo = lo, bin_hi = hi,
    n = as.integer(cnt),
    pct = round(100 * cnt / max(length(v), 1), 3),
    n_total = length(v),
    stringsAsFactors = FALSE)
}


# ---- Cluster bootstrap for a clustered proportion ---------------------------
# Events are clustered within episodes, so a naive binomial interval is too
# narrow. Resamples EPISODES and reports both the interval and the design
# effect -- how much wider the honest interval is than the naive one.
#
# Shared by 05_titration.R and 06_unit_variation.R. `n_resamples` and
# `episodes` are REQUIRED arguments with no default: both are properties of the
# calling analysis, not of this function, and a silent default for either is
# exactly how one caller's universe leaked into another's (lessons.md #26).
boot_group <- function(d, group_col, seed, outcome, episodes, n_resamples) {
  n_ep <- length(episodes)
  stopifnot("an event belongs to no episode in the resampling universe" =
              all(d$encounter_block %in% episodes),
            "n_resamples must be a positive whole number" =
              is.numeric(n_resamples) && n_resamples > 0)
  # Per (episode, group) counts, so a resample is a weighted sum rather than a
  # re-tabulation of every event.
  g <- factor(d[[group_col]])
  cell <- data.frame(ei = match(d$encounter_block, episodes),
                     gi = as.integer(g), n = 1L, k = as.integer(d[[outcome]]))
  cell <- aggregate(cbind(n, k) ~ ei + gi, data = cell, FUN = sum)
  G <- nlevels(g)

  point_n <- numeric(G); point_k <- numeric(G)
  s <- rowsum(as.matrix(cell[, c("k", "n")]), cell$gi, reorder = TRUE)
  idx <- as.integer(rownames(s))
  point_k[idx] <- s[, "k"]; point_n[idx] <- s[, "n"]

  # Re-seeded immediately before the draw rather than relying on a seed set at
  # the top of the script: any RNG consumed in between would silently change
  # which intervals a reader sees.
  set.seed(seed)
  reps <- matrix(NA_real_, nrow = n_resamples, ncol = G)
  overall <- numeric(n_resamples)
  for (b in seq_len(n_resamples)) {
    mult <- tabulate(sample.int(n_ep, n_ep, replace = TRUE), nbins = n_ep)
    w <- mult[cell$ei]
    num <- numeric(G); den <- numeric(G)
    sb <- rowsum(cbind(cell$k * w, cell$n * w), cell$gi, reorder = TRUE)
    ib <- as.integer(rownames(sb))
    num[ib] <- sb[, 1]; den[ib] <- sb[, 2]
    reps[b, ] <- ifelse(den > 0, 100 * num / den, NA_real_)
    overall[b] <- 100 * sum(num) / sum(den)
  }

  out <- data.frame(
    group = levels(g), n_events = point_n, n_paired = point_k,
    pct_paired = 100 * point_k / point_n,
    ci_lo = apply(reps, 2, quantile, probs = 0.025, na.rm = TRUE),
    ci_hi = apply(reps, 2, quantile, probs = 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE)
  out$n_episodes <- as.vector(tapply(d$encounter_block, g,
                                     function(x) length(unique(x)))[out$group])
  # PER-GROUP design effect. Read by SINGLE-GROUP callers -- 05_titration.R
  # passes g = "all", so this column is the design effect of the overall rate
  # and is what that script ships.
  #
  # MULTI-GROUP callers deliberately discard it: 06_unit_variation.R overwrites
  # the column with the overall attribute below, because a per-unit deff would
  # only be needed to combine a NAMED unit across sites, and
  # docs/cohort_and_outputs.md section 8.7 forbids that -- care_setting level
  # sets differ by site by construction, so the centre pools the between-unit
  # spread or an ICC, never the named units. Do not "fix" that overwrite by
  # shipping both; the per-unit column would have no permitted consumer.
  naive <- out$pct_paired * (100 - out$pct_paired) / out$n_events
  out$deff <- apply(reps, 2, stats::var, na.rm = TRUE) / naive

  p0 <- 100 * sum(point_k) / sum(point_n)
  attr(out, "overall_pct") <- p0
  attr(out, "deff") <- stats::var(overall) / (p0 * (100 - p0) / sum(point_n))
  out
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
