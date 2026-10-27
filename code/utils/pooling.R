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

# ---- Between-group dispersion, in TWO families that must not be conflated -----
#
# Both describe "do these groups differ by more than sampling error", and they
# share most of their algebra, which is exactly why they are prefixed apart. This
# project has already shipped one defect of precisely this shape -- a
# within-episode design effect applied under a citation for between-unit
# over-dispersion -- and the fix cost a figure rebuild. Read the family names.
#
# FUNNEL FAMILY (no prefix: q_stat, phi, tau2, share_not_sampling)
#   Variances are the NULL binomial p0(100-p0)/n, on the PERCENTAGE scale, which
#   is the scale the funnel envelope is drawn on. It describes the limits as
#   DRAWN. Spiegelhalter, Qual Saf Health Care 2005;14:347-51.
#   Consequence worth knowing: with null variances the DL weights are
#   proportional to n, so a weighted mean is p0 BY CONSTRUCTION -- which is why
#   no random-effects estimate belongs in this family. It would be the
#   fixed-effect estimate wearing another name.
#
# META-ANALYTIC FAMILY (ma_ prefix)
#   Variances are each group's OWN p_i(100-p_i)/n_i, which is the standard
#   meta-analytic input and the only one under which I-squared, tau-squared and a
#   random-effects estimate mean what the literature means by them. Use these
#   when reporting heterogeneity, never the funnel family.
#
# GROUPING-AGNOSTIC BY DESIGN. Takes vectors, not a frame, so the same function
# serves units-within-a-site (06_unit_variation.R) and sites-within-the-
# consortium (pool_sites.R). k is the number of groups either way.

winsorise <- function(z, trim) {
  if (trim <= 0) return(z)
  q <- stats::quantile(z, c(trim, 1 - trim), names = FALSE)
  pmin(pmax(z, q[1]), q[2])
}

# I-squared and its test-based interval, via H. Higgins & Thompson,
# Stat Med 2002;21:1539-58. I2 = (H^2 - 1)/H^2 = (Q - df)/Q, so the point
# estimate is computed the direct way and only the INTERVAL needs H.
.i2_interval <- function(q_stat, df) {
  if (!is.finite(q_stat) || df < 1) return(c(NA_real_, NA_real_))
  h <- sqrt(max(q_stat, df) / df)              # H is bounded below at 1
  se_ln_h <- if (q_stat > df) {
    (log(q_stat) - log(df)) / (2 * (sqrt(2 * q_stat) - sqrt(2 * df - 1)))
  } else if (df > 1) {
    sqrt((1 / (2 * (df - 1))) * (1 - 1 / (3 * (df - 1)^2)))
  } else {
    NA_real_                                    # k = 2: no closed form here
  }
  if (!is.finite(se_ln_h)) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(0.975)
  h_lo <- max(1, exp(log(h) - z * se_ln_h))
  h_hi <- exp(log(h) + z * se_ln_h)
  c((h_lo^2 - 1) / h_lo^2, (h_hi^2 - 1) / h_hi^2)
}

overdispersion <- function(pct, n, p0, trim, min_k = 5L) {
  k <- length(pct)
  na <- NA_real_
  empty <- data.frame(
    k = k, q_stat = na, q_df = NA_integer_, q_p = na, phi = na,
    phi_winsorised = na, phi_crit_95 = na, tau2 = na, mean_binomial_var = na,
    share_not_sampling = na, ma_k = NA_integer_, ma_q = na,
    ma_i2 = na, ma_i2_lo = na,
    ma_i2_hi = na, ma_h2 = na, ma_tau2 = na, ma_re_pct = na, ma_re_lo = na,
    ma_re_hi = na, ma_dl_lo = na, ma_dl_hi = na,
    reportable = FALSE, stringsAsFactors = FALSE)
  if (k < 2L) return(empty)                    # df = 0; nothing is estimable

  df <- k - 1L

  # ---- funnel family: null variances, percentage scale --------------------
  v <- p0 * (100 - p0) / n
  z <- (pct - p0) / sqrt(v)
  q_stat <- sum(z^2)
  w <- 1 / v
  tau2 <- max(0, (q_stat - df) / (sum(w) - sum(w^2) / sum(w)))

  # ---- meta-analytic family: each group's OWN variance --------------------
  # Zero-variance groups (a rate at exactly 0 or 100) would carry infinite
  # weight. They are dropped from this family and the drop is visible in ma_k.
  v_i <- pct * (100 - pct) / n
  ok <- is.finite(v_i) & v_i > 0
  if (sum(ok) < 2L) {
    out <- empty
    out$q_stat <- q_stat; out$q_df <- df
    out$q_p <- stats::pchisq(q_stat, df, lower.tail = FALSE)
    out$phi <- q_stat / k
    out$phi_winsorised <- sum(winsorise(z, trim)^2) / k
    out$phi_crit_95 <- stats::qchisq(0.95, df) / k
    out$tau2 <- tau2; out$mean_binomial_var <- mean(v)
    out$share_not_sampling <- tau2 / (tau2 + mean(v))
    return(out)
  }
  y <- pct[ok]; vi <- v_i[ok]; k_ma <- length(y); df_ma <- k_ma - 1L
  wf <- 1 / vi
  theta_fe <- sum(wf * y) / sum(wf)
  ma_q <- sum(wf * (y - theta_fe)^2)
  ma_tau2 <- max(0, (ma_q - df_ma) / (sum(wf) - sum(wf^2) / sum(wf)))

  wr <- 1 / (vi + ma_tau2)                     # random-effects weights
  theta_re <- sum(wr * y) / sum(wr)
  # DerSimonian-Laird interval: normal, variance 1/sum(w*)
  se_dl <- sqrt(1 / sum(wr))
  zc <- stats::qnorm(0.975)
  # Hartung-Knapp-Sidik-Jonkman: t reference and a variance correction. This is
  # the method recommended for SMALL k, which is the whole reason it is here --
  # at k = 3 the t multiplier is 4.3 against the normal's 1.96, and that
  # difference is the honesty.
  se_hk <- sqrt(sum(wr * (y - theta_re)^2) / (df_ma * sum(wr)))
  tc <- stats::qt(0.975, df_ma)

  i2 <- max(0, (ma_q - df_ma) / ma_q)
  i2_ci <- .i2_interval(ma_q, df_ma)

  data.frame(
    k = k,
    q_stat = q_stat,
    q_df = df,
    q_p = stats::pchisq(q_stat, df, lower.tail = FALSE),
    phi = q_stat / k,
    phi_winsorised = sum(winsorise(z, trim)^2) / k,
    # On phi's own scale. phi carries divisor k while Q is chi-square on k-1 df,
    # so the threshold must carry the same divisor or the comparison is off by
    # k/(k-1) -- 20% at six units.
    phi_crit_95 = stats::qchisq(0.95, df) / k,
    tau2 = tau2,
    mean_binomial_var = mean(v),
    # The share of between-group spread that is NOT sampling error. This is what
    # the funnel's volume axis exists to ask -- "are the small units just
    # noisy?" -- answered as a number rather than left to the eye. RELATED TO
    # BUT NOT EQUAL TO ma_i2: this uses mean(v), where I-squared uses Higgins'
    # typical within-study variance. Measured gap on the 2026-10-26 UCMC file:
    # 0.2 to 1.7 percentage points. Do not report this as I-squared.
    share_not_sampling = tau2 / (tau2 + mean(v)),
    # Groups with zero variance (a rate at exactly 0 or 100) carry infinite
    # weight and are dropped from this family only. ma_k < k says so, and the
    # I-squared identity must be checked against ma_k - 1, not k - 1.
    ma_k = k_ma,
    ma_q = ma_q,
    ma_i2 = i2,
    ma_i2_lo = i2_ci[1],
    ma_i2_hi = i2_ci[2],
    ma_h2 = ma_q / df_ma,
    ma_tau2 = ma_tau2,
    ma_re_pct = theta_re,
    ma_re_lo = theta_re - tc * se_hk,
    ma_re_hi = theta_re + tc * se_hk,
    ma_dl_lo = theta_re - zc * se_dl,
    ma_dl_hi = theta_re + zc * se_dl,
    # COMPUTED ALWAYS, REPORTABLE ABOVE A THRESHOLD. tau-squared needs roughly
    # 5-10 groups to be stable and random-effects intervals are
    # anti-conservative below that even with Hartung-Knapp. The numbers are
    # never hidden -- only their USE is governed -- so they can be watched
    # stabilising as sites arrive, and the day k reaches the threshold nothing
    # is built: this flag flips.
    reportable = k >= min_k,
    stringsAsFactors = FALSE)
}


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
# the overall row only. Config: protocol.json disclosure.full_levels_strata_scope
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
# to within one bin width. Config: protocol.json pooling.histograms.
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
  # A proportion at exactly 0% or 100% makes the binomial variance zero, so the
  # design effect is 0/0 or x/0 -- NaN or Inf, both of which would ship as a
  # number and both of which abort the `deff >= 1` guard downstream with
  # "missing value where TRUE/FALSE needed". NA is the honest value: the
  # quantity is undefined, not large. Inert at UCMC, where every unit pairs at
  # least once, and live at any site with a unit that never did.
  naive <- out$pct_paired * (100 - out$pct_paired) / out$n_events
  out$deff <- ifelse(naive > 0,
                     apply(reps, 2, stats::var, na.rm = TRUE) / naive,
                     NA_real_)

  p0 <- 100 * sum(point_k) / sum(point_n)
  attr(out, "overall_pct") <- p0
  naive0 <- p0 * (100 - p0) / sum(point_n)
  attr(out, "deff") <- if (naive0 > 0) stats::var(overall) / naive0 else NA_real_
  out
}


# Collapse category levels for display through a map declared in protocol.json.
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
