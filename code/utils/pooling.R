# Federated-pooling exports, shared by the descriptive and landmark scripts.
#
# A median cannot be pooled across sites; a mean can, exactly, which is why sum
# and sum_sq travel rather than mean and sd:
#     mean_pooled = sum(sum) / sum(n)
#     var_pooled  = (sum(sum_sq) - sum(sum)^2 / sum(n)) / (sum(n) - 1)
# Medians travel alongside because the doses are zero-inflated: pool the means,
# report the medians. Cells below reporting.small_cell_min_den are suppressed.

# ---- Between-group dispersion, in TWO families that must not be conflated -----
#
# Both ask "do these groups differ by more than sampling error" and share most of
# their algebra, so read the family name before using either.
#
# FUNNEL FAMILY (no prefix) uses the NULL binomial variance p0(100-p0)/n on the
#   percentage scale, the scale the envelope is drawn on (Spiegelhalter, Qual
#   Saf Health Care 2005;14:347-51). Its DL weights go as n, so a weighted mean
#   is p0 by construction: no random-effects estimate belongs here.
# META-ANALYTIC FAMILY (ma_) uses each group's OWN p_i(100-p_i)/n_i, the only
#   input under which I-squared, tau-squared and random effects mean what the
#   literature means. Use these when reporting heterogeneity.
#
# Takes vectors, not a frame: the groups are units within a site, or sites.

winsorise <- function(z, trim) {
  if (trim <= 0) return(z)
  q <- stats::quantile(z, c(trim, 1 - trim), names = FALSE)
  pmin(pmax(z, q[1]), q[2])
}

# I-squared and its test-based interval, via H. Higgins & Thompson, Stat Med
# 2002;21:1539-58. I2 = (Q - df)/Q direct; only the INTERVAL needs H.
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
  # Zero-variance groups carry infinite weight; dropped here, visible in ma_k.
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
  # Hartung-Knapp-Sidik-Jonkman: t reference and a variance correction, for SMALL
  # k -- at k = 3 the multiplier is 4.3 against the normal's 1.96.
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
    # On phi's own scale: phi divides by k, Q is chi-square on k-1 df.
    phi_crit_95 = stats::qchisq(0.95, df) / k,
    tau2 = tau2,
    mean_binomial_var = mean(v),
    # The share of between-group spread that is NOT sampling error. Related to but
    # NOT equal to ma_i2, which uses Higgins' typical within-study variance.
    share_not_sampling = tau2 / (tau2 + mean(v)),
    # ma_k < k means groups were dropped; check the I2 identity against ma_k - 1.
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
    # Computed always, reportable only above a threshold: tau-squared needs 5-10
    # groups to be stable, so its USE is governed, never its value.
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

# Level counts per stratum, `strata` as long as `v`, plus an "overall" row.
# `strata_scope` (protocol.json disclosure.full_levels_strata_scope): "all" emits
# overall plus each named stratum, "overall" only the overall row -- a full and
# a collapsed view over the same strata make suppressed cells recoverable.
pool_cat <- function(variable, v, strata, min_cell, strata_scope = "all") {
  # as.character FIRST: c("overall", <factor>) coerces the factor to its INTEGER
  # CODES, so every named stratum matches nothing. Table 1's `grp` is a factor.
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
    # A row of counts sums to its PUBLISHED denominator, so one suppressed cell is
    # denominator minus the rest; suppress the next-smallest too.
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
# A median has no closed form across sites, and the mean is no substitute for a
# zero-inflated variable. Counts on FIXED bins do pool: sum them, then read any
# quantile to within one bin width (config: protocol.json pooling.histograms).
# `spec` is one entry; returns one row per bin, lowest first, numeric edges.
pool_hist <- function(variable, v, spec) {
  v <- v[!is.na(v)]
  from <- as.numeric(spec$from)
  to   <- as.numeric(spec$to)
  by   <- as.numeric(spec$by)
  zero_bin <- isTRUE(spec$zero_bin)
  open_top <- isTRUE(spec$open_top)
  # round(): seq() accumulates float error at 0.05 steps, and the bins are the
  # federation contract, so two sites' edges must compare equal.
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
    # Half-open [a, b), per windows.interval_convention -- except the final
    # interval of a BOUNDED variable, which must close or its maximum falls out.
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

  # A value that lands in no bin is silently lost exposure.
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
# narrow: resample EPISODES, and report the design effect beside the interval.
# `n_resamples` and `episodes` are REQUIRED, being properties of the calling
# analysis -- a silent default leaks one caller's universe into another's.
boot_group <- function(d, group_col, seed, outcome, episodes, n_resamples) {
  n_ep <- length(episodes)
  stopifnot("an event belongs to no episode in the resampling universe" =
              all(d$encounter_block %in% episodes),
            "n_resamples must be a positive whole number" =
              is.numeric(n_resamples) && n_resamples > 0)
  # Per (episode, group) counts, so a resample is a weighted sum.
  g <- factor(d[[group_col]])
  cell <- data.frame(ei = match(d$encounter_block, episodes),
                     gi = as.integer(g), n = 1L, k = as.integer(d[[outcome]]))
  cell <- aggregate(cbind(n, k) ~ ei + gi, data = cell, FUN = sum)
  G <- nlevels(g)

  point_n <- numeric(G); point_k <- numeric(G)
  s <- rowsum(as.matrix(cell[, c("k", "n")]), cell$gi, reorder = TRUE)
  idx <- as.integer(rownames(s))
  point_k[idx] <- s[, "k"]; point_n[idx] <- s[, "n"]

  # Re-seeded immediately before the draw: RNG consumed since the top of the
  # script would otherwise change which intervals a reader sees.
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
  # PER-GROUP design effect, read by SINGLE-GROUP callers (05_titration.R passes
  # g = "all"). Multi-group callers overwrite it with the overall attribute: a
  # per-unit deff would only serve pooling a NAMED unit across sites, which
  # docs/cohort_and_outputs.md section 8.7 forbids. A rate at exactly 0% or 100%
  # makes the naive variance zero, so NA ships rather than a NaN or Inf that
  # would abort the `deff >= 1` guard downstream.
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
  # A NAMED VECTOR, not list subsetting: `unlist(map[v])` drops NULLs for
  # unmatched values, so ifelse recycles a short vector and misaligns rows.
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
# per class. Layout follows Chen et al., BMC Infect Dis 2026;26:950 Table 1.
# stats: one row per model with an `ng` column; sizes: a named list, one
# element per model, each a vector of class counts.
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
# Key builder, local so these functions carry no hidden globals: a generalised
# function that reads an enclosing MIN_K or LABEL_COLS gives a silently wrong
# answer when called from anywhere else (lessons.md #26).
meta_key <- function(d, keys) {
  do.call(paste, c(lapply(keys, function(k) as.character(d[[k]])), sep = "\r"))
}

# Inverse-variance fixed-effect meta-analysis of the exposure coefficients.
# ON THE LOG-ODDS, never the odds ratio: exp() does not commute with a weighted
# mean. Only `poolable` rows combine -- ICU fixed effects are blanked at the site
# and thresholds are not shipped, so the pooled estimand is the within-site,
# ICU-adjusted association. See pooling_spec.json coefficient_meta._why.
pool_coefficient_meta <- function(d, keys, spec, rel, min_k,
                                  label_cols = character(0),
                                  not_poolable = character(0),
                                  key_fn = meta_key) {
  bcol <- spec[["coefficient"]]; scol <- spec[["standard_error"]]
  stopifnot("coefficient_meta needs coefficient and standard_error declared" =
              !is.null(bcol) && !is.null(scol))
  if (!all(c(bcol, scol) %in% names(d))) return(NULL)

  ok <- rep(TRUE, nrow(d))
  if ("poolable" %in% names(d)) ok <- toupper(as.character(d$poolable)) == "TRUE"
  b  <- suppressWarnings(as.numeric(d[[bcol]]))
  se <- suppressWarnings(as.numeric(d[[scol]]))
  ok <- ok & is.finite(b) & is.finite(se) & se > 0
  if (!any(ok)) return(NULL)
  dd <- d[ok, , drop = FALSE]; b <- b[ok]; se <- se[ok]

  k <- key_fn(dd, keys)
  split_by <- split(seq_along(b), k)
  rows <- lapply(names(split_by), function(g) {
    i <- split_by[[g]]
    w <- 1 / se[i]^2
    sw <- sum(w)
    est <- sum(w * b[i]) / sw
    sep <- sqrt(1 / sw)
    # Cochran's Q and DerSimonian-Laird tau^2. DESCRIPTIVE ONLY below
    # heterogeneity_min_k: Q has almost no power at small k, so the fixed-effect
    # estimate stays primary and the role is stated in its own column rather than
    # left to a reader comparing k against a number in another file.
    Q <- sum(w * (b[i] - est)^2)
    df <- length(i) - 1L
    I2 <- if (df > 0 && Q > df) 100 * (Q - df) / Q else 0
    C <- sw - sum(w^2) / sw
    tau2 <- if (df > 0 && C > 0) max(0, (Q - df) / C) else 0
    kd <- dd[i[1], intersect(c(keys, label_cols), names(dd)), drop = FALSE]
    data.frame(
      kd, site = "POOLED",
      estimate = round(est, 6), std_error = round(sep, 6),
      z = round(est / sep, 4),
      p_value = signif(2 * stats::pnorm(-abs(est / sep)), 4),
      odds_ratio = round(exp(est), 4),
      ci_lo = round(est - 1.959964 * sep, 6),
      ci_hi = round(est + 1.959964 * sep, 6),
      or_ci_lo = round(exp(est - 1.959964 * sep), 4),
      or_ci_hi = round(exp(est + 1.959964 * sep), 4),
      n_admissions = if ("n_admissions" %in% names(dd))
        sum(suppressWarnings(as.numeric(dd$n_admissions[i])), na.rm = TRUE) else NA_real_,
      k_sites = length(i), Q = round(Q, 4), Q_df = df,
      Q_p = if (df > 0) signif(stats::pchisq(Q, df, lower.tail = FALSE), 4) else NA_real_,
      I2_pct = round(I2, 1), tau2 = round(tau2, 6),
      heterogeneity_role = if (length(i) < min_k) "descriptive_only" else "reportable",
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)

  blanked <- setdiff(names(dd), c(names(out), "site"))
  numeric_blanked <- blanked[vapply(blanked, function(nm) is.numeric(dd[[nm]]), logical(1))]
  undeclared <- setdiff(numeric_blanked, c(not_poolable, unlist(spec[["not_poolable"]])))
  if (length(undeclared))
    stop(rel, ": numeric column(s) ", paste(undeclared, collapse = ", "),
         " would be blanked on the POOLED row but are not declared in ",
         "not_poolable_columns.", call. = FALSE)
  for (nm in blanked) out[[nm]] <- NA

  # The covariance pool must use the SAME weights, so they travel with the result
  # rather than being recomputed and silently drifting.
  attr(out, "weights") <- stats::setNames(
    lapply(names(split_by), function(g) {
      i <- split_by[[g]]
      stats::setNames(1 / se[i]^2, as.character(dd$site[i]))
    }),
    vapply(names(split_by), function(g) {
      i <- split_by[[g]][1]
      paste(dd$model[i], dd$term[i], sep = "\r")
    }, character(1)))
  cat(sprintf("  %s: pooled %d term(s) over %d site(s)%s\n", rel, nrow(out),
              length(unique(dd$site)),
              if (length(unique(dd$site)) < min_k)
                sprintf("; below heterogeneity_min_k = %d, tau^2 descriptive only", min_k)
              else ""))
  out
}

# The 2x2 between the two exposures, pooled so the centre can test them jointly.
# Under site independence, with the SAME weights the estimates used:
#   Cov(b1_bar, b2_bar) = sum_i(w1_i * w2_i * Cov_i) / (sum w1 * sum w2)
# The diagonal must reproduce the pooled variances, which is asserted below.
pool_covariance_meta <- function(d, spec, rel, coef_pooled) {
  need <- c("model", "term_row", "term_col", "covariance")
  if (!all(need %in% names(d))) return(NULL)
  if (is.null(coef_pooled) || !nrow(coef_pooled)) return(NULL)
  wmap <- attr(coef_pooled, "weights")
  if (is.null(wmap)) return(NULL)
  keyed <- paste(coef_pooled$model, coef_pooled$term, sep = "\r")

  rows <- list()
  for (m in unique(d$model)) {
    dm <- d[d$model == m, , drop = FALSE]
    for (a in unique(dm$term_row)) for (bb in unique(dm$term_col)) {
      sel <- dm[dm$term_row == a & dm$term_col == bb, , drop = FALSE]
      wa <- wmap[[paste(m, a, sep = "\r")]]; wb <- wmap[[paste(m, bb, sep = "\r")]]
      if (is.null(wa) || is.null(wb)) next
      s_use <- intersect(intersect(names(wa), names(wb)), as.character(sel$site))
      if (!length(s_use)) next
      num <- sum(vapply(s_use, function(s) {
        cv <- suppressWarnings(as.numeric(sel$covariance[as.character(sel$site) == s][1]))
        if (!is.finite(cv)) 0 else wa[[s]] * wb[[s]] * cv
      }, numeric(1)))
      den <- sum(wa[s_use]) * sum(wb[s_use])
      rows[[length(rows) + 1L]] <- data.frame(
        model = m, term_row = a, term_col = bb, site = "POOLED",
        covariance = round(num / den, 12), k_sites = length(s_use),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  diagr <- out[out$term_row == out$term_col, , drop = FALSE]
  for (i in seq_len(nrow(diagr))) {
    j <- which(keyed == paste(diagr$model[i], diagr$term_row[i], sep = "\r"))
    if (!length(j)) next
    stopifnot("the pooled covariance diagonal disagrees with the pooled standard error" =
                abs(sqrt(diagr$covariance[i]) - coef_pooled$std_error[j[1]]) < 1e-5)
  }
  for (nm in setdiff(names(d), c(names(out), "site"))) out[[nm]] <- NA
  cat(sprintf("  %s: pooled %d covariance cell(s)\n", rel, nrow(out)))
  out
}

# ---- Aggregate-data meta-analysis of medians (McGrath 2019) -------------------
# references/Statistics in Medicine - 2018 - McGrath - One-sample aggregate data
# meta-analysis of medians.pdf, section 2.2. Hand-rolled rather than taking
# metamedian + Hmisc: this is fifteen lines of arithmetic against a definition in
# print, and the pinned stack is the reproducibility guarantee.

# Weighted quantile: McGrath's step definition with linear interpolation between
# consecutive order statistics, on MIDPOINT plotting positions -- the cumulative
# weight BEFORE x_i plus half of x_i's own weight.
#
# The midpoint is load-bearing, not a refinement. Interpolating on the plain
# cumulative weight puts three equally weighted studies at 1.5 instead of the
# middle value: cw = (1/3, 2/3, 1), q = 0.5 lands a third of the way between the
# first and second. With midpoints the positions are (1/6, 1/2, 5/6), q = 0.5
# sits exactly on the second, and two equally weighted studies interpolate to
# their midpoint -- both matching the unweighted median they reduce to.
weighted_quantile <- function(x, w, q) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (!length(x)) return(NA_real_)
  if (length(x) == 1L) return(x[1])
  o <- order(x); x <- x[o]; w <- w[o] / sum(w[o])
  pos <- cumsum(w) - w / 2
  if (q <= pos[1]) return(x[1])
  if (q >= pos[length(pos)]) return(x[length(x)])
  j <- which(pos >= q)[1]
  t <- (q - pos[j - 1]) / (pos[j] - pos[j - 1])
  x[j - 1] + t * (x[j] - x[j - 1])
}

# Weighted median of medians, weights proportional to each study's n.
# THE INTERVAL DEGENERATES AT SMALL k AND THAT IS NOT A BUG: the limits are the
# 1/2 -+ min(1/2, z/(2*sqrt(k))) weighted quantiles, and z/(2*sqrt(k)) exceeds
# 1/2 until k >= 5 -- so at k = 2, 3 or 4 the interval IS (essentially) the full
# range of the study medians. Report it as a range at that k rather than as an
# interval that narrows with evidence. McGrath's own simulations used k = 15 or
# 50, so three sites is outside the range the method was evaluated on.
wm_of_medians <- function(medians, n, conf = 0.95) {
  ok <- is.finite(medians) & is.finite(n) & n > 0
  m <- medians[ok]; nn <- n[ok]
  k <- length(m)
  if (!k) return(list(estimate = NA_real_, lo = NA_real_, hi = NA_real_,
                      k = 0L, min = NA_real_, max = NA_real_, degenerate = NA))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  half <- min(0.5, z / (2 * sqrt(k)))
  lo <- weighted_quantile(m, nn, 0.5 - half)
  hi <- weighted_quantile(m, nn, 0.5 + half)
  list(estimate = weighted_quantile(m, nn, 0.5), lo = lo, hi = hi,
       k = k, min = min(m), max = max(m),
       # MEASURED, not predicted from k: with midpoint positions the interval
       # stays pinned to the range for longer than the half < 1/2 condition
       # suggests, because the requested quantile can still fall outside the
       # first and last plotting positions.
       degenerate = isTRUE(all.equal(lo, min(m))) &&
                    isTRUE(all.equal(hi, max(m))))
}
