# Table 1 row builders, shared so one generator owns the table.
#
# Moved out of 02_descriptive_cohort.R on 2026-10-28 when 07 became the sole
# generator: 07 holds the exposures and the fitted model frame, and it reads
# trajectory_long window 0 unrestricted, so it can describe the whole cohort too.
#
# EVERY SCOPE IS PASSED IN, never read from an enclosing frame. These used to
# close over `grp`, `STRATA`, `POOL` and `MIN_CELL`; a shared helper that reads
# its caller's globals gives a silently wrong answer from a second caller
# (lessons.md #26).
#
# Each builder returns BOTH halves at once -- the formatted display row and the
# tidy poolable rows -- because they must describe the same numbers. Computing
# them in two places is how a table and its pooling export drift apart.
#
# Needs: pool_row() and pool_cat() from utils/pooling.R.

t1_fmt_iqr <- function(v, digits = 1) {
  v <- v[!is.na(v)]
  if (!length(v)) return("")
  q <- unname(stats::quantile(v, c(0.25, 0.5, 0.75)))
  sprintf("%.*f (%.*f, %.*f)", digits, q[2], digits, q[1], digits, q[3])
}

t1_fmt_pct <- function(k, n) {
  if (!is.finite(k) || n <= 0) return("")
  sprintf("%s (%.1f%%)", format(k, big.mark = ","), 100 * k / n)
}

t1_header <- function(label, n) sprintf("%s (N = %s)", label, format(n, big.mark = ","))

# values: a NAMED list of numeric vectors, one per scope, in display order. A
# scope that does not apply to this row is simply absent -- that is how VFD-28
# and the two exposures stay out of the Overall column instead of being computed
# over a denominator that does not apply to them.
t1_continuous <- function(label, values, unit = NA_character_, digits = 1,
                          min_cell = 10, scope_order = names(values)) {
  stopifnot("t1_continuous needs a named list of scopes" =
              !is.null(names(values)) && all(nzchar(names(values))))
  pool <- do.call(rbind, lapply(names(values), function(sc)
    pool_row("baseline", label, unit, sc, NA, NA, values[[sc]], min_cell)))
  disp <- data.frame(characteristic = sprintf("__%s__, median (IQR)", label),
                     stringsAsFactors = FALSE)
  for (sc in scope_order)
    disp[[sc]] <- if (is.null(values[[sc]])) "" else t1_fmt_iqr(values[[sc]], digits)
  list(display = disp, pool = pool, pool_cat = NULL)
}

# levels_shown: the levels that get a display row, in order. A level present in
# the data but not listed still SHIPS in the tidy export -- it is only the
# display that is curated, so a pooled count is never lost to a display choice.
t1_categorical <- function(label, values, collapse = NULL, min_cell = 10,
                           scope_order = names(values), levels_shown = NULL) {
  stopifnot("t1_categorical needs a named list of scopes" =
              !is.null(names(values)) && all(nzchar(names(values))))
  raw <- lapply(values, function(v) ifelse(is.na(v), "Missing", as.character(v)))
  shown <- if (is.null(collapse)) raw else lapply(raw, collapse_levels, map = collapse)

  # pool_cat is built per scope and stamped with it, rather than letting pool_cat
  # derive strata: the scopes here are NESTED, so a single "strata" vector cannot
  # express them -- one admission belongs to several.
  # NAMING MATCHES 02's: the BASE label carries the raw mCIDE levels and
  # "<label> (collapsed)" carries the display collapse. Inverting it silently
  # breaks every consumer keyed on the variable name, tests included.
  tidy_for <- function(vals, nm) do.call(rbind, lapply(names(vals), function(sc) {
    r <- pool_cat(nm, vals[[sc]], rep("x", length(vals[[sc]])), min_cell,
                  strata_scope = "overall")
    r$stratum <- sc
    r
  }))
  pc <- tidy_for(raw, label)
  if (!is.null(collapse)) {
    pcc <- tidy_for(shown, paste(label, "(collapsed)"))
    # Cross-view suppression. pool_cat suppresses within one view, but a full and
    # a collapsed view of the same scope are two equations over the same cells:
    # a published collapsed level whose sources are all published but one gives
    # that one by subtraction (Other = 4017 exposed a suppressed level as 67, and
    # the row total then exposed Missing). So a collapsed level is withheld
    # whenever any raw level feeding it is suppressed.
    for (sc in unique(pc$stratum)) {
      src <- pc[pc$variable == label & pc$stratum == sc &
                  pc$n_suppressed_small_cell == 1L, "level"]
      hit <- unique(collapse_levels(src, map = collapse))
      i <- pcc$stratum == sc & pcc$level %in% hit
      pcc$n[i] <- NA_integer_
      pcc$pct[i] <- NA_real_
      pcc$n_suppressed_small_cell[i] <- 1L
    }
    pc <- rbind(pc, pcc)
  }

  lv <- if (!is.null(levels_shown)) levels_shown else
    sort(unique(unlist(lapply(shown, unique))))
  head <- data.frame(characteristic = sprintf("__%s__, n (%%)", label),
                     stringsAsFactors = FALSE)
  for (sc in scope_order) head[[sc]] <- ""
  rows <- lapply(lv, function(l) {
    r <- data.frame(characteristic = sprintf("    %s", l), stringsAsFactors = FALSE)
    for (sc in scope_order) {
      x <- shown[[sc]]
      r[[sc]] <- if (is.null(x)) "" else t1_fmt_pct(sum(x == l), length(x))
    }
    r
  })
  list(display = do.call(rbind, c(list(head), rows)), pool = NULL, pool_cat = pc)
}

# Assemble: stack the display rows, rename the scope columns to their labels with
# each column's own N, and return the three pieces the caller ships.
t1_assemble <- function(parts, scope_order, scope_labels, scope_n) {
  disp <- do.call(rbind, lapply(parts, function(p) p$display))
  names(disp)[match(scope_order, names(disp))] <-
    vapply(scope_order, function(sc) t1_header(scope_labels[[sc]], scope_n[[sc]]),
           character(1))
  list(display = disp,
       continuous = do.call(rbind, Filter(Negate(is.null), lapply(parts, function(p) p$pool))),
       categorical = do.call(rbind, Filter(Negate(is.null), lapply(parts, function(p) p$pool_cat))))
}
