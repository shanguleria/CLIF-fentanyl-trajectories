# Purpose : Combine N site bundles into one set of pooled numbers. Run at the
#           coordinating centre, never at a site.
# Inputs  : one directory per site, each a copy of that site's output/final_no_phi/
# Outputs : output/pooled/ : one CSV per input file, plus heterogeneity.csv
#
#   Rscript code/pool_sites.R multisite/ucmc multisite/rush multisite/nu
#
# Output is long with a `site` column; site == "POOLED" is the combined row.
# Percentages are recomputed from summed counts, never averaged. Unit-keyed rows,
# missingness and disagreeing protocols are not pooled. The per-file method is
# iterated from config/pooling_spec.json, so a new export is never silently skipped.

suppressMessages({
  library(here)
  library(jsonlite)
})
source(here("code", "utils", "paths.R"))
source(here("code", "utils", "pooling.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  cat("usage: Rscript code/pool_sites.R <site_dir> [<site_dir> ...]\n\n",
      "Each <site_dir> is a copy of that site's output/final_no_phi/.\n",
      "Run code/check_provenance.py over the same directories first.\n", sep = "")
  quit(status = 1)
}
# `multisite/*` globs the README in there too, so skip a non-directory argument.
skipped <- args[!dir.exists(args)]
if (length(skipped))
  cat(sprintf("  skipping non-directory argument(s): %s\n",
              paste(basename(skipped), collapse = ", ")))
args <- args[dir.exists(args)]
if (!length(args))
  stop("no directories among the arguments. Each must be a copy of a site's ",
       "output/final_no_phi/ -- see multisite/README.md.", call. = FALSE)

SPEC <- fromJSON(here("config", "pooling_spec.json"), simplifyVector = FALSE)
FILES <- SPEC$files
MIN_K <- as.integer(SPEC$heterogeneity_min_k)
LABEL_COLS <- unlist(SPEC$label_columns)
NOT_POOLABLE <- unlist(SPEC$not_poolable_columns)
OUT <- here("output", "pooled")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
# Cleared every run: the summary step re-reads these files by name, so one left from
# a previous pool is rendered as this run's number. Safe only before any write.
stale <- list.files(OUT, full.names = TRUE)
if (length(stale)) {
  unlink(stale, recursive = TRUE)
  message(sprintf("  cleared %d file(s) from a previous pool", length(stale)))
}


# ---- 1. Read the bundles -----------------------------------------------------

read_bundle <- function(path) {
  # A site that zips output/ sends final_no_phi's parent; descend, do not refuse.
  if (!file.exists(file.path(path, "manifest.json")) &&
      file.exists(file.path(path, "final_no_phi", "manifest.json"))) {
    cat(sprintf("  NOTE %s: descending into final_no_phi/\n", basename(path)))
    path <- file.path(path, "final_no_phi")
  }
  prov_f <- file.path(path, "manifest.json")
  if (!file.exists(prov_f))
    stop(path, " has no manifest.json -- is it really an output/final_no_phi/ copy?",
         call. = FALSE)
  prov <- fromJSON(prov_f, simplifyVector = FALSE)
  csvs <- list.files(path, pattern = "\\.csv$", recursive = TRUE)
  tabs <- lapply(file.path(path, csvs), function(f)
    utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE))
  names(tabs) <- csvs
  list(site = prov$site_name, path = path, prov = prov, tabs = tabs)
}

bundles <- lapply(args, read_bundle)
sites <- vapply(bundles, function(b) b$site, character(1))
if (anyDuplicated(sites))
  stop("two bundles report the same site_name: ", paste(sites, collapse = ", "),
       ". Pooling a site with itself would double its weight.", call. = FALSE)
K <- length(bundles)
cat(sprintf("\n%d site(s): %s\n", K, paste(sites, collapse = ", ")))


# ---- 2. Refuse before pooling ------------------------------------------------
# Protocol divergence is a hard stop. Code divergence is reported, not refused:
# a site may legitimately extend _to_mcg_hr for a locally charted dose unit.

one_of <- function(f) unique(vapply(bundles, function(b) {
  v <- f(b$prov); if (is.null(v)) "" else paste(unlist(v), collapse = "|")
}, character(1)))

bad <- list(
  protocol_digests = one_of(function(p) p$protocol_digests),
  definition_versions = one_of(function(p) p$definition_versions),
  clif_version = one_of(function(p) p$clif_version))
diverged <- names(bad)[vapply(bad, length, integer(1)) > 1L]
if (length(diverged)) {
  for (d in diverged)
    cat(sprintf("  PROTOCOL DIVERGENCE  %s:\n%s\n", d,
                paste("     ", bad[[d]], collapse = "\n")))
  stop("sites did not run the same protocol, so a pooled number would be ",
       "meaningless. Resolve the divergence, or re-run the divergent site. ",
       "Run code/check_provenance.py for a per-file breakdown.", call. = FALSE)
}
cat("  protocol OK  all sites ran the same protocol and definitions\n")

code_digests <- unique(vapply(bundles, function(b) b$prov$code_digest %||% "",
                              character(1)))
if (length(code_digests) > 1L)
  cat(sprintf("  NOTE code divergence across sites (%s). Reported, not refused;\n",
              paste(substr(code_digests, 1, 8), collapse = ", ")),
      "       run check_provenance.py to see WHICH files differ.\n", sep = "")

# An unclaimed file is never pooled; a phantom entry describes a dead export.
present <- sort(unique(unlist(lapply(bundles, function(b) names(b$tabs)))))
claimed <- names(FILES)
unclaimed <- setdiff(present, claimed)
phantom <- setdiff(claimed, present)
if (length(unclaimed))
  stop("these files are in a bundle and claimed by no entry in ",
       "config/pooling_spec.json: ", paste(unclaimed, collapse = ", "),
       ". Add an entry -- a shipped export that nothing pools is the silent ",
       "no-op the consumption contract exists to catch.", call. = FALSE)
if (length(phantom))
  cat(sprintf("  NOTE declared but absent from every bundle: %s\n",
              paste(phantom, collapse = ", ")))
cat(sprintf("  files OK     %d CSV(s), all claimed by the pooling spec\n",
            length(present)))


# ---- 3. The methods ----------------------------------------------------------

# One row per key per site plus a POOLED row. Real names; anonymised at figure time.
stack_sites <- function(rel, keys) {
  fr <- Filter(Negate(is.null), lapply(bundles, function(b) {
    d <- b$tabs[[rel]]
    if (is.null(d)) return(NULL)
    cbind(site = b$site, d, stringsAsFactors = FALSE)
  }))
  if (!length(fr)) return(NULL)
  # Union the columns rather than requiring a match: baseline_characteristics.csv
  # embeds each site's own N in its column names, so a plain rbind fails.
  all_cols <- unique(unlist(lapply(fr, names)))
  fr <- lapply(fr, function(d) {
    for (nm in setdiff(all_cols, names(d))) d[[nm]] <- NA
    d[, all_cols, drop = FALSE]
  })
  do.call(rbind, fr)
}

key_of <- function(d, keys) {
  if (is.null(keys) || !length(keys)) return(rep("", nrow(d)))
  do.call(paste, c(lapply(keys, function(k) as.character(d[[k]])), sep = "\r"))
}

# Sum the summables, carry the labels, recompute the declared percentages, leave
# everything else NA. label_columns, not_poolable_columns and the per-file
# `percents` map are declared in config/pooling_spec.json, not inferred.
pool_sum_counts <- function(d, keys, spec, rel) {
  # Spec keys are read with [[ ]], never $: $ partial-matches on lists, and once
  # returned denominator_within for denominator.
  nums <- unlist(spec[["numerators"]]); den <- spec[["denominator"]]
  pcts <- spec[["percents"]]; within <- unlist(spec[["denominator_within"]])
  sum_cols <- intersect(c(nums, den), names(d))
  if (!length(sum_cols)) return(NULL)
  k <- key_of(d, keys)
  first <- !duplicated(k)

  # na.rm = TRUE would turn an all-NA group into a false 0; strobe's first row has
  # no n_before by construction.
  sum_keep_na <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
  agg <- lapply(sum_cols, function(cc)
    tapply(suppressWarnings(as.numeric(d[[cc]])), k, sum_keep_na))
  names(agg) <- sum_cols
  out <- as.data.frame(agg, stringsAsFactors = FALSE)

  kd <- d[first, intersect(c(keys, LABEL_COLS), names(d)), drop = FALSE]
  rownames(kd) <- k[first]
  out <- cbind(kd[rownames(out), , drop = FALSE], out)
  out$site <- "POOLED"

  # Percentages are recomputed from the summed counts, into their own column name.
  for (pc in names(pcts)) {
    num <- pcts[[pc]]
    if (!(num %in% names(out))) next
    denom <- if (length(within) && !identical(within, "__first_window__")) {
      g <- key_of(out, intersect(within, names(out)))
      as.vector(tapply(out[[num]], g, sum, na.rm = TRUE)[g])
    } else if (identical(unname(within), "__first_window__")) {
      # retention's denominator: this file's window 0 count, fixed at every window.
      rep(out[[num]][which.min(suppressWarnings(as.numeric(out$window_idx)))], nrow(out))
    } else if (!is.null(den) && den %in% names(out)) {
      out[[den]]
    } else NULL
    recomputed <- if (is.null(denom)) NA_real_ else 100 * out[[num]] / denom
    # If every site shipped NA for this percentage in a key group the blank was
    # deliberate (dose_summary receivers_only), so the pooled row blanks too.
    if (pc %in% names(d)) {
      site_all_na <- tapply(is.na(suppressWarnings(as.numeric(d[[pc]]))),
                            key_of(d, keys), all)
      blank <- site_all_na[rownames(out)]
      recomputed[!is.na(blank) & blank] <- NA_real_
    }
    # Round to the precision the sites shipped for this column, read off their own
    # values rather than fixed (1 decimal in some files, 3 in the histograms).
    dec <- max(c(0L, vapply(as.character(d[[pc]]), function(v) {
      v <- sub("^.*\\.", "", v[grepl("\\.", v)])
      if (!length(v)) 0L else nchar(sub("0+$", "", v))
    }, integer(1))), na.rm = TRUE)
    out[[pc]] <- round(recomputed, dec)
  }
  # DECLARED weighted means. A mean does not pool from means -- but mean * n IS
  # the sum, so where the file ships the weight beside the mean the pooled value
  # is exact rather than an average of averages. Declared per file in
  # weighted_means as {column: weight_column}; never inferred.
  for (mc in names(spec[["weighted_means"]])) {
    wc <- spec[["weighted_means"]][[mc]]
    if (!all(c(mc, wc) %in% names(d))) next
    x <- suppressWarnings(as.numeric(d[[mc]]))
    w <- suppressWarnings(as.numeric(d[[wc]]))
    use <- is.finite(x) & is.finite(w) & w > 0
    num <- tapply(ifelse(use, x * w, NA_real_), k, sum_keep_na)
    den <- tapply(ifelse(use, w, NA_real_), k, sum_keep_na)
    dec <- max(c(0L, vapply(as.character(d[[mc]]), function(v) {
      v <- sub("^.*\\.", "", v[grepl("\\.", v)])
      if (!length(v)) 0L else nchar(sub("0+$", "", v))
    }, integer(1))), na.rm = TRUE)
    out[[mc]] <- round(as.vector(num[rownames(out)] / den[rownames(out)]), dec)
  }

  # Numeric columns neither summed, labels nor recomputed percents become NA.
  blanked <- setdiff(names(d), c(names(out), "site"))
  numeric_blanked <- blanked[vapply(blanked, function(nm)
    is.numeric(d[[nm]]), logical(1))]
  # Global not_poolable list, unioned with this file's own.
  undeclared <- setdiff(numeric_blanked,
                        c(NOT_POOLABLE, unlist(spec[["not_poolable"]])))
  if (length(undeclared))
    stop(rel, ": numeric column(s) ", paste(undeclared, collapse = ", "),
         " are being blanked in the POOLED row but are not declared in ",
         "not_poolable_columns. Either add them to this file's numerators so ",
         "they are summed, or declare them unpoolable with a reason. A column ",
         "that silently becomes NA is the gap this guard exists to find.",
         call. = FALSE)
  for (nm in blanked) out[[nm]] <- NA
  out
}

pool_sufficient <- function(d, keys, stats) {
  k <- key_of(d, keys)
  n <- tapply(as.numeric(d$n), k, sum, na.rm = TRUE)
  s <- tapply(as.numeric(d$sum), k, sum, na.rm = TRUE)
  ss <- tapply(as.numeric(d$sum_sq), k, sum, na.rm = TRUE)
  mean_p <- s / n
  # Exact from the sufficient statistics; assumes nothing about equal variances.
  var_p <- (ss - s^2 / n) / (n - 1)
  var_p[!is.finite(var_p) | var_p < 0] <- NA_real_
  first <- !duplicated(k)
  kd <- d[first, intersect(c(keys, LABEL_COLS), names(d)), drop = FALSE]
  rownames(kd) <- k[first]
  out <- cbind(kd[names(n), , drop = FALSE],
               data.frame(site = "POOLED", n = as.vector(n), sum = as.vector(s),
                          sum_sq = as.vector(ss), mean = as.vector(mean_p),
                          sd = sqrt(as.vector(var_p)), stringsAsFactors = FALSE))
  out
}

pool_bins <- function(d, keys, rel) {
  # Grouping comes from the declared keys, not a hardcoded column name: bins sum to
  # a total within whatever else the file is keyed on (`variable`, `exposure`).
  group_keys <- setdiff(keys, "bin")
  stopifnot("a sum_bins file must declare a `bin` key" = "bin" %in% keys)
  grp_of <- function(x) if (!length(group_keys)) rep("__all__", nrow(x))
                        else key_of(x, group_keys)

  # Summing bins is valid only if they partition identically, so check the edges.
  edge_cols <- intersect(c(group_keys, "bin", "bin_lo", "bin_hi"), names(d))
  edges <- unique(vapply(split(d, d$site), function(s)
    paste(do.call(paste, c(unname(as.list(s[, edge_cols, drop = FALSE])),
                           list(sep = "\u001f"))), collapse = ";"), character(1)))
  if (length(edges) > 1L)
    stop(rel, ": histogram bin edges differ between sites, so the bins do not ",
         "partition the data identically and summing them would be wrong.",
         call. = FALSE)
  k <- key_of(d, keys)
  n <- tapply(as.numeric(d$n), k, sum, na.rm = TRUE)
  first <- !duplicated(k)
  kd <- d[first, intersect(c(keys, LABEL_COLS), names(d)), drop = FALSE]
  rownames(kd) <- k[first]
  out <- cbind(kd[names(n), , drop = FALSE],
               data.frame(site = "POOLED", n = as.vector(n), stringsAsFactors = FALSE))
  g <- grp_of(out)
  tot <- tapply(out$n, g, sum, na.rm = TRUE)
  out$n_total <- as.vector(tot[as.character(g)])
  # Same precision rule as the counts pooler: match what the sites shipped.
  dec <- max(c(0L, vapply(as.character(d$pct), function(v) {
    v <- sub("^.*\\.", "", v[grepl("\\.", v)])
    if (!length(v)) 0L else nchar(sub("0+$", "", v))
  }, integer(1))), na.rm = TRUE)
  out$pct <- round(100 * out$n / out$n_total, dec)
  out
}


# ---- 4. Pool ------------------------------------------------------------------

written <- character(0)
write_out <- function(rel, d) {
  dup <- unique(names(d)[duplicated(names(d))])
  if (length(dup))
    stop(rel, ": pooled frame has duplicated column name(s) ",
         paste(dup, collapse = ", "), ". A CSV hides this -- the reader gets a ",
         "silently renamed column -- so it is refused here.", call. = FALSE)
  f <- file.path(OUT, gsub("/", "__", rel))
  utils::write.csv(d, f, row.names = FALSE)
  written <<- c(written, basename(f))
}

cat("\nPooling\n")
COEF_POOLS <- list()
for (rel in claimed) {
  spec <- FILES[[rel]]
  method <- spec[["method"]]
  keys <- unlist(spec[["keys"]])
  d <- stack_sites(rel, keys)
  if (is.null(d) || !nrow(d)) next

  pooled <- switch(
    method,
    sum_counts = pool_sum_counts(d, keys, spec, rel),
    sufficient_stats = pool_sufficient(d, keys, unlist(spec[["stats"]])),
    sum_bins = pool_bins(d, keys, rel),
    count_matrix = pool_sum_counts(d, keys, spec, rel),
    coefficient_meta = {
      r <- pool_coefficient_meta(d, keys, spec, rel, min_k = MIN_K,
                                 label_cols = LABEL_COLS,
                                 not_poolable = NOT_POOLABLE, key_fn = key_of)
      # The covariance pool needs these weights; stash under the file name the
      # covariance entry declares in paired_with, so the link is by declaration
      # and not by position in the loop.
      # `<-`, NOT `<<-`: this loop body runs in the global environment, and `<<-`
      # begins its search in the PARENT of the current scope -- so it would skip
      # COEF_POOLS entirely and fail with "object not found".
      COEF_POOLS[[rel]] <- r
      r
    },
    covariance_meta = {
      paired <- spec[["paired_with"]]
      stopifnot("covariance_meta needs paired_with declared" = !is.null(paired))
      cp <- COEF_POOLS[[paired]]
      if (is.null(cp)) {
        stop(rel, ": its paired_with file ", paired, " has not been pooled yet, so ",
             "the inverse-variance weights are unavailable. The coefficient entry ",
             "must come FIRST in pooling_spec.json's files object -- that ordering ",
             "is load-bearing and this is the guard that says so.", call. = FALSE)
      }
      pool_covariance_meta(d, spec, rel, cp)
    },
    # Everything else falls through to NULL: per-site rows, no POOLED row --
    # correct for side_by_side, per_site_only and heterogeneity_source.
    NULL)

  # POOLED rows come back in whatever order the key split produced, which is
  # ALPHABETICAL and therefore wrong for any key that is not: the ordinal bands
  # came out (0,7] | (14,21] | (21,25] | (7,14] | 0, putting (14,21] before
  # (7,14]. Reorder to the first-appearance order of the site rows, which is the
  # order the sites declared. Applies to every method, not just sum_counts.
  if (!is.null(pooled) && all(keys %in% names(pooled))) {
    want <- unique(key_of(d, keys))
    pooled <- pooled[order(match(key_of(pooled, keys), want)), , drop = FALSE]
  }

  out <- if (is.null(pooled)) d else {
    common <- union(names(d), names(pooled))
    for (nm in setdiff(common, names(d))) d[[nm]] <- NA
    for (nm in setdiff(common, names(pooled))) pooled[[nm]] <- NA
    rbind(d[, common, drop = FALSE], pooled[, common, drop = FALSE])
  }
  write_out(rel, out)
  cat(sprintf("  %-52s %-18s %4d row(s) x %d site(s)%s\n",
              rel, method, nrow(d), K,
              if (is.null(pooled)) "  [per-site only]" else " + POOLED"))
}


# ---- The readable summary ----------------------------------------------------
# Rendered from config/pooling_spec.json rather than retyped, so prose cannot drift.

md_table <- function(d, cols) {
  hdr <- paste0("| ", paste(cols, collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", length(cols)), collapse = "|"), "|")
  # Thousands separators on whole counts.
  fmt <- function(v) {
    x <- suppressWarnings(as.numeric(v))
    ifelse(is.na(x), ifelse(is.na(v) | v == "", "--", as.character(v)),
           ifelse(x == round(x) & abs(x) >= 1000,
                  formatC(x, format = "d", big.mark = ","), as.character(v)))
  }
  rows <- apply(d[, cols, drop = FALSE], 1, function(r) paste0("| ", paste(fmt(r), collapse = " | "), " |"))
  c(hdr, sep, rows)
}

write_summary <- function() {
  # Every percent-bearing entry must carry a description of what it counts.
  gaps <- Filter(function(rel) {
    e <- FILES[[rel]]
    !is.null(e[["percents"]]) && is.null(e[["_numerator"]])
  }, claimed)
  if (length(gaps))
    stop("these entries declare `percents` but no `_numerator` in ",
         "config/pooling_spec.json, so their pooled value would ship with no ",
         "statement of what it counts: ", paste(gaps, collapse = ", "),
         call. = FALSE)

  have <- Filter(function(rel) !is.null(FILES[[rel]][["_numerator"]]) &&
                   file.exists(file.path(OUT, gsub("/", "__", rel))), claimed)
  ord <- vapply(have, function(rel) {
    o <- FILES[[rel]][["_order"]]; if (is.null(o)) 999 else as.numeric(o)
  }, numeric(1))
  have <- have[order(ord)]

  L <- c(
    "# Pooled summary",
    "",
    sprintf("%d site(s): **%s**. Generated %s from pooling spec %s.",
            K, paste(sites, collapse = ", "),
            format(Sys.time(), "%Y-%m-%d %H:%M %Z"), SPEC$definition_version),
    "",
    "Every percentage below is **sum of numerators over sum of denominators**,",
    "never an average of the sites' percentages. The site exports carry both",
    "counts precisely so that no percentage is reconstructed from a rounded one.",
    "",
    "**Real site names appear here.** Anonymise before anything audience-facing.",
    "")

  for (rel in have) {
    e <- FILES[[rel]]
    d <- utils::read.csv(file.path(OUT, gsub("/", "__", rel)),
                         stringsAsFactors = FALSE, check.names = FALSE)
    keys <- intersect(unlist(e[["keys"]]), names(d))
    L <- c(L, sprintf("## %s", e[["_title"]]), "",
           sprintf("`%s` &mdash; grain: **%s**", rel, e[["_grain"]]), "",
           sprintf("- **Numerator:** %s", e[["_numerator"]]),
           sprintf("- **Denominator:** %s", e[["_denominator"]]), "")
    if (!is.null(e[["_caveat"]]))
      L <- c(L, sprintf("> %s", e[["_caveat"]]), "")

    # A coefficient file has no numerator, so the percents loop below writes
    # nothing for it -- and the pooled odds ratio is the headline result of the
    # whole study. Emit the exposure rows explicitly.
    if (identical(e[["method"]], "coefficient_meta") && "is_exposure" %in% names(d)) {
      ex <- d[toupper(as.character(d$is_exposure)) == "TRUE" &
                d$model == "primary", , drop = FALSE]
      if (nrow(ex)) {
        ex$`OR (95% CI)` <- sprintf("%s (%s to %s)", ex$odds_ratio,
                                    ex$or_ci_lo, ex$or_ci_hi)
        ex <- ex[order(match(ex$term, unique(ex$term)),
                       match(ex$site, c(sites, "POOLED"))), , drop = FALSE]
        L <- c(L, "**Primary model.** Pooled by inverse variance on the log odds;",
               "sensitivity models are in the CSV.", "",
               md_table(ex, c("term", "site", "OR (95% CI)", "p_value",
                              "n_admissions", "k_sites", "I2_pct",
                              "heterogeneity_role")), "")
      }
    }

    for (pc in names(e[["percents"]])) {
      num <- e[["percents"]][[pc]]
      if (!all(c(pc, num) %in% names(d))) next
      den <- e[["denominator"]]
      show <- c(keys, "site", num,
                if (!is.null(den) && den %in% names(d)) den, pc)
      nrow_per_site <- nrow(d) / max(1L, length(unique(d$site)))
      if (nrow_per_site > 24) {
        # Too many rows to read. Ship the pooled rows and point at the CSV.
        sub <- d[d$site == "POOLED", , drop = FALSE]
        L <- c(L, sprintf("`%s` has %d rows per site &mdash; POOLED rows only here; full table in the CSV.",
                          pc, round(nrow_per_site)), "",
               md_table(utils::head(sub, 12), show), "")
      } else {
        # Keep the order the site wrote: sorting scrambles the titration sequence.
        kk <- if (length(keys)) do.call(paste, c(lapply(keys, function(k)
                 as.character(d[[k]])), sep = "\r")) else rep("", nrow(d))
        o <- d[order(match(kk, unique(kk)),
                     match(d$site, c(sites, "POOLED"))), , drop = FALSE]
        L <- c(L, md_table(o, show), "")
      }
    }
  }
  writeLines(L, file.path(OUT, "pooled_summary.md"))
  written <<- c(written, "pooled_summary.md")
  cat(sprintf("  pooled_summary.md  (%d measure(s) described)\n", length(have)))
}


# ---- 5. Between-site heterogeneity --------------------------------------------
# The same overdispersion() the funnel uses, grouped by site instead of by unit.

cat("\nBetween-site heterogeneity\n")
het <- list()
ua <- stack_sites("06_unit_variation/unit_adherence.csv",
                  c("metric", "unit_key", "unit", "series"))
if (!is.null(ua) && nrow(ua)) {
  # One point per site: its own overall rate for the metric and series.
  agg <- aggregate(cbind(n_events, n_paired) ~ site + metric + series,
                   data = ua[ua$unit_key == "care_setting", ], FUN = sum)
  agg$pct <- 100 * agg$n_paired / agg$n_events
  for (m in unique(agg$metric)) for (s in unique(agg$series[agg$metric == m])) {
    g <- agg[agg$metric == m & agg$series == s, ]
    if (!nrow(g)) next
    p0 <- 100 * sum(g$n_paired) / sum(g$n_events)
    r <- overdispersion(g$pct, g$n_events, p0, trim = 0, min_k = MIN_K)
    het[[length(het) + 1L]] <- cbind(
      data.frame(grouping = "site", metric = m, series = s,
                 overall_pct = p0, stringsAsFactors = FALSE), r)
  }
}
if (length(het)) {
  het_tbl <- do.call(rbind, het)
  rownames(het_tbl) <- NULL
  utils::write.csv(het_tbl, file.path(OUT, "heterogeneity.csv"), row.names = FALSE)
  written <- c(written, "heterogeneity.csv")
  for (i in seq_len(nrow(het_tbl))) {
    r <- het_tbl[i, ]
    cat(sprintf("  %-12s %-16s k=%d  I2=%s  tau2=%s  %s\n", r$metric, r$series, r$k,
                if (is.na(r$ma_i2)) "NA" else sprintf("%.3f", r$ma_i2),
                if (is.na(r$ma_tau2)) "NA" else sprintf("%.2f", r$ma_tau2),
                if (isTRUE(r$reportable)) "reportable"
                else sprintf("NOT reportable (k < %d)", MIN_K)))
  }
} else {
  cat("  no unit_adherence rows, so nothing to compare between sites\n")
}


cat("\nSummary\n")
write_summary()


# ---- 6. Provenance -------------------------------------------------------------

prov <- list(
  pooled_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  n_sites = K,
  sites = sites,
  pooling_spec_version = SPEC$definition_version,
  protocol_digests = bundles[[1]]$prov$protocol_digests,
  definition_versions = bundles[[1]]$prov$definition_versions,
  code_digests = as.list(setNames(
    lapply(bundles, function(b) b$prov$code_digest), sites)),
  code_versions = as.list(setNames(
    lapply(bundles, function(b) b$prov$code_version), sites)),
  heterogeneity_min_k = MIN_K,
  files_written = written)
write(toJSON(prov, auto_unbox = TRUE, pretty = TRUE),
      file.path(OUT, "pooled_provenance.json"))

cat(sprintf("\n%d file(s) written to output/pooled/\n", length(written) + 1L))
cat("Review before sharing. Site names are REAL in these files; anonymise at ",
    "figure time.\n", sep = "")
