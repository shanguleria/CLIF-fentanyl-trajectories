# Purpose : Combine N site bundles into one set of pooled numbers. Run at the
#           COORDINATING CENTRE, never at a site.
# Author  : Shan Guleria
# Created : 2026-10-26
# Inputs  : one directory per site, each a copy of that site's
#           output/final_no_phi/ exactly as uploaded
# Outputs : output/pooled/ : one CSV per input file, plus heterogeneity.csv
#
#   Rscript code/pool_sites.R multisite/ucmc multisite/rush multisite/nu
#
# HOW TO READ THE OUTPUT. Every pooled file is LONG with a `site` column. Real
# site names appear in the per-site rows; the row where site == "POOLED" is the
# combined estimate. That is the same shape the CRRT work used, so a reader who
# knows that convention already knows this one.
#
# WHAT THIS SCRIPT WILL NOT DO, and why each matters:
#   - It will not average a percentage. Every proportion is recomputed from
#     summed numerators and denominators. The site exports carry both precisely
#     so that no percentage is ever reconstructed from a rounded one.
#   - It will not pool a row keyed by a unit name (section 8.7). care_setting
#     level sets differ by site by construction and location_name is site-local
#     free text, so those rows pass through per-site and feed the between-site
#     dispersion statistic instead.
#   - It will not pool missingness. That file is the comparability diagnostic --
#     what you read to decide whether pooling is honest at all.
#   - It will not pool across disagreeing protocols. That is a hard stop.
#
# The file-by-file method comes from config/pooling_spec.json, which this script
# ITERATES. A hardcoded list here would mean a newly shipped export is silently
# never pooled.

suppressMessages({
  library(here)
  library(jsonlite)
})
source(here("code", "utils", "paths.R"))
source(here("code", "utils", "pooling.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a   # used from step 2 onward

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  cat("usage: Rscript code/pool_sites.R <site_dir> [<site_dir> ...]\n\n",
      "Each <site_dir> is a copy of that site's output/final_no_phi/.\n",
      "Run code/check_provenance.py over the same directories first.\n", sep = "")
  quit(status = 1)
}
# `multisite/*` is the documented invocation and it globs the README in there too,
# so a non-directory argument is expected rather than exceptional. Skipped with a
# note, the way check_provenance.py already tolerates one. Erroring here meant
# the documented command failed on first use.
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
# CLEARED EVERY RUN, and it has to be. This script owns every file in here and
# writes each one from the bundles, so anything left behind is from a previous
# pool -- a different k, a different protocol, or a path that has since moved.
# The summary step re-READS these files by name (`file.exists` at step 6), so a
# stale twin is not inert: it is read back and rendered as though it were this
# run's number, with nothing to distinguish it.
#
# MEASURED 2026-10-27: relocating four method checks into
# 05_titration/diagnostics/ left 05_titration__window_sensitivity.csv and three
# siblings in here, holding a two-site v0.7.1 pool, beside the new
# 05_titration__diagnostics__*.csv from a one-site v0.8.0 pool. Clearing is the
# only fix that does not require a retired-path list to be kept accurate by hand.
# Safe at this point and nowhere later: nothing has been written yet.
stale <- list.files(OUT, full.names = TRUE)
if (length(stale)) {
  unlink(stale, recursive = TRUE)
  message(sprintf("  cleared %d file(s) from a previous pool", length(stale)))
}


# ---- 1. Read the bundles -----------------------------------------------------

read_bundle <- function(path) {
  # A site that zips its output/ directory sends the PARENT of final_no_phi, not
  # its contents -- which is a reasonable reading of "send output/final_no_phi/"
  # and is what Rush did. Descend one level rather than refusing: the alternative
  # is a confusing error on a bundle that is completely fine.
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
# Protocol divergence is a HARD STOP: bundles built to different cohort rules
# produce a number that means nothing. Code divergence is reported, not refused
# -- 01_build_cohort.py instructs a site to extend _to_mcg_hr for a locally
# charted dose unit, so some divergence is legitimate (SG, 2026-10-02).

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

# Every CSV present must be claimed by exactly one manifest entry, and every
# claimed file must be present. Both directions matter: an unclaimed file is
# silently never pooled, and a phantom entry means the spec describes an export
# that no longer exists.
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

# Long frame: one row per key per site, plus a POOLED row. Real site names in the
# per-site rows; anonymisation happens at figure time, not here, because the
# centre needs the real names to sanity-check a site's own numbers.
stack_sites <- function(rel, keys) {
  fr <- Filter(Negate(is.null), lapply(bundles, function(b) {
    d <- b$tabs[[rel]]
    if (is.null(d)) return(NULL)
    cbind(site = b$site, d, stringsAsFactors = FALSE)
  }))
  if (!length(fr)) return(NULL)
  # UNION the columns rather than requiring them to match. Nearly every file has
  # an identical header at every site, but baseline_characteristics.csv embeds
  # each site's own N in its column NAMES -- "Overall (N = 14,897)" against
  # "Overall (N = 10,627)" -- so a plain rbind dies with "names do not match
  # previous names". It is a per-site display table, so the honest result is each
  # site's rows under its own columns with NA elsewhere, which is lossless and
  # still reads side by side.
  #
  # Found at k = 2, on the first real second bundle. The idempotence test could
  # not have caught it: with one site there is nothing to disagree with.
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

# Sum the summables, carry the labels, recompute the declared percentages, and
# leave everything else explicitly NA. The three categories are DECLARED in
# config/pooling_spec.json rather than inferred: label_columns, not_poolable_columns
# and a per-file `percents` map. Every one of those lists exists because the
# idempotence test caught the engine getting it wrong.
pool_sum_counts <- function(d, keys, spec, rel) {
  # NOTE every spec key below is read with [[ ]], never $. R's $ on a list does
  # PREFIX MATCHING, so spec$denominator silently returned denominator_within on
  # the prevalence files -- which put "window_idx" into the summed columns and
  # produced a second window_idx column that cbind quietly renamed window_idx.1.
  # Found 2026-10-26 by the idempotence check.
  nums <- unlist(spec[["numerators"]]); den <- spec[["denominator"]]
  pcts <- spec[["percents"]]; within <- unlist(spec[["denominator_within"]])
  sum_cols <- intersect(c(nums, den), names(d))
  if (!length(sum_cols)) return(NULL)
  k <- key_of(d, keys)
  first <- !duplicated(k)

  # na.rm = TRUE would turn an all-NA group into 0. strobe's first row has no
  # n_before by construction, and 0 excluded-from-nothing is not the same claim.
  sum_keep_na <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
  agg <- lapply(sum_cols, function(cc)
    tapply(suppressWarnings(as.numeric(d[[cc]])), k, sum_keep_na))
  names(agg) <- sum_cols
  out <- as.data.frame(agg, stringsAsFactors = FALSE)

  kd <- d[first, intersect(c(keys, LABEL_COLS), names(d)), drop = FALSE]
  rownames(kd) <- k[first]
  out <- cbind(kd[rownames(out), , drop = FALSE], out)
  out$site <- "POOLED"

  # Percentages are RECOMPUTED from the summed counts and written back to their
  # own column name, so a reader reads pct_paired and not some pct_n_paired twin.
  for (pc in names(pcts)) {
    num <- pcts[[pc]]
    if (!(num %in% names(out))) next
    denom <- if (length(within) && !identical(within, "__first_window__")) {
      g <- key_of(out, intersect(within, names(out)))
      as.vector(tapply(out[[num]], g, sum, na.rm = TRUE)[g])
    } else if (identical(unname(within), "__first_window__")) {
      # retention's denominator is the anchor cohort, fixed at every window --
      # which is this file's own window 0 count.
      rep(out[[num]][which.min(suppressWarnings(as.numeric(out$window_idx)))], nrow(out))
    } else if (!is.null(den) && den %in% names(out)) {
      out[[den]]
    } else NULL
    recomputed <- if (is.null(denom)) NA_real_ else 100 * out[[num]] / denom
    # If EVERY site shipped NA for this percentage in a key group, the site code
    # blanked it on purpose and the pooled row must blank it too. dose_summary's
    # receivers_only rows are the case: the subset is already filtered to
    # dose > 0, so a 100 there carries no information, and section 8.5 removed it
    # on 2026-10-05. Recomputing reintroduced exactly that defect -- caught by
    # the idempotence check, which is why it runs over every column.
    if (pc %in% names(d)) {
      site_all_na <- tapply(is.na(suppressWarnings(as.numeric(d[[pc]]))),
                            key_of(d, keys), all)
      blank <- site_all_na[rownames(out)]
      recomputed[!is.na(blank) & blank] <- NA_real_
    }
    # Round to the precision the SITES shipped for this column, read off their
    # own values rather than fixed: coadministration ships 1 decimal and the
    # histograms ship 3, and a pooled row reading 81.222814864283% next to a
    # site row reading 81.2% is the kind of thing that reaches a manuscript.
    dec <- max(c(0L, vapply(as.character(d[[pc]]), function(v) {
      v <- sub("^.*\\.", "", v[grepl("\\.", v)])
      if (!length(v)) 0L else nchar(sub("0+$", "", v))
    }, integer(1))), na.rm = TRUE)
    out[[pc]] <- round(recomputed, dec)
  }
  # Anything numeric that is neither summed nor a label nor a recomputed percent
  # is a quantity a sum cannot produce. Say NA rather than carry a site's value.
  blanked <- setdiff(names(d), c(names(out), "site"))
  numeric_blanked <- blanked[vapply(blanked, function(nm)
    is.numeric(d[[nm]]), logical(1))]
  # Global list, unioned with this file's own. Poolability is a property of what
  # the file ships: mean pools where sum_sq ships and cannot where it does not.
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
  # Exact, and assumes nothing about equal variances. This is the whole reason
  # sum_sq ships rather than only mean and sd.
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
  # THE GROUPING COMES FROM THE DECLARED KEYS, not from a hardcoded column name.
  # Bins sum to a total WITHIN whatever the file is keyed on besides `bin`:
  # `variable` for the pooling histograms, `exposure` for 07's exposure
  # distribution. This hardcoded `variable` and crashed pool_sites.R outright on
  # the first bundle carrying 07 -- tapply(out$n, out$variable, ...) with a NULL
  # second argument -- after 36 of the files had already been written. The
  # bin-edge fingerprint below had the same bug more quietly: a NULL column drops
  # out of paste() silently, so the check still ran but compared less than it
  # claimed to.
  group_keys <- setdiff(keys, "bin")
  stopifnot("a sum_bins file must declare a `bin` key" = "bin" %in% keys)
  grp_of <- function(x) if (!length(group_keys)) rep("__all__", nrow(x))
                        else key_of(x, group_keys)

  # Summing bins is valid ONLY if they partition identically. They are
  # config-declared and protocol agreement is already asserted, so this can
  # only fail if a site edited the protocol -- but the exactness claim rests on
  # it, so it is checked rather than assumed.
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
    # Everything else falls through to NULL, which ships the per-site rows with
    # NO POOLED ROW. That is correct for side_by_side, per_site_only and
    # heterogeneity_source, which are per-site by design.
    #
    # coefficient_meta IS NOT IMPLEMENTED AND FALLS THROUGH HERE ON PURPOSE.
    # 07_outcome_model.R ships regression coefficients for inverse-variance
    # meta-analysis, which is a different operation from summing numerators and
    # denominators -- it needs the standard errors, a weighting scheme and a
    # decision about metafor that is deliberately deferred (protocol.json
    # model.outcome_model._pooling_why). Declaring the method in
    # pooling_spec.json means the file is CLAIMED, so a bundle carrying it is not
    # refused as unclaimed; landing here means the coefficients pass through
    # per-site and a reader sees every site's estimate with no combined row,
    # which is honest rather than wrong. Implement it before quoting a pooled
    # odds ratio.
    NULL)

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
# A pooled percentage with no statement of what it counts is the thing this
# document exists to prevent. Every measure gets its grain, its numerator, its
# denominator and its caveat, rendered from config/pooling_spec.json rather than
# retyped here -- so the prose cannot drift from the method it describes.

md_table <- function(d, cols) {
  hdr <- paste0("| ", paste(cols, collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", length(cols)), collapse = "|"), "|")
  # Counts get thousands separators: this document exists to be READ, and
  # 62998 against 62,998 is the difference between scanning and parsing.
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
  # Every percent-bearing entry must carry a description, or a number ships with
  # nothing saying what it counts.
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
        # Keep the ORDER THE SITE WROTE -- initiation, uptitration,
        # downtitration, any increase is a logical sequence that alphabetical
        # sorting scrambles into "any increase" first. Group sites under each
        # key, keys in first-appearance order.
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
# The same overdispersion() the funnel uses, grouped by SITE instead of by unit.
# Same algebra, different grouping -- which is the whole point: if units inside
# one hospital vary this much, the question is whether hospitals do too.

cat("\nBetween-site heterogeneity\n")
het <- list()
ua <- stack_sites("06_unit_variation/unit_adherence.csv",
                  c("metric", "unit_key", "unit", "series"))
if (!is.null(ua) && nrow(ua)) {
  # One point per SITE: that site's own overall rate for the metric and series.
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
