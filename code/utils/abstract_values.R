# Regenerate the VALUE table in docs/abstract_source_map.md from output/pooled/.
#
# The map used to carry hand-copied numbers beside their addresses, and it went
# stale within hours twice: the addresses stayed right while every value moved.
# A number that lives in a pooled CSV must not also live in prose, so the prose
# half of the document is kept by hand and the table between the two markers is
# written here, by the same run that produced the numbers.
#
# Called at the end of pool_sites.R. If docs/abstract_source_map.md is absent --
# a site that cloned the repo and never wrote one -- this does nothing rather
# than inventing a file nobody asked for.

AV_BEGIN <- "<!-- BEGIN GENERATED VALUES -- written by code/utils/abstract_values.R -->"
AV_END   <- "<!-- END GENERATED VALUES -->"

# Read a pooled file and pull one cell, naming the selector in the failure so a
# changed schema says which claim it broke rather than returning NA quietly.
av_cell <- function(pool, rel, where, col, label) {
  f <- file.path(pool, gsub("/", "__", rel))
  if (!file.exists(f)) return(list(value = NA_character_, note = "file absent"))
  d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  keep <- rep(TRUE, nrow(d))
  for (k in names(where)) {
    if (!(k %in% names(d))) return(list(value = NA_character_,
                                        note = sprintf("no column %s", k)))
    keep <- keep & as.character(d[[k]]) == where[[k]]
  }
  if (!any(keep) || !(col %in% names(d)))
    return(list(value = NA_character_, note = "no matching row"))
  list(value = as.character(d[[col]][which(keep)[1]]), note = NA_character_)
}

av_row <- function(label, value, rel, where, col, derived = FALSE,
                   how = NA_character_) {
  sel <- paste(sprintf("`%s=%s`", names(where), unlist(where)), collapse = ", ")
  sprintf("| %s | **%s** | %s | %s | %s |", label, value,
          if (derived) "_derived_" else sprintf("`%s`", rel),
          if (derived) how else if (nzchar(sel)) sel else "-",
          if (derived) "-" else sprintf("`%s`", col))
}

# Counts are read off this table by a human, so they wear separators. Anything
# that is already a percentage or a ratio is left exactly as the CSV wrote it --
# re-rounding a published number here would be a second source of truth.
av_count <- function(x) {
  n <- suppressWarnings(as.numeric(x))
  if (is.na(n)) return(as.character(x))
  formatC(n, format = "d", big.mark = ",")
}

write_abstract_values <- function(pool, doc = file.path(dirname(pool), "..", "docs",
                                                        "abstract_source_map.md"),
                                  sites = character(0)) {
  doc <- normalizePath(doc, mustWork = FALSE)
  if (!file.exists(doc)) {
    cat("  abstract_source_map.md absent -- skipped\n"); return(invisible(NULL))
  }
  P <- "POOLED"
  g <- function(rel, where, col) av_cell(pool, rel, where, col, "")$value
  num <- function(x) suppressWarnings(as.numeric(x))

  # The claims the abstract makes, each with the selector that resolves it.
  rows <- c()
  add <- function(...) rows <<- c(rows, av_row(...))

  add("Analytic cohort", av_count(g("07_outcome_model/strobe.csv",
      list(site = P, branch = "trunk", step = "ANALYTIC COHORT -- trajectory_long"), "n_after")),
      "07_outcome_model/strobe.csv",
      list(site = P, branch = "trunk", step = "ANALYTIC COHORT -- trajectory_long"), "n_after")
  add("Admissions with >= 1 escalation", av_count(g("07_outcome_model/strobe.csv",
      list(site = P, branch = "practice_patterns"), "n_after")),
      "07_outcome_model/strobe.csv", list(site = P, branch = "practice_patterns"), "n_after")
  add("Analysed in the model", av_count(g("07_outcome_model/strobe.csv",
      list(site = P, branch = "model", step = "ANALYSED IN THE MODEL"), "n_after")),
      "07_outcome_model/strobe.csv",
      list(site = P, branch = "model", step = "ANALYSED IN THE MODEL"), "n_after")

  ind <- list(site = P, cohort = "all", event_type = "increase", class = "either")
  add("Escalations", av_count(g("05_titration/indication.csv", ind, "n_events")),
      "05_titration/indication.csv", ind, "n_events")
  add("Score charted within 1h, %", g("05_titration/indication.csv", ind, "pct_assessed"),
      "05_titration/indication.csv", ind, "pct_assessed")
  add("Met a qualifying threshold, %", g("05_titration/indication.csv", ind, "pct_indicated_all"),
      "05_titration/indication.csv", ind, "pct_indicated_all")

  for (grp in c("any increase", "downtitration")) {
    w <- list(site = P, cohort = "all", group = grp)
    add(sprintf("Bolus-paired, %s, %%", grp), g("05_titration/coadministration.csv", w, "pct_paired"),
        "05_titration/coadministration.csv", w, "pct_paired")
  }

  add("VFD-28 = 0, %", g("07_outcome_model/vfd_ordinal.csv", list(site = P, level = "0"), "pct"),
      "07_outcome_model/vfd_ordinal.csv", list(site = P, level = "0"), "pct")

  cf <- file.path(pool, "07_outcome_model__model_coefficients.csv")
  or_pind <- NA_character_
  if (file.exists(cf)) {
    d <- utils::read.csv(cf, stringsAsFactors = FALSE)
    for (tm in c("Pindicated_per10pp", "Pbolus_per10pp")) {
      w <- list(site = P, model = "primary", term = tm)
      r <- d[d$site == P & d$model == "primary" & d$term == tm, , drop = FALSE]
      if (!nrow(r)) next
      if (tm == "Pindicated_per10pp") or_pind <- r$odds_ratio[1]
      add(sprintf("%s, pooled OR (95%% CI)", tm),
          sprintf("%s (%s-%s)", r$odds_ratio[1], r$or_ci_lo[1], r$or_ci_hi[1]),
          "07_outcome_model/model_coefficients.csv", w, "odds_ratio, or_ci_lo, or_ci_hi")
      add(sprintf("%s, p", tm), r$p_value[1],
          "07_outcome_model/model_coefficients.csv", w, "p_value")
      add(sprintf("%s, I-squared %%", tm), r$I2_pct[1],
          "07_outcome_model/model_coefficients.csv", w, "I2_pct")
    }
  }
  # The only derived figure the abstract uses, with its arithmetic stated.
  grid <- NA_real_
  pr <- file.path(pool, "07_outcome_model__adjusted_prediction.csv")
  if (file.exists(pr)) {
    d <- utils::read.csv(pr, stringsAsFactors = FALSE)
    grid <- max(num(d$exposure_pct), na.rm = TRUE)
  }
  if (!is.na(num(or_pind)) && is.finite(grid)) {
    k <- grid / 10
    add(sprintf("Fold-change across the 0-%d%% range", as.integer(grid)),
        sprintf("%.2f-fold", num(or_pind)^k),
        "", list(), "", derived = TRUE,
        how = sprintf("`%s ^ %g` -- the OR is per 10 pp and the grid spans %d",
                      or_pind, k, as.integer(grid)))
  }

  # Per-unit ranges are over PER-SITE rows: location_name is site-local, so there
  # is no pooled per-unit row to take a range from.
  ua <- file.path(pool, "06_unit_variation__unit_adherence.csv")
  if (file.exists(ua)) {
    d <- utils::read.csv(ua, stringsAsFactors = FALSE)
    d <- d[d$unit_key == "location_name" & d$site != P & is.finite(num(d$pct_paired)), ]
    for (m in c("indication", "coadmin")) {
      v <- num(d$pct_paired[d$metric == m &
                              d$series %in% c("any escalation", "any increase")])
      if (length(v))
        add(sprintf("%s across %d units, %%", m, length(v)),
            sprintf("%.1f to %.1f", min(v), max(v)),
            "06_unit_variation/unit_adherence.csv",
            list(unit_key = "location_name", metric = m, site = "per-site rows"),
            "min/max of pct_paired")
    }
  }

  body <- c(
    AV_BEGIN,
    "",
    sprintf("_Regenerated %s by `code/pool_sites.R` over %d site(s): %s._",
            format(Sys.time(), "%Y-%m-%d %H:%M %Z"), length(sites),
            paste(sites, collapse = ", ")),
    "",
    "**Do not hand-edit this block.** Every value is read from `output/pooled/`",
    "by the run that produced it; the prose outside the markers is maintained by",
    "hand. A number that lives in a pooled CSV must not also live in prose.",
    "",
    "| Claim | Value | Pooled file | Row selector | Column |",
    "|---|---|---|---|---|",
    rows,
    "",
    AV_END)

  txt <- readLines(doc, warn = FALSE)
  i <- which(txt == AV_BEGIN); j <- which(txt == AV_END)
  if (length(i) == 1 && length(j) == 1 && j > i) {
    txt <- c(txt[seq_len(i - 1)], body, txt[(j + 1):length(txt)])
  } else {
    txt <- c(txt, "", "---", "", "## Current pooled values", "", body)
  }
  writeLines(txt, doc)
  cat(sprintf("  docs/abstract_source_map.md  (%d value(s) refreshed)\n", length(rows)))
}
