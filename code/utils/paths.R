# Where this run reads and writes, and the provenance stamped on what it shares.
#
# One site, one data source, one output tree -- the standard CLIF layout. The R
# half is a transcription of code/utils/paths.py; the two must agree, because
# Phase 0 (Python) writes what Phases 1-6 (R) read.

PHI_LABEL <- paste(
  "# intermediate_phi",
  "",
  "Patient-level intermediates: one row per patient, or per patient-window.",
  "",
  "**This directory never leaves the site.** It is not part of any export bundle",
  "and must not be committed, copied to a shared drive, or read into an analysis",
  "transcript. Only `output/final_no_phi/` is shareable.",
  sep = "\n")

# Create and return the directories this pipeline may write to. Per-script
# subfolders of out_final come from phase_dir(), not from here.
site_dirs <- function() {
  root <- here::here()
  d <- list(
    out_phi     = file.path(root, "output", "intermediate_phi"),
    out_final   = file.path(root, "output", "final_no_phi"),
    logs        = file.path(root, "logs")
  )
  for (p in d) dir.create(p, recursive = TRUE, showWarnings = FALSE)

  for (phi in c(d$out_phi)) {
    label <- file.path(phi, "README.md")
    if (!file.exists(label)) writeLines(PHI_LABEL, label)
  }

  d
}

# A per-script subfolder of the shareable output tree. Mirrors phase_dir() in
# paths.py; the two must agree. The folder says which script produced the file,
# the phaseN_ prefix says which phase it belongs to.
phase_dir <- function(dirs, name) {
  d <- file.path(dirs$out_final, name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}


# Delete this script's own outputs before it starts, so a crash leaves nothing
# rather than a stale file that still looks valid. Mirrors clear_owned_outputs()
# in paths.py; the two must agree.
clear_owned_outputs <- function(dirs, owned, retired = character(0)) {
  n <- 0L
  for (key in names(owned)) {
    for (name in owned[[key]]) {
      f <- file.path(dirs[[key]], name)
      if (file.exists(f)) { unlink(f); n <- n + 1L }
    }
  }
  # A relocated output leaves a stale twin the owned list no longer names.
  root <- here::here()
  for (rel in retired) {
    f <- file.path(root, rel)
    if (file.exists(f)) { unlink(f); n <- n + 1L }
  }
  n
}


# The block stamped onto every shareable output.
provenance <- function(config) {
  # `git describe` exits non-zero outside a checkout AND in a repo with no commits
  # yet. system(intern = TRUE) signals that as a WARNING, not an error, so
  # tryCatch(error=) alone does not catch it -- it returns character(0) and prints
  # a warning. Handle both, and treat an empty result as unknown.
  git_sha <- tryCatch(
    suppressWarnings(
      system("git describe --always --dirty", intern = TRUE, ignore.stderr = TRUE)),
    error = function(e) character(0)
  )
  list(
    site_name       = config$site_name,
    clif_version    = config$clif_version,     # CLIF SPEC version
    dataset_version = if (is.null(config$dataset_version)) "" else config$dataset_version,
    code_version    = if (length(git_sha)) git_sha[1] else "unknown",
    generated       = format(Sys.time(), tz = config$timezone, usetz = TRUE)
  )
}


# Fail loudly if the tables on disk were not produced by this code and config.
# Mirrors require_manifest() in paths.py; the two must agree.
require_manifest <- function(dirs, root) {
  f <- file.path(dirs$out_final, "phase0_manifest.json")
  if (!file.exists(f)) {
    stop("phase0_manifest.json is absent, so Phase 0 either never completed or ",
         "was cleared. Re-run code/01_build_cohort.py; do not read the parquet ",
         "files that may still be on disk.", call. = FALSE)
  }
  m <- jsonlite::fromJSON(f)
  now <- vapply(c("config.json", "covariates.json", "outlier_config.json"),
                function(n) {
                  p <- file.path(root, "config", n)
                  # tools::sha256sum is base R; it matches Python's
                  # hashlib.sha256 of the same bytes, verified 2026-09-07.
                  if (file.exists(p)) substr(tools::sha256sum(p), 1, 16)
                  else NA_character_
                }, character(1))
  was <- unlist(m$config_digests)
  drift <- names(was)[!is.na(now[names(was)]) & now[names(was)] != was]
  if (length(drift)) {
    stop("config has changed since Phase 0 ran: ", paste(drift, collapse = ", "),
         ". Re-run code/01_build_cohort.py rather than analysing stale tables.",
         call. = FALSE)
  }
  m
}
