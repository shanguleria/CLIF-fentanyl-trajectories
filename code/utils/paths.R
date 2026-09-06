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

# Create and return the four directories this pipeline may write to.
site_dirs <- function() {
  root <- here::here()
  d <- list(
    out_phi   = file.path(root, "output", "intermediate_phi"),
    out_final = file.path(root, "output", "final_no_phi"),
    logs      = file.path(root, "logs")
  )
  for (p in d) dir.create(p, recursive = TRUE, showWarnings = FALSE)

  for (phi in c(d$out_phi)) {
    label <- file.path(phi, "README.md")
    if (!file.exists(label)) writeLines(PHI_LABEL, label)
  }

  d
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
