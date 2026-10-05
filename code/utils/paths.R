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
  # And a renumbered SCRIPT leaves the whole directory behind. Unlinking the
  # files emptied it but left the folder standing in the shareable tree, which
  # is a smaller version of the same problem: a reader finds an `03_states/`
  # that nothing writes and cannot tell whether it matters. Remove a retired
  # directory once it is empty -- never one that still holds anything, since
  # that would delete a file nobody declared.
  # NEVER a directory this run writes to. 02 retires old files that live in its
  # OWN phase directory, so after clearing them the folder is legitimately empty
  # and the rule above would delete the directory out from under the script --
  # which it did, taking all twelve of 02's outputs with it before the first
  # ggsave failed. A retired path may sit in a live directory; only the dead
  # ones may be removed.
  live <- normalizePath(unlist(dirs), mustWork = FALSE)
  for (d in unique(dirname(file.path(root, retired)))) {
    if (normalizePath(d, mustWork = FALSE) %in% live) next
    if (dir.exists(d) && !length(list.files(d, all.files = TRUE, no.. = TRUE))) {
      unlink(d, recursive = TRUE)
    }
  }
  n
}


# Mirrors paths.py. CODE_GLOBS, PROTOCOL_FILES and SITE_CONFIG must match it
# exactly or the two languages stamp different digests for the same tree;
# tests/test_paths.py asserts they agree.
CODE_GLOBS <- c("code/*.py", "code/*.R", "code/utils/*.py", "code/utils/*.R")
PROTOCOL_FILES <- c("covariates.json", "outlier_config.json")
SITE_CONFIG <- "config.json"

# sha256 of a STRING. tools::sha256sum only takes file paths, and `digest` is not
# in the renv library, so the string goes via a tempfile. Base R only, by design.
.sha_string <- function(s) {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  # writeChar without the trailing NUL, so the bytes match Python's .encode()
  con <- file(tmp, "wb"); writeChar(s, con, eos = NULL); close(con)
  substr(tools::sha256sum(tmp), 1, 16)
}

.sha_file <- function(p) substr(tools::sha256sum(p), 1, 16)

# Content fingerprint of the code that produced a run: list(overall, per_file).
# Content-based rather than git-based -- see the docstring in paths.py.
code_digests <- function(root = here::here()) {
  files <- unique(unlist(lapply(CODE_GLOBS,
                                function(g) Sys.glob(file.path(root, g)))))
  rel <- substring(files, nchar(root) + 2L)
  # method = "radix" is C-locale BYTE order, which is what Python's sorted()
  # gives. R's DEFAULT sort uses locale collation and ordered paths.py before
  # paths.R (collating case-insensitively, 'p' < 'r') where Python puts paths.R
  # first ('R' = 0x52 < 'p' = 0x70). Identical tree, different overall digest --
  # measured 2026-10-02, and invisible in the per-file digests, which matched.
  o <- order(rel, method = "radix")
  files <- files[o]
  rel <- rel[o]
  per <- setNames(vapply(files, .sha_file, character(1)), rel)
  # The overall digest hashes the (path, digest) PAIRS, so a rename moves it too.
  joined <- paste(sprintf("%s %s", rel, unname(per)), collapse = "\n")
  list(overall = .sha_string(joined), per_file = as.list(per))
}

# The protocol version each config declares. One field per file, not one
# overall: they are bumped independently.
definition_versions <- function(root = here::here()) {
  out <- list()
  for (n in PROTOCOL_FILES) {
    p <- file.path(root, "config", n)
    if (file.exists(p)) {
      v <- jsonlite::fromJSON(p)$definition_version
      out[[n]] <- if (is.null(v)) "" else v
    }
  }
  out
}

# SHA-256 of every file that governs what the outputs contain. Used by
# require_manifest() for WITHIN-site staleness; provenance() splits the same
# digests into protocol-vs-site for the CROSS-site check.
config_digests <- function(root = here::here()) {
  out <- list()
  for (n in c(SITE_CONFIG, PROTOCOL_FILES)) {
    p <- file.path(root, "config", n)
    # tools::sha256sum is base R; it matches Python's hashlib.sha256 of the same
    # bytes, verified 2026-09-07.
    if (file.exists(p)) out[[n]] <- .sha_file(p)
  }
  out
}

# The block stamped onto every shareable output. `code_digest` is the
# authoritative answer to "did two sites run the same code?"; `code_version` is
# the human-readable git label and may be "unknown".
provenance <- function(config) {
  root <- here::here()
  # `git describe` exits non-zero outside a checkout AND in a repo with no commits
  # yet. system(intern = TRUE) signals that as a WARNING, not an error, so
  # tryCatch(error=) alone does not catch it -- it returns character(0) and prints
  # a warning. Handle both, and treat an empty result as unknown.
  git_sha <- tryCatch(
    suppressWarnings(
      system("git describe --always --dirty", intern = TRUE, ignore.stderr = TRUE)),
    error = function(e) character(0)
  )
  cd <- code_digests(root)
  cfg <- config_digests(root)
  list(
    site_name       = config$site_name,
    clif_version    = config$clif_version,     # CLIF SPEC version
    dataset_version = if (is.null(config$dataset_version)) "" else config$dataset_version,
    definition_versions = definition_versions(root),
    code_version    = if (length(git_sha)) git_sha[1] else "unknown",
    code_digest     = cd$overall,
    protocol_digests = cfg[PROTOCOL_FILES],
    site_config_digest = if (is.null(cfg[[SITE_CONFIG]])) "" else cfg[[SITE_CONFIG]],
    file_digests    = cd$per_file,
    generated       = format(Sys.time(), tz = config$timezone, usetz = TRUE)
  )
}


# Fail loudly if the tables on disk were not produced by this code and config.
# Mirrors require_manifest() in paths.py; the two must agree.
require_manifest <- function(dirs, root) {
  f <- file.path(dirs$out_final, "manifest.json")
  if (!file.exists(f)) {
    stop("manifest.json is absent, so Phase 0 either never completed or ",
         "was cleared. Re-run code/01_build_cohort.py; do not read the parquet ",
         "files that may still be on disk.", call. = FALSE)
  }
  m <- jsonlite::fromJSON(f)
  was <- unlist(m$config_digests)
  # REQUIRED, not optional: a manifest with no digests cannot be checked against
  # the config on disk, and treating an absent block as empty made this guard
  # pass vacuously -- the same hole as a test that cannot fail.
  if (!length(was)) {
    stop("manifest.json carries no config_digests, so the tables on disk ",
         "cannot be checked against the current config. ",
         "Re-run code/01_build_cohort.py.", call. = FALSE)
  }
  now <- unlist(config_digests(root))
  drift <- names(was)[!is.na(now[names(was)]) & now[names(was)] != was]
  if (length(drift)) {
    stop("config has changed since Phase 0 ran: ", paste(drift, collapse = ", "),
         ". Re-run code/01_build_cohort.py rather than analysing stale tables.",
         call. = FALSE)
  }
  # Code drift is RECORDED, never refused (SG, 2026-10-02): a site may have to
  # edit _to_mcg_hr to run at all. The digest in each provenance.json is what
  # the coordinating centre compares.
  if (!is.null(m$code_digest)) {
    # unname() IS LOAD-BEARING. .sha_string() returns a NAMED character, so
    # identical() compared the names attribute too and was FALSE for equal
    # digests -- the NOTE fired on every run of every R script, printing the
    # same hash twice ("424e4460bc02333e -> 424e4460bc02333e"). A guard that
    # can never say "no drift" cannot report drift either. Measured and fixed
    # 2026-10-05. The Python twin at paths.py:243 compares plain strings and
    # was always correct, so the two mirrors had silently diverged.
    nowc <- unname(code_digests(root)$overall)
    if (!identical(as.character(m$code_digest), nowc)) {
      cat(sprintf("  NOTE code changed since Phase 0 ran (%s -> %s); these outputs mix two code versions\n",
                  m$code_digest, nowc))
    }
  }
  m
}
