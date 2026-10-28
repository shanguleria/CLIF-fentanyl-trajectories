# Where this run reads and writes, and the provenance stamped on what it shares.
# Transcription of code/utils/paths.py; the two must agree.

PHI_LABEL <- paste(
  "# intermediate_phi",
  "",
  "Patient-level intermediates: one row per patient, or per patient-window.",
  "",
  "**This directory never leaves the site.** It is not part of any export bundle",
  "and must not be committed, copied to a shared drive, or read into an analysis",
  "transcript. Only `output/final_no_phi/` is shareable.",
  sep = "\n")

# Directories this pipeline may write to. Per-script subfolders: phase_dir().
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

# A per-script subfolder of the shareable output tree. Mirrors paths.py.
phase_dir <- function(dirs, name) {
  d <- file.path(dirs$out_final, name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}


# Delete this script's own outputs first, so a crash leaves nothing stale.
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
  # A renumbered script leaves its directory behind. Remove a retired directory
  # only once empty, and never one this run writes to (02 retires into its own).
  live <- normalizePath(unlist(dirs), mustWork = FALSE)
  for (d in unique(dirname(file.path(root, retired)))) {
    if (normalizePath(d, mustWork = FALSE) %in% live) next
    if (dir.exists(d) && !length(list.files(d, all.files = TRUE, no.. = TRUE))) {
      unlink(d, recursive = TRUE)
    }
  }
  n
}


# Must match paths.py exactly or the two languages stamp different digests for
# the same tree; tests/test_paths.py asserts they agree.
CODE_GLOBS <- c("code/*.py", "code/*.R", "code/utils/*.py", "code/utils/*.R")
PROTOCOL_FILES <- c("protocol.json", "outlier_config.json")
SITE_CONFIG <- "config.json"

# sha256 of a STRING, via a tempfile: tools::sha256sum takes paths only, and
# `digest` is not in the renv library. Base R by design.
.sha_string <- function(s) {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  # writeChar without the trailing NUL, so the bytes match Python's .encode()
  con <- file(tmp, "wb"); writeChar(s, con, eos = NULL); close(con)
  # unname() as in .sha_file(); the stray name here is the tempfile path.
  unname(substr(tools::sha256sum(tmp), 1, 16))
}

# unname() is LOAD-BEARING: sha256sum() returns a NAMED character and unlist()
# folds that name into the key, which once disabled the config-drift guard.
.sha_file <- function(p) unname(substr(tools::sha256sum(p), 1, 16))

# Content fingerprint of the code that produced a run: list(overall, per_file).
code_digests <- function(root = here::here()) {
  files <- unique(unlist(lapply(CODE_GLOBS,
                                function(g) Sys.glob(file.path(root, g)))))
  rel <- substring(files, nchar(root) + 2L)
  # radix = C-locale BYTE order, matching Python's sorted(); R's default sort
  # collates by locale and would order the same tree differently.
  o <- order(rel, method = "radix")
  files <- files[o]
  rel <- rel[o]
  per <- setNames(vapply(files, .sha_file, character(1)), rel)
  # The overall digest hashes the (path, digest) PAIRS, so a rename moves it too.
  joined <- paste(sprintf("%s %s", rel, unname(per)), collapse = "\n")
  list(overall = .sha_string(joined), per_file = as.list(per))
}

# The protocol version each config declares, one field per file.
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

# SHA-256 of every file that governs what the outputs contain: require_manifest()
# checks it WITHIN site, provenance() splits it protocol-vs-site for ACROSS.
config_digests <- function(root = here::here()) {
  out <- list()
  for (n in c(SITE_CONFIG, PROTOCOL_FILES)) {
    p <- file.path(root, "config", n)
    # base R, and matches Python's hashlib.sha256 of the same bytes.
    if (file.exists(p)) out[[n]] <- .sha_file(p)
  }
  out
}

# The block stamped onto every shareable output. `code_digest` answers "did two
# sites run the same code?"; `code_version` is the git label, maybe "unknown".
provenance <- function(config) {
  root <- here::here()
  # `git describe` fails outside a checkout and in a repo with no commits, which
  # system(intern = TRUE) reports as a warning returning character(0), not an
  # error. Treat empty as unknown.
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
require_manifest <- function(dirs, root) {
  f <- file.path(dirs$out_final, "manifest.json")
  if (!file.exists(f)) {
    stop("manifest.json is absent, so Phase 0 either never completed or ",
         "was cleared. Re-run code/01_build_cohort.py; do not read the parquet ",
         "files that may still be on disk.", call. = FALSE)
  }
  m <- jsonlite::fromJSON(f)
  was <- unlist(m$config_digests)
  # REQUIRED, not optional: treating an absent digest block as empty made this
  # guard pass vacuously.
  if (!length(was)) {
    stop("manifest.json carries no config_digests, so the tables on disk ",
         "cannot be checked against the current config. ",
         "Re-run code/01_build_cohort.py.", call. = FALSE)
  }
  # vapply rather than unlist(), so one malformed element cannot reshape the
  # names of all of them and turn every lookup into an unnoticed NA.
  nowl <- config_digests(root)
  now <- vapply(nowl, function(x) unname(as.character(x))[1], character(1))
  # A manifest key with no current counterpart is NOT agreement: the config that
  # governed Phase 0 is gone under that name.
  missing_now <- setdiff(names(was), names(now))
  if (length(missing_now)) {
    stop("manifest.json records config file(s) that are no longer present: ",
         paste(missing_now, collapse = ", "),
         ". The tables on disk were built under a config this repo no longer ",
         "has, so they cannot be checked. Re-run code/01_build_cohort.py.",
         call. = FALSE)
  }
  drift <- names(was)[now[names(was)] != was]
  if (length(drift)) {
    stop("config has changed since Phase 0 ran: ", paste(drift, collapse = ", "),
         ". Re-run code/01_build_cohort.py rather than analysing stale tables.",
         call. = FALSE)
  }
  # Code drift is RECORDED, never refused: a site may have to edit _to_mcg_hr to
  # run at all. The centre compares the digest in each provenance.json.
  if (!is.null(m$code_digest)) {
    # unname() is LOAD-BEARING: identical() would otherwise compare the names
    # attribute too and report drift on equal digests.
    nowc <- unname(code_digests(root)$overall)
    if (!identical(as.character(m$code_digest), nowc)) {
      cat(sprintf("  NOTE code changed since Phase 0 ran (%s -> %s); these outputs mix two code versions\n",
                  m$code_digest, nowc))
    }
  }
  m
}
