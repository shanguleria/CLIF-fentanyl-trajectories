# Checks on require_manifest()'s config-drift guard -- the R side, which the
# Python test suite cannot reach.
#
# WHY THIS EXISTS. The guard was INERT from the day it was written until
# 2026-10-28. tools::sha256sum() returns a NAMED character and substr() keeps
# the name, so every digest carried the file path as a name; unlist()ing them
# folded that into the keys ("protocol.json" -> "protocol.json./Users/.../
# protocol.json"), every lookup by the clean name returned NA, and the
# !is.na(...) clause read NA as "no drift". A changed protocol.json never
# stopped anything. Found by changing protocol.json and watching 07 run.
#
# It is the SECOND time this exact attribute has cost this file a guard -- the
# code_digest NOTE had the same cause on 2026-10-05 and was fixed at the call
# site rather than at the source. Hence a test rather than a third comment.
#
# Run standalone:  Rscript tests/test_require_manifest.R
suppressPackageStartupMessages({library(here); library(jsonlite)})
source(here("code", "utils", "paths.R"))
root <- here()

pass <- 0L; fail <- 0L
check <- function(name, ok, detail = "") {
  if (isTRUE(ok)) { pass <<- pass + 1L; cat(sprintf("  PASS  %s\n", name)) }
  else { fail <<- fail + 1L; cat(sprintf("  FAIL  %s\n        %s\n", name, detail)) }
}
attempt <- function(dir) {
  tryCatch({require_manifest(list(out_final = dir), root); list(ok = TRUE, msg = "")},
           error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
}
truth <- vapply(config_digests(root), function(x) unname(as.character(x))[1],
                character(1))

# A digest must be a bare 16-hex string with NO name attribute. This is the
# root-cause check: everything below is downstream of it.
check("digests carry no stray name attribute",
      is.null(names(.sha_file(file.path(root, "config", "protocol.json")))),
      "tools::sha256sum's name survived -- unlist() will corrupt every key")
check("config_digests keys are the bare file names",
      identical(sort(names(unlist(config_digests(root)))),
                sort(names(truth))),
      paste("got:", paste(names(unlist(config_digests(root))), collapse = ", ")))

base <- tempfile("final_"); dir.create(base, recursive = TRUE)
write_manifest <- function(cd) {
  m <- list(code_version = "test", code_digest = unname(code_digests(root)$overall),
            generated = "test", config_digests = as.list(cd))
  write_json(m, file.path(base, "manifest.json"), auto_unbox = TRUE)
}

write_manifest(truth)
r <- attempt(base)
check("matching digests pass", r$ok, r$msg)

drifted <- truth; drifted[["protocol.json"]] <- "0000000000000000"
write_manifest(drifted)
r <- attempt(base)
check("a changed protocol.json STOPS the run",
      !r$ok && grepl("config has changed", r$msg, fixed = TRUE), r$msg)
check("the stop names the file that drifted",
      !r$ok && grepl("protocol.json", r$msg, fixed = TRUE), r$msg)

drifted2 <- truth; drifted2[["config.json"]] <- "0000000000000000"
write_manifest(drifted2)
r <- attempt(base)
check("a changed site config.json STOPS the run",
      !r$ok && grepl("config has changed", r$msg, fixed = TRUE), r$msg)

gone <- c(truth, a_retired_config.json = "deadbeefdeadbeef")
write_manifest(gone)
r <- attempt(base)
check("a manifest key with no file on disk STOPS the run",
      !r$ok && grepl("no longer present", r$msg, fixed = TRUE), r$msg)

write_manifest(list())
r <- attempt(base)
check("a manifest with no config_digests STOPS the run",
      !r$ok && grepl("no config_digests", r$msg, fixed = TRUE), r$msg)

unlink(file.path(base, "manifest.json"))
r <- attempt(base)
check("an absent manifest STOPS the run",
      !r$ok && grepl("manifest.json is absent", r$msg, fixed = TRUE), r$msg)

cat(sprintf("\n%d/%d passed\n", pass, pass + fail))
if (fail) quit(status = 1)
