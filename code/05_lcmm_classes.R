# ==============================================================================
# 05_lcmm_classes.R  --  Phase 4 -- lcmm classes, ARI vs Phase 3
#
# Purpose : Latent class mixed model on the same data; compare partitions against Phase 3 by adjusted Rand index.
# Author  : Shan Guleria
# Created : 2026-09-05
# Inputs  : output/intermediate_phi/landmark_cohort.parquet
# Outputs : output/intermediate_phi/lcmm_fits.rds; ARI comparison
#
# Spec: docs/design_notes.md section 10.
# ==============================================================================

# Run this in a FRESH R session (RStudio: Cmd+Shift+F10).
# Do not use rm(list = ls()) -- it does not unload packages or reset options,
# so it only gives the appearance of a clean slate.


# ---- 1. Packages -------------------------------------------------------------

pkgs <- c("here", "jsonlite", "lcmm", "arrow")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

source(here("code", "utils", "paths.R"))


# ---- 2. Config (never setwd(); here() anchors to the .Rproj) -----------------

config <- fromJSON(here("config", "config.json"), simplifyVector = FALSE)
set.seed(config$model$seed)


# ---- 3. Paths and provenance -------------------------------------------------
# One site, one output tree. site_dirs() creates them and labels the PHI ones.

dirs <- site_dirs()
prov <- provenance(config)

message(sprintf("[05_lcmm_classes] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))


# ---- 4. TODO: Phase 4 -- lcmm classes, ARI vs Phase 3 ----
# read from  : dirs$data_phi
# PHI out    : dirs$out_phi
# aggregate  : dirs$out_final   (stamp `prov` onto anything shareable)

stop("05_lcmm_classes not yet implemented -- see docs/design_notes.md")


# ---- 5. Provenance -----------------------------------------------------------
# Which package versions produced these numbers?

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/05_lcmm_classes.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "05_lcmm_classes_sessioninfo.txt")
)
