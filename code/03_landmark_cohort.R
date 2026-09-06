# ==============================================================================
# 03_landmark_cohort.R  --  Phase 2 -- apply landmark T, retention reporting
#
# Purpose : Restrict to patients alive and ventilated at landmark T; report retention and failed-extubation counts.
# Author  : Shan Guleria
# Created : 2026-09-05
# Inputs  : output/intermediate_phi/trajectory_long.parquet
# Outputs : output/intermediate_phi/landmark_cohort.parquet; retention table
#
# Spec: docs/design_notes.md section 10.
# ==============================================================================

# Run this in a FRESH R session (RStudio: Cmd+Shift+F10).
# Do not use rm(list = ls()) -- it does not unload packages or reset options,
# so it only gives the appearance of a clean slate.


# ---- 1. Packages -------------------------------------------------------------

pkgs <- c("here", "jsonlite", "ggplot2", "arrow")
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

message(sprintf("[03_landmark_cohort] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))


# ---- 4. TODO: Phase 2 -- apply landmark T, retention reporting ----
# read from  : dirs$out_phi
# PHI out    : dirs$out_phi
# aggregate  : dirs$out_final   (stamp `prov` onto anything shareable)

stop("03_landmark_cohort not yet implemented -- see docs/design_notes.md")


# ---- 5. Provenance -----------------------------------------------------------
# Which package versions produced these numbers?

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/03_landmark_cohort.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "03_landmark_cohort_sessioninfo.txt")
)
