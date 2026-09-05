# ==============================================================================
# 04_gbmt_classes.R  --  Phase 3 -- gbmt classes on combined dose
#
# Purpose : Group-based trajectory model on total fentanyl dose; sweep ng and select by the full criterion conjunction (design notes section 7).
# Author  : Shan Guleria
# Created : 2026-09-05
# Inputs  : output/intermediate_phi/landmark_cohort.parquet
# Outputs : output/intermediate_phi/gbmt_fits.rds; IC table, scree plot
#
# Spec: docs/design_notes.md section 10.
# ==============================================================================

# Run this in a FRESH R session (RStudio: Cmd+Shift+F10).
# Do not use rm(list = ls()) -- it does not unload packages or reset options,
# so it only gives the appearance of a clean slate.


# ---- 1. Packages -------------------------------------------------------------

pkgs <- c("here", "jsonlite", "gbmt", "arrow")
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

message(sprintf("[04_gbmt_classes] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))


# ---- 4. TODO: Phase 3 -- gbmt classes on combined dose ----
# read from  : dirs$data_phi
# PHI out    : dirs$out_phi
# aggregate  : dirs$out_final   (stamp `prov` onto anything shareable)

stop("04_gbmt_classes not yet implemented -- see docs/design_notes.md")


# ---- 5. Provenance -----------------------------------------------------------
# Which package versions produced these numbers?

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/04_gbmt_classes.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "04_gbmt_classes_sessioninfo.txt")
)
