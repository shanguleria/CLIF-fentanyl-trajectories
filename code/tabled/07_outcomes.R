# ==============================================================================
# 07_outcomes.R  --  Phase 6 -- competing-risks outcome models
#
# Purpose : Trajectory class as predictor of successful extubation (competing with death and tracheostomy) and of in-hospital death (competing with discharge alive), followed from the landmark.
# Author  : Shan Guleria
# Created : 2026-09-05
# Inputs  : class assignments + output/intermediate_phi/time_to_event.parquet
# Outputs : output/final_no_phi/ : CIF curves, cause-specific and Fine-Gray models
#
# ==============================================================================

# Run this in a FRESH R session (RStudio: Cmd+Shift+F10).
# Do not use rm(list = ls()) -- it does not unload packages or reset options,
# so it only gives the appearance of a clean slate.


# ---- 1. Packages -------------------------------------------------------------

pkgs <- c("here", "jsonlite", "survival", "cmprsk", "arrow")
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

message(sprintf("[07_outcomes] site=%s  clif=%s  data=%s",
                config$site_name, config$clif_version, config$data_directory))


# ---- 4. TODO: Phase 6 -- competing-risks outcome models ----
# read from  : dirs$out_phi
# PHI out    : dirs$out_phi
# aggregate  : dirs$out_final   (stamp `prov` onto anything shareable)

stop("07_outcomes not yet implemented")


# ---- 5. Provenance -----------------------------------------------------------
# Which package versions produced these numbers?

writeLines(
  c(paste("Run at:", format(Sys.time(), tz = config$timezone, usetz = TRUE)),
    paste("Script :", "code/07_outcomes.R"),
    "",
    capture.output(sessionInfo())),
  here("logs", "07_outcomes_sessioninfo.txt")
)
