# ==============================================================================
# dependencies.R  --  declare packages so renv can see them
#
# NEVER sourced. The scripts load packages through a loop,
#     for (p in pkgs) library(p, character.only = TRUE)
# which hides the names from renv's static scanner; without the literal
# library() calls below, renv::snapshot() records almost nothing. Keep in sync
# with the `pkgs` vectors in code/*.R.
# ==============================================================================

library(here)       # project-relative paths
library(jsonlite)   # config.json
library(arrow)      # parquet handoff from the Python phase
library(ggplot2)    # figures
library(ggalluvial)  # state-transition alluvial (Phase 1)
library(gbmt)       # group-based trajectory models (Phases 3, 5)
library(lcmm)       # latent class mixed models (Phase 4)
library(survival)   # competing risks (Phase 6)
library(cmprsk)     # cumulative incidence, Gray's test (Phase 6)

# 07_outcome_model.R: clm() has no base-R substitute. metafor is deliberately
# absent -- the pooling is coordinating-centre code, so no site installs it.
library(ordinal)
