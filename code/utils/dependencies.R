# ==============================================================================
# dependencies.R  --  declare packages so renv can see them
#
# This file is NEVER sourced. It exists only to be read by renv's dependency
# scanner.
#
# Why it is needed: our scripts load packages through a loop,
#     for (p in pkgs) library(p, character.only = TRUE)
# which keeps `pkgs` as the single source of truth but makes the package names
# invisible to static analysis. renv scans source text for literal library()
# and pkg:: calls, so without this file `renv::snapshot()` records almost
# nothing and `renv::restore()` rebuilds a broken library on another machine.
#
# Keep this list in sync with the `pkgs` vectors in code/*.R.
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
