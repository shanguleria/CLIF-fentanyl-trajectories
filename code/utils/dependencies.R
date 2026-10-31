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
library(arrow)      # parquet handoff from the Python cohort build
library(ggplot2)    # figures
library(ggalluvial)  # state-transition alluvial (04_delivery_states.R)
library(gbmt)       # group-based trajectory models (code/tabled/)
library(lcmm)       # latent class mixed models (code/tabled/)
# RETIRED 2026-10-30: `survival` and `cmprsk` were declared here for "competing
# risks (Phase 6)" and NO SCRIPT IN THE REPO EVER CALLED EITHER. A declared
# dependency nothing consumes reads as policy -- the next person assumes the
# study estimates cumulative incidence with cmprsk, and it does not.
#
# The 28-day liberation figure (02_descriptive_cohort.R section 11b) is the
# place they would have gone, and it does not need them: discharge alive is
# carried forward to the horizon, exactly as ventilator_free_days() does it, so
# EVERY block has a known state at every point, there is no censoring, and the
# Aalen-Johansen cumulative incidence collapses to an exact proportion. Calling
# survfit() on uncensored data to borrow the name of a method would be
# decoration. If a future outcome does need a censored estimator, re-declare
# here and name the script.

# 07_outcome_model.R: clm() has no base-R substitute. metafor is deliberately
# absent -- the pooling is coordinating-centre code, so no site installs it.
library(ordinal)
