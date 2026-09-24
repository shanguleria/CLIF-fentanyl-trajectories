# Shared figure look, in one place.
#
# The palette line and the theme were copy-pasted into four scripts, which is
# tolerable while one script draws figures and a problem once two do. Sourcing
# this is what keeps a figure drawn in 03_delivery_states.R looking like one
# drawn in 02_descriptive_cohort.R.
#
# Needs: ggplot2 attached by the caller.

INK      <- "#0b0b0b"
MUTED    <- "#898781"
GRIDLINE <- "#e1e0d9"
PAPER    <- "#fcfcfb"

# Kept as lowercase aliases too: the scripts were written against these names and
# renaming every reference buys nothing.
ink <- INK; muted <- MUTED; gridline <- GRIDLINE

# The house theme. Horizontal gridlines only -- every figure here has a
# continuous y and a time or categorical x, so vertical rules add ink without
# adding information.
house <- function(p) {
  p + theme_minimal(base_size = 12) +
    theme(plot.title = element_text(colour = INK, face = "bold"),
          plot.subtitle = element_text(colour = MUTED, margin = margin(b = 10)),
          legend.position = "top",
          legend.text = element_text(colour = MUTED, size = 9),
          axis.title = element_text(colour = MUTED),
          axis.text = element_text(colour = MUTED),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = GRIDLINE, linewidth = 0.4),
          plot.background = element_rect(fill = PAPER, colour = NA),
          strip.text = element_text(colour = INK, face = "bold"))
}

# ---- State palettes ----------------------------------------------------------
# ONE palette per state definition, shared by every figure that draws it. The
# prevalence plot, the alluvial and the per-patient raster must agree on colour
# or a reader cannot move between them, and a per-figure palette is how they
# stop agreeing. Levels come from code/utils/states.R.
#
# Family in hue, level in lightness: the four fentanyl states are one blue ramp
# so "more fentanyl" reads as "darker", and the three terminal states are
# separate hues so leaving the cohort never looks like a dose.
STATE_COLOURS <- c(
  "no fentanyl"        = "#d9e2ec",
  "bolus only"         = "#9fb3c8",
  "continuous only"    = "#4a8bd8",
  "continuous + bolus" = "#14427e",
  "extubated"          = "#7fb069",
  "discharged alive"   = "#b8b2a7",
  "died"               = "#5c5552")

# The intensity-band definition shares the three terminal colours, so the two
# state figures stay comparable where they describe the same thing.
DOSE_STATE_COLOURS <- c(
  "zero"               = "#d9e2ec",
  "low"                = "#9fb3c8",
  "medium"             = "#4a8bd8",
  "high"               = "#14427e",
  "extubated"          = STATE_COLOURS[["extubated"]],
  "discharged alive"   = STATE_COLOURS[["discharged alive"]],
  "died"               = STATE_COLOURS[["died"]])

# Fail loudly rather than silently dropping a level to grey: a state with no
# colour is a state the reader cannot see.
state_palette <- function(levels) {
  pal <- if (all(levels %in% names(STATE_COLOURS))) STATE_COLOURS
         else if (all(levels %in% names(DOSE_STATE_COLOURS))) DOSE_STATE_COLOURS
         else NULL
  if (is.null(pal)) {
    stop("no palette covers these state levels: ",
         paste(setdiff(levels, union(names(STATE_COLOURS),
                                     names(DOSE_STATE_COLOURS))), collapse = ", "),
         ". Add them to code/utils/figures.R rather than letting ggplot grey them.",
         call. = FALSE)
  }
  pal[levels]
}
