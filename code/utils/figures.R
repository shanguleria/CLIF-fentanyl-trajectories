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
# stop agreeing -- these values lived in 03_delivery_states.R until 2026-09-24
# and had already drifted from the copy here on six of seven states.
# Levels come from code/utils/states.R.
#
# CATEGORICAL, not a lightness ramp: delivery routes are not ordered, so
# "continuous + bolus" is a different colour from "continuous only" rather than
# a darker one. The three terminal states are separate hues, so leaving the
# cohort never looks like a dose.
#
# Checked, not eyeballed. Every pair clears both separation floors under
# ALL-pairs adjacency, which is what the raster needs -- in it any two states
# can share an edge, unlike the area and alluvial where the stacking order
# fixes the neighbours:
#   worst simulated-CVD pair  extubated <-> continuous + bolus  dE 13.7 (floor 8)
#   worst normal-vision pair  died      <-> discharged alive    dE 15.6 (floor 15)
# `extubated` was #7fb069 and `discharged alive` #a8c8ee until 2026-09-24; those
# collided at dE 5.8 with `bolus only` under deuteranopia, and at dE 9.1 with
# `no fentanyl` for readers with NORMAL colour vision. Re-step with the
# validator, not by eye, if these ever move again.
STATE_COLOURS <- c(
  "no fentanyl"        = "#d8d6cf",
  "continuous only"    = "#14427e",
  "bolus only"         = "#eb6834",
  "continuous + bolus" = "#4a8bd8",
  "extubated"          = "#09b3a6",
  "discharged alive"   = "#86667b",
  "died"               = "#8c2f18")

# The intensity-band definition. The four bands ARE ordered, so they get a
# sequential ramp rather than categorical hues, and the three terminal states
# are inherited from above so the two state figures agree about what leaving
# the cohort looks like. Labels come from covariates.json via the caller --
# never restated here, or a consortium change to the band names would silently
# grey every band.
dose_palette <- function(labels) {
  stopifnot("dose_palette() expects the four band labels" = length(labels) == 4)
  setNames(c("#e8e6df", "#a8c8ee", "#4a8bd8", "#14427e",
             STATE_COLOURS[["extubated"]],
             STATE_COLOURS[["discharged alive"]],
             STATE_COLOURS[["died"]]),
           c(labels, "extubated", "discharged alive", "died"))
}

# ---- Numbers-at-risk table ---------------------------------------------------
# A survival-curve risk table, for a stacked area. A reader cannot take a
# percentage off a band to better than about five points, and cannot take the
# denominator off it at all -- which matters most here, where the same bands are
# drawn twice on denominators that differ by more than half the cohort by 72h.
#
# Rows are in STACKING order read top to bottom, so a row sits under the band it
# describes: ggplot's position_stack() puts factor level 1 at the TOP, so level
# order IS top-to-bottom order, and a discrete y axis has to be reversed to
# match. The denominator sits at the foot, where a risk table puts it.
#
# The panel shares the main plot's x scale and breaks so the columns land under
# the ticks. End columns are anchored rather than centred (hjust 0 at the first
# hour, 1 at the last); centring them would push half of each label outside the
# panel, and widening the x expansion to make room would leave the area chart
# floating off its own axis.
risk_table_panel <- function(prevalence, windows, denom_label, x_title = NULL,
                             text_size = 2.3) {
  stopifnot("risk_table_panel needs a factor state, in stacking order" =
              is.factor(prevalence$state),
            "windows must carry start_hr and label" =
              all(c("start_hr", "label") %in% names(windows)))
  lv  <- levels(prevalence$state)
  tab <- prevalence[prevalence$window_start_hr %in% windows$start_hr, ]
  den <- aggregate(n ~ window_start_hr, data = tab, FUN = sum)
  num <- function(x) formatC(x, format = "d", big.mark = ",")
  HEADER <- "hour"

  cells <- rbind(
    # The header names the hour each column reports. The figure no longer
    # carries a subtitle to explain the windowing, so the table says it itself.
    data.frame(window_start_hr = windows$start_hr, row = HEADER,
               label = windows$label),
    data.frame(window_start_hr = tab$window_start_hr,
               row   = as.character(tab$state),
               label = sprintf("%s (%.1f%%)", num(tab$n), tab$pct)),
    data.frame(window_start_hr = den$window_start_hr,
               row   = denom_label, label = num(den$n)))

  # Reversed, because a discrete y axis counts UP from the bottom while
  # position_stack() puts factor level 1 at the TOP of the bands.
  rows <- c(HEADER, lv, denom_label)
  cells$row   <- factor(cells$row, levels = rev(rows))
  cells$hjust <- ifelse(cells$window_start_hr == min(windows$start_hr), 0,
                 ifelse(cells$window_start_hr == max(windows$start_hr), 1, 0.5))
  cells$face  <- ifelse(cells$row %in% c(HEADER, denom_label), "bold", "plain")

  # Rules BETWEEN rows, not through them. house() sets panel.grid.major.y
  # explicitly, and `panel.grid = element_blank()` does NOT override an
  # explicitly-set child element -- the y gridlines survived and struck through
  # every line of text. Each grid component has to be blanked by name.
  p <- house(
    ggplot(cells, aes(window_start_hr, row)) +
      geom_hline(yintercept = seq_len(length(rows) - 1L) + 0.5,
                 colour = GRIDLINE, linewidth = 0.3) +
      geom_text(aes(label = label, hjust = hjust, fontface = face),
                size = text_size, colour = INK) +
      scale_x_continuous(breaks = windows$start_hr) +
      labs(x = x_title, y = NULL)) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          axis.text.x        = element_blank(),
          axis.ticks         = element_blank(),
          axis.text.y        = element_text(colour = INK, hjust = 0,
                                            size = rel(0.72)),
          plot.margin        = margin(t = 2, r = 5, b = 2, l = 5))
  attr(p, "n_rows") <- length(rows)
  p
}

# Stack a plot over its risk table, sharing one x axis. Panel widths are matched
# so the table columns line up under the bands -- the two panels have different
# y-label widths, so without this the table drifts sideways. The table's panel
# is given an absolute height; left as the default `null` unit it would claim
# half the figure.
stack_with_risk_table <- function(p_main, p_table, row_height = 0.20) {
  # ggplotGrob() measures text, which needs an open graphics device. Under
  # Rscript there is none, so R opens the DEFAULT one and leaves an Rplots.pdf
  # in the repo root on every pipeline run. Measure on a null device instead.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  g1 <- ggplotGrob(p_main)
  g2 <- ggplotGrob(p_table)
  w  <- grid::unit.pmax(g1$widths, g2$widths)
  g1$widths <- w
  g2$widths <- w
  panel <- g2$layout$t[g2$layout$name == "panel"]
  g2$heights[panel] <- grid::unit(row_height * attr(p_table, "n_rows"), "in")
  rbind(g1, g2, size = "first")
}
