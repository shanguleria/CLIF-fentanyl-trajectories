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
#
# LEGEND AT THE FOOT, in a box (SG, 2026-09-24). A legend on top is awkward and
# nonstandard; on the right it would steal width from an x axis that is always
# time, and on the exemplar it would collide with the three right-hand scales.
house <- function(p) {
  p + theme_minimal(base_size = 12) +
    theme(plot.title = element_text(colour = INK, face = "bold"),
          plot.subtitle = element_text(colour = MUTED, margin = margin(b = 10)),
          legend.position = "bottom",
          legend.background = element_rect(fill = PAPER, colour = INK,
                                           linewidth = 0.3),
          legend.margin = margin(t = 5, r = 8, b = 5, l = 8),
          legend.box.margin = margin(t = 4),
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

# ---- Stacking panels on one shared x axis -------------------------------------
# Several plots drawn one above the other, sharing an x axis: an area chart over
# its risk table, or the exemplar's four series. Two things have to be forced.
#
# WIDTHS. Each panel has its own y-axis labels, of different widths, so their
# plot panels would otherwise start at different x positions and the shared axis
# would be a lie. unit.pmax() across every gtable column aligns them.
#
# HEIGHTS. A ggplot panel is `1null`, meaning "share the space equally". Stacking
# four of those gives four equal panels, and stacking a risk table under a chart
# gives the table half the figure. Any panel passed an absolute height gets it;
# those left NA stay `null` and divide what is left.
#
# `heights` is in inches, one per plot, NA to leave a panel flexible.
stack_panels <- function(plots, heights = NULL) {
  stopifnot("stack_panels needs at least two plots" = length(plots) >= 2)
  if (is.null(heights)) heights <- rep(NA_real_, length(plots))
  stopifnot("one height per plot" = length(heights) == length(plots))

  # ggplotGrob() measures text, which needs an open graphics device. Under
  # Rscript there is none, so R opens the DEFAULT one and leaves an Rplots.pdf
  # in the repo root on every pipeline run. Measure on a null device instead.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  gs <- lapply(plots, ggplotGrob)

  # A bottom legend belongs under the WHOLE stack, not between the first plot
  # and whatever is stacked beneath it. On the prevalence figures that would
  # push the risk table away from the axis ticks its columns are aligned to,
  # which is the entire reason the panels share a width in the first place.
  legend <- NULL
  gs <- lapply(gs, function(g) {
    i <- which(g$layout$name == "guide-box-bottom")
    if (length(i) == 1L && !inherits(g$grobs[[i]], "zeroGrob")) {
      if (is.null(legend)) legend <<- g$grobs[[i]]
      g$grobs[[i]] <- grid::nullGrob()
      g$heights[g$layout$t[i]] <- grid::unit(0, "cm")
    }
    g
  })

  w  <- Reduce(grid::unit.pmax, lapply(gs, function(g) g$widths))
  gs <- lapply(seq_along(gs), function(i) {
    g <- gs[[i]]
    g$widths <- w
    if (!is.na(heights[i])) {
      panel <- g$layout$t[g$layout$name == "panel"]
      g$heights[panel] <- grid::unit(heights[i], "in")
    }
    g
  })
  out <- Reduce(function(a, b) rbind(a, b, size = "first"), gs)
  if (!is.null(legend)) {
    out <- gtable::gtable_add_rows(out, grid::grobHeight(legend))
    out <- gtable::gtable_add_grob(out, legend, t = nrow(out), b = nrow(out),
                                   l = 1, r = ncol(out),
                                   name = "guide-box-bottom")
  }
  out
}

# The risk-table case: the chart stays flexible, the table gets exactly the
# height its rows need.
stack_with_risk_table <- function(p_main, p_table, row_height = 0.20) {
  stack_panels(list(p_main, p_table),
               c(NA, row_height * attr(p_table, "n_rows")))
}

# ---- Captions ----------------------------------------------------------------
# Figures here are journal panels: no title, no subtitle. What would have sat on
# the panel is registered instead and written out as captions.md beside the
# figures, so the n, the denominator, the seed and the selection rules survive in
# a form a manuscript can lift. Held in an environment rather than a list so a
# draw function can register from inside its own scope without `<<-` reaching
# through whatever happens to enclose it.
.captions <- new.env(parent = emptyenv())

register_caption <- function(file, title, note) {
  assign(file, list(title = title, note = note), envir = .captions)
  invisible(NULL)
}

# `figures` is the authority on what must exist -- pass the .png entries of the
# script's OWNED list. A panel whose provenance lives only in the code is a panel
# nobody can caption later, so a missing entry is an error, not a warning.
write_captions <- function(path, script, figures, prov) {
  missing <- setdiff(figures, ls(.captions))
  stopifnot("every figure must register a caption" = length(missing) == 0)
  get_cap <- function(f) get(f, envir = .captions)
  writeLines(c(
    sprintf("# Figure captions -- %s", script),
    "",
    sprintf("Site %s | code %s | generated %s", prov$site_name, prov$code_version,
            prov$generated),
    "",
    "Figures are drawn journal-style, with no title or subtitle on the panel.",
    "These are the captions; edit for house style, but do not restate the numbers",
    "from memory -- they are written here by the run that drew the figures.",
    "",
    unlist(lapply(figures, function(f) c(
      sprintf("## `%s`", f), "",
      sprintf("**%s.** %s", get_cap(f)$title, get_cap(f)$note), "")))),
    path)
  invisible(length(figures))
}

# ---- A third y axis ----------------------------------------------------------
# ggplot gives a panel one secondary axis and no more, so a Baker-style panel
# with three scales -- one left, two right -- cannot be built from scales alone.
# This borrows a right-hand axis from a second plot drawn on the SAME panel range
# and hangs it outside the first plot's existing right axis.
#
# The donor's y scale must span the identical range, or the borrowed ticks will
# point at the wrong heights while looking perfectly plausible. `donor` is
# expected to be a bare plot carrying nothing but that scale.
add_third_axis <- function(p_main, donor, title = NULL, title_size = 11) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  g  <- ggplotGrob(p_main)
  gd <- ggplotGrob(donor)

  ax <- gd$grobs[[which(gd$layout$name == "axis-r")]]
  pos <- g$layout[g$layout$name == "axis-r", ]
  # ggplot lays out an axis-r SLOT whether or not a secondary axis exists, and
  # fills it with a zeroGrob when it does not -- so the slot's presence proves
  # nothing. Check the grob itself, or a third axis silently hangs off an empty
  # second one and the reader gets two scales where three are labelled.
  stopifnot(
    "the main plot has no right axis to hang a third scale beside" =
      nrow(pos) == 1,
    "the main plot's right axis is empty: give it a sec.axis first" =
      !inherits(g$grobs[[which(g$layout$name == "axis-r")]], "zeroGrob"),
    "the donor plot has no right axis to borrow" = !inherits(ax, "zeroGrob"))

  # Insert OUTSIDE the existing right axis AND its title, so each scale is
  # followed by its own label. Hanging it between the second axis and the second
  # axis's title puts one axis's name next to the other axis's ticks, which is
  # worse than having no title at all.
  ttl <- g$layout[g$layout$name %in% c("ylab-r", "axis-r-title", "ylab-r-title"), ]
  at <- if (nrow(ttl)) max(ttl$r) else max(pos$r)

  g <- gtable::gtable_add_cols(g, grid::unit(3, "mm"), pos = at)
  g <- gtable::gtable_add_cols(g, gd$widths[gd$layout$r[gd$layout$name == "axis-r"]],
                               pos = at + 1)
  g <- gtable::gtable_add_grob(g, ax, t = pos$t, b = pos$b, l = at + 2,
                               name = "axis-r-third")
  if (!is.null(title)) {
    lab <- grid::textGrob(title, rot = -90,
                          gp = grid::gpar(col = MUTED, fontsize = title_size))
    g <- gtable::gtable_add_cols(g, grid::unit(1.2, "lines"), pos = at + 2)
    g <- gtable::gtable_add_grob(g, lab, t = pos$t, b = pos$b, l = at + 3,
                                 name = "ylab-r-third")
  }
  g
}
