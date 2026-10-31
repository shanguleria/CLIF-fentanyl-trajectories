# Shared figure look, in one place, so every script's figures agree.
# Needs: ggplot2 attached by the caller.

INK      <- "#0b0b0b"
# MUTED is for MARKS that must recede (a reference rule, an axis spine); it sets
# no text -- every label reads INK, because these figures are for print.
MUTED    <- "#898781"
GRIDLINE <- "#e1e0d9"
PAPER    <- "#fcfcfb"

# Lowercase aliases too: the scripts were written against these names.
ink <- INK; muted <- MUTED; gridline <- GRIDLINE

# The house theme. Horizontal gridlines only; legend at the foot in a box, since
# a right-hand legend would steal width from an x axis that is always time.
house <- function(p) {
  p + theme_minimal(base_size = 12) +
    theme(plot.title = element_text(colour = INK, face = "bold"),
          plot.subtitle = element_text(colour = INK, margin = margin(b = 10)),
          legend.position = "bottom",
          legend.background = element_rect(fill = PAPER, colour = INK,
                                           linewidth = 0.3),
          legend.margin = margin(t = 5, r = 8, b = 5, l = 8),
          legend.box.margin = margin(t = 4),
          legend.text = element_text(colour = INK, size = 9),
          axis.title = element_text(colour = INK),
          axis.text = element_text(colour = INK),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = GRIDLINE, linewidth = 0.4),
          plot.background = element_rect(fill = PAPER, colour = NA),
          strip.text = element_text(colour = INK, face = "bold"))
}

# ---- State palettes ----------------------------------------------------------
# ONE palette per state definition, shared by every figure that draws it, or the
# prevalence plot, the alluvial and the raster stop agreeing. Levels come from
# code/utils/states.R. Categorical, not a lightness ramp: routes are not
# ordered. Separation measured with the validator under ALL-pairs adjacency,
# which is what the raster needs; re-step there, not by eye, if these move.
STATE_COLOURS <- c(
  "no fentanyl"        = "#d8d6cf",
  "continuous only"    = "#14427e",
  "bolus only"         = "#eb6834",
  "continuous + bolus" = "#4a8bd8",
  "extubated"          = "#09b3a6",
  "discharged alive"   = "#86667b",
  "died"               = "#8c2f18")

# The three titration series in 06_unit_variation.R. Categorical -- "any
# increase" is the union of the other two. Deliberately NOT the state hues:
# those fail the per-colour lightness and contrast checks against the #fcfcfb
# ground, which does not show in a filled area but does in 2mm dots and thin
# bars. Stepped to the nearest passing triple with the validator, not by eye.
SERIES_COLOURS <- c(
  "initiation"   = "#234e90",
  "uptitration"  = "#ed6531",
  "any increase" = "#39a390")

# The intensity-band definition. The four bands ARE ordered, so a sequential
# ramp; the terminal three are inherited so the two state figures agree. Labels
# come from protocol.json via the caller -- restating them here would silently
# grey a renamed band.
dose_palette <- function(labels) {
  stopifnot("dose_palette() expects the four band labels" = length(labels) == 4)
  setNames(c("#e8e6df", "#a8c8ee", "#4a8bd8", "#14427e",
             STATE_COLOURS[["extubated"]],
             STATE_COLOURS[["discharged alive"]],
             STATE_COLOURS[["died"]]),
           c(labels, "extubated", "discharged alive", "died"))
}

# The 28-day liberation / survival figure's four states. Categorical, not a ramp:
# the states are not ordered. These are the STATE_COLOURS hues stepped to the
# nearest passing quadruple under the validator at --pairs all (a stacked band
# can sit next to any other once a reader compares across the legend), on the
# #fcfcfb ground. The originals FAIL: #14427e is below the lightness band at
# L 0.384 and #09b3a6 gives only 2.55:1 against the paper. Worst measured pair
# here is #e4682f vs #0e9a90, dE 12.6 protan / 17.0 normal vision.
# Re-step with the validator, not by eye, if these move.
LIBERATION_COLOURS <- c(
  "alive, on IMV"                    = "#35649f",
  "alive, free of IMV"               = "#0e9a90",
  "discharged alive, charted on IMV" = "#e4682f",
  "died"                             = "#93331a")

# ---- Numbers-at-risk table ---------------------------------------------------
# A survival-curve risk table, for a stacked area: a reader cannot take a
# percentage off a band to better than about five points, nor the denominator
# off it at all, and here the same bands are drawn twice on denominators that
# differ by half the cohort. Rows run in STACKING order top to bottom --
# position_stack() puts factor level 1 at the TOP, so the discrete y axis is
# reversed to match. End columns are anchored, not centred, to stay in panel.
risk_table_panel <- function(prevalence, windows, denom_label, x_title = NULL,
                             text_size = 2.3, header = "hour", limits = NULL,
                             theme_fn = house) {
  stopifnot("risk_table_panel needs a factor state, in stacking order" =
              is.factor(prevalence$state),
            "windows must carry start_hr and label" =
              all(c("start_hr", "label") %in% names(windows)))
  lv  <- levels(prevalence$state)
  tab <- prevalence[prevalence$window_start_hr %in% windows$start_hr, ]
  den <- aggregate(n ~ window_start_hr, data = tab, FUN = sum)
  num <- function(x) formatC(x, format = "d", big.mark = ",")
  # What the column labels COUNT IN. A day-scaled caller passes "day"; leaving
  # this hardcoded would head a column of days with the word hour.
  HEADER <- header

  cells <- rbind(
    # The header names the hour each column reports; there is no subtitle.
    data.frame(window_start_hr = windows$start_hr, row = HEADER,
               label = windows$label),
    data.frame(window_start_hr = tab$window_start_hr,
               row   = as.character(tab$state),
               label = sprintf("%s (%.1f%%)", num(tab$n), tab$pct)),
    data.frame(window_start_hr = den$window_start_hr,
               row   = denom_label, label = num(den$n)))

  # Reversed: a discrete y axis counts up, position_stack() puts level 1 on top.
  rows <- c(HEADER, lv, denom_label)
  cells$row   <- factor(cells$row, levels = rev(rows))
  # Anchored at the PANEL edges, which are the limits when given and the first
  # and last column otherwise. Anchoring a column that is not at the edge pushes
  # its text at the neighbour instead of away from the border.
  edge <- if (is.null(limits)) range(windows$start_hr) else limits
  cells$hjust <- ifelse(cells$window_start_hr == edge[1], 0,
                 ifelse(cells$window_start_hr == edge[2], 1, 0.5))
  cells$face  <- ifelse(cells$row %in% c(HEADER, denom_label), "bold", "plain")

  # `panel.grid = element_blank()` does NOT override house()'s explicitly-set
  # panel.grid.major.y, which struck through the text; blank each by name.
  p <- theme_fn(
    ggplot(cells, aes(window_start_hr, row)) +
      geom_hline(yintercept = seq_len(length(rows) - 1L) + 0.5,
                 colour = GRIDLINE, linewidth = 0.3) +
      geom_text(aes(label = label, hjust = hjust, fontface = face),
                size = text_size, colour = INK) +
      # `limits` lets the table report FEWER columns than the chart has ticks
      # without drifting off them: without it the table's x range is its own
      # data's, and two adjacent ticks become two misaligned axes.
      scale_x_continuous(breaks = windows$start_hr, limits = limits) +
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
# Plots drawn one above the other on a shared x axis. WIDTHS: y-label widths
# differ, so unit.pmax() across every gtable column aligns the panels, or the
# shared axis is a lie. HEIGHTS: a ggplot panel is `1null` (equal shares),
# which would give a risk table half the figure; a panel passed an absolute
# height gets it. `heights` is in inches, one per plot, NA to stay flexible.
stack_panels <- function(plots, heights = NULL) {
  stopifnot("stack_panels needs at least two plots" = length(plots) >= 2)
  if (is.null(heights)) heights <- rep(NA_real_, length(plots))
  stopifnot("one height per plot" = length(heights) == length(plots))

  # ggplotGrob() measures text, which needs an open device; under Rscript R opens
  # the DEFAULT one and leaves an Rplots.pdf behind. Measure on a null device.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  gs <- lapply(plots, ggplotGrob)

  # A bottom legend belongs under the WHOLE stack: between panels it would push
  # the risk table off the axis ticks its columns are aligned to.
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

# The risk-table case: the chart stays flexible, the table gets its rows' height.
stack_with_risk_table <- function(p_main, p_table, row_height = 0.20) {
  stack_panels(list(p_main, p_table),
               c(NA, row_height * attr(p_table, "n_rows")))
}

# ---- The 28-day liberation / survival figure ---------------------------------
# After MODE (Semler et al., NEJM Evid 2024) Figure 3: patients alive against
# liberation from invasive ventilation, over the outcome horizon. Two additions
# that figure did not need, because a trial ascertains 28-day vital status and a
# retrospective CLIF extract does not:
#
#   - a STATE-OCCUPANCY panel beneath, so "discharged alive while still charted
#     on IMV" is a visible band rather than an invisible term inside the gap
#     between the two curves. That group is scored VFD = 0 by rule, and a reader
#     has to be able to see how much of the zero mass it is.
#   - the cumulative curve is the first QUALIFYING liberation (>= the protocol's
#     successful_extubation_hours), not any extubation.
#
# Shared by the site script and the pooled one so the two cannot drift. `occ` is
# the shipped CSV, read back rather than recomputed; everything drawn here comes
# off it.
#
# NO LEGEND ON THE TOP PANEL. stack_panels() lifts the FIRST bottom legend it
# finds and nulls the rest, so a second legend would be silently dropped; the
# two curves are direct-labelled instead, which is what <= 4 series want anyway.
# `theme_fn` is the dress, not the drawing: the pooled script wears house_ats()
# and every figure beside this one in output/pooled/figures/ already does.
liberation_panels <- function(occ, denom, landmark_d, horizon_d,
                              breaks_d = c(0, 3, 7, 14, 21, 28),
                              theme_fn = house) {
  stopifnot(
    "occ must carry hour, measure, n_admissions" =
      all(c("hour", "measure", "n_admissions") %in% names(occ)),
    "occ must carry the four occupancy states and the cumulative liberation" =
      all(c(names(LIBERATION_COLOURS), "cum_liberated") %in% occ$measure))

  brk <- breaks_d[breaks_d <= horizon_d] * 24
  x_scale <- scale_x_continuous(
    breaks = brk, labels = brk / 24, limits = c(0, horizon_d * 24),
    expand = expansion(mult = c(0.01, 0.03)))
  landmark <- geom_vline(xintercept = landmark_d * 24, colour = MUTED,
                         linewidth = 0.4, linetype = "22")
  pct <- function(n) 100 * n / denom

  # ---- top: alive, and cumulative incidence of first qualifying liberation ---
  wide <- stats::reshape(occ[, c("hour", "measure", "n_admissions")],
                         idvar = "hour", timevar = "measure",
                         direction = "wide")
  names(wide) <- sub("^n_admissions\\.", "", names(wide))
  # Alive is the complement of the dead BAND, not a separate count: one number
  # for the cohort's mortality, so the panel and the area cannot disagree.
  curves <- rbind(
    data.frame(hour = wide$hour, series = "Alive",
               value = pct(denom - wide$died)),
    # NOT the band below: that one is the share free of IMV AT t, this is the
    # share ever liberated BY t. At day 28 they differ by everyone who was
    # liberated and then died or was reintubated.
    data.frame(hour = wide$hour, series = "Liberated from IMV",
               value = pct(wide$cum_liberated)))
  curves$series <- factor(curves$series,
                          levels = c("Alive", "Liberated from IMV"))
  ends <- do.call(rbind, lapply(split(curves, curves$series), function(d)
    d[which.max(d$hour), ]))

  p_top <- theme_fn(
    ggplot(curves, aes(hour, value, colour = series, linetype = series)) +
      landmark +
      geom_step(linewidth = 0.7) +
      geom_text(data = ends, aes(label = series), hjust = 1, vjust = -0.8,
                size = 3.0, colour = INK, show.legend = FALSE) +
      scale_colour_manual(
        values = c("Alive" = INK,
                   "Liberated from IMV" =
                     unname(LIBERATION_COLOURS[["alive, free of IMV"]])),
        guide = "none") +
      scale_linetype_manual(values = c("Alive" = "solid",
                                       "Liberated from IMV" = "21"),
                            guide = "none") +
      x_scale +
      scale_y_continuous(limits = c(0, 104), breaks = seq(0, 100, 25),
                         labels = function(v) paste0(v, "%")) +
      labs(x = NULL, y = "Patients"))

  # ---- bottom: state occupancy ---------------------------------------------
  # Reversed levels: position_stack() puts level 1 at the TOP, and `free of IMV`
  # belongs on the baseline where its area can be read against the axis.
  states <- occ[occ$measure %in% names(LIBERATION_COLOURS), ]
  states$measure <- factor(states$measure, levels = rev(names(LIBERATION_COLOURS)))
  p_mid <- theme_fn(
    ggplot(states, aes(hour, pct(n_admissions), fill = measure)) +
      geom_area(position = "stack", colour = NA) +
      landmark +
      scale_fill_manual(values = LIBERATION_COLOURS, name = NULL,
                        breaks = names(LIBERATION_COLOURS)) +
      x_scale +
      # NO HARD ylim. The four shares are computed independently and summed in
      # floating point, so a full stack lands on 100.00000000000003 as often as
      # on 100 -- and ggplot DROPS data outside the scale limits rather than
      # clipping it, silently losing the top of the area.
      scale_y_continuous(breaks = seq(0, 100, 25),
                         labels = function(v) paste0(v, "%"),
                         expand = expansion(mult = c(0, 0.02))) +
      labs(x = NULL, y = "Patients") +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE)))

  # ---- the risk table, on the same x -----------------------------------------
  # Day 0 is dropped from the TABLE only: it is 100% ventilated by
  # construction, and at 3 days from the next column on a 28-day axis the two
  # cells overprint.
  tab_brk <- brk[brk > 0]
  tab <- states[states$hour %in% tab_brk, ]
  names(tab)[names(tab) == "hour"] <- "window_start_hr"
  names(tab)[names(tab) == "measure"] <- "state"
  tab$state <- factor(as.character(tab$state), levels = names(LIBERATION_COLOURS))
  tab$n <- tab$n_admissions
  tab$pct <- pct(tab$n_admissions)
  p_tab <- risk_table_panel(
    tab, data.frame(start_hr = tab_brk, label = as.character(tab_brk / 24)),
    denom_label = "all", header = "day", limits = c(0, horizon_d * 24),
    theme_fn = theme_fn,
    x_title = "Days from the start of invasive mechanical ventilation")

  # The top panel is two lines and wants a fixed, modest height; the occupancy
  # area carries the finding and takes whatever the page has left. Pinning BOTH
  # charts in inches (as a first draft did) overruns the page, and ggsave then
  # crops the top of the figure and the legend without warning.
  stack_panels(list(p_top, p_mid, p_tab),
               c(2.0, NA, 0.20 * attr(p_tab, "n_rows")))
}

# ---- A full-width panel over a side-by-side pair ------------------------------
# patchwork is not in renv and the pooled script is the only caller, so the
# composition is done on gtables directly: each child keeps its own axes and the
# outer table only has to place them. This cannot be rbind()/cbind() the way
# stack_panels() does -- a side-by-side row has twice the gtable columns of the
# panel above it, so there is nothing to align -- hence an outer table with one
# column and one row per band.
#
# The x-axis TITLE and the LEGEND are lifted out of the children and drawn once
# under the whole figure: two 3.4in panels would otherwise each repeat a title
# wider than the panel, and a legend per panel says the same thing three times.
compose_panels <- function(top = NULL, bottom, heights, x_title = NULL,
                           x_title_gp = grid::gpar(col = INK, fontsize = 12)) {
  stopifnot("bottom must be a non-empty list of plots" =
              is.list(bottom) && length(bottom) > 0,
            # One height per BAND, not per plot: the side-by-side pair shares a
            # row, and the x title and legend size themselves off their grobs.
            "heights is the top band (if any) then the bottom row" =
              length(heights) == 1L + !is.null(top))

  # ggplotGrob() measures text, which needs an open device; under Rscript R opens
  # the DEFAULT one and leaves an Rplots.pdf behind. Measure on a null device.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  legend <- NULL
  as_grob <- function(p) {
    g <- ggplotGrob(p)
    i <- which(g$layout$name == "guide-box-bottom")
    if (length(i) == 1L && !inherits(g$grobs[[i]], "zeroGrob")) {
      if (is.null(legend)) legend <<- g$grobs[[i]]
      g$grobs[[i]] <- grid::nullGrob()
      g$heights[g$layout$t[i]] <- grid::unit(0, "cm")
    }
    g
  }

  gb <- lapply(bottom, as_grob)
  # Equal gtable widths across the pair, so the two PANELS line up rather than
  # the two y-axis labels: a wider tick label would otherwise shift one panel.
  w  <- Reduce(grid::unit.pmax, lapply(gb, function(g) g$widths))
  gb <- lapply(gb, function(g) { g$widths <- w; g })
  row_bottom <- gtable::gtable(widths  = grid::unit(rep(1, length(gb)), "null"),
                               heights = grid::unit(1, "null"))
  for (j in seq_along(gb))
    row_bottom <- gtable::gtable_add_grob(row_bottom, gb[[j]], t = 1, l = j,
                                          name = paste0("bottom", j))

  bands <- list()
  add <- function(g, h) bands[[length(bands) + 1L]] <<- list(g = g, h = h)
  if (!is.null(top)) add(as_grob(top), grid::unit(heights[1], "in"))
  add(row_bottom, grid::unit(heights[length(heights)], "in"))
  if (!is.null(x_title)) {
    lab <- grid::textGrob(x_title, gp = x_title_gp)
    add(lab, grid::grobHeight(lab) + grid::unit(8, "pt"))
  }
  if (!is.null(legend)) add(legend, grid::grobHeight(legend))

  out <- gtable::gtable(widths  = grid::unit(1, "null"),
                        heights = do.call(grid::unit.c, lapply(bands, `[[`, "h")))
  for (i in seq_along(bands))
    out <- gtable::gtable_add_grob(out, bands[[i]]$g, t = i, l = 1,
                                   name = paste0("band", i))
  out
}

# ---- Captions ----------------------------------------------------------------
# Journal panels carry no title or subtitle, so what would have sat there is
# registered here and written out as captions.md beside the figures. An
# environment, not a list, so a draw function can register without `<<-`.
.captions <- new.env(parent = emptyenv())

register_caption <- function(file, title, note) {
  assign(file, list(title = title, note = note), envir = .captions)
  invisible(NULL)
}

# Pass the .png entries of the script's OWNED list; a missing caption is an error.
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
# ggplot gives a panel one secondary axis, so three scales -- one left, two
# right -- cannot come from scales alone. This borrows a right-hand axis from a
# donor drawn on the SAME panel range and hangs it outside the existing right
# axis. The donor must span the identical y range or the ticks plausibly lie.
add_third_axis <- function(p_main, donor, title = NULL, title_size = 11) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  g  <- ggplotGrob(p_main)
  gd <- ggplotGrob(donor)

  ax <- gd$grobs[[which(gd$layout$name == "axis-r")]]
  pos <- g$layout[g$layout$name == "axis-r", ]
  # ggplot lays out an axis-r SLOT whether or not a secondary axis exists, filling
  # it with a zeroGrob, so check the grob rather than the slot.
  stopifnot(
    "the main plot has no right axis to hang a third scale beside" =
      nrow(pos) == 1,
    "the main plot's right axis is empty: give it a sec.axis first" =
      !inherits(g$grobs[[which(g$layout$name == "axis-r")]], "zeroGrob"),
    "the donor plot has no right axis to borrow" = !inherits(ax, "zeroGrob"))

  # OUTSIDE the existing right axis AND its title, so each scale keeps its label.
  ttl <- g$layout[g$layout$name %in% c("ylab-r", "axis-r-title", "ylab-r-title"), ]
  at <- if (nrow(ttl)) max(ttl$r) else max(pos$r)

  g <- gtable::gtable_add_cols(g, grid::unit(3, "mm"), pos = at)
  g <- gtable::gtable_add_cols(g, gd$widths[gd$layout$r[gd$layout$name == "axis-r"]],
                               pos = at + 1)
  g <- gtable::gtable_add_grob(g, ax, t = pos$t, b = pos$b, l = at + 2,
                               name = "axis-r-third")
  if (!is.null(title)) {
    lab <- grid::textGrob(title, rot = -90,
                          gp = grid::gpar(col = INK, fontsize = title_size))
    g <- gtable::gtable_add_cols(g, grid::unit(1.2, "lines"), pos = at + 2)
    g <- gtable::gtable_add_grob(g, lab, t = pos$t, b = pos$b, l = at + 3,
                                 name = "ylab-r-third")
  }
  g
}

# ---- Adjusted vs unadjusted ---------------------------------------------------
# The prediction figure (07_outcome_model.R) reuses the validated
# initiation/uptitration hues. One entry, because it now draws one series.
PREDICTION_COLOURS <- c("adjusted" = "#234e90")

# ---- ATS presentation palette -------------------------------------------------
# For the POOLED (audience-facing) figures. Brand values supplied 2026-10-28 and
# VALIDATED rather than applied by eye -- dataviz scripts/validate_palette.js,
# light surface #fcfcfb, --pairs all because a caterpillar sorts by value so
# adjacency is arbitrary:
#
#   #0C0F36  chroma 0.075, L 0.196 -> measures as an INK, not a hue. 17.99:1
#            against the paper (near-black is 19.17:1), so it carries type,
#            spines and reference rules at no legibility cost.
#   #1977F9  kept as categorical 1.
#   #0FADEE  normal-vision dE 14.8 against #1977F9, BELOW the 15 floor: two blues
#            a full-colour reader cannot separate. That is the one failure
#            secondary encoding does not excuse, so this is FILL-ONLY (ribbons,
#            low alpha) and never a categorical neighbour of #1977F9.
#   #F1A41F  L 0.775 and 2.03:1 against the paper -> darkened to #C77F12
#            (3.16:1, inside the band). SG approved the substitution.
#   #863A68  kept as categorical 3.
#
# The surviving trio passes every check: worst all-pairs dE 22.6 (deutan), 25.8
# (normal), all >= 3:1 contrast. Re-run the validator if any value moves.
#
# THREE HUES DOES NOT COVER EVERY POOLED FIGURE. The unit caterpillars and
# funnels colour by care-setting display group (7 at k=2) and the delivery-state
# prevalence draws 7 states. Cycling three hues over seven categories is
# forbidden, so those two keep their existing qualitative ramps and only their
# TYPOGRAPHY becomes ATS. Supply more brand hues, or collapse the groups, if they
# should go on-brand.
ATS_INK    <- "#0C0F36"
ATS_SERIES <- c("#1977F9", "#C77F12", "#863A68")
ATS_FILL   <- "#0FADEE"
ATS_FONT   <- "Arial"

# The house theme in ATS dress: same geometry and the same legend-at-the-foot
# rule, with brand ink and Arial. Type still reads near-black-dark, never grey.
house_ats <- function(p) {
  house(p) +
    theme(text = element_text(family = ATS_FONT, colour = ATS_INK),
          plot.title = element_text(family = ATS_FONT, colour = ATS_INK,
                                    face = "bold"),
          plot.subtitle = element_text(family = ATS_FONT, colour = ATS_INK),
          axis.title = element_text(family = ATS_FONT, colour = ATS_INK),
          axis.text = element_text(family = ATS_FONT, colour = ATS_INK),
          legend.text = element_text(family = ATS_FONT, colour = ATS_INK,
                                     size = 9),
          strip.text = element_text(family = ATS_FONT, colour = ATS_INK,
                                    face = "bold"),
          legend.background = element_rect(fill = PAPER, colour = ATS_INK,
                                           linewidth = 0.3))
}
