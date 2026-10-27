# Purpose : Figures for the pooled multi-site analysis. Run at the COORDINATING
#           CENTRE over output/pooled/, after code/pool_sites.R.
# Author  : Shan Guleria
# Created : 2026-10-27
# Inputs  : output/pooled/*.csv
# Outputs : output/pooled/figures/ : F2B, F3, F4 (+ F1 from pool_strobe.py)
#
#   Rscript code/pool_figures.R
#
# ANONYMISED BY DEFAULT. Sites become "Site 1, Site 2, ..." by the declared rule
# in pooling_spec.json (alphabetical on site_name) and units become "S1-U1" and
# so on, numbered ALPHABETICALLY BY THEIR REAL NAME -- never by their rate, which
# would leak the ordering into the label. Real names stay in the CSVs, which SG
# reads; these panels are what an audience sees. A unit code is not PHI, but it
# is site-identifying, and Rush's names run ~45 characters and are unreadable on
# a 160mm panel besides.
#
# JOURNAL STYLE: no title, no subtitle. The n, the denominator, the site count,
# the suppression count and the reference-line definition go in captions.md,
# which is manuscript text. Legend at the foot in a box. Type in near-black.

suppressMessages({
  library(here); library(jsonlite); library(ggplot2)
})
source(here("code", "utils", "figures.R"))

POOL <- here("output", "pooled")
FIG  <- file.path(POOL, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(POOL, "pooled_provenance.json")))
  stop("run code/pool_sites.R first -- ", POOL, " has no pooled_provenance.json",
       call. = FALSE)

SPEC <- fromJSON(here("config", "pooling_spec.json"), simplifyVector = FALSE)
PROV <- fromJSON(file.path(POOL, "pooled_provenance.json"), simplifyVector = FALSE)
SITES <- sort(unlist(PROV$sites))          # the declared rule: alphabetical
K <- length(SITES)
CS_GROUP <- SPEC$care_setting_display_groups

rd <- function(f) utils::read.csv(file.path(POOL, f), stringsAsFactors = FALSE,
                                  check.names = FALSE)

# ---- anonymisation -----------------------------------------------------------
site_label <- setNames(sprintf("Site %d", seq_along(SITES)), SITES)
cs_label <- function(v) {
  out <- vapply(as.character(v), function(x) {
    g <- CS_GROUP[[x]]; if (is.null(g)) x else g
  }, character(1), USE.NAMES = FALSE)
  out
}


# ---- F2B. Pooled delivery-state prevalence ------------------------------------
# Pools EXACTLY: n per (window, state) over a denominator fixed at every window,
# so summing both and recomputing is the whole method. This is the panel that
# replaces an alluvial, which cannot pool -- a ribbon is one patient's path and
# the per-episode sequences never leave a site.

sp <- rd("04_states__state_prevalence.csv")
sp <- sp[sp$site == "POOLED", ]
sp$state <- factor(sp$state, levels = names(STATE_COLOURS))
stopifnot("a pooled prevalence state is outside the declared palette" =
            !any(is.na(sp$state)))
# Rows sum to 100% within a window -- verified here rather than assumed, because
# it is the property that makes the denominator self-contained.
tot <- tapply(sp$pct, sp$window_start_hr, sum)
stopifnot("pooled prevalence does not sum to 100% within a window" =
            all(abs(tot - 100) < 0.5))

p_prev <- house(
  ggplot(sp, aes(window_start_hr, pct, fill = state)) +
    geom_area(colour = PAPER, linewidth = 0.2) +
    scale_fill_manual(values = STATE_COLOURS, name = NULL, drop = FALSE) +
    scale_x_continuous(breaks = sort(unique(sp$window_start_hr))[c(TRUE, FALSE)],
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(x = "Hours since the first IMV episode",
         y = "% of ventilation episodes"))
n_pooled <- max(tapply(sp$n, sp$window_start_hr, sum))
register_caption("F2B_state_prevalence_pooled.png",
  "Delivery-state prevalence over the first 72 hours, pooled",
  sprintf(paste0("%s ventilation episodes across %d sites. Numerator: episodes ",
                 "in that state in the window. Denominator: all anchor episodes, ",
                 "FIXED at every window, with terminal states absorbing -- so the ",
                 "bands sum to 100%% throughout and liberation and death are ",
                 "visible as the fentanyl states drain. Pools exactly by summing ",
                 "both counts. An alluvial of the same data CANNOT be pooled: a ",
                 "ribbon is one patient's path and the per-episode sequences do ",
                 "not leave a site."),
          format(n_pooled, big.mark = ","), K))
ggsave(file.path(FIG, "F2B_state_prevalence_pooled.png"), p_prev,
       width = 8.4, height = 5.0, dpi = 200)
cat("  F2B_state_prevalence_pooled.png\n")


# ---- F3 / F4. Unit variation, both metrics ------------------------------------

ua <- rd("06_unit_variation__unit_adherence.csv")
ua <- ua[ua$unit_key == "location_name" & ua$site != "POOLED", ]

METRIC <- list(
  indication = list(series = "any escalation", panel = "A",
    x = "% of escalations with a qualifying pain or agitation score within 1 h",
    # A y-axis title has the PANEL HEIGHT to fit in, not its width. The long
    # form truncated to "...within 30 m" on the funnel; the full definition is
    # in the caption either way.
    y = "% with a qualifying score"),
  coadmin = list(series = "any increase", panel = "B",
    x = "% of rate increases with a fentanyl bolus within 30 min",
    y = "% paired with a bolus"))

prep <- function(metric) {
  spec <- METRIC[[metric]]
  d <- ua[ua$metric == metric & ua$series == spec$series, ]
  n_all <- nrow(d)
  # A suppressed unit cannot be drawn. It is COUNTED in the caption rather than
  # vanishing: a panel that silently drops a unit misstates how many there are.
  d <- d[!is.na(d$pct_paired) & !is.na(d$n_paired), ]
  list(d = d, n_suppressed = n_all - nrow(d), spec = spec)
}

# The reference line is the POOLED rate: summed numerator over summed
# denominator across every unit drawn. It must share the denominator of the
# points it is drawn against -- these are attributed events only, so the overall
# rate from coadministration.csv (which includes unattributed events) is a
# DIFFERENT quantity and would sit in the wrong place.
pooled_rate <- function(d) 100 * sum(d$n_paired) / sum(d$n_events)

draw_pair <- function(kind) {
  plots <- list(); caps <- list()
  for (metric in names(METRIC)) {
    pr <- prep(metric); d <- pr$d; spec <- pr$spec
    if (!nrow(d)) next
    d$grp <- cs_label(d$care_setting)
    d$site_lab <- site_label[d$site]
    # Unit labels numbered alphabetically BY REAL NAME within site, so the label
    # says nothing about the rate.
    d <- d[order(d$site, d$unit), ]
    d$unit_lab <- mapply(function(s, i) sprintf("S%d-U%d", match(s, SITES), i),
                         d$site, ave(seq_len(nrow(d)), d$site, FUN = seq_along))
    ref <- pooled_rate(d)
    pal <- grDevices::hcl.colors(length(unique(d$grp)), "Dark 3")
    names(pal) <- sort(unique(d$grp))

    if (kind == "caterpillar") {
      d <- d[order(d$pct_paired), ]
      d$unit_lab <- factor(d$unit_lab, levels = d$unit_lab)
      p <- house(
        ggplot(d, aes(pct_paired, unit_lab, colour = grp)) +
          geom_vline(xintercept = ref, colour = INK, linetype = "22",
                     linewidth = 0.45) +
          geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y",
                        width = 0, linewidth = 0.5) +
          geom_point(size = 2.1) +
          scale_colour_manual(values = pal, name = NULL) +
          guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
          annotate("text", x = ref, y = Inf, label = sprintf(" pooled %.1f%%", ref),
                   hjust = 0, vjust = 1.4, colour = INK, size = 3.1) +
          scale_x_continuous(expand = expansion(mult = c(0.06, 0.10))) +
          labs(x = spec$x, y = NULL))
      h <- max(3.2, 1.3 + 0.34 * nrow(d))
    } else {
      grid_n <- seq(max(50, min(d$n_events) * 0.6), max(d$n_events) * 1.08,
                    length.out = 220)
      se <- sqrt(ref * (100 - ref) / grid_n)
      band <- rbind(
        data.frame(n = grid_n, lo = ref - 1.96 * se, hi = ref + 1.96 * se,
                   band = "95%"),
        data.frame(n = grid_n, lo = ref - 3.09 * se, hi = ref + 3.09 * se,
                   band = "99.8%"))
      p <- house(
        ggplot() +
          geom_ribbon(data = band, aes(n, ymin = lo, ymax = hi, group = band),
                      fill = GRIDLINE, alpha = 0.55) +
          geom_hline(yintercept = ref, colour = INK, linetype = "22",
                     linewidth = 0.45) +
          geom_point(data = d, aes(n_events, pct_paired, colour = grp),
                     size = 2.2) +
          scale_colour_manual(values = pal, name = NULL) +
          guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
          annotate("text", x = max(d$n_events), y = ref,
                   label = sprintf("pooled %.1f%%", ref), hjust = 1.02,
                   vjust = -0.6, colour = INK, size = 3.1) +
          scale_x_continuous(labels = function(x) format(x, big.mark = ",")) +
          scale_y_continuous(expand = expansion(mult = c(0.08, 0.10))) +
          labs(x = "Events contributed by the unit (volume)", y = spec$y))
      h <- 5.0
    }
    plots[[metric]] <- list(p = p, h = h)
    caps[[metric]] <- list(
      n_units = nrow(d), n_sup = pr$n_suppressed, ref = ref,
      lo = min(d$pct_paired), hi = max(d$pct_paired),
      n_ev = sum(d$n_events), n_pa = sum(d$n_paired))
  }
  list(plots = plots, caps = caps)
}

for (kind in c("caterpillar", "funnel")) {
  res <- draw_pair(kind)
  fnum <- if (kind == "caterpillar") "F3" else "F4"
  for (metric in names(res$plots)) {
    spec <- METRIC[[metric]]; cc <- res$caps[[metric]]
    fn <- sprintf("%s%s_%s_%s.png", fnum, spec$panel, metric, kind)
    ggsave(file.path(FIG, fn), res$plots[[metric]]$p,
           width = 8.4, height = res$plots[[metric]]$h, dpi = 200)
    register_caption(fn,
      sprintf("%s %s. %s by unit, %d units across %d sites",
              fnum, spec$panel,
              if (metric == "indication") "Documented indication" else
                "Bolus co-administration",
              cc$n_units, K),
      sprintf(paste0("Each point is one ICU or ED, ANONYMISED: sites are ",
                     "numbered alphabetically by name and units alphabetically ",
                     "within site, so no label encodes a rate. Colour is the ",
                     "unit's care setting, grouped for display (only 2 of 10 raw ",
                     "mCIDE settings occur at both sites, so the raw value would ",
                     "be a site label). Observed range %.1f%% to %.1f%%. The ",
                     "dashed reference is the POOLED rate, %s/%s = %.1f%%, summed ",
                     "numerator over summed denominator across the units drawn -- ",
                     "a description, NOT a target, and no evidence establishes ",
                     "that a higher rate is better care. %s%s"),
              cc$lo, cc$hi, format(cc$n_pa, big.mark = ","),
              format(cc$n_ev, big.mark = ","), cc$ref,
              if (cc$n_sup > 0) sprintf("%d unit(s) suppressed for small cell size and not drawn. ", cc$n_sup) else "",
              if (kind == "funnel")
                paste0("Bands are UNADJUSTED binomial limits at 95% and 99.8% ",
                       "around the pooled rate: where a unit of this volume would ",
                       "land IF it were truly at that rate. They are not ",
                       "corrected for over-dispersion, which is measured ",
                       "separately and is large.") else
                paste0("Bars are each unit's own bootstrap interval, which is ",
                       "WITHIN-unit precision and covers sampling only. ",
                       "Overlap is an informal read, not a pairwise test.")))
    cat(sprintf("  %s\n", fn))
  }
}

# Only the figures THIS script drew. F1 is rendered by the Python STROBE step
# and registers its own caption there; write_captions asserts that every file it
# is handed has one, which is the guard that caught the mismatch.
MINE <- grep("^(F2B|F3|F4)", list.files(FIG, pattern = "\\.png$"), value = TRUE)
write_captions(file.path(FIG, "captions.md"), "pool_figures.R", MINE,
               list(site_name = sprintf("%d pooled", K),   # write_captions prefixes "Site "
                    clif_version = PROV$definition_versions$`protocol.json`,
                    code_version = PROV$code_versions[[1]],
                    generated = format(Sys.time(), "%Y-%m-%d %H:%M %Z")))
cat(sprintf("\n%d figure(s) in output/pooled/figures/\n",
            length(list.files(FIG, pattern = "\\.png$"))))
