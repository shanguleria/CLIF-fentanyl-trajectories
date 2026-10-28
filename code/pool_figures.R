# Purpose : Figures for the pooled multi-site analysis. Run at the coordinating
#           centre over output/pooled/, after code/pool_sites.R.
# Inputs  : output/pooled/*.csv
# Outputs : output/pooled/figures/ : F2B, F3, F4 (+ F1 from pool_strobe.py)
#
#   Rscript code/pool_figures.R
#
# Anonymised by default: sites become "Site 1, ..." and units "S1-U1", numbered
# alphabetically by real name so no label encodes a rate (see pooling_spec.json).
#
# Journal style: no title or subtitle -- n, denominators, counts and the reference
# line definition go in captions.md. Legend at the foot in a box, type near-black.

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
# Pools exactly: n per (window, state) over a denominator fixed at every window.
# Replaces an alluvial, which cannot pool -- a ribbon is one patient's path.

sp <- rd("04_states__state_prevalence.csv")
sp <- sp[sp$site == "POOLED", ]
sp$state <- factor(sp$state, levels = names(STATE_COLOURS))
stopifnot("a pooled prevalence state is outside the declared palette" =
            !any(is.na(sp$state)))
# Rows must sum to 100% within a window: checked rather than assumed.
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
    # Short form -- the long one truncates. Full definition is in the caption.
    y = "% with a qualifying score"),
  coadmin = list(series = "any increase", panel = "B",
    x = "% of rate increases with a fentanyl bolus within 30 min",
    y = "% paired with a bolus"))

prep <- function(metric) {
  spec <- METRIC[[metric]]
  d <- ua[ua$metric == metric & ua$series == spec$series, ]
  n_all <- nrow(d)
  # A suppressed unit cannot be drawn, so it is counted in the caption instead.
  d <- d[!is.na(d$pct_paired) & !is.na(d$n_paired), ]
  list(d = d, n_suppressed = n_all - nrow(d), spec = spec)
}

# Reference line is the pooled rate over the units drawn, and must share their
# denominator: attributed events only, unlike coadministration.csv's overall rate.
pooled_rate <- function(d) 100 * sum(d$n_paired) / sum(d$n_events)

draw_pair <- function(kind) {
  plots <- list(); caps <- list()
  for (metric in names(METRIC)) {
    pr <- prep(metric); d <- pr$d; spec <- pr$spec
    if (!nrow(d)) next
    d$grp <- cs_label(d$care_setting)
    d$site_lab <- site_label[d$site]
    # Numbered alphabetically by real name within site, so no label encodes a rate.
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


# ---- F5. Forest plot of the exposure odds ratios ------------------------------
# Per-site estimates and the inverse-variance pool, on the PRIMARY model. The
# sensitivities are in the CSV; six models on one panel is unreadable.
cf_f <- "07_outcome_model__model_coefficients.csv"
if (file.exists(file.path(POOL, cf_f))) {
  cf <- rd(cf_f)
  cf <- cf[cf$model == "primary" &
             toupper(as.character(cf$is_exposure)) == "TRUE", , drop = FALSE]
  cf <- cf[is.finite(suppressWarnings(as.numeric(cf$odds_ratio))), , drop = FALSE]

  if (nrow(cf)) {
    for (nm in c("odds_ratio", "or_ci_lo", "or_ci_hi", "std_error"))
      cf[[nm]] <- suppressWarnings(as.numeric(cf[[nm]]))
    cf$is_pool <- cf$site == "POOLED"
    # n on the row label, the way the per-unit figures do it: a reader cannot
    # take a denominator off a square.
    nfmt <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")
    cf$label <- ifelse(cf$is_pool,
                       sprintf("Pooled (n = %s)", nfmt(cf$n_admissions)),
                       sprintf("%s (n = %s)", site_label[cf$site], nfmt(cf$n_admissions)))
    # POOLED last, then sites in the declared alphabetical order, and the y axis
    # reversed so Site 1 reads at the top.
    lab_of <- function(s) unique(cf$label[cf$site == s])
    lev <- c(lab_of("POOLED"), rev(unlist(lapply(SITES, lab_of))))
    cf$label <- factor(cf$label, levels = lev)
    # Point area proportional to inverse-variance weight, the forest convention.
    cf$wt <- ifelse(cf$is_pool, NA, 1 / cf$std_error^2)
    cf$wt_rel <- cf$wt / max(cf$wt, na.rm = TRUE)

    # "Pbolus_per10pp" -> "P(bolus)": the unit is on the x axis, so repeating
    # per10pp in every strip adds ink without information.
    cf$term_label <- sub("^P", "P(", sub("_per[0-9]+pp$", ")", cf$term))
    het <- cf[cf$is_pool, c("term", "term_label", "k_sites", "Q", "Q_df", "Q_p", "I2_pct",
                            "tau2", "heterogeneity_role")]
    rng <- range(c(cf$or_ci_lo, cf$or_ci_hi), na.rm = TRUE)

    p_f5 <- ggplot(cf, aes(x = odds_ratio, y = label)) +
      geom_vline(xintercept = 1, colour = MUTED, linewidth = 0.4) +
      geom_errorbar(aes(xmin = or_ci_lo, xmax = or_ci_hi), orientation = "y",
                    width = 0.18, linewidth = 0.5, colour = INK) +
      geom_point(aes(size = wt_rel, shape = is_pool), colour = INK, fill = INK) +
      scale_shape_manual(values = c(`FALSE` = 22, `TRUE` = 23), guide = "none") +
      scale_size_area(max_size = 4.2, guide = "none", na.value = 3.6) +
      scale_x_continuous(trans = "log10",
                         breaks = scales::breaks_log(5),
                         labels = function(x) formatC(x, format = "fg", digits = 3)) +
      facet_wrap(~ term_label, ncol = 1, scales = "free_y") +
      labs(x = sprintf("Odds of a better VFD-28 band, per %d percentage points (log scale)",
                       round(100 * fromJSON(here("config", "protocol.json"),
                                            simplifyVector = FALSE)$model$outcome_model$exposure_scale_per)),
           y = NULL)
    p_f5 <- house(p_f5) + theme(panel.grid.major.y = element_blank(),
                                panel.grid.major.x = element_line(colour = GRIDLINE,
                                                                  linewidth = 0.4))
    ggsave(file.path(FIG, "F5_exposure_forest_pooled.png"), p_f5,
           width = 7.0, height = 1.4 + nrow(het) * (0.34 * length(lev) + 0.55),
           dpi = 300,
           bg = PAPER)
    register_caption(
      "F5_exposure_forest_pooled.png",
      "Does documenting an indication for a fentanyl increase track ventilator-free days, across sites?",
      sprintf(paste0(
        "Per-admission titration exposures against VFD-28, one panel per ",
        "exposure, from the PRIMARY proportional-odds model at each of %d ",
        "site(s). SQUARES ARE SITES, area proportional to inverse-variance ",
        "weight; the DIAMOND is the pooled estimate. Bars are 95%% intervals. ",
        "An odds ratio above 1 is higher odds of landing in a BETTER VFD-28 ",
        "band, read against the band ordering and not against intuition. ",
        "POOLED BY INVERSE VARIANCE ON THE LOG-ODDS, which is the scale that ",
        "pools -- odds ratios are never averaged. The pooled estimand is the ",
        "WITHIN-SITE, ICU-adjusted association: unit fixed effects are ",
        "site-local and are not shipped, so this is not a between-site ",
        "contrast. %s SITES ARE ANONYMISED as Site 1..n by alphabetical ",
        "site_name, the declared and reproducible rule; real names appear only ",
        "on a title page. Sensitivity models are in ",
        "07_outcome_model__model_coefficients.csv and are not drawn here. AN ",
        "ASSOCIATION, NOT AN EFFECT -- confounding by indication is ",
        "unaddressed."),
        as.integer(max(het$k_sites, na.rm = TRUE)),
        paste(sprintf("%s: Q = %.2f on %d df (p = %s), I-squared = %.0f%%, tau-squared = %.4f%s.",
                      het$term_label, het$Q, as.integer(het$Q_df), format(het$Q_p),
                      het$I2_pct, het$tau2,
                      ifelse(het$heterogeneity_role == "descriptive_only",
                             sprintf(" -- DESCRIPTIVE ONLY, k is below heterogeneity_min_k = %d, so the fixed-effect estimate is primary",
                                     as.integer(SPEC$heterogeneity_min_k)), "")),
              collapse = " ")))
    cat("  F5_exposure_forest_pooled.png\n")
  }
}

# Only the figures this script drew; F1 registers its own caption in the Python
# STROBE step.
MINE <- grep("^(F2B|F3|F4|F5)", list.files(FIG, pattern = "\\.png$"), value = TRUE)
write_captions(file.path(FIG, "captions.md"), "pool_figures.R", MINE,
               list(site_name = sprintf("%d pooled", K),   # write_captions prefixes "Site "
                    clif_version = PROV$definition_versions$`protocol.json`,
                    code_version = PROV$code_versions[[1]],
                    generated = format(Sys.time(), "%Y-%m-%d %H:%M %Z")))
cat(sprintf("\n%d figure(s) in output/pooled/figures/\n",
            length(list.files(FIG, pattern = "\\.png$"))))
