# Purpose : Figures for the pooled multi-site analysis. Run at the coordinating
#           centre over output/pooled/, after code/pool_sites.R.
# Inputs  : output/pooled/*.csv
# Outputs : output/pooled/figures/ : named descriptively, NOT F-numbered --
#           manuscript numbering changes and the abstract takes one figure,
#           so numbering here only invited a collision
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
# Every PNG here is drawn by this script, so clear them first: a renamed
# figure otherwise leaves its old name behind to be read as a current one.
old_pngs <- list.files(FIG, pattern = "\\.png$", full.names = TRUE)
if (length(old_pngs)) {
  invisible(file.remove(old_pngs))
  cat(sprintf("  cleared %d figure(s) from a previous run\n", length(old_pngs)))
}
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

# ---- exposure panel titles ----------------------------------------------------
# Declared in pooling_spec.json, not derived from the column name: "P(indicated)"
# is working shorthand that a 400-word abstract never defines. The KEY ORDER is
# the display order -- the primary, significant exposure is met before the null.
EXP_LAB <- SPEC$exposure_display_labels
stopifnot("pooling_spec.json must declare exposure_display_labels" =
            !is.null(EXP_LAB) && length(EXP_LAB) > 0)
EXP_ORDER <- unlist(EXP_LAB, use.names = FALSE)
# Keyed on the UNSCALED exposure names, so a term like Pindicated_per10pp resolves
# by prefix and the labels cannot drift from protocol.json's exposures list.
exp_label <- function(term) {
  bare <- sub("_per[0-9]+pp$", "", as.character(term))
  miss <- setdiff(unique(bare), names(EXP_LAB))
  stopifnot("an exposure has no declared display label" = !length(miss))
  factor(unlist(EXP_LAB[bare], use.names = FALSE), levels = EXP_ORDER)
}


# ---- anonymisation -----------------------------------------------------------
site_label <- setNames(sprintf("Site %d", seq_along(SITES)), SITES)

# ONE site -> hue map for the whole script, so a site is the same colour in
# every figure that colours by site and one legend decodes all of them. Three
# validated brand hues (figures.R) and cycling three over four is the one thing
# that palette comment forbids, so past three sites this is NULL and the figures
# that need it fall back to a form that does not colour by site.
SITE_HUE <- if (K <= length(ATS_SERIES))
  setNames(ATS_SERIES[seq_len(K)], unname(site_label[SITES])) else NULL
cs_label <- function(v) {
  out <- vapply(as.character(v), function(x) {
    g <- CS_GROUP[[x]]; if (is.null(g)) x else g
  }, character(1), USE.NAMES = FALSE)
  out
}


# ---- Pooled delivery-state prevalence ------------------------------------
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

p_prev <- house_ats(
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
register_caption("state_prevalence_pooled.png",
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
ggsave(file.path(FIG, "state_prevalence_pooled.png"), p_prev,
       width = 8.4, height = 5.0, dpi = 200)
cat("  state_prevalence_pooled.png\n")


# ---- Unit variation, both metrics: caterpillar and funnel ------------------------------------

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
    # SEVEN care-setting display groups at k = 2, so this keeps a qualitative
    # ramp: cycling the three validated ATS hues over seven categories is
    # forbidden, and the ATS palette has no fourth categorical slot. Only the
    # TYPOGRAPHY here is ATS. Collapse the groups or supply more brand hues if
    # this should go on-brand -- and re-run the validator if so.
    pal <- grDevices::hcl.colors(length(unique(d$grp)), "Dark 3")
    names(pal) <- sort(unique(d$grp))

    if (kind == "caterpillar") {
      d <- d[order(d$pct_paired), ]
      d$unit_lab <- factor(d$unit_lab, levels = d$unit_lab)
      p <- house_ats(
        ggplot(d, aes(pct_paired, unit_lab, colour = grp)) +
          geom_vline(xintercept = ref, colour = ATS_INK, linetype = "22",
                     linewidth = 0.45) +
          geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y",
                        width = 0, linewidth = 0.5) +
          geom_point(size = 2.1) +
          scale_colour_manual(values = pal, name = NULL) +
          guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
          annotate("text", x = ref, y = Inf, label = sprintf(" pooled %.1f%%", ref),
                   hjust = 0, vjust = 1.4, colour = ATS_INK, size = 3.1) +
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
      p <- house_ats(
        ggplot() +
          geom_ribbon(data = band, aes(n, ymin = lo, ymax = hi, group = band),
                      fill = GRIDLINE, alpha = 0.55) +
          geom_hline(yintercept = ref, colour = ATS_INK, linetype = "22",
                     linewidth = 0.45) +
          geom_point(data = d, aes(n_events, pct_paired, colour = grp),
                     size = 2.2) +
          scale_colour_manual(values = pal, name = NULL) +
          guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
          annotate("text", x = max(d$n_events), y = ref,
                   label = sprintf("pooled %.1f%%", ref), hjust = 1.02,
                   vjust = -0.6, colour = ATS_INK, size = 3.1) +
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
  for (metric in names(res$plots)) {
    spec <- METRIC[[metric]]; cc <- res$caps[[metric]]
    fn <- sprintf("%s_%s_pooled.png", metric, kind)
    ggsave(file.path(FIG, fn), res$plots[[metric]]$p,
           width = 8.4, height = res$plots[[metric]]$h, dpi = 200)
    register_caption(fn,
      sprintf("%s by unit, %d units across %d sites",
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


# ---- Forest plot of the exposure odds ratios ------------------------------
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
    cf$term_label <- exp_label(cf$term)
    het <- cf[cf$is_pool, c("term", "term_label", "k_sites", "Q", "Q_df", "Q_p", "I2_pct",
                            "tau2", "heterogeneity_role")]
    # Caption order must follow PANEL order, or a reader matching one to the
    # other reads the panels top-down and the stats in the opposite sequence.
    het <- het[order(match(as.character(het$term_label), EXP_ORDER)), , drop = FALSE]
    rng <- range(c(cf$or_ci_lo, cf$or_ci_hi), na.rm = TRUE)
    # The forest's own colour key: a site, or the pool, which is not a site.
    cf$hue <- ifelse(cf$is_pool, "Pooled", site_label[cf$site])

    FOREST_X <- sprintf("Odds of a better VFD-28 band, per %d percentage points (log scale)",
                        round(100 * fromJSON(here("config", "protocol.json"),
                                             simplifyVector = FALSE)$model$outcome_model$exposure_scale_per))
    # `hues` NULL draws the conventional monochrome forest, which is what the
    # STANDALONE figure stays. The composite passes SITE_HUE so its single
    # legend at the foot decodes every panel, including this one -- a reader
    # told "Site 1 is blue" and shown an all-ink forest has been misled.
    draw_forest <- function(hues = NULL) {
      p <- ggplot(cf, aes(x = odds_ratio, y = label)) +
        geom_vline(xintercept = 1, colour = MUTED, linewidth = 0.4)
      p <- if (is.null(hues))
        p + geom_errorbar(aes(xmin = or_ci_lo, xmax = or_ci_hi), orientation = "y",
                          width = 0.18, linewidth = 0.5, colour = ATS_INK) +
            geom_point(aes(size = wt_rel, shape = is_pool),
                       colour = ATS_INK, fill = ATS_INK)
      else
        p + geom_errorbar(aes(xmin = or_ci_lo, xmax = or_ci_hi, colour = hue),
                          orientation = "y", width = 0.18, linewidth = 0.5) +
            geom_point(aes(size = wt_rel, shape = is_pool, colour = hue, fill = hue)) +
            # The pooled diamond stays ink: it is the pool, not a fourth site,
            # and giving it a hue would put it in the site legend.
            scale_colour_manual(values = c(hues, Pooled = ATS_INK), guide = "none") +
            scale_fill_manual(values = c(hues, Pooled = ATS_INK), guide = "none")
      p <- p +
        scale_shape_manual(values = c(`FALSE` = 22, `TRUE` = 23), guide = "none") +
        scale_size_area(max_size = 4.2, guide = "none", na.value = 3.6) +
        scale_x_continuous(trans = "log10",
                           breaks = scales::breaks_log(5),
                           labels = function(x) formatC(x, format = "fg", digits = 3)) +
        facet_wrap(~ term_label, ncol = 1, scales = "free_y") +
        labs(x = FOREST_X, y = NULL)
      house_ats(p) +
        theme(panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_line(colour = GRIDLINE, linewidth = 0.4))
    }
    p_f5     <- draw_forest()
    p_f5_hue <- draw_forest(SITE_HUE)
    ggsave(file.path(FIG, "exposure_forest_pooled.png"), p_f5,
           width = 7.0, height = 1.4 + nrow(het) * (0.34 * length(lev) + 0.55),
           dpi = 300,
           bg = PAPER)
    register_caption(
      "exposure_forest_pooled.png",
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
    cat("  exposure_forest_pooled.png\n")
  }
}


# ---- Pooled distribution of the two exposures -----------------------------
# Bars are the POOLED bins -- sum of counts over sum of denominators, never an
# average of the sites' percentages. Site outlines sit behind them so between-site
# spread is visible without a second figure.
ed_f <- "07_outcome_model__exposure_distribution.csv"
if (file.exists(file.path(POOL, ed_f))) {
  ed <- rd(ed_f)
  for (nm in c("bin_lo", "bin_hi", "n", "n_total", "pct")) ed[[nm]] <- as.numeric(ed[[nm]])
  ed$term_label <- exp_label(ed$exposure)
  pooled_ed <- ed[ed$site == "POOLED", , drop = FALSE]
  sites_ed  <- ed[ed$site != "POOLED", , drop = FALSE]
  sites_ed$label <- site_label[sites_ed$site]

  if (nrow(pooled_ed)) {
    p_f6 <- ggplot() +
      geom_rect(data = pooled_ed,
                aes(xmin = bin_lo, xmax = bin_hi, ymin = 0, ymax = pct),
                fill = ATS_FILL, alpha = 0.55,
                colour = NA) +
      geom_step(data = sites_ed, aes(x = bin_lo, y = pct, group = label),
                colour = ATS_INK, linewidth = 0.35, direction = "hv") +
      facet_wrap(~ term_label, ncol = 1) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                         labels = function(x) sprintf("%.0f%%", 100 * x)) +
      labs(x = "Share of an admission's fentanyl increases meeting the condition",
           y = "% of admissions")
    p_f6 <- house_ats(p_f6) + theme(legend.position = "none")
    ggsave(file.path(FIG, "exposure_histogram_pooled.png"), p_f6,
           width = 7.0, height = 5.4, dpi = 300, bg = PAPER)
    register_caption(
      "exposure_histogram_pooled.png",
      "How often is a fentanyl increase paired with a bolus, or justified by a score, across sites?",
      sprintf(paste0(
        "Distribution of the two per-admission exposures over the first 72h, ",
        "each the share of INCREASE events (initiation union up-titration) ",
        "meeting its condition. FILLED BARS ARE THE POOLED DISTRIBUTION -- sum ",
        "of counts over sum of denominators across %d site(s), never an average ",
        "of the sites' percentages. THIN STEPS ARE THE INDIVIDUAL SITES, drawn ",
        "so between-site spread is visible here rather than in a second figure; ",
        "they are unlabelled because the pooled shape is the subject. Bin width ",
        "%.2f on the share scale. Pooled denominators: %s. BOTH DISTRIBUTIONS ",
        "ARE EXPECTED TO BE DOMINATED BY ZERO; the shared denominator is small, ",
        "so a proportion on two or three events can only take a few values and ",
        "the comb is the small-denominator limitation, not a rendering ",
        "artifact. Sites are anonymised by the declared rule."),
        K, unique(pooled_ed$bin_hi - pooled_ed$bin_lo)[1],
        paste(sprintf("%s n = %s", unique(pooled_ed$term_label),
                      formatC(vapply(unique(pooled_ed$term_label), function(t)
                        max(pooled_ed$n_total[pooled_ed$term_label == t]), numeric(1)),
                        format = "d", big.mark = ",")), collapse = "; ")))
    cat("  exposure_histogram_pooled.png\n")
  }
}


# ---- Adjusted predicted VFD-28, sites on one set of axes -------------------
# NOT POOLED, and deliberately so: the curve is standardised over each site's own
# covariate distribution, so three curves are not three estimates of one quantity
# (pooling_spec.json adjusted_prediction._why). Overlaid on shared axes instead
# of one column per site -- the ranges are close enough that a common scale
# reads, and it frees the width for the three-part composite below.
ap_f <- "07_outcome_model__adjusted_prediction.csv"
if (file.exists(file.path(POOL, ap_f))) {
  ap <- rd(ap_f)
  ap <- ap[ap$site != "POOLED", , drop = FALSE]
  for (nm in grep("mean_vfd_28|pct_vfd_zero|exposure_pct", names(ap), value = TRUE))
    ap[[nm]] <- as.numeric(ap[[nm]])
  # Three validated brand hues (figures.R), so the overlay holds to three sites.
  # A fourth must BRING a fourth validated hue -- cycling three over four is the
  # one thing the palette comment forbids -- so past that the figure goes back to
  # one column per site and the composite is skipped rather than drawn wrong.
  AP_OVERLAY <- nrow(ap) > 0 && !is.null(SITE_HUE)
  AP_XLAB <- "Increases with a documented indication, % of an admission's increases"
  if (nrow(ap)) {
    ap$label <- factor(site_label[ap$site], levels = unname(site_label[SITES]))
    ap$label <- droplevels(ap$label)
    long <- rbind(
      data.frame(label = ap$label, x = ap$exposure_pct, y = ap$mean_vfd_28,
                 lo = ap$mean_vfd_28_lo, hi = ap$mean_vfd_28_hi,
                 metric = "Adjusted mean VFD-28, days", stringsAsFactors = FALSE),
      data.frame(label = ap$label, x = ap$exposure_pct, y = ap$pct_vfd_zero,
                 lo = ap$pct_vfd_zero_lo, hi = ap$pct_vfd_zero_hi,
                 metric = "Adjusted % with VFD-28 = 0", stringsAsFactors = FALSE))
    long$metric <- factor(long$metric,
                          levels = c("Adjusted mean VFD-28, days",
                                     "Adjusted % with VFD-28 = 0"))
    ap_n <- vapply(levels(ap$label), function(l)
      max(as.numeric(ap$n_admissions[ap$label == l])), numeric(1))
    # n on the KEY, the way the forest puts it on the row label: a reader cannot
    # take a denominator off a line, and there is no per-site strip left to hold it.
    SITE_KEY <- sprintf("%s (n = %s)", levels(ap$label),
                        formatC(ap_n, format = "d", big.mark = ","))

    # One metric, every site on it. The x title is dropped here and drawn ONCE
    # by compose_panels(): side by side, each panel is narrower than that title.
    pred_panel <- function(m) {
      d <- long[long$metric == m, , drop = FALSE]
      p <- ggplot(d, aes(x = x, colour = label, fill = label)) +
        geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
        geom_line(aes(y = y), linewidth = 0.7) +
        geom_point(aes(y = y), size = 1.7) +
        scale_colour_manual(values = SITE_HUE, labels = SITE_KEY, name = NULL) +
        scale_fill_manual(values = SITE_HUE, guide = "none") +
        scale_x_continuous(breaks = sort(unique(long$x))) +
        labs(x = NULL, y = m)
      house_ats(p) +
        # The y title would otherwise sit on the tick labels, and the stripped
        # legend row leaves its box spacing behind as a gap above the shared
        # x title. Both are only visible once the panels are composed.
        theme(axis.title.y = element_text(margin = margin(r = 6)),
              legend.box.spacing = unit(0, "pt"),
              plot.margin = margin(t = 4, r = 8, b = 0, l = 4))
    }
    # Tagged at COMPOSITION, not here: the same two panels are A/B standalone and
    # B/C under the forest, so the tag cannot live in the panel that is reused.
    tag_as <- function(p, tag)
      p + labs(tag = tag) +
        theme(plot.tag = element_text(family = ATS_FONT, colour = ATS_INK,
                                      face = "bold", size = 13),
              plot.tag.position = "topleft")

    p_mean <- pred_panel("Adjusted mean VFD-28, days")
    p_zero <- pred_panel("Adjusted % with VFD-28 = 0")
  }

  if (AP_OVERLAY) {
    g_f7 <- compose_panels(
      bottom = list(tag_as(p_mean, "A"), tag_as(p_zero, "B")),
      heights = 3.3, x_title = AP_XLAB,
      x_title_gp = grid::gpar(col = ATS_INK, fontfamily = ATS_FONT, fontsize = 12))
    ggsave(file.path(FIG, "adjusted_prediction_by_site.png"), g_f7,
           width = 7.2, height = 4.3, dpi = 300, bg = PAPER)
  } else if (nrow(ap)) {
    # Fallback past three sites: facets, which is what this figure was before a
    # shared-axis overlay became legible at k = 3.
    p_f7 <- ggplot(long, aes(x = x)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = ATS_FILL, alpha = 0.16) +
      geom_line(aes(y = y), colour = ATS_SERIES[1], linewidth = 0.7) +
      geom_point(aes(y = y), colour = ATS_SERIES[1], size = 1.7) +
      facet_grid(metric ~ label, scales = "free_y", switch = "y") +
      scale_x_continuous(breaks = sort(unique(long$x))) +
      labs(x = AP_XLAB, y = NULL)
    p_f7 <- house_ats(p_f7) +
      theme(legend.position = "none", strip.placement = "outside",
            panel.spacing.x = unit(0.9, "lines"))
    ggsave(file.path(FIG, "adjusted_prediction_by_site.png"), p_f7,
           width = max(6.6, 2.4 + 2.1 * nlevels(ap$label)), height = 5.6,
           dpi = 300, bg = PAPER)
  }

  if (nrow(ap)) {
    register_caption(
      "adjusted_prediction_by_site.png",
      "Does documenting an indication for a fentanyl increase track ventilator-free days, at each site?",
      sprintf(paste0(
        "%d site(s), %s. A: adjusted mean VFD-28 in days. B: the adjusted share ",
        "of admissions with VFD-28 exactly zero. Each curve is that site's ",
        "PRIMARY proportional-odds model STANDARDISED over ITS OWN observed ",
        "covariate distribution, with the band-to-days mapping it measured ",
        "locally. DELIBERATELY NOT POOLED: standardising over different ",
        "covariate distributions means %d curves are not %d estimates of one ",
        "quantity, so combining them would average over %d different ",
        "populations while looking like a meta-analysis. SHARED AXES ARE A ",
        "VISUAL CONVENIENCE, NOT A POOLING CLAIM -- the vertical offset between ",
        "sites is a difference in case mix and in local band-to-days mapping, ",
        "and is NOT a between-site comparison. What DOES pool is the exposure ",
        "coefficient, by inverse variance, in exposure_forest_pooled.png. A ",
        "pooled curve would have to be drawn FROM that pooled coefficient ",
        "against a chosen reference covariate distribution, which is a protocol ",
        "decision and is not done here. THE TWO PANELS HAVE DIFFERENT UNITS AND ",
        "SCALES and are not comparable to each other. Bands are bootstrap ",
        "percentile intervals at each site's declared seed; they overlap, and ",
        "overlap between sites is not a test. Sites are anonymised by the ",
        "declared rule. AN ASSOCIATION, NOT AN EFFECT -- confounding by ",
        "indication is unaddressed."),
        nlevels(ap$label),
        if (AP_OVERLAY) "colour-coded on shared axes" else
          "one column per site -- the brand palette has no further validated hue, so the overlay is not drawn",
        nlevels(ap$label), nlevels(ap$label), nlevels(ap$label)))
    cat("  adjusted_prediction_by_site.png\n")
  }
}


# ---- Three-part composite: the forest over the two adjusted-prediction panels --
# The manuscript figure. B and C are the SAME plot objects the standalone
# prediction figure saved, and A comes from the same draw_forest() the standalone
# forest used, so a change to either cannot leave the composite showing the
# previous version of itself. A differs only in taking the site hues.
if (exists("p_f5_hue") && exists("p_mean") && isTRUE(AP_OVERLAY)) {
  # A is sized by the forest's OWN row formula, so adding a site grows the band
  # rather than squeezing it; the pair below keeps a fixed height, and the last
  # inch is the shared x title and the legend, which size themselves.
  H_A  <- 0.45 + nrow(het) * (0.34 * length(lev) + 0.55)
  H_BC <- 3.1
  g_abc <- compose_panels(
    top    = tag_as(p_f5_hue, "A"),
    bottom = list(tag_as(p_mean, "B"), tag_as(p_zero, "C")),
    heights = c(H_A, H_BC), x_title = AP_XLAB,
    x_title_gp = grid::gpar(col = ATS_INK, fontfamily = ATS_FONT, fontsize = 12))
  ggsave(file.path(FIG, "outcome_model_composite.png"), g_abc,
         width = 7.2, height = H_A + H_BC + 1.0, dpi = 300, bg = PAPER)
  register_caption(
    "outcome_model_composite.png",
    "Does documenting an indication for a fentanyl increase track ventilator-free days?",
    paste0(
      "A: per-admission titration exposures against VFD-28, one sub-panel per ",
      "exposure, from the PRIMARY proportional-odds model at each site. ",
      "SQUARES ARE SITES, area proportional to inverse-variance weight; the ",
      "DIAMOND is the pooled estimate, POOLED BY INVERSE VARIANCE ON THE ",
      "LOG-ODDS. An odds ratio above 1 is higher odds of landing in a BETTER ",
      "VFD-28 band, read against the band ordering and not against intuition. ",
      "B and C: adjusted mean VFD-28 and the adjusted share with VFD-28 exactly ",
      "zero, against the documented-indication exposure, each site standardised ",
      "over ITS OWN covariate distribution and so DELIBERATELY NOT POOLED -- ",
      "the shared axes are a visual convenience and the vertical offset between ",
      "sites is case mix, not a between-site contrast. The pooled estimand in A ",
      "is the WITHIN-SITE, ICU-adjusted association: unit fixed effects are ",
      "site-local and are not shipped. Sites carry the same colour in B and C ",
      "and in A, and are named once in the legend at the foot, which carries ",
      "each site's n; A repeats n on its row labels and draws the POOLED ",
      "estimate in ink, because the pool is not a site. ",
      "Full per-panel notes, including the heterogeneity statistics for A and ",
      "the bootstrap seeds for B and C, are in the captions for ",
      "exposure_forest_pooled.png and adjusted_prediction_by_site.png, which ",
      "draw the same numbers. Sites are anonymised by the declared rule. AN ",
      "ASSOCIATION, NOT AN EFFECT -- confounding by indication is unaddressed."))
  cat("  outcome_model_composite.png\n")
}


# ---- Per-unit exposure caterpillars, mirroring the unit-variation pattern -------------------------
# The site-side figure draws violins from one row per admission; that stays
# PHI-side, so exposure_by_unit.csv ships n / n_at_zero / median / q1 / q3 / mean
# / sd instead. A caterpillar is what those columns support.
#
# THE POINT IS THE MEAN, not the median the unit-variation caterpillars would
# half of each exposure is exactly zero, so most units have median = q1 = 0 and a
# caterpillar of medians is a column on the floor. The mean also POOLS EXACTLY
# from n and mean, so the reference line is the same quantity as the points.
eu_f <- "07_outcome_model__exposure_by_unit.csv"
if (file.exists(file.path(POOL, eu_f))) {
  eu <- rd(eu_f)
  eu <- eu[eu$site != "POOLED", , drop = FALSE]
  for (nm in c("n_admissions", "n_at_zero", "median", "q1", "q3", "mean", "sd"))
    eu[[nm]] <- suppressWarnings(as.numeric(eu[[nm]]))

  for (ex in sort(unique(eu$exposure))) {
    d <- eu[eu$exposure == ex, , drop = FALSE]
    n_all <- nrow(d)
    # A suppressed unit cannot be drawn; it is counted in the caption instead.
    d <- d[is.finite(d$mean) & is.finite(d$n_admissions), , drop = FALSE]
    if (!nrow(d)) next
    n_suppressed <- n_all - nrow(d)

    d$grp <- tolower(ifelse(is.na(d$location_category) |
                              d$location_category == "", "unattributed",
                            d$location_category))
    # Numbered alphabetically by real name within site, so no label encodes a rate.
    d <- d[order(d$site, d$location_name), ]
    d$unit_lab <- sprintf("S%d-U%d (n = %s)",
                          match(d$site, SITES),
                          ave(seq_len(nrow(d)), d$site, FUN = seq_along),
                          formatC(d$n_admissions, format = "d", big.mark = ","))
    # Pooled mean over the units DRAWN, sharing their denominator.
    ref <- 100 * sum(d$mean * d$n_admissions) / sum(d$n_admissions)
    zero_pct <- 100 * d$n_at_zero / d$n_admissions

    # A 95% INTERVAL FOR THE MEAN, from the shipped sd and n -- not the IQR.
    # The point is the mean (the median is 0 at most units), and an interquartile
    # range brackets the MEDIAN: pairing the two let the bar miss its own point
    # wherever q3 < mean -- 8 of 44 units on Pbolus -- and collapse to an
    # invisible zero-length bar wherever q1 = q3 = 0.
    # CLAMPED to [0, 1]: these are means of per-admission SHARES, so the normal
    # approximation can put a bound past the limit at small n, and did.
    d$se <- d$sd / sqrt(d$n_admissions)
    d$ci_lo <- pmax(0, d$mean - 1.959964 * d$se)
    d$ci_hi <- pmin(1, d$mean + 1.959964 * d$se)
    stopifnot("a unit's interval must contain its own point estimate" =
                all(is.na(d$se) | (d$ci_lo <= d$mean & d$mean <= d$ci_hi)))
    n_clamped <- sum(!is.na(d$se) & (d$mean - 1.959964 * d$se < 0 |
                                       d$mean + 1.959964 * d$se > 1))

    d <- d[order(d$mean), ]
    d$unit_lab <- factor(d$unit_lab, levels = d$unit_lab)
    # Three groups here (ed / icu / unattributed), so the validated ATS trio
    # fits without cycling. Assigned alphabetically so a group keeps its colour
    # when a site with a different mix joins -- colour follows the entity, never
    # its rank or the row count.
    grps <- sort(unique(d$grp))
    stopifnot("more care categories than the validated palette has hues" =
                length(grps) <= length(ATS_SERIES))
    pal <- stats::setNames(ATS_SERIES[seq_along(grps)], grps)

    fn <- sprintf("exposure_by_unit_%s.png", tolower(ex))
    p <- house_ats(
      ggplot(d, aes(100 * mean, unit_lab, colour = grp)) +
        geom_vline(xintercept = ref, colour = ATS_INK, linetype = "22",
                   linewidth = 0.45) +
        geom_errorbar(aes(xmin = 100 * ci_lo, xmax = 100 * ci_hi),
                      orientation = "y", width = 0, linewidth = 0.5) +
        geom_point(size = 2.1) +
        scale_colour_manual(values = pal, name = NULL) +
        guides(colour = guide_legend(nrow = 1)) +
        annotate("text", x = ref, y = Inf,
                 label = sprintf(" pooled mean %.1f%%", ref),
                 hjust = 0, vjust = 1.4, colour = ATS_INK, size = 3.1) +
        scale_x_continuous(expand = expansion(mult = c(0.06, 0.10)),
                           labels = function(x) sprintf("%.0f%%", x)) +
        labs(x = sprintf("%s, %% of an admission's fentanyl increases",
                         if (ex == "Pindicated")
                           "Increases with a documented indication" else
                           "Increases paired with a fentanyl bolus"),
             y = NULL))
    ggsave(file.path(FIG, fn), p, width = 7.2,
           height = max(3.4, 1.5 + 0.34 * nrow(d)), dpi = 300, bg = PAPER)
    register_caption(
      fn,
      sprintf("Does %s differ between ICUs, across sites?",
              if (ex == "Pindicated") "documented indication" else
                "bolus pairing"),
      sprintf(paste0(
        "One row per ICU, %d unit(s) across %d site(s), ordered by the unit's ",
        "MEAN exposure. POINTS ARE MEANS, not medians: about half of this ",
        "exposure is exactly zero, so most units have median = Q1 = 0 and a ",
        "caterpillar of medians is a column on the floor. The mean also pools ",
        "exactly from n and mean, so the dashed reference line -- the pooled ",
        "mean over the units drawn, %.1f%% -- is the same quantity as the ",
        "points. BARS ARE A 95%% INTERVAL FOR THAT MEAN, from the unit's shipped ",
        "standard deviation and n, so a bar always contains its own point. It is ",
        "PRECISION, not spread. The interquartile range still ships in the CSV ",
        "and is deliberately NOT drawn: an IQR brackets the MEDIAN, so it misses ",
        "the mean wherever the third quartile falls below it -- which happens ",
        "wherever most admissions pair no bolus at all. %s Units are coloured by the ",
        "mCIDE location_category, which means the same thing at every site ",
        "where location_name is site-local free text. THE SHARE AT EXACTLY ZERO ",
        "is what a violin would have shown and is not drawable from a summary; ",
        "across the units drawn it runs %.0f%% to %.0f%% and ships per unit as ",
        "n_at_zero in 07_outcome_model__exposure_by_unit.csv. %s Labels are ",
        "S<site>-U<unit>, numbered alphabetically by real name within site so no ",
        "label encodes a rate; n is the unit's admissions. NOT POOLED ACROSS ",
        "SITES -- location_name is site-local, so every row is one site's unit."),
        nrow(d), length(unique(d$site)), ref,
        if (n_clamped)
          sprintf(paste0("%d interval(s) are CLAMPED to the [0, 100%%] bound, ",
                         "because a normal approximation on a mean of shares ",
                         "can reach past it at small n."), n_clamped) else
          "No interval needed clamping to the [0, 100%] bound.",
        min(zero_pct, na.rm = TRUE), max(zero_pct, na.rm = TRUE),
        if (n_suppressed)
          sprintf("%d unit(s) are absent: below the site's small-cell floor, so their statistics were blanked at the site and only their counts ship.", n_suppressed)
        else "Every unit shipped drawable statistics."))
    cat(sprintf("  %s\n", fn))
  }
}

# ---- Pooled liberation and survival to the horizon ---------------------------
# The SAME liberation_panels() the site script draws, on the POOLED row of the
# same CSV, so the pooled figure and the site figures cannot drift apart. The
# measures are counts over a fixed grid, so the pooled stack is still a
# partition -- summing them is the whole reason this file pools as sum_counts.
ls_f <- "02_descriptive__liberation_survival.csv"
if (file.exists(file.path(POOL, ls_f))) {
  ls <- rd(ls_f)
  ls <- ls[ls$site == "POOLED", , drop = FALSE]
  for (nm in c("hour", "day", "n_admissions", "denominator"))
    ls[[nm]] <- as.numeric(ls[[nm]])

  if (nrow(ls)) {
    denom <- max(ls$denominator)
    stopifnot("a pooled grid point reports a different denominator" =
                all(ls$denominator == denom))
    horizon_d <- max(ls$day)
    # The landmark is the protocol's, not read off the curve: a pooled file
    # whose first liberation happened to fall elsewhere would move the rule.
    COVP <- fromJSON(here("config", "protocol.json"), simplifyVector = FALSE)
    p_ls <- liberation_panels(ls, denom = denom,
                              landmark_d = COVP$cohort$landmark_hours / 24,
                              horizon_d = horizon_d, theme_fn = house_ats)
    fn <- "liberation_survival.png"
    ggsave(file.path(FIG, fn), p_ls, width = 7.2, height = 8.4, dpi = 300,
           bg = PAPER)
    dead_end <- ls$n_admissions[ls$measure == "died" & ls$hour == max(ls$hour)]
    lost_end <- ls$n_admissions[ls$measure == "discharged alive, charted on IMV" &
                                  ls$hour == max(ls$hour)]
    register_caption(fn,
      "Liberation from invasive ventilation and survival through the horizon, pooled",
      sprintf(paste0(
        "All %s admissions ventilated at the landmark across %d site(s), ",
        "followed to day %g from the start of invasive ventilation; the ",
        "dashed rule marks the landmark. Counts are summed across sites over ",
        "an identical grid, never averaged over site percentages. UPPER ",
        "PANEL, after MODE (Semler et al., NEJM Evid 2024) Figure 3: the share ",
        "alive, and the cumulative incidence of a first successful ",
        "liberation. LOWER PANEL: where the cohort is, hour by hour. At day ",
        "%g, %s (%.1f%%) had died and %s (%.1f%%) were discharged alive while ",
        "respiratory_support was still charting invasive ventilation -- both ",
        "score VFD-28 = 0, for different reasons. THE ORANGE BAND IS PARTLY A ",
        "CHARTING ARTEFACT and its size differs between sites with how ",
        "completely each populates device_category; read it against each ",
        "site's 01_cohort/diagnostics/waterfall_carryforward.csv before ",
        "reading a between-site difference as clinical. Mortality is ",
        "in-hospital and so a lower bound."),
        format(denom, big.mark = ","), K, horizon_d, horizon_d,
        format(dead_end, big.mark = ","), 100 * dead_end / denom,
        format(lost_end, big.mark = ","), 100 * lost_end / denom))
    cat(sprintf("  %s\n", fn))
  }
}

# FIG is cleared at the start and everything in it is drawn above, so no name
# pattern is needed -- and none can go stale the next time a figure is renamed.
MINE <- list.files(FIG, pattern = "\\.png$")
write_captions(file.path(FIG, "captions.md"), "pool_figures.R", MINE,
               list(site_name = sprintf("%d pooled", K),   # write_captions prefixes "Site "
                    clif_version = PROV$definition_versions$`protocol.json`,
                    code_version = PROV$code_versions[[1]],
                    generated = format(Sys.time(), "%Y-%m-%d %H:%M %Z")))
cat(sprintf("\n%d figure(s) in output/pooled/figures/\n",
            length(list.files(FIG, pattern = "\\.png$"))))
