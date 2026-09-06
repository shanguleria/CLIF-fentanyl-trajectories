# ==============================================================================
# 04_composition_bias_demo.R
#
# Purpose : Illustrate why a whole-cohort dosing curve needs a balanced-panel
#           overlay (design notes section 10, Phase 1). Simulates a cohort in
#           which NO patient's dose changes over time, yet the cohort mean rises
#           steeply -- purely because higher-dose patients stay ventilated longer.
#           Also demonstrates NESTED balanced panels (>=24h, >=48h, >=72h), whose
#           vertical separation displays the selection gradient directly.
# Author  : Shan Guleria
# Created : 2026-09-04
# Inputs  : none (simulated)
# Outputs : output/final_no_phi/validation/composition_bias_demo.png
# ==============================================================================

pkgs <- c("here", "ggplot2")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}
set.seed(20260904)

# ---- Simulate ----------------------------------------------------------------
# Each patient has a CONSTANT dose. Sicker patients get more fentanyl AND stay
# ventilated longer -- the only link between dose and time.

n <- 400
baseline <- pmax(0, rnorm(n, 1.2, 0.6))                    # mcg/kg/hr, constant
vent_hours <- pmin(72, 6 + rexp(n, rate = 1 / (10 * baseline + 5)))

hours <- seq(3, 72, by = 3)
long <- do.call(rbind, lapply(seq_len(n), function(i) {
  h <- hours[hours <= vent_hours[i]]
  if (!length(h)) return(NULL)
  data.frame(id = i, hour = h, dose = baseline[i])          # dose NEVER changes
}))

# ---- Curves ------------------------------------------------------------------

all_comers <- aggregate(dose ~ hour, long, mean)
all_comers$n <- as.integer(table(long$hour))
all_comers$curve <- "All at-risk patients"

# Nested balanced panels. Each panel is plotted ONLY over the window it is
# balanced for -- a >=24h panel extended past hour 24 would no longer be
# balanced, which is the whole point of the construction.
thresholds <- c(24, 48, 72)
panels <- do.call(rbind, lapply(thresholds, function(k) {
  ids <- which(vent_hours >= k)
  sub <- long[long$id %in% ids & long$hour <= k, ]
  out <- aggregate(dose ~ hour, sub, mean)
  out$n <- length(ids)
  out$curve <- sprintf("Balanced panel, ventilated >=%dh (n=%d)", k, length(ids))
  out
}))

plot_df <- rbind(all_comers[, c("hour","dose","n","curve")], panels)

cat("Truth: every patient dose is CONSTANT. No within-patient change exists.\n\n")
cmp <- merge(all_comers[, c("hour","dose","n")],
             panels[panels$curve == unique(panels$curve)[3], c("hour","dose")],
             by = "hour", suffixes = c("_all", "_bal72"))
print(round(cmp[cmp$hour %in% c(3, 24, 48, 72), ], 2), row.names = FALSE)

cat(sprintf("\nAll-comers mean rises %.0f%% from h3 to h72.\n",
            100 * (tail(all_comers$dose, 1) / all_comers$dose[1] - 1)))

cat("\nNested panels -- each FLAT within itself; the gradient is between them:\n")
for (k in thresholds) {
  pk <- panels[grepl(sprintf(">=%dh", k), panels$curve), ]
  cat(sprintf("  >=%2dh  n=%3d  level %.2f mcg/kg/hr  (within-panel change %+.1f%%)\n",
              k, pk$n[1], mean(pk$dose),
              100 * (tail(pk$dose,1) / pk$dose[1] - 1)))
}

# ---- Plot --------------------------------------------------------------------
# The three panels are ordinal (a duration threshold), so they get a SEQUENTIAL
# single-hue ramp, light -> dark. All-comers is a different kind of thing, so it
# gets a contrasting categorical hue.

ink <- "#0b0b0b"; muted <- "#898781"; gridline <- "#e1e0d9"
lv <- c("All at-risk patients", unique(panels$curve))
plot_df$curve <- factor(plot_df$curve, levels = lv)
pal <- setNames(c("#eb6834", "#a8c8ee", "#4a8bd8", "#14427e"), lv)

ends <- do.call(rbind, lapply(split(plot_df, plot_df$curve),
                              function(g) g[which.max(g$hour), ]))
ends$lab <- c("all at-risk", ">=24h", ">=48h", ">=72h")
# all-comers and the >=72h panel both terminate at hour 72 at the same value,
# so their end labels collide -- push them apart vertically.
ends$dy <- c(-0.035, 0, 0, 0.035)

p <- ggplot(plot_df, aes(hour, dose, colour = curve)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  geom_text(data = ends, aes(label = lab, y = dose + dy), hjust = -0.15,
            size = 3.2, colour = ink, show.legend = FALSE) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(breaks = seq(0, 72, 12), limits = c(0, 82)) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(title    = "Nested balanced panels separate composition from real change",
       subtitle = "Simulated: doses are constant; higher-dose patients stay ventilated longer",
       x = "Hours since intubation", y = "Mean dose (mcg/kg/hr)", colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", colour = ink),
        plot.subtitle = element_text(colour = muted, margin = margin(b = 10)),
        legend.position = "top", legend.text = element_text(colour = muted, size = 9),
        axis.title = element_text(colour = muted),
        axis.text = element_text(colour = muted),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = gridline, linewidth = 0.4),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA))

out <- here("output", "final_no_phi", "validation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out, "composition_bias_demo.png"), p, width = 7.5, height = 4.8, dpi = 200)
cat("\nWrote:", file.path(out, "composition_bias_demo.png"), "\n")
