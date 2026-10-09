# Figure 4: the temporal block resolved by timesteps, sample size and spatial coverage
#   (a) temporal Shapley contribution against T, one line per sample size n (results/expanded_T6_temporal_by_T_n.csv)
#   (b) share of the total gain carried by the temporal block against n, one line per T >= 3 (same file)
#   (c) temporal Shapley against the median distance from a held out cell to its nearest site, by design, at T = 5
#       (results/temporal_coverage_T5.csv for the datasets, results/expanded_T8_temporal_by_coverage.csv for the binned means)
# Points and intervals that share an x position are offset horizontally (deterministic dodge); the dataset level
# points in (c) are jittered with a fixed seed so that the figure is reproducible.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figstyle.R"))
t6 <- read.csv(res("expanded_T6_temporal_by_T_n.csv"))
cov <- read.csv(res("temporal_coverage_T5.csv")); t8 <- read.csv(res("expanded_T8_temporal_by_coverage.csv"))
ns <- sort(unique(t6$n)); Ts <- c(3, 5, 10); TEMP <- COMP_COL[["Temporal"]]
RAMP <- c("30" = "#F2A66F", "50" = "#E8833F", "100" = "#D55E00", "200" = "#A84A00", "500" = "#6E3000")   # one hue, light to dark
grid_y <- theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))
ylab_shap <- "Temporal Shapley contribution (ΔAUC)"

# ---- (a) by T and n
a <- t6[t6$quantity == "Shapley", ]; a$n <- factor(a$n, levels = ns)
pda <- position_dodge(width = 0.9)
pa <- ggplot(a, aes(n_time, value, colour = n, group = n)) + grid_y +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0, linewidth = 0.45, position = pda, alpha = 0.9) +
  geom_line(linewidth = 0.6, position = pda) +
  geom_point(size = 2.1, position = pda, shape = 21, fill = "white", stroke = 1) +
  annotate("text", x = 1.6, y = 0.0045, label = "placebo", size = 2.7, colour = "grey35", hjust = 0) +
  scale_colour_manual(values = RAMP, name = "n =") +
  scale_x_continuous(breaks = c(1, 3, 5, 10)) + scale_y_continuous(breaks = seq(0, 0.06, 0.01)) +
  labs(x = "Number of timesteps", y = ylab_shap) +
  guides(colour = guide_legend(nrow = 1, title.position = "left", title.vjust = 0.5, override.aes = list(linewidth = 1.2)))

# ---- (b) share of the total gain by n and T
b <- t6[t6$quantity == "Share of total gain (%)", ]; b$n_time <- factor(b$n_time, levels = Ts, labels = paste0("T = ", Ts))
pdl <- position_dodge(width = 0.16)
pb <- ggplot(b, aes(n, value, shape = n_time, linetype = n_time)) + grid_y +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0, linewidth = 0.45, colour = TEMP, linetype = "solid", position = pdl, alpha = 0.9) +
  geom_line(linewidth = 0.6, colour = TEMP, position = pdl) +
  geom_point(size = 2.1, colour = TEMP, fill = "white", stroke = 1, position = pdl) +
  scale_x_log10(breaks = ns) + scale_y_continuous(limits = c(10, 42), breaks = seq(10, 40, 10)) +
  scale_shape_manual(values = c(21, 22, 24)) + scale_linetype_manual(values = c("solid", "22", "11")) +
  labs(x = "Sites per timestep", y = "Temporal share of total gain (%)", shape = NULL, linetype = NULL)

# ---- (c) coverage curve: binned means at the bin centres, dodged by design; dataset points jittered
DES <- c("random", "stratified", "clustered"); cov$design <- factor(cov$design, levels = DES); t8$design <- factor(t8$design, levels = DES)
edges <- as.numeric(sub("\\((.*),.*", "\\1", t8$bin)); uppers <- as.numeric(sub(".*,(.*)\\]", "\\1", t8$bin))
t8$x_mid <- sqrt(edges * uppers)                       # geometric bin centre on the log axis
pdc <- position_dodge(width = 0.09)
set.seed(20261009)
pc <- ggplot() + grid_y +
  geom_point(data = cov, aes(nn_median_km, phi_TEMP, shape = design), colour = TEMP, alpha = 0.14, size = 0.9,
             position = position_jitter(width = 0.025, height = 0, seed = 1)) +
  geom_errorbar(data = t8, aes(x_mid, ymin = ci_lo, ymax = ci_hi, group = design), width = 0, linewidth = 0.45, colour = TEMP, position = pdc) +
  geom_line(data = t8, aes(x_mid, value, linetype = design, group = design), colour = TEMP, linewidth = 0.6, position = pdc) +
  geom_point(data = t8, aes(x_mid, value, shape = design), colour = TEMP, fill = "white", stroke = 1, size = 2.1, position = pdc) +
  scale_x_log10(breaks = c(5, 10, 20, 40, 80, 160)) + coord_cartesian(ylim = c(-0.01, 0.1)) +
  scale_linetype_manual(values = c("11", "22", "solid")) + scale_shape_manual(values = c(24, 22, 21)) +
  labs(x = "Median nearest site distance (km)", y = ylab_shap, linetype = NULL, shape = NULL)

fig4 <- (pa + pb + pc) + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 7.5),
        legend.key.width = unit(16, "pt"))
save_png(fig4, "main/Figure4.png", 11.5, 4.4)
