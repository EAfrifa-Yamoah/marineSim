# Figure 4: the temporal block resolved by timesteps, sample size and design
#   (a) temporal Shapley contribution against T, one line per sample size n (results/expanded_T6_temporal_by_T_n.csv)
#   (b) share of the total gain carried by the temporal block against n, one line per T >= 3 (same file)
#   (c) temporal Shapley against the median distance from a held out cell to its nearest site, by design, at T = 5
#       (results/temporal_coverage_T5.csv for the datasets, results/expanded_T8_temporal_by_coverage.csv for the binned means)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figstyle.R"))
t6 <- read.csv(res("expanded_T6_temporal_by_T_n.csv"))
cov <- read.csv(res("temporal_coverage_T5.csv")); t8 <- read.csv(res("expanded_T8_temporal_by_coverage.csv"))
ns <- sort(unique(t6$n)); Ts <- c(3, 5, 10)
RAMP <- c("30" = "#F2A66F", "50" = "#E8833F", "100" = "#D55E00", "200" = "#A84A00", "500" = "#6E3000")   # one hue, light to dark
pdg <- position_dodge(width = 0.9)      # deterministic horizontal offsets so that intervals do not overlap
pdl <- position_dodge(width = 0.16)     # the same on the log10 n axis

a <- t6[t6$quantity == "Shapley", ]; a$n <- factor(a$n, levels = ns)
pa <- ggplot(a, aes(n_time, value, colour = n, group = n)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.25, linewidth = 0.35, position = pdg) +
  geom_line(linewidth = 0.7, position = pdg) + geom_point(size = 1.8, position = pdg) +
  annotate("text", x = 1.6, y = 0.004, label = "placebo", size = 2.6, colour = "grey35", hjust = 0) +
  scale_colour_manual(values = RAMP, name = "n =") + scale_x_continuous(breaks = c(1, 3, 5, 10)) +
  labs(x = "Number of timesteps", y = "Temporal Shapley contribution (ΔAUC)") +
  guides(colour = guide_legend(nrow = 1, title.position = "left", title.vjust = 0.5))

b <- t6[t6$quantity == "Share of total gain (%)", ]; b$n_time <- factor(b$n_time, levels = Ts, labels = paste0("T = ", Ts))
pb <- ggplot(b, aes(n, value, linetype = n_time, shape = n_time)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.04, linewidth = 0.35, colour = COMP_COL[["Temporal"]], linetype = "solid", position = pdl) +
  geom_line(linewidth = 0.7, colour = COMP_COL[["Temporal"]], position = pdl) + geom_point(size = 1.9, colour = COMP_COL[["Temporal"]], position = pdl) +
  scale_x_log10(breaks = ns) + scale_linetype_manual(values = c("solid", "22", "11")) + scale_shape_manual(values = c(16, 15, 17)) +
  labs(x = "Sites per timestep", y = "Temporal share of total gain (%)", linetype = NULL, shape = NULL)

DES <- c("random", "stratified", "clustered"); cov$design <- factor(cov$design, levels = DES); t8$design <- factor(t8$design, levels = DES)
pc <- ggplot() +
  geom_point(data = cov, aes(nn_median_km, phi_TEMP, shape = design), colour = COMP_COL[["Temporal"]], alpha = 0.16, size = 1) +
  geom_errorbar(data = t8, aes(x_km, ymin = ci_lo, ymax = ci_hi, group = design), width = 0.03, linewidth = 0.35, colour = COMP_COL[["Temporal"]]) +
  geom_line(data = t8, aes(x_km, value, linetype = design), colour = COMP_COL[["Temporal"]], linewidth = 0.7) +
  geom_point(data = t8, aes(x_km, value, shape = design), colour = COMP_COL[["Temporal"]], size = 1.9) +
  scale_x_log10(breaks = c(5, 10, 20, 40, 80, 160)) + scale_linetype_manual(values = c("11", "22", "solid")) + scale_shape_manual(values = c(17, 15, 16)) +
  coord_cartesian(ylim = c(-0.01, 0.1)) +
  labs(x = "Median nearest site distance (km)", y = "Temporal Shapley contribution (\u0394AUC)", linetype = NULL, shape = NULL)

fig4 <- (pa + pb + pc) + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 7.5))
save_png(fig4, "main/Figure4.png", 11, 4.3)
