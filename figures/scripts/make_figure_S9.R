# Figure S9: the pooled temporal contribution by T under three estimators (exact Shapley; increment when the
# block is added last; the same increment with the same timestep neighbour control present). Formerly main text Figure 3.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figstyle.R"))
t2 <- read.csv(res("expanded_T2_temporal_by_T.csv"))
t2$quantity <- factor(t2$quantity, levels = c("Shapley", "Added last", "Added last, NN0 controlled"))
figS9 <- ggplot(t2, aes(n_time, value, shape = quantity, linetype = quantity)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.25, linewidth = 0.4, colour = COMP_COL[["Temporal"]], linetype = "solid") +
  geom_line(colour = COMP_COL[["Temporal"]], linewidth = 0.7) +
  geom_point(colour = COMP_COL[["Temporal"]], fill = COMP_COL[["Temporal"]], size = 2.2) +
  annotate("text", x = 1.2, y = 0.0055, label = "placebo", size = 2.6, colour = "grey35", hjust = 0) +
  scale_shape_manual(values = c(16, 15, 17)) + scale_linetype_manual(values = c("solid", "22", "11")) +
  scale_x_continuous(breaks = c(1, 3, 5, 10)) +
  labs(x = "Number of timesteps", y = "Temporal contribution (ΔAUC)", shape = NULL, linetype = NULL)
save_png(figS9, "si/FigureS9.png", 6.4, 3.9)

