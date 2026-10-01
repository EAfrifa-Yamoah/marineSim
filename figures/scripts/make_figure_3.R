# Figure 3: the five method core benchmark from results/core_five_methods.csv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figstyle.R"))
c0 <- read.csv(res("core_five_methods.csv")); c0$field <- paste(c0$world_cfg, c0$realisation)
M <- c(GLM = "GLM", GAM = "GAM", RF = "Standard RF", SRF = "Spatial RF", STRF = "stRF"); ns <- sort(unique(c0$n))
long <- do.call(rbind, lapply(names(M), function(k) data.frame(field = c0$field, design = c0$design, n = c0$n, method = M[[k]], auc = c0[[k]])))
long$method <- factor(long$method, levels = M)

# (a) pooled over designs: mean AUC and 2.5 to 97.5 percentile band of field realisation means
fm <- aggregate(auc ~ field + n + method, long, mean)
a <- do.call(rbind, lapply(split(fm, list(fm$n, fm$method)), function(s)
  data.frame(n = s$n[1], method = s$method[1], mean = mean(s$auc), lo = quantile(s$auc, 0.025), hi = quantile(s$auc, 0.975))))
pa <- ggplot(a, aes(n, mean, colour = method, fill = method)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_hline(yintercept = 0.5, linetype = "13", linewidth = 0.4) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.9) +
  scale_x_log10(breaks = ns) + scale_colour_manual(values = METHOD_COL) + scale_fill_manual(values = METHOD_COL) +
  labs(x = "Sites per timestep", y = "Mean test AUC", colour = NULL, fill = NULL)

# (b) change in mean AUC relative to random sampling, by method
dm <- aggregate(auc ~ design + n + method, long, mean)
b <- merge(dm[dm$design != "random", ], dm[dm$design == "random", c("n", "method", "auc")], by = c("n", "method"), suffixes = c("", "_random"))
b$diff <- b$auc - b$auc_random
b$design <- factor(b$design, levels = c("clustered", "stratified"), labels = c("Clustered vs random", "Stratified vs random"))
pb <- ggplot(b, aes(n, diff, colour = method, linetype = design, shape = design)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.9) +
  scale_x_log10(breaks = ns) + scale_colour_manual(values = METHOD_COL, guide = "none") +
  scale_linetype_manual(values = c("solid", "22")) + scale_shape_manual(values = c(16, 15)) +
  guides(linetype = guide_legend(override.aes = list(colour = "black")), shape = "legend") +
  labs(x = "Sites per timestep", y = "Change in mean AUC relative to random", linetype = NULL, shape = NULL)
fig3 <- (pa + pb) + plot_layout(widths = c(1.15, 1), guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")
save_png(fig3, "main/Figure3.png", 10.4, 4.4)
