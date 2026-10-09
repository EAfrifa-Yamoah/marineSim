# Figure 2: the five method core benchmark from results/core_five_methods.csv
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
  theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_hline(yintercept = 0.5, linetype = "13", linewidth = 0.4) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.9) +
  scale_x_log10(breaks = ns) + scale_colour_manual(values = METHOD_COL) + scale_fill_manual(values = METHOD_COL) +
  labs(x = "Sites per timestep", y = "Mean test AUC", colour = NULL, fill = NULL)

# (b, c) change in mean AUC relative to random sampling, by method: clustered and stratified on their own scales
dm <- aggregate(auc ~ design + n + method, long, mean)
b <- merge(dm[dm$design != "random", ], dm[dm$design == "random", c("n", "method", "auc")], by = c("n", "method"), suffixes = c("", "_random"))
b$diff <- b$auc - b$auc_random
b$design <- factor(b$design, levels = c("clustered", "stratified"), labels = c("Clustered vs random", "Stratified vs random"))
pdl <- position_dodge(width = 0.12)      # horizontal offsets on the log10 axis so that points do not overlap
diff_panel <- function(dsg, ylab) {
  ggplot(b[b$design == dsg, ], aes(n, diff, colour = method)) +
    geom_line(linewidth = 0.7, position = pdl) + geom_point(size = 1.9, position = pdl) +
    scale_x_log10(breaks = ns) + scale_colour_manual(values = METHOD_COL, guide = "none") +
    labs(x = "Sites per timestep", y = ylab) +
    theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))
}
pb <- diff_panel("Clustered vs random", "Clustered minus random (\u0394AUC)") + scale_y_continuous(breaks = seq(-0.20, -0.05, 0.05))
pc <- diff_panel("Stratified vs random", "Stratified minus random (\u0394AUC)")
fig2 <- (pa + pb + pc) + plot_layout(widths = c(1.3, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")
save_png(fig2, "main/Figure2.png", 12.5, 4.4)
