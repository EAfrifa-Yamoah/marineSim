# Figure S8: stRF versus sdmTMB from results/si/tableS6_full.csv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figstyle.R"))
d <- read.csv(res("si", "tableS6_full.csv"))
MOD <- c(stRF = "stRF", `stRF (pooled rows)` = "stRF_pooled", `sdmTMB spatial` = "sdmTMB_spatial",
         `sdmTMB spatiotemporal AR(1)` = "sdmTMB_spatiotemporal", GLM = "GLM")
MCOL <- c(stRF = OI[["vermilion"]], `stRF (pooled rows)` = OI[["orange"]], `sdmTMB spatial` = OI[["blue"]],
          `sdmTMB spatiotemporal AR(1)` = OI[["green"]], GLM = OI[["grey"]])
TCOL <- c(stRF = "t_stRF", `sdmTMB spatial` = "t_spatial", `sdmTMB spatiotemporal AR(1)` = "t_st")
acc <- do.call(rbind, lapply(names(MOD), function(m) { a <- aggregate(d[[MOD[[m]]]], list(n = d$n), mean, na.rm = TRUE); data.frame(n = a$n, model = m, auc = a$x) }))
tim <- do.call(rbind, lapply(names(TCOL), function(m) { a <- aggregate(d[[TCOL[[m]]]], list(n = d$n), median, na.rm = TRUE); data.frame(n = a$n, model = m, t = a$x) }))
acc$model <- factor(acc$model, levels = names(MOD)); tim$model <- factor(tim$model, levels = names(MOD))
ns <- sort(unique(d$n))
pa <- ggplot(acc, aes(n, auc, colour = model)) + geom_line(linewidth = 0.7) + geom_point(size = 1.9) +
  scale_x_log10(breaks = ns) + scale_colour_manual(values = MCOL, drop = FALSE) + labs(x = "Sites per timestep (n)", y = "Mean test AUC", colour = NULL)
pb <- ggplot(tim, aes(n, t, colour = model)) + geom_line(linewidth = 0.7) + geom_point(size = 1.9) +
  scale_x_log10(breaks = ns) + scale_y_log10() + scale_colour_manual(values = MCOL, drop = FALSE) + labs(x = "Sites per timestep (n)", y = "Median fit time (s)", colour = NULL)
figS8 <- (pa + pb) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom") & guides(colour = guide_legend(nrow = 2))
save_png(figS8, "si/FigureS8.png", 9.2, 4.3)
