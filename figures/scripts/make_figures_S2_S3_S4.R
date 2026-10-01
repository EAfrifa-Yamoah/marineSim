# Figures S2, S3, S4 from results/si/ (worlds and designs; tuning sensitivity; calibration)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figstyle.R"))
G <- 50; CELL <- 8
cell_xy <- function(cell) data.frame(x = ((cell - 1) %/% G + 0.5) * CELL, y = ((cell - 1) %% G + 0.5) * CELL)

# ------------------------------------------------------------------ Figure S2
w <- read.csv(res("si", "figS2_worlds.csv")); des <- read.csv(res("si", "figS2_designs.csv")); rep <- read.csv(res("si", "figS2_repworld.csv"))
w <- cbind(w, cell_xy(w$cell)); rep <- cbind(rep, cell_xy(rep$cell))
w$panel <- factor(sprintf("(%s) %s, %d km", letters[(!w$stationary) * 3 + match(w$range_km, c(30, 80, 200))],
                          ifelse(w$stationary, "stationary", "non stationary"), w$range_km))
p_world <- ggplot(w, aes(x, y, fill = p)) + geom_raster() + facet_wrap(~panel, ncol = 3) + coord_equal(expand = FALSE) +
  scale_fill_viridis_c(limits = c(0, 1), name = "Occurrence\nprobability") +
  scale_x_continuous(breaks = c(0, 200, 400)) + scale_y_continuous(breaks = c(0, 200, 400)) +
  labs(x = NULL, y = NULL) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, face = "bold", size = 8),
        legend.position = "right", legend.title = element_text(size = 8), axis.text = element_text(size = 7))
des$row <- ifelse(grepl("^design_", des$panel), "design", "n")
des$panel <- factor(des$panel, levels = c("design_random", "design_clustered", "design_stratified", "n_30", "n_100", "n_500"),
                    labels = c("(g) random, n = 100", "(h) clustered, n = 100", "(i) stratified, n = 100", "(j) random, n = 30", "(k) random, n = 100", "(l) random, n = 500"))
des$obs <- factor(ifelse(des$y_obs == 1, "presence", "absence"), levels = c("presence", "absence"))
bg <- merge(rep[, c("x", "y", "p")], data.frame(panel = levels(des$panel)))
p_des <- ggplot() + geom_raster(data = bg, aes(x, y, alpha = p), fill = "black") + scale_alpha_continuous(range = c(0, 0.6), guide = "none") +
  geom_point(data = des, aes(x, y, fill = obs), shape = 21, size = 1.5, stroke = 0.25) +
  scale_fill_manual(values = c(presence = OI[["vermilion"]], absence = OI[["skyblue"]]), name = "Sampled cells") +
  facet_wrap(~panel, ncol = 3) + coord_equal(expand = FALSE, xlim = c(0, 400), ylim = c(0, 400)) +
  scale_x_continuous(breaks = c(0, 200, 400)) + scale_y_continuous(breaks = c(0, 200, 400)) + labs(x = "km", y = NULL) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, face = "bold", size = 8),
        legend.position = "right", legend.title = element_text(size = 8), axis.text = element_text(size = 7))
figS2 <- p_world / p_des + plot_layout(heights = c(2, 2))
save_png(figS2, "si/FigureS2.png", 8.4, 10.8)

# ------------------------------------------------------------------ Figure S3
t <- read.csv(res("si", "figS3_tuning.csv")); agg <- aggregate(auc ~ method + min_node_size + mtry_rule, t, mean)
agg$mtry_rule <- factor(agg$mtry_rule, levels = c("sqrt", "half", "p8"), labels = c("√p", "0.5p", "0.8p"))
agg <- agg[order(agg$min_node_size, agg$mtry_rule), ]
agg$setting <- factor(paste0("node ", agg$min_node_size, "\nmtry ", agg$mtry_rule), levels = unique(paste0("node ", agg$min_node_size, "\nmtry ", agg$mtry_rule)))
agg$method <- factor(agg$method, levels = c("Standard RF", "Spatial RF", "stRF"))
agg$xi <- as.integer(agg$setting)
figS3 <- ggplot(agg, aes(xi, auc, colour = method, group = method)) +
  scale_x_continuous(breaks = seq_along(levels(agg$setting)), labels = levels(agg$setting)) +
  geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.4, linetype = "22") +
  annotate("text", x = 1.15, y = max(agg$auc) + 0.004, label = "default", size = 2.5, colour = "grey40", hjust = 0) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) + scale_colour_manual(values = METHOD_COL) +
  labs(x = "Hyperparameter setting (minimum node size, features per split)", y = "Mean test AUC", colour = NULL) +
  theme(axis.text.x = element_text(size = 7))
save_png(figS3, "si/FigureS3.png", 7.2, 4.0)

# ------------------------------------------------------------------ Figure S4
cal <- read.csv(res("si", "figS4_calibration.csv")); cal$bin <- pmin(pmax(floor(cal$p_hat * 10), 0), 9)
rel <- aggregate(cbind(p_hat, y, n = 1) ~ method + bin, cal, sum); rel$p_hat <- rel$p_hat / rel$n; rel$y <- rel$y / rel$n; rel <- rel[rel$n >= 30, ]
rel$method <- factor(rel$method, levels = names(METHOD_COL))
pa <- ggplot(rel, aes(p_hat, y, colour = method)) + geom_abline(linetype = "22", colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) + scale_colour_manual(values = METHOD_COL) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + labs(x = "Predicted probability", y = "Observed frequency", colour = NULL)
tab5 <- read.csv(res("si", "tableS5_brier.csv")); tab5$method <- factor(tab5$method, levels = rev(tab5$method))
pb <- ggplot(tab5, aes(brier, method, fill = method)) + geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
  geom_errorbarh(aes(xmin = brier - sd, xmax = brier + sd), height = 0.25, linewidth = 0.4) +
  scale_fill_manual(values = METHOD_COL, guide = "none") + labs(x = "Brier score (lower is better)", y = NULL)
figS4 <- (pa + pb) + plot_layout(widths = c(1, 1.1), guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")
save_png(figS4, "si/FigureS4.png", 8.6, 4.2)
