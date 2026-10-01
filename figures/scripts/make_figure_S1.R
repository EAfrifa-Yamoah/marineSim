# Figure S1: conceptual diagram of the stRF framework (no data input)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "figstyle.R"))
C <- c(C = OI[["skyblue"]], E = OI[["blue"]], L = OI[["green"]], T = OI[["vermilion"]])
boxes <- rbind(
  data.frame(x = 6, y = 86, w = 40, h = 10, title = "Survey data", sub = "n sites, revisited over T timesteps", fc = "#f4f4f4", ec = "grey35"),
  data.frame(x = 54, y = 86, w = 40, h = 10, title = "Environmental covariates", sub = "measured predictors z(s, t)", fc = "#f4f4f4", ec = "grey35"),
  data.frame(x = 5, y = 55, w = 21.5, h = 15, title = "Coordinates (C)", sub = "projected position", fc = "#dff0fa", ec = C[["C"]]),
  data.frame(x = 29, y = 55, w = 21.5, h = 15, title = "Distance fields (E)", sub = "distance to K anchors", fc = "#d6e6f4", ec = C[["E"]]),
  data.frame(x = 53, y = 55, w = 21.5, h = 15, title = "Spatial lags (L)", sub = "means at three radii", fc = "#d9f0e8", ec = C[["L"]]),
  data.frame(x = 77, y = 55, w = 21.5, h = 15, title = "Temporal features (T)", sub = "previous step + site history", fc = "#fbe4d7", ec = C[["T"]]),
  data.frame(x = 22, y = 31, w = 56, h = 14, title = "Random forest ensemble", sub = "150 trees; any coalition S of the four blocks", fc = "#f4f4f4", ec = "grey35"),
  data.frame(x = 2, y = 8, w = 30, h = 12, title = "Prediction surface", sub = "occurrence or abundance", fc = "#f4f4f4", ec = "grey35"),
  data.frame(x = 35, y = 8, w = 30, h = 12, title = "Shapley decomposition", sub = "over {C, E, L, T} + algorithmic term", fc = "#f4f4f4", ec = "grey35"),
  data.frame(x = 68, y = 8, w = 30, h = 12, title = "Area of applicability", sub = "flags extrapolation in space and time", fc = "#f4f4f4", ec = "grey35"))
arrows <- data.frame(x = c(26, 74, 50, 50, 35, 65), y = c(86, 86, 52, 26, 26, 26), xend = c(26, 74, 50, 50, 17, 83), yend = c(79, 79, 45, 20, 20, 20))
figS1 <- ggplot() + coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
  annotate("rect", xmin = 3, xmax = 97, ymin = 52, ymax = 79, fill = NA, colour = "grey40", linetype = "22", linewidth = 0.4) +
  annotate("text", x = 50, y = 76, label = "bold('Augmented feature map  ') * Phi[S](s, t) == group('[', list(z(s, t), 'blocks in S'), ']')", parse = TRUE, size = 3.4) +
  geom_rect(data = boxes, aes(xmin = x, xmax = x + w, ymin = y, ymax = y + h), fill = boxes$fc, colour = boxes$ec, linewidth = 0.4) +
  geom_text(data = boxes, aes(x + w / 2, y + h * 0.62, label = title), size = 3.3, fontface = "bold") +
  geom_text(data = boxes, aes(x + w / 2, y + h * 0.27, label = sub), size = 2.8, colour = "grey25") +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend), colour = "grey30", linewidth = 0.4,
               arrow = arrow(length = unit(5, "pt"), type = "closed")) +
  annotate("text", x = 50, y = 29.2, label = "RF = covariates only     Spatial RF = {C}     stRF = {C, E, L, T}", size = 2.8, colour = "grey25", vjust = 1) +
  annotate("text", x = 50, y = 2.5, label = "Controls: single timestep placebo (T constant)  ·  same timestep neighbour (NN0) for the temporal block", size = 2.6, colour = "grey35") +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
save_png(figS1, "si/FigureS1.png", 9.2, 7.0)
