# figstyle.R — shared publication style for every figure (sourced by the make_*.R scripts)
# Okabe Ito palette; fixed meaning per colour across the article; no titles inside
# figures; one shared legend per figure.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(grid)})
invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "C.UTF-8")))   # Greek letters in axis labels

HERE    <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = FALSE)
if (is.na(HERE) || !nzchar(HERE)) HERE <- getwd()
RESULTS <- Sys.getenv("MS_RESULTS", file.path(HERE, "..", "..", "results"))
FIGDIR  <- Sys.getenv("MS_FIGURES", file.path(HERE, ".."))
res <- function(...) file.path(RESULTS, ...)

OI <- c(grey = "#999999", orange = "#E69F00", skyblue = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", vermilion = "#D55E00", purple = "#CC79A7")

# Methods (Figures 3, 6, S3, S4)
METHOD_COL <- c("GLM" = OI[["grey"]], "GAM" = OI[["green"]], "Standard RF" = OI[["orange"]],
                "Spatial RF" = OI[["skyblue"]], "stRF" = OI[["vermilion"]])
# Feature blocks and the algorithmic term (Figures 2, 4, 5, S1, 6b/c)
COMP_COL <- c("Algorithmic" = OI[["grey"]], "Coordinates" = OI[["skyblue"]], "Distance fields" = OI[["blue"]],
              "Spatial lags" = OI[["green"]], "Temporal" = OI[["vermilion"]])
# Marine parks (Figures 7, 8, S5)
REG   <- c("JBMP", "MMP", "CSMC", "SIMP", "NCMP")
REG_L <- c(JBMP = "Jurien Bay", MMP = "Marmion", CSMC = "Cockburn Sound", SIMP = "Shoalwater", NCMP = "Ngari Capes")
REG_COL <- c(JBMP = "#0072B2", MMP = "#E69F00", CSMC = "#009E73", SIMP = "#CC79A7", NCMP = "#5D3A9B")

theme_ms <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "DejaVu Sans") +
    theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
          axis.text = element_text(colour = "black", size = base_size - 1),
          legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),
          legend.text = element_text(size = base_size - 1), legend.key.size = unit(10, "pt"),
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(4, "pt"),
          plot.tag = element_text(face = "bold", size = base_size + 1),
          plot.margin = margin(4, 6, 2, 4), panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA))
}
theme_set(theme_ms())
# an in-panel tag at the top left, as in the article's figures
tag_in <- function(lab, x = 0.02, y = 0.97, hjust = 0, size = 3.4)
  annotation_custom(grid::textGrob(lab, x = unit(x, "npc"), y = unit(y, "npc"), hjust = hjust, vjust = 1,
                                   gp = gpar(fontface = "bold", fontsize = size * .pt / 1.2)))
save_png <- function(p, name, width, height, dpi = 300) {
  out <- file.path(FIGDIR, name); dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  ggsave(out, p, width = width, height = height, dpi = dpi, units = "in", bg = "white",
         device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else grDevices::png)
  cat(name, "\n")
}
