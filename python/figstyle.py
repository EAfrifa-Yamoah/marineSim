"""
figstyle.py  - shared publication style (Okabe-Ito, 300 dpi, clean axes).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Okabe-Ito colourblind-safe palette
OI = dict(black="#000000", orange="#E69F00", skyblue="#56B4E9",
          green="#009E73", yellow="#F0E442", blue="#0072B2",
          vermillion="#D55E00", purple="#CC79A7", grey="#999999")

METHOD_COLOR = {"GLM": OI["grey"], "GAM": OI["orange"], "RF": OI["skyblue"],
                "SpatialRF": OI["blue"], "stRF": OI["vermillion"]}
METHOD_MARKER = {"GLM": "o", "GAM": "s", "RF": "^", "SpatialRF": "D", "stRF": "P"}
COMP_COLOR = {"ml": OI["skyblue"], "spatial": OI["blue"],
              "temporal": OI["green"], "total": OI["black"]}


def setup():
    rcParams.update({
        "figure.dpi": 120, "savefig.dpi": 300,
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 11, "axes.labelsize": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.9, "xtick.direction": "out", "ytick.direction": "out",
        "legend.frameon": False, "legend.fontsize": 9,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
        "figure.facecolor": "white", "axes.facecolor": "white",
    })


def save(fig, path):
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", path)
