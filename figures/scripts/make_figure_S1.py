"""Figure S1: conceptual diagram of the stRF framework (no data input)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
R = os.environ.get("MS_RESULTS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")) + "/"
FIG = os.environ.get("MS_FIGURES", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")) + "/"

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

try:
    import figstyle as fs; fs.setup()
except Exception: pass
COL={"C":"#56B4E9","E":"#0072B2","L":"#009E73","T":"#D55E00","grey":"#e6e6e6","alg":"#999999"}
fig,ax=plt.subplots(figsize=(9.2,7.0)); ax.set_xlim(0,100); ax.set_ylim(0,100); ax.axis("off")
def box(x,y,w,h,title,sub,fc="#f4f4f4",ec="0.35",tfs=9.5,sfs=8):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.4,rounding_size=1.5",fc=fc,ec=ec,lw=1)); ax.text(x+w/2,y+h*0.62,title,ha="center",va="center",fontsize=tfs,fontweight="bold"); ax.text(x+w/2,y+h*0.27,sub,ha="center",va="center",fontsize=sfs,color="0.25")
def arrow(x1,y1,x2,y2): ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=12,lw=1,color="0.3"))
# inputs
box(6,86,40,10,"Survey data","n sites, revisited over T timesteps"); box(54,86,40,10,"Environmental covariates","measured predictors z(s, t)")
arrow(26,86,26,79); arrow(74,86,74,79)
# feature map container
ax.add_patch(FancyBboxPatch((3,52),94,27,boxstyle="round,pad=0.4,rounding_size=2",fc="none",ec="0.4",lw=1,ls="--"))
ax.text(50,76,"Augmented feature map  Φ_S(s, t) = [ z(s, t), blocks in S ]",ha="center",va="center",fontsize=9.5,fontweight="bold")
bw=21.5; xs=[5,29,53,77]
box(xs[0],55,bw,15,"Coordinates (C)","projected position","#dff0fa",COL["C"]); box(xs[1],55,bw,15,"Distance fields (E)","distance to K anchors","#d6e6f4",COL["E"])
box(xs[2],55,bw,15,"Spatial lags (L)","means at three radii","#d9f0e8",COL["L"]); box(xs[3],55,bw,15,"Temporal features (T)","previous step + site history","#fbe4d7",COL["T"])
arrow(50,52,50,45)
# forest with coalitions
box(22,31,56,14,"Random forest ensemble","150 trees; any coalition S of the four blocks","#f4f4f4","0.35")
ax.text(50,29.2,"RF = covariates only     Spatial RF = {C}     stRF = {C, E, L, T}",ha="center",va="top",fontsize=8,color="0.25")
arrow(50,26,50,20); arrow(35,26,17,20); arrow(65,26,83,20)
# outputs
box(2,8,30,12,"Prediction surface","occurrence or abundance"); box(35,8,30,12,"Shapley decomposition","over {C, E, L, T} + algorithmic term"); box(68,8,30,12,"Area of applicability","flags extrapolation in space and time")
# controls note
ax.text(50,2.5,"Controls: single timestep placebo (T constant)  ·  same timestep neighbour (NN0) for the temporal block",ha="center",va="center",fontsize=7.6,color="0.35")
os.makedirs(FIG + "si", exist_ok=True); fig.savefig(FIG + "si/FigureS1.png", dpi=300, bbox_inches="tight"); print("Figure S1 redrawn")
