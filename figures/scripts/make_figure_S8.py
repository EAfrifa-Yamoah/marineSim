"""Figure S8: stRF versus sdmTMB from results/si/tableS6_full.csv."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
R = os.environ.get("MS_RESULTS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")) + "/"
FIG = os.environ.get("MS_FIGURES", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")) + "/"

import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

try:
    import figstyle as fs; fs.setup()
except Exception: pass
d=pd.read_csv(R + "si/tableS6_full.csv")
M={"stRF":"#D55E00","stRF (pooled rows)":"#E69F00","sdmTMB spatial":"#0072B2","sdmTMB spatiotemporal AR(1)":"#009E73","GLM":"#999999"}
col={"stRF":"stRF","stRF (pooled rows)":"stRF_pooled","sdmTMB spatial":"sdmTMB_spatial","sdmTMB spatiotemporal AR(1)":"sdmTMB_spatiotemporal","GLM":"GLM"}
tcol={"stRF":"t_stRF","sdmTMB spatial":"t_spatial","sdmTMB spatiotemporal AR(1)":"t_st"}
fig,ax=plt.subplots(1,2,figsize=(9.2,3.9)); fig.subplots_adjust(left=0.07,right=0.985,top=0.94,bottom=0.30,wspace=0.26)
g=d.groupby("n")
for m,c in M.items():
    mu=g[col[m]].mean(); ax[0].plot(mu.index,mu.values,"-o",color=c,ms=4.5,lw=1.6)
ax[0].set_xscale("log"); ax[0].set_xticks([30,50,100,200,500]); ax[0].set_xticklabels([30,50,100,200,500]); ax[0].set_xlabel("Sites per timestep (n)"); ax[0].set_ylabel("Mean test AUC")
ax[0].text(0.02,0.97,"(a)",transform=ax[0].transAxes,va="top",fontweight="bold")
for m in tcol:
    mt=g[tcol[m]].median(); ax[1].plot(mt.index,mt.values,"-o",color=M[m],ms=4.5,lw=1.6)
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xticks([30,50,100,200,500]); ax[1].set_xticklabels([30,50,100,200,500]); ax[1].set_xlabel("Sites per timestep (n)"); ax[1].set_ylabel("Median fit time (s)")
ax[1].text(0.02,0.97,"(b)",transform=ax[1].transAxes,va="top",fontweight="bold")
fig.legend(handles=[Line2D([],[],color=c,marker="o",ms=4.5,lw=1.6,label=m) for m,c in M.items()],loc="lower center",ncol=3,frameon=False,fontsize=8.5,bbox_to_anchor=(0.5,0.0))
os.makedirs(FIG + "si", exist_ok=True); fig.savefig(FIG + "si/FigureS8.png", dpi=300); print("FigureS8 written")
