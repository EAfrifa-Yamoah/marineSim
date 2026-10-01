"""Figures 2, 4 and 5: feature block Shapley decomposition from results/expanded_T*.csv."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
R = os.environ.get("MS_RESULTS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")) + "/"
FIG = os.environ.get("MS_FIGURES", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")) + "/"

import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D; from matplotlib.patches import Patch
try:
    import figstyle as fs; fs.setup()
except Exception: fs=None
O=R; F=FIG+"main/"; os.makedirs(F,exist_ok=True)
COMP=["Algorithmic flexibility","Coordinates","Distance fields","Spatial lags","Temporal"]
COL={"Algorithmic flexibility":"#999999","Coordinates":"#56B4E9","Distance fields":"#0072B2","Spatial lags":"#009E73","Temporal":"#D55E00"}
SHORT={"Algorithmic flexibility":"Algorithmic","Coordinates":"Coordinates","Distance fields":"Distance fields","Spatial lags":"Spatial lags","Temporal":"Temporal"}
t1=pd.read_csv(O+"expanded_T1_decomposition_by_design.csv"); t2=pd.read_csv(O+"expanded_T2_temporal_by_T.csv"); t3=pd.read_csv(O+"expanded_T3_detection.csv"); t4=pd.read_csv(O+"expanded_T4_abundance.csv"); t5=pd.read_csv(O+"expanded_T5_NN0_by_design.csv")
def bars(ax,tab,groups,ylabel):
    xs=np.arange(len(groups)); w=0.16
    for i,c in enumerate(COMP):
        v=[tab[(tab.subset==g)&(tab.component==c)].value.iloc[0] for g in groups]; lo=[v[j]-tab[(tab.subset==g)&(tab.component==c)].ci_lo.iloc[0] for j,g in enumerate(groups)]; hi=[tab[(tab.subset==g)&(tab.component==c)].ci_hi.iloc[0]-v[j] for j,g in enumerate(groups)]
        ax.bar(xs+(i-2)*w,v,w,color=COL[c],yerr=[lo,hi],capsize=2,error_kw=dict(lw=0.7,capthick=0.7))
    ax.axhline(0,color="0.5",lw=0.6); ax.set_xticks(xs); ax.set_ylabel(ylabel)
leg=[Patch(fc=COL[c],label=SHORT[c]) for c in COMP]
# ---- Figure 2: (a) by design, (b) NN0 control by design
fig,ax=plt.subplots(1,2,figsize=(10,4.0),gridspec_kw=dict(width_ratios=[1.35,1])); fig.subplots_adjust(left=0.07,right=0.985,top=0.94,bottom=0.26,wspace=0.28)
g=["pooled (all designs)","random","stratified","clustered"]; bars(ax[0],t1,g,"Shapley contribution (ΔAUC)"); ax[0].set_xticklabels(["pooled","random","stratified","clustered"]); ax[0].text(0.01,0.97,"(a)",transform=ax[0].transAxes,va="top",fontweight="bold")
d=["random","stratified","clustered"]; xs=np.arange(3); w=0.34
for i,(q,lab,hatch) in enumerate([("Temporal marginal, no control","no control",""),("Temporal marginal, NN0 controlled","same time neighbour control","///")]):
    s=t5[t5.quantity==q].set_index("design").loc[d]; ax[1].bar(xs+(i-0.5)*w,s.value,w,color=COL["Temporal"],hatch=hatch,edgecolor="black",lw=0.5,yerr=[s.value-s.ci_lo,s.ci_hi-s.value],capsize=2,error_kw=dict(lw=0.7))
ax[1].set_xticks(xs); ax[1].set_xticklabels(d); ax[1].set_ylabel("Temporal increment when added last (ΔAUC)"); ax[1].text(0.01,0.97,"(b)",transform=ax[1].transAxes,va="top",fontweight="bold")
fig.legend(handles=leg+[Patch(fc=COL["Temporal"],ec="black",label="temporal, no control (b)"),Patch(fc=COL["Temporal"],ec="black",hatch="///",label="temporal, NN0 controlled (b)")],loc="lower center",ncol=4,frameon=False,fontsize=8,bbox_to_anchor=(0.5,0.0))
fig.savefig(F+"Figure2.png",dpi=300); plt.close(fig)
# ---- Figure 4: temporal by T
fig,ax=plt.subplots(figsize=(6.4,3.9)); fig.subplots_adjust(left=0.12,right=0.98,top=0.95,bottom=0.30)
st={"Shapley":dict(marker="o",ls="-"),"Added last":dict(marker="s",ls="--"),"Added last, NN0 controlled":dict(marker="^",ls=":")}
for q,s_ in st.items():
    dd=t2[t2.quantity==q].sort_values("n_time"); ax.errorbar(dd.n_time,dd.value,yerr=[dd.value-dd.ci_lo,dd.ci_hi-dd.value],color=COL["Temporal"],capsize=2,lw=1.5,ms=5,**s_)
ax.axhline(0,color="0.5",lw=0.6); ax.set_xticks([1,3,5,10]); ax.set_xlabel("Number of timesteps"); ax.set_ylabel("Temporal contribution (ΔAUC)")
ax.text(1,0.004,"placebo",fontsize=7.5,color="0.35",ha="center")
fig.legend(handles=[Line2D([],[],color=COL["Temporal"],lw=1.5,ms=5,label=q,**s_) for q,s_ in st.items()],loc="lower center",ncol=3,frameon=False,fontsize=8,bbox_to_anchor=(0.5,0.0))
fig.savefig(F+"Figure4.png",dpi=300); plt.close(fig)
# ---- Figure 5: (a) detection, (b) abundance shares
fig,ax=plt.subplots(1,2,figsize=(10,4.0)); fig.subplots_adjust(left=0.07,right=0.985,top=0.94,bottom=0.26,wspace=0.26)
gd=sorted(t3.subset.unique(),key=lambda s:-float(s.split()[-1])); bars(ax[0],t3,gd,"Shapley contribution (ΔAUC)"); ax[0].set_xticklabels([f"detection {float(s.split()[-1]):.1f}" for s in gd]); ax[0].text(0.01,0.97,"(a)",transform=ax[0].transAxes,va="top",fontweight="bold")
bp=t1[t1.subset=="pooled (all designs)"]; ab=t4[t4.subset=="abundance pooled"]; xs=np.arange(5); w=0.36
ax[1].bar(xs-w/2,[bp[bp.component==c].share_pct.iloc[0] for c in COMP],w,color=[COL[c] for c in COMP],edgecolor="black",lw=0.5)
ax[1].bar(xs+w/2,[ab[ab.component==c].share_pct.iloc[0] for c in COMP],w,color=[COL[c] for c in COMP],edgecolor="black",lw=0.5,hatch="///")
ax[1].set_xticks(xs); ax[1].set_xticklabels([SHORT[c] for c in COMP],fontsize=7.5,rotation=12); ax[1].set_ylabel("Share of total gain (%)"); ax[1].text(0.01,0.97,"(b)",transform=ax[1].transAxes,va="top",fontweight="bold")
fig.legend(handles=leg+[Patch(fc="white",ec="black",label="occurrence, AUC (b)"),Patch(fc="white",ec="black",hatch="///",label="abundance, Spearman ρ (b)")],loc="lower center",ncol=4,frameon=False,fontsize=8,bbox_to_anchor=(0.5,0.0))
fig.savefig(F+"Figure5.png",dpi=300); plt.close(fig)
