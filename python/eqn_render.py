"""Render the stRF display equations to tight PNGs via matplotlib mathtext."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

F = os.path.join(os.path.dirname(__file__), "..", "figures")
DPI = 220

EQS = {
    "eq1_features":
        r"$\Phi(\mathbf{s},t)=[\,\mathbf{z}(\mathbf{s},t),\;"
        r"L_1(\mathbf{s},t),\dots,L_M(\mathbf{s},t),\;"
        r"\tau_1(\mathbf{s},t),\;\tau_H(\mathbf{s},t),\;"
        r"e_1(\mathbf{s}),\dots,e_K(\mathbf{s})\,]^{\top}$",
    "eq2_lag":
        r"$L_m(\mathbf{s},t)=\frac{1}{|N_m(\mathbf{s})|}\sum_{j\in N_m(\mathbf{s})}y_j,"
        r"\qquad N_m(\mathbf{s})=\{\,j:0<\|\mathbf{s}-\mathbf{s}_j\|\leq r_m\,\}$",
    "eq3_temporal":
        r"$\tau_1(\mathbf{s},t)=y(\mathbf{s},\,t-1),\qquad"
        r"\tau_H(\mathbf{s},t)=\frac{1}{t-1}\sum_{u=1}^{t-1}y(\mathbf{s},u)$",
    "eq4_edf":
        r"$e_k(\mathbf{s})=\|\mathbf{s}-\mathbf{a}_k\|,\qquad k=1,\dots,K$",
    "eq5_forest":
        r"$\hat{f}(\Phi)=\frac{1}{B}\sum_{b=1}^{B}T_b(\Phi;\Theta_b),"
        r"\qquad T_b(\Phi;\Theta_b)=\sum_{\ell}\bar{y}_{b\ell}\,"
        r"\mathbf{1}\{\Phi\in R_{b\ell}\}$",
    "eq6_split":
        r"$(q^{*},c^{*})=\mathrm{arg\,min}_{q,\,c}\,[\,\mathcal{L}(R_{\mathrm{L}})"
        r"+\mathcal{L}(R_{\mathrm{R}})\,],\quad q\in\mathcal{M}_{\mathrm{try}}"
        r"\subset\{1,\dots,p{+}M{+}2{+}K\}$",
    "eq7_additive":
        r"$\mathrm{E}\!\left[\,Y\mid\mathbf{s},t\,\right]\approx\; "
        r"g(\mathbf{z}(\mathbf{s},t))+\varphi_{S}(\mathbf{s})+\varphi_{T}(\mathbf{s},t)$",
    "eq8_decomp":
        r"$A_{\mathrm{stRF}}-A_{\mathrm{GLM}}=(A_{\mathrm{RF}}-A_{\mathrm{GLM}})"
        r"+(A_{\mathrm{sRF}}-A_{\mathrm{RF}})+(A_{\mathrm{stRF}}-A_{\mathrm{sRF}})$",
    "eq9_di":
        r"$\mathrm{DI}(\mathbf{s})=\frac{\min_{i}\,d(\tilde{\Phi}(\mathbf{s}),"
        r"\tilde{\Phi}(\mathbf{s}_i))}{\bar{d}},\qquad"
        r"\bar{d}=\frac{1}{n(n-1)}\sum_{i\neq j}d(\tilde{\Phi}(\mathbf{s}_i),"
        r"\tilde{\Phi}(\mathbf{s}_j))$",
    "eq10_diw":
        r"$d(\tilde{\Phi}(\mathbf{s}),\tilde{\Phi}(\mathbf{s}'))="
        r"\sqrt{\sum_{q}w_q\left(\frac{\Phi_q(\mathbf{s})-\mu_q}{\sigma_q}-"
        r"\frac{\Phi_q(\mathbf{s}')-\mu_q}{\sigma_q}\right)^{2}},\qquad w_q=\sqrt{I_q}$",
}


def render():
    ok = []
    for name, tex in EQS.items():
        try:
            fig = plt.figure(figsize=(0.1, 0.1))
            fig.text(0.0, 0.0, tex, fontsize=15, color="black")
            path = f"{F}/{name}.png"
            fig.savefig(path, dpi=DPI, bbox_inches="tight", pad_inches=0.06,
                        facecolor="white", transparent=False)
            plt.close(fig)
            w, h = Image.open(path).size
            ok.append((name, w, h))
        except Exception as e:
            print(f"FAIL {name}: {str(e)[:90]}")
    for name, w, h in ok:
        print(f"{name}: {w}x{h}")
    return ok


if __name__ == "__main__":
    render()
    print("DONE", )
