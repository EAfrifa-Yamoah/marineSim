$$\Phi(\mathbf{s},t) = \left[\, \mathbf{z}(\mathbf{s},t),\ L_1(\mathbf{s},t),\dots,L_M(\mathbf{s},t),\ \tau_1(\mathbf{s},t),\ \tau_H(\mathbf{s},t),\ e_1(\mathbf{s}),\dots,e_K(\mathbf{s}) \,\right]^{\top}$$

$$L_m(\mathbf{s},t) = \frac{1}{|N_m(\mathbf{s})|} \sum_{j \in N_m(\mathbf{s})} y_j, \qquad N_m(\mathbf{s}) = \{\, j : 0 < \|\mathbf{s} - \mathbf{s}_j\| \le r_m \,\}$$

$$\tau_1(\mathbf{s},t) = y(\mathbf{s},\, t-1), \qquad \tau_H(\mathbf{s},t) = \frac{1}{t-1} \sum_{u=1}^{t-1} y(\mathbf{s}, u)$$

$$e_k(\mathbf{s}) = \|\mathbf{s} - \mathbf{a}_k\|, \qquad \{\mathbf{a}_k\}_{k=1}^{K} = k\text{-means centroids of } \{\mathbf{s}_i\}$$

$$\hat{f}(\Phi) = \frac{1}{B} \sum_{b=1}^{B} T_b(\Phi; \Theta_b), \qquad T_b(\Phi; \Theta_b) = \sum_{\ell} \bar{y}_{b\ell}\, \mathbf{1}\{\Phi \in R_{b\ell}\}$$

$$(q^{*}, c^{*}) = \operatorname*{arg\,min}_{q,\, c} \left[\, \mathcal{L}(R_{\mathrm{L}}) + \mathcal{L}(R_{\mathrm{R}}) \,\right], \qquad q \in \mathcal{M}_{\mathrm{try}} \subset \{1, \dots, p + M + 2 + K\}$$

$$\mathbb{E}\!\left[\, Y \mid \mathbf{s}, t \,\right] \approx g(\mathbf{z}(\mathbf{s},t)) + \varphi_S(\mathbf{s}) + \varphi_T(\mathbf{s}, t)$$

$$A_{\mathrm{stRF}} - A_{\mathrm{GLM}} = \underbrace{(A_{\mathrm{RF}} - A_{\mathrm{GLM}})}_{\text{ML}} + \underbrace{(A_{\mathrm{sRF}} - A_{\mathrm{RF}})}_{\text{spatial}} + \underbrace{(A_{\mathrm{stRF}} - A_{\mathrm{sRF}})}_{\text{temporal}}$$

$$\mathrm{DI}(\mathbf{s}) = \frac{\min_{i}\, d\!\left(\tilde{\Phi}(\mathbf{s}), \tilde{\Phi}(\mathbf{s}_i)\right)}{\bar{d}}, \qquad \bar{d} = \frac{1}{n(n-1)} \sum_{i \neq j} d\!\left(\tilde{\Phi}(\mathbf{s}_i), \tilde{\Phi}(\mathbf{s}_j)\right)$$

$$d\!\left(\tilde{\Phi}(\mathbf{s}), \tilde{\Phi}(\mathbf{s}')\right) = \sqrt{\sum_{q} w_q \left( \frac{\Phi_q(\mathbf{s}) - \mu_q}{\sigma_q} - \frac{\Phi_q(\mathbf{s}') - \mu_q}{\sigma_q} \right)^{2}}, \qquad w_q = \sqrt{I_q}$$
