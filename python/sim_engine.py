"""
sim_engine.py
=============
Core engine for the spatio-temporal random forest (stRF) simulation benchmark.

Implements, in NumPy/scikit-learn (a faithful translation of the marineSim R package):
  * Matern Gaussian random field generator (FFT spectral method, smoothness nu = 3/2)
  * Ground-truth marine ecosystem generator (static + AR(1) dynamic drivers,
    optional spatial non stationarity, delta zero inflation)
  * Three sampling designs (random / clustered / stratified) on a revisited network
  * Imperfect detection
  * Feature engineering: multi scale spatial lags, temporal lag, Euclidean distance fields
  * Five methods: GLM, GAM, RF, Spatial RF, stRF  (+ a Matern GP geostatistical comparator)

All spatial units are kilometres on a square domain. Author: E. Afrifa-Yamoah et al.
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftfreq
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.cluster import KMeans
from sklearn.metrics import roc_auc_score
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
import warnings
warnings.filterwarnings("ignore")

# ----------------------------------------------------------------------------- #
#  Domain configuration
# ----------------------------------------------------------------------------- #
GRID = 50              # grid cells per side
EXTENT_KM = 400.0      # domain side length (km)
CELL_KM = EXTENT_KM / GRID   # 8 km cells
RADII_KM = (25.0, 75.0, 150.0)   # multi scale spatial lag radii
N_ANCHORS = 15         # Euclidean distance field anchors
RF_TREES = 150         # trees per forest
N_TEST = 1500          # held out test cells for evaluation


# ----------------------------------------------------------------------------- #
#  Matern Gaussian random field via FFT spectral synthesis
# ----------------------------------------------------------------------------- #
def grf_matern(shape, range_cells, rng, nu=1.5):
    """Return a standardised (mean 0, var 1) Matern GRF on a regular grid.

    Uses spectral synthesis: white noise multiplied by sqrt of the Matern spectral
    density in the Fourier domain. Simulated on a padded grid to suppress wrap-around.
    """
    ny, nx = shape
    pad = 2  # pad factor to limit periodicity artefacts
    My, Mx = ny * pad, nx * pad
    ky = fftfreq(My)[:, None]
    kx = fftfreq(Mx)[None, :]
    k2 = (kx**2 + ky**2)
    alpha = np.sqrt(2.0 * nu) / max(range_cells, 1e-6)   # inverse range scale
    # 2-D Matern spectral density  ~ (alpha^2 + |2 pi k|^2)^{-(nu + d/2)}, d = 2
    spec = (alpha**2 + (2 * np.pi) ** 2 * k2) ** (-(nu + 1.0))
    spec[0, 0] = 0.0   # remove DC -> zero mean
    amp = np.sqrt(spec)
    noise = rng.standard_normal((My, Mx)) + 1j * rng.standard_normal((My, Mx))
    field = np.real(ifft2(fft2(noise) * amp))
    field = field[:ny, :nx]
    field = (field - field.mean()) / (field.std() + 1e-12)
    return field


# ----------------------------------------------------------------------------- #
#  Ground-truth world
# ----------------------------------------------------------------------------- #
class World:
    """A synthetic marine ecosystem (one of the 6 ground-truth worlds)."""

    def __init__(self, stationary, range_km, n_time=5, seed=0,
                 zero_inflation=0.25, sigma_re=2.55, sigma_persist=1.30,
                 persist_range_km=26.0):
        self.stationary = stationary
        self.range_km = range_km
        self.n_time = n_time
        rng = np.random.default_rng(seed)
        self.rng = rng

        g = GRID
        # coordinate grids (km, cell centres)
        xs = (np.arange(g) + 0.5) * CELL_KM
        self.xx, self.yy = np.meshgrid(xs, xs)          # (g, g)

        # static drivers with distinct autocorrelation ranges (km -> cells).
        # These are deliberately short-to-moderate range (rough/patchy): they
        # represent measured local covariates (substrate, rugosity, local
        # bathymetry, local thermal stress). They are rougher than the broad
        # unmeasured spatial process (spatial_re) below, so tree ensembles
        # cannot use the covariate vector as a smooth proxy for location.
        self.depth = grf_matern((g, g), 45 / CELL_KM, rng)
        self.subst = grf_matern((g, g), 35 / CELL_KM, rng)
        self.wave = grf_matern((g, g), 60 / CELL_KM, rng)
        sst0 = grf_matern((g, g), 45 / CELL_KM, rng)

        # AR(1) dynamic SST over time (rho = 0.6)
        rho = 0.6
        ssts = [sst0]
        for _ in range(1, n_time):
            innov = grf_matern((g, g), 45 / CELL_KM, rng)
            ssts.append(rho * ssts[-1] + np.sqrt(1 - rho**2) * innov)
        self.sst = np.stack(ssts, axis=0)               # (n_time, g, g)

        # response-autocorrelation field at the world's range (the 'range' axis)
        self.spatial_re = sigma_re * grf_matern((g, g), range_km / CELL_KM, rng)

        # temporally-persistent, spatially-rough site effect:
        # constant across time, short correlation range -> coordinate/lag features
        # cannot resolve it, but a site's own temporal lag can. This is the
        # mechanism underpinning the conditional value of temporal features.
        self.persist = sigma_persist * grf_matern((g, g), persist_range_km / CELL_KM, rng)

        # spatially varying coefficients (non-stationary worlds only)
        if stationary:
            self.b_depth = np.full((g, g), -0.55)
            self.b_sst = np.full((g, g), -0.45)
        else:
            self.b_depth = -0.55 + 0.20 * grf_matern((g, g), 250 / CELL_KM, rng)
            self.b_sst = -0.45 + 0.20 * grf_matern((g, g), 250 / CELL_KM, rng)

        # fixed linear-predictor terms (modest, mostly linear)
        self.b0 = 0.10
        self.b_depth2 = -0.22      # mild quadratic in depth (niche optimum)
        self.b_subst = 0.45
        self.b_wave = -0.30
        self.b_int = 0.22          # depth x wave interaction: the controlled
                                   # source of pure-ML uplift (linear GLM omits it)

        # structural (delta) zeros: spatially clustered suitable-habitat mask
        hab = grf_matern((g, g), 120 / CELL_KM, rng)
        thr = np.quantile(hab, zero_inflation)
        self.struct = (hab > thr).astype(float)         # 1 = potentially suitable

        # realised occurrence probability and labels for every cell & time
        self.p = np.zeros((n_time, g, g))
        self.y = np.zeros((n_time, g, g), dtype=int)
        for t in range(n_time):
            eta = (self.b0
                   + self.b_depth * self.depth + self.b_depth2 * self.depth**2
                   + self.b_subst * self.subst + self.b_wave * self.wave
                   + self.b_int * self.depth * self.wave
                   + self.b_sst * self.sst[t] + self.spatial_re + self.persist)
            p = self.struct * (1.0 / (1.0 + np.exp(-eta)))
            self.p[t] = p
            self.y[t] = (rng.random((g, g)) < p).astype(int)

    # --- helpers -----------------------------------------------------------
    def covariates(self, idx, t):
        """Static + dynamic covariates for flattened cell indices at time t."""
        d = self.depth.ravel()[idx]
        return np.column_stack([
            d, d**2,
            self.subst.ravel()[idx],
            self.wave.ravel()[idx],
            self.sst[t].ravel()[idx],
        ])

    def coords(self, idx):
        return np.column_stack([self.xx.ravel()[idx], self.yy.ravel()[idx]])


# ----------------------------------------------------------------------------- #
#  Sampling designs (revisited monitoring network: same cells across timesteps)
# ----------------------------------------------------------------------------- #
def sample_sites(n, design, rng):
    """Return n flattened cell indices under the chosen design."""
    g = GRID
    if design == "random":
        return rng.choice(g * g, size=n, replace=False)

    if design == "clustered":
        centres = rng.uniform(0.15, 0.85, size=(3, 2)) * g
        per = np.full(3, n // 3)
        per[: n - per.sum()] += 1
        rows, cols = [], []
        for c, k in zip(centres, per):
            rr = np.clip(rng.normal(c[0], g * 0.06, size=k), 0, g - 1).astype(int)
            cc = np.clip(rng.normal(c[1], g * 0.06, size=k), 0, g - 1).astype(int)
            rows.append(rr); cols.append(cc)
        rows = np.concatenate(rows); cols = np.concatenate(cols)
        flat = rows * g + cols
        flat = np.unique(flat)
        # top up if collisions removed some
        while len(flat) < n:
            extra = rng.choice(g * g, size=n - len(flat), replace=False)
            flat = np.unique(np.concatenate([flat, extra]))
        return flat[:n]

    if design == "stratified":
        s = int(np.ceil(np.sqrt(n)))
        edges = np.linspace(0, g, s + 1).astype(int)
        cells = []
        for i in range(s):
            for j in range(s):
                r = rng.integers(edges[i], max(edges[i] + 1, edges[i + 1]))
                c = rng.integers(edges[j], max(edges[j] + 1, edges[j + 1]))
                cells.append(r * g + c)
        cells = np.unique(np.array(cells))
        rng.shuffle(cells)
        return cells[:n]

    raise ValueError(design)


# ----------------------------------------------------------------------------- #
#  Feature engineering
# ----------------------------------------------------------------------------- #
def spatial_lags(train_xy, train_y, query_xy, radii_km, is_train):
    """Mean neighbour response within concentric rings (km). Self excluded for train."""
    D = cdist(query_xy, train_xy)                # (nq, ntr) km
    gmean = train_y.mean()
    out = np.empty((query_xy.shape[0], len(radii_km)))
    for k, r in enumerate(radii_km):
        mask = D <= r
        if is_train:
            np.fill_diagonal(mask, False)        # exclude self
        cnt = mask.sum(1)
        s = (mask * train_y[None, :]).sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            m = np.where(cnt > 0, s / np.maximum(cnt, 1), gmean)
        out[:, k] = m
    return out


def temporal_lag(prev_xy, prev_y, query_xy, gmean):
    """Nearest prior-timestep response via NN lookup; global mean if none."""
    if prev_xy is None or len(prev_xy) == 0:
        return np.full(query_xy.shape[0], gmean)
    tree = cKDTree(prev_xy)
    _, idx = tree.query(query_xy, k=1)
    return prev_y[idx]


def edf(train_xy, query_xy, anchors):
    return cdist(query_xy, anchors)


# ----------------------------------------------------------------------------- #
#  Fit + evaluate one scenario replicate
# ----------------------------------------------------------------------------- #
def fit_eval(world, n, design, detection, n_time, rng, methods=None,
             include_geo=False, rf_params=None, return_proba=False):
    """Draw one replicate; fit methods; return dict of test-set AUCs."""
    if methods is None:
        methods = ["GLM", "GAM", "RF", "SpatialRF", "stRF"]
    g = GRID
    t_eval = n_time - 1                       # evaluate at last timestep

    # revisited sampling network
    site_idx = sample_sites(n, design, rng)
    site_xy = world.coords(site_idx)

    # observed labels across timesteps (imperfect detection -> false negatives)
    obs_y = {}
    for t in range(n_time):
        yt = world.y[t].ravel()[site_idx].copy()
        if detection < 1.0:
            keep = rng.random(len(yt)) < detection
            yt = yt * keep
        obs_y[t] = yt

    # ------- held-out test cells at evaluation timestep (perfect detection) ----
    pool = np.setdiff1d(np.arange(g * g), site_idx, assume_unique=False)
    test_idx = rng.choice(pool, size=min(N_TEST, len(pool)), replace=False)
    test_xy = world.coords(test_idx)
    y_test = world.y[t_eval].ravel()[test_idx]
    if y_test.sum() == 0 or y_test.sum() == len(y_test):
        return None                            # degenerate; skip

    # ------- design matrices at evaluation timestep ---------------------------
    Xc_tr = world.covariates(site_idx, t_eval)
    Xc_te = world.covariates(test_idx, t_eval)
    y_tr = obs_y[t_eval]
    if y_tr.sum() == 0 or y_tr.sum() == len(y_tr):
        return None

    res = {}

    def _rf():
        kw = dict(n_estimators=RF_TREES, n_jobs=1, random_state=int(rng.integers(1e9)))
        if rf_params:
            kw.update(rf_params)
        return RandomForestClassifier(**kw)

    def safe_auc(p):
        try:
            return roc_auc_score(y_test, p)
        except Exception:
            return np.nan

    # GLM  (non-spatial baseline: covariates only)
    if "GLM" in methods:
        m = LogisticRegression(max_iter=400, C=1.0)
        m.fit(Xc_tr, y_tr)
        _p = m.predict_proba(Xc_te)[:, 1]; res["GLM"] = safe_auc(_p)
        if return_proba: res["GLM_p"] = _p

    # GAM  (penalised splines + spatial tensor on coords)
    if "GAM" in methods:
        try:
            import os as _os, sys as _sys, contextlib as _ctx
            from pygam import LogisticGAM, s, te
            xy_tr = world.coords(site_idx); xy_te = test_xy
            Xg_tr = np.column_stack([xy_tr, Xc_tr[:, [0, 2, 3, 4]]])
            Xg_te = np.column_stack([xy_te, Xc_te[:, [0, 2, 3, 4]]])
            gam = LogisticGAM(te(0, 1) + s(2) + s(3) + s(4) + s(5),
                              max_iter=25, verbose=False)
            with open(_os.devnull, "w") as _dn, _ctx.redirect_stdout(_dn), \
                    _ctx.redirect_stderr(_dn):
                gam.fit(Xg_tr, y_tr)
                proba = gam.predict_proba(Xg_te)
            res["GAM"] = safe_auc(proba)
            if return_proba: res["GAM_p"] = proba
        except Exception:
            res["GAM"] = np.nan

    # RF (standard, covariates only)
    if "RF" in methods:
        rf = _rf()
        rf.fit(Xc_tr, y_tr)
        _p = rf.predict_proba(Xc_te)[:, 1]; res["RF"] = safe_auc(_p)
        if return_proba: res["RF_p"] = _p

    # Spatial RF (covariates + coordinates)
    if "SpatialRF" in methods:
        Xs_tr = np.column_stack([Xc_tr, world.coords(site_idx)])
        Xs_te = np.column_stack([Xc_te, test_xy])
        rf = _rf()
        rf.fit(Xs_tr, y_tr)
        _p = rf.predict_proba(Xs_te)[:, 1]; res["SpatialRF"] = safe_auc(_p)
        if return_proba: res["SpatialRF_p"] = _p

    # stRF (covariates + multi-scale spatial lags + temporal lag + EDF)
    if "stRF" in methods:
        # spatial lags from observed labels at the evaluation timestep
        lag_tr = spatial_lags(site_xy, y_tr.astype(float), site_xy,
                              RADII_KM, is_train=True)
        lag_te = spatial_lags(site_xy, y_tr.astype(float), test_xy,
                              RADII_KM, is_train=False)
        # temporal lag from previous timestep (t-1) AND mean over all prior
        # timesteps (a site-history feature whose reliability grows with the
        # number of timesteps available -> the timesteps gradient in Fig 3).
        gmean = y_tr.mean()
        if t_eval >= 1:
            prev_y = obs_y[t_eval - 1].astype(float)
            tl_tr = temporal_lag(site_xy, prev_y, site_xy, gmean)
            tl_te = temporal_lag(site_xy, prev_y, test_xy, gmean)
            # site-history mean across every prior timestep at the SAME network site
            hist = np.stack([obs_y[tt].astype(float) for tt in range(t_eval)], axis=0)
            hist_site = hist.mean(0)                     # length n (per site)
            # map history mean onto query points by nearest training site
            tree = cKDTree(site_xy)
            _, nn_tr = tree.query(site_xy, k=1)
            _, nn_te = tree.query(test_xy, k=1)
            hm_tr = hist_site[nn_tr]
            hm_te = hist_site[nn_te]
        else:
            tl_tr = np.full(len(site_xy), gmean); tl_te = np.full(len(test_xy), gmean)
            hm_tr = np.full(len(site_xy), gmean); hm_te = np.full(len(test_xy), gmean)
        # EDF anchors via spatial stratification (k-means on training coords)
        k = min(N_ANCHORS, len(site_xy))
        anchors = KMeans(n_clusters=k, n_init=3,
                         random_state=0).fit(site_xy).cluster_centers_
        edf_tr = edf(site_xy, site_xy, anchors)
        edf_te = edf(site_xy, test_xy, anchors)

        Xt_tr = np.column_stack([Xc_tr, lag_tr, tl_tr[:, None], hm_tr[:, None], edf_tr])
        Xt_te = np.column_stack([Xc_te, lag_te, tl_te[:, None], hm_te[:, None], edf_te])
        rf = _rf()
        rf.fit(Xt_tr, y_tr)
        _p = rf.predict_proba(Xt_te)[:, 1]; res["stRF"] = safe_auc(_p)
        if return_proba: res["stRF_p"] = _p

    # Geostatistical comparator (Matern GP classifier; sdmTMB-class)
    if include_geo:
        import time as _t
        from sklearn.gaussian_process import GaussianProcessClassifier
        from sklearn.gaussian_process.kernels import Matern, ConstantKernel, WhiteKernel
        try:
            xy_tr = world.coords(site_idx)
            kern = (ConstantKernel(1.0)
                    * Matern(length_scale=50.0, nu=1.5)
                    + WhiteKernel(1e-2))
            Xgp_tr = np.column_stack([xy_tr, Xc_tr])
            Xgp_te = np.column_stack([test_xy, Xc_te])
            gp = GaussianProcessClassifier(kernel=kern, n_restarts_optimizer=0,
                                           max_iter_predict=100,
                                           random_state=0)
            t0 = _t.time(); gp.fit(Xgp_tr, y_tr); ft = _t.time() - t0
            res["GEO"] = safe_auc(gp.predict_proba(Xgp_te)[:, 1])
            res["GEO_time"] = ft
            res["GEO_conv"] = 1
        except Exception:
            res["GEO"] = np.nan; res["GEO_conv"] = 0; res["GEO_time"] = np.nan

    if return_proba:
        res["_ytest"] = y_test
    return res


# ----------------------------------------------------------------------------- #
#  World registry
# ----------------------------------------------------------------------------- #
def make_worlds(base_seed=20260601):
    """Return the 6 ground-truth worlds: {stationary, nonstationary} x {30,80,200} km."""
    worlds = {}
    wid = 0
    for stat in (True, False):
        for rng_km in (30, 80, 200):
            name = f"{'stat' if stat else 'nonstat'}_r{rng_km}"
            worlds[name] = World(stationary=stat, range_km=rng_km,
                                 n_time=5, seed=base_seed + wid)
            wid += 1
    return worlds


if __name__ == "__main__":
    # quick self-test / calibration on one world
    import time
    w = World(stationary=False, range_km=80, n_time=5, seed=1)
    print("mean occurrence prob (t0):", round(float(w.p[0].mean()), 3),
          "| prevalence:", round(float(w.y[0].mean()), 3))
    rng = np.random.default_rng(7)
    t0 = time.time()
    r = fit_eval(w, n=200, design="random", detection=1.0, n_time=5, rng=rng)
    print("AUCs:", {k: round(v, 3) for k, v in r.items()})
    print("one replicate (5 methods):", round(time.time() - t0, 2), "s")
