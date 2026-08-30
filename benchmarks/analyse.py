"""Score AR(2) |lambda| against JTK_CYCLE, ARSER, Lomb-Scargle and RAIN on the simulations.

Produces:
  results/lambda_fits.csv      per-series AR(2) fit, bootstrap CI and oscillation support
  results/detection.csv        AUC / sensitivity / specificity per method, per n
  results/recovery.csv         |lambda| bias, RMSE and bootstrap CI coverage
  results/sensitivity_by_scenario.csv
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from par2.core import fit_ar2  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SIM, RES = os.path.join(HERE, "sim"), os.path.join(HERE, "results")
NS = [12, 18, 24, 48]
SAMPLING_H = 2.0
PERIOD_WINDOW = (20.0, 28.0)   # the window JTK/ARS/LS were given
N_BOOT = 200


def bh(p):
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    q = np.full(p.shape, np.nan)
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(pv)
    adj = ranked * m / (np.arange(m) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


def auc(score, positive):
    """AUC with higher score = more likely positive; ties handled by rank averaging."""
    score, positive = np.asarray(score, float), np.asarray(positive, bool)
    ok = ~np.isnan(score)
    score, positive = score[ok], positive[ok]
    order = np.argsort(score)
    ranks = np.empty(len(score), float)
    ranks[order] = np.arange(1, len(score) + 1)
    # average ranks within ties
    df = pd.DataFrame({"s": score, "r": ranks})
    ranks = df.groupby("s")["r"].transform("mean").to_numpy()
    n1, n0 = positive.sum(), (~positive).sum()
    if n1 == 0 or n0 == 0:
        return np.nan
    return (ranks[positive].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def oscillation_support(x, n_boot=N_BOOT, seed=0):
    """Fraction of residual-bootstrap refits with complex roots and period in window.

    This turns the AR(2) fit into a rhythmicity score comparable with the p-value
    based methods: it asks how stable the 'damped oscillator with circadian period'
    conclusion is under resampling, rather than reading it off a single fit.
    """
    x = np.asarray(x, float)
    x = x - x.mean()
    n = len(x)
    base = fit_ar2(x)
    phi1, phi2 = base["phi1"], base["phi2"]
    resid = x[2:] - (phi1 * x[1:-1] + phi2 * x[:-2])
    resid = resid - resid.mean()
    rng = np.random.default_rng(seed)
    hits = 0
    lam_draws = []
    for _ in range(n_boot):
        eps = rng.choice(resid, n - 2, replace=True)
        xb = np.empty(n)
        xb[:2] = x[:2]
        for t in range(2, n):
            xb[t] = phi1 * xb[t - 1] + phi2 * xb[t - 2] + eps[t - 2]
        try:
            f = fit_ar2(xb)
        except Exception:
            continue
        lam_draws.append(f["eigenvalue"])
        per = f["eigenperiod"]
        if f["root_type"] == "Complex" and per is not None:
            per_h = per * SAMPLING_H
            if PERIOD_WINDOW[0] <= per_h <= PERIOD_WINDOW[1]:
                hits += 1
    ci = (np.percentile(lam_draws, [2.5, 97.5]) if lam_draws else [np.nan, np.nan])
    return base, hits / max(n_boot, 1), ci


def fit_all():
    rows = []
    for n in NS:
        dat = pd.read_csv(os.path.join(SIM, f"sim_n{n}.csv"))
        ids = dat["series_id"].to_numpy()
        mat = dat.drop(columns=["series_id"]).to_numpy(float)
        for i, sid in enumerate(ids):
            base, support, ci = oscillation_support(mat[i], seed=i)
            rows.append(dict(n=n, series_id=sid, eigenvalue=base["eigenvalue"],
                             root_type=base["root_type"], r2=base["r2"],
                             eigenperiod_h=(base["eigenperiod"] * SAMPLING_H
                                            if base["eigenperiod"] else np.nan),
                             osc_support=support, lam_lo=ci[0], lam_hi=ci[1]))
        print(f"  fitted n={n}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RES, "lambda_fits.csv"), index=False)
    return df


def load_incumbents(n):
    d = os.path.join(RES, f"meta_n{n}")
    out = {}
    jtk = pd.read_csv(os.path.join(d, f"JTKresult_sim_n{n}.csv"))
    out["JTK"] = jtk.set_index("CycID")["ADJ.P"]
    ls = pd.read_csv(os.path.join(d, f"LSresult_sim_n{n}.csv"))
    out["LS"] = ls.set_index("CycID")["p"] if "p" in ls.columns else ls.set_index("CycID")["ADJ.P"]
    ars_path = os.path.join(d, f"ARSresult_sim_n{n}.csv")
    if os.path.exists(ars_path):
        ars = pd.read_csv(ars_path)
        col = "pvalue" if "pvalue" in ars.columns else ars.columns[-1]
        out["ARS"] = ars.set_index("CycID")[col]
    rain = pd.read_csv(os.path.join(RES, f"rain_n{n}.csv")).set_index("series_id")["rain_p"]
    out["RAIN"] = rain
    return out


def main():
    truth = pd.read_csv(os.path.join(SIM, "truth.csv"))
    fits_path = os.path.join(RES, "lambda_fits.csv")
    fits = pd.read_csv(fits_path) if os.path.exists(fits_path) else fit_all()

    det_rows, sens_rows, rec_rows = [], [], []
    for n in NS:
        t = truth[truth.n == n].set_index("id")
        f = fits[fits.n == n].set_index("series_id")
        inc = load_incumbents(n)

        scores = {m: -np.log10(np.clip(s.reindex(t.index).to_numpy(float), 1e-300, None))
                  for m, s in inc.items()}
        scores["PAR2_support"] = f["osc_support"].reindex(t.index).to_numpy(float)
        scores["PAR2_lambda"] = f["eigenvalue"].reindex(t.index).to_numpy(float)

        pos = t.rhythmic.to_numpy(bool)
        for m, sc in scores.items():
            det_rows.append(dict(n=n, method=m, auc=auc(sc, pos)))

        # sensitivity / specificity per scenario at a matched 5% FDR-style threshold
        for m, s in list(inc.items()):
            q = bh(s.reindex(t.index).to_numpy(float))
            called = q < 0.05
            for scen, grp in t.groupby("scenario"):
                idx = t.index.get_indexer(grp.index)
                sens_rows.append(dict(n=n, method=m, scenario=scen,
                                      rhythmic=int(grp.rhythmic.iloc[0]),
                                      called_rate=float(np.nanmean(called[idx]))))
        # PAR(2): call rhythmic when >=95% of bootstrap refits are circadian oscillators
        called = scores["PAR2_support"] >= 0.95
        for scen, grp in t.groupby("scenario"):
            idx = t.index.get_indexer(grp.index)
            sens_rows.append(dict(n=n, method="PAR2", scenario=scen,
                                  rhythmic=int(grp.rhythmic.iloc[0]),
                                  called_rate=float(np.nanmean(called[idx]))))

        # |lambda| recovery on the AR(2) scenarios
        ar2 = t[t.true_lambda.notna()]
        for lam, grp in ar2.groupby("true_lambda"):
            sub = f.reindex(grp.index)
            est = sub.eigenvalue.to_numpy(float)
            cov = ((sub.lam_lo <= lam) & (sub.lam_hi >= lam)).mean()
            rec_rows.append(dict(n=n, true_lambda=lam, mean_est=est.mean(),
                                 bias=est.mean() - lam, rmse=np.sqrt(np.mean((est - lam) ** 2)),
                                 ci_coverage=float(cov),
                                 mean_ci_width=float((sub.lam_hi - sub.lam_lo).mean())))

    pd.DataFrame(det_rows).to_csv(os.path.join(RES, "detection.csv"), index=False)
    pd.DataFrame(sens_rows).to_csv(os.path.join(RES, "sensitivity_by_scenario.csv"), index=False)
    pd.DataFrame(rec_rows).to_csv(os.path.join(RES, "recovery.csv"), index=False)

    print("\n== AUC, rhythmic vs non-rhythmic ==")
    print(pd.DataFrame(det_rows).pivot(index="method", columns="n", values="auc").round(3))
    print("\n== |lambda| recovery ==")
    print(pd.DataFrame(rec_rows).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
