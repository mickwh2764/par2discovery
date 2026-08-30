"""Compare methods at a matched null: threshold each score so the white-noise
false-positive rate is 5%, then report the call rate in every other scenario.

This is the only fair comparison across a p-value method and a model-fit score,
and it makes the red-noise column directly interpretable as a false-positive rate
against autocorrelated-but-arrhythmic data.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SIM, RES = os.path.join(HERE, "sim"), os.path.join(HERE, "results")
NS = [12, 18, 24, 48]
TARGET_FPR = 0.05

from analyse import load_incumbents  # noqa: E402


def main():
    truth = pd.read_csv(os.path.join(SIM, "truth.csv"))
    fits = pd.read_csv(os.path.join(RES, "lambda_fits.csv"))
    rows = []
    for n in NS:
        t = truth[truth.n == n].set_index("id")
        f = fits[fits.n == n].set_index("series_id")
        scores = {m: -np.log10(np.clip(s.reindex(t.index).to_numpy(float), 1e-300, None))
                  for m, s in load_incumbents(n).items()}
        scores["PAR2_support"] = f["osc_support"].reindex(t.index).to_numpy(float)
        scores["PAR2_lambda"] = f["eigenvalue"].reindex(t.index).to_numpy(float)

        white = (t.scenario == "white").to_numpy()
        for m, sc in scores.items():
            thr = np.nanquantile(sc[white], 1 - TARGET_FPR)
            called = sc > thr
            for scen, grp in t.groupby("scenario"):
                idx = t.index.get_indexer(grp.index)
                rows.append(dict(n=n, method=m, scenario=scen,
                                 rate=float(np.nanmean(called[idx]))))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "calibrated_power.csv"), index=False)
    order = ["ar2_0.300", "ar2_0.500", "ar2_0.618", "ar2_0.750", "ar2_0.900",
             "cosinor", "damped", "ultradian", "red", "white"]
    for n in NS:
        print(f"\n=== n = {n} (threshold calibrated to 5% FPR on white noise) ===")
        p = out[out.n == n].pivot(index="scenario", columns="method", values="rate")
        print(p.reindex(order).round(3).to_string())


if __name__ == "__main__":
    main()
