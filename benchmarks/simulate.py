"""Ground-truth simulation for benchmarking AR(2) |lambda| against rhythm-detection methods.

Writes, for each (scenario, n_timepoints) cell, a CSV in MetaCycle input format
(first column = series id, remaining columns = timepoints) plus a truth table.

Scenarios
---------
ar2_<lam>   damped AR(2) oscillator, 24 h eigenperiod, known |lambda|  -> rhythmic
cosinor     sustained 24 h cosine + white noise                       -> rhythmic
damped      cosine with exponentially decaying amplitude              -> rhythmic
ultradian   sustained 12 h cosine + white noise                       -> rhythmic (non-24 h)
white       white noise                                               -> not rhythmic
red         AR(1) red noise, phi = 0.8 (no oscillation)                -> not rhythmic
"""
import os
import numpy as np
import pandas as pd

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sim")
os.makedirs(OUT, exist_ok=True)

N_REP = 200
SAMPLING_H = 2.0          # sampling interval in hours
PERIOD_H = 24.0
NOISE_SD = 1.0
AMPLITUDE = 2.0           # signal amplitude for deterministic-signal scenarios
TRUE_LAMBDAS = [0.30, 0.50, 0.618, 0.75, 0.90]
N_TIMEPOINTS = [12, 18, 24, 48]


def ar2_from_lambda_period(lam: float, period_points: float):
    """AR(2) coefficients with complex roots of modulus lam and given period."""
    omega = 2 * np.pi / period_points
    phi1 = 2 * lam * np.cos(omega)
    phi2 = -(lam ** 2)
    return phi1, phi2


def sim_ar2(rng, n, lam, burn=200):
    period_points = PERIOD_H / SAMPLING_H
    phi1, phi2 = ar2_from_lambda_period(lam, period_points)
    x = np.zeros(n + burn)
    eps = rng.normal(0, NOISE_SD, n + burn)
    for t in range(2, n + burn):
        x[t] = phi1 * x[t - 1] + phi2 * x[t - 2] + eps[t]
    return x[burn:]


def sim_cosinor(rng, n, period_h=PERIOD_H, decay=0.0):
    t = np.arange(n) * SAMPLING_H
    phase = rng.uniform(0, 2 * np.pi)
    env = np.exp(-decay * t)
    return AMPLITUDE * env * np.cos(2 * np.pi * t / period_h + phase) + rng.normal(0, NOISE_SD, n)


def sim_white(rng, n):
    return rng.normal(0, NOISE_SD, n)


def sim_red(rng, n, phi=0.8, burn=200):
    x = np.zeros(n + burn)
    eps = rng.normal(0, NOISE_SD, n + burn)
    for t in range(1, n + burn):
        x[t] = phi * x[t - 1] + eps[t]
    return x[burn:]


def build():
    rng = np.random.default_rng(20260828)
    truth_rows = []
    for n in N_TIMEPOINTS:
        rows, ids = [], []
        for lam in TRUE_LAMBDAS:
            for r in range(N_REP):
                sid = f"ar2_{lam:.3f}_{r}"
                rows.append(sim_ar2(rng, n, lam))
                ids.append(sid)
                truth_rows.append(dict(n=n, id=sid, scenario=f"ar2_{lam:.3f}",
                                       rhythmic=1, true_lambda=lam, period_h=PERIOD_H))
        for scenario, fn, rhythmic, period in [
            ("cosinor", lambda: sim_cosinor(rng, n), 1, PERIOD_H),
            ("damped", lambda: sim_cosinor(rng, n, decay=0.02), 1, PERIOD_H),
            ("ultradian", lambda: sim_cosinor(rng, n, period_h=12.0), 1, 12.0),
            ("white", lambda: sim_white(rng, n), 0, np.nan),
            ("red", lambda: sim_red(rng, n), 0, np.nan),
        ]:
            for r in range(N_REP):
                sid = f"{scenario}_{r}"
                rows.append(fn())
                ids.append(sid)
                truth_rows.append(dict(n=n, id=sid, scenario=scenario, rhythmic=rhythmic,
                                       true_lambda=np.nan, period_h=period))

        times = np.arange(n) * SAMPLING_H
        df = pd.DataFrame(rows, columns=[f"CT{t:g}" for t in times])
        df.insert(0, "series_id", ids)
        df.to_csv(os.path.join(OUT, f"sim_n{n}.csv"), index=False)
        print(f"n={n}: {df.shape[0]} series x {n} timepoints")

    pd.DataFrame(truth_rows).to_csv(os.path.join(OUT, "truth.csv"), index=False)


if __name__ == "__main__":
    build()
