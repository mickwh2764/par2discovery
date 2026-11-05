"""Known-answer check of the timescale tau under measurement noise.

Simulates series whose true decay time constant is known, fits them with
par2.fit_ar2(..., sampling_interval=dt) and compares tau with the truth.
Writes results/noise_bias.csv. Run from the repository root:

    python benchmarks/noise_bias.py
"""

import csv
import math
import os
import sys
import warnings

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from par2.core import fit_ar2  # noqa: E402

PERIOD = 24.0
OMEGA = 2 * math.pi / PERIOD
REPS = 100
rng = np.random.default_rng(20261105)


def damped_wave(tau, dt, duration, noise):
    """Deterministic damped cosine plus white measurement noise (explant-like)."""
    t = np.arange(0, duration, dt)
    return np.exp(-t / tau) * np.cos(OMEGA * t) + rng.normal(0, noise, t.size)


def noisy_oscillator(tau, dt, duration, noise):
    """Noise-driven damped oscillator, unit variance, plus white measurement noise."""
    n = int(duration / dt)
    r = math.exp(-dt / tau)
    c, s = math.cos(OMEGA * dt), math.sin(OMEGA * dt)
    rot = r * np.array([[c, -s], [s, c]])
    q = math.sqrt(1 - r * r)
    z = rng.normal(size=2)
    x = np.empty(n)
    for i in range(n):
        z = rot @ z + q * rng.normal(size=2)
        x[i] = z[0]
    return x / x.std() + rng.normal(0, noise, n)


SCENARIOS = [
    ("damped wave, 10 min, no noise", damped_wave, 1 / 6, lambda tau: min(3 * tau, 168), 0.0, [24, 48, 72, 120]),
    ("damped wave, 10 min, 5% noise", damped_wave, 1 / 6, lambda tau: min(3 * tau, 168), 0.05, [24, 48, 72, 120]),
    ("damped wave, 10 min, 20% noise", damped_wave, 1 / 6, lambda tau: min(3 * tau, 168), 0.20, [24, 48, 72, 120]),
    ("noisy oscillator, 2 h x 48 h, 30% noise", noisy_oscillator, 2.0, lambda tau: 48, 0.30, [6, 12, 24, 48, 96]),
]


def main():
    rows = []
    for name, generate, dt, duration_for, noise, taus in SCENARIOS:
        truth, estimates = [], []
        for tau in taus:
            duration = duration_for(tau)
            est = []
            for _ in range(REPS):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    est.append(fit_ar2(generate(tau, dt, duration, noise), sampling_interval=dt)["timescale"])
            est = np.array(est, dtype=float)
            q25, med, q75 = np.percentile(est, [25, 50, 75])
            rows.append([name, dt, duration, noise, tau, round(med, 2), round(q25, 2), round(q75, 2), round(med / tau, 3)])
            truth += [tau] * REPS
            estimates += list(est)
        rho = spearmanr(truth, estimates)[0]
        for row in rows[-len(taus):]:
            row.append(round(rho, 3))
        print(f"{name}: rank correlation {rho:.2f}")
    out = os.path.join(os.path.dirname(__file__), "results", "noise_bias.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["scenario", "interval_h", "duration_h", "noise_sd", "true_tau_h",
                    "median_tau_h", "q25_tau_h", "q75_tau_h", "ratio_to_truth", "spearman_rho"])
        w.writerows(rows)
    print("wrote", out)


if __name__ == "__main__":
    main()
