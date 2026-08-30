"""Benchmark figures: calibrated power by scenario, and |lambda| estimator behaviour."""
import os
import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

METHODS = ["PAR2_support", "PAR2_lambda", "JTK", "LS", "ARS", "RAIN"]
COLOURS = {"PAR2_support": "#0e7490", "PAR2_lambda": "#67e8f9", "JTK": "#b45309",
           "LS": "#a16207", "ARS": "#7c2d12", "RAIN": "#4d7c0f"}
LAMS = [0.300, 0.500, 0.618, 0.750, 0.900]


def fig_power():
    d = pd.read_csv(os.path.join(RES, "calibrated_power.csv"))
    ns = sorted(d.n.unique())
    fig, axes = plt.subplots(1, len(ns), figsize=(4 * len(ns), 3.6), sharey=True)
    for ax, n in zip(axes, ns):
        sub = d[d.n == n]
        for m in METHODS:
            s = sub[sub.method == m].set_index("scenario")["rate"]
            if s.empty:
                continue
            ys = [s.get(f"ar2_{lam:.3f}") for lam in LAMS]
            ax.plot(LAMS, ys, marker="o", ms=4, color=COLOURS[m], label=m)
        ax.axhline(0.05, ls=":", c="grey", lw=1)
        ax.set_title(f"n = {n} timepoints")
        ax.set_xlabel("true |λ| of simulated AR(2) oscillator")
        ax.set_ylim(0, 1.02)
    axes[0].set_ylabel("detection rate at 5% FPR")
    axes[-1].legend(fontsize=7, loc="lower right")
    fig.suptitle("Power to detect a noise-driven 24 h AR(2) oscillator", y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "power_ar2.png"), dpi=160, bbox_inches="tight")


def fig_specificity():
    d = pd.read_csv(os.path.join(RES, "calibrated_power.csv"))
    scen = ["cosinor", "damped", "ultradian", "red"]
    labels = ["sustained\ncosine", "amplitude-decaying\ncosine", "12 h ultradian\n(non-24 h)",
              "AR(1) red noise\n(arrhythmic)"]
    sub = d[(d.n == 24)]
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    width = 0.14
    present = [m for m in METHODS if m in set(sub.method)]
    for i, m in enumerate(present):
        s = sub[sub.method == m].set_index("scenario")["rate"]
        xs = [j + i * width for j in range(len(scen))]
        ax.bar(xs, [s.get(k, 0) for k in scen], width, color=COLOURS[m], label=m)
    ax.set_xticks([j + width * (len(present) - 1) / 2 for j in range(len(scen))])
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("call rate at 5% FPR")
    ax.set_title("n = 24: where the methods differ (last panel is a false-positive rate)", fontsize=10)
    ax.legend(fontsize=7, ncol=3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "specificity_n24.png"), dpi=160, bbox_inches="tight")


def fig_recovery():
    r = pd.read_csv(os.path.join(RES, "recovery.csv"))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.4))
    for n, grp in r.groupby("n"):
        a1.plot(grp.true_lambda, grp.mean_est, marker="o", ms=4, label=f"n = {n}")
        a2.plot(grp.true_lambda, grp.ci_coverage, marker="o", ms=4, label=f"n = {n}")
    a1.plot([0.25, 0.95], [0.25, 0.95], ls="--", c="grey", lw=1, label="unbiased")
    a1.set_xlabel("true |λ|")
    a1.set_ylabel("mean estimated |λ|")
    a1.set_title("Upward bias concentrates at low persistence", fontsize=10)
    a1.legend(fontsize=7)
    a2.axhline(0.95, ls="--", c="grey", lw=1)
    a2.set_xlabel("true |λ|")
    a2.set_ylabel("bootstrap 95% CI coverage")
    a2.set_title("CI under-covers when |λ| is low", fontsize=10)
    a2.set_ylim(0.6, 1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "estimator_behaviour.png"), dpi=160, bbox_inches="tight")


if __name__ == "__main__":
    fig_power()
    fig_specificity()
    fig_recovery()
    print("figures written to", FIG)
