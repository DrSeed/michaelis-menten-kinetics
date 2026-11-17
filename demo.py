# Self-contained demo: 500 noisy kinetics experiments, three ways to estimate Km.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from kinetics import mm_rate, fit_lineweaver_burk, fit_eadie_hofstee, fit_nonlinear

VMAX, KM = 100.0, 2.5
S = np.array([0.25, 0.5, 1, 2, 4, 8, 16, 32])


def experiment(rng, sd=3.0):
    # Constant measurement noise (e.g. absorbance read-out), 3% of Vmax.
    return np.clip(mm_rate(S, VMAX, KM) + rng.normal(0, sd, len(S)), 0.5, None)


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    rng = np.random.default_rng(11)
    est = {"Lineweaver-Burk": [], "Eadie-Hofstee": [], "Nonlinear least squares": []}
    for _ in range(500):
        v = experiment(rng)
        est["Lineweaver-Burk"].append(fit_lineweaver_burk(S, v)[1])
        est["Eadie-Hofstee"].append(fit_eadie_hofstee(S, v)[1])
        est["Nonlinear least squares"].append(fit_nonlinear(S, v, p0=fit_eadie_hofstee(S, v))[1])
    rows = []
    for name, k in est.items():
        k = np.array(k)
        rows.append({"method": name, "true_km": KM, "median_km": round(float(np.median(k)), 3),
                     "bias_pct": round(100 * (np.mean(k) - KM) / KM, 1),
                     "rmse": round(float(np.sqrt(np.mean((k - KM) ** 2))), 3),
                     "iqr": round(float(np.subtract(*np.percentile(k, [75, 25]))), 3)})
    table = pd.DataFrame(rows)
    table.to_csv("results/summary.csv", index=False)

    v = experiment(np.random.default_rng(5))
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    grid = np.linspace(0.01, 34, 300)
    axes[0].scatter(S, v, color="k", zorder=3, label="one experiment")
    for name, fit, c in [("Lineweaver-Burk", fit_lineweaver_burk, "#c0504d"),
                         ("Nonlinear", lambda s, y: fit_nonlinear(s, y)[:2], "#2a7f62")]:
        vm, km = fit(S, v)
        axes[0].plot(grid, mm_rate(grid, vm, km), color=c, label=f"{name}: Km={km:.2f}")
    axes[0].plot(grid, mm_rate(grid, VMAX, KM), "k--", lw=1, label=f"truth: Km={KM}")
    axes[0].set_xlabel("[S] (mM)")
    axes[0].set_ylabel("v")
    axes[0].set_title("Same data, different answers")
    axes[0].legend(fontsize=8)

    reps = np.array([experiment(np.random.default_rng(100 + i)) for i in range(30)])
    axes[1].errorbar(1 / S, (1 / reps).mean(axis=0), yerr=(1 / reps).std(axis=0), fmt="o", color="#c0504d")
    axes[1].set_xlabel("1/[S]")
    axes[1].set_ylabel("1/v")
    axes[1].set_title("Lineweaver-Burk: the same noise explodes at low [S]")

    axes[2].boxplot([est[n] for n in est], showfliers=False)
    axes[2].set_xticks([1, 2, 3], ["Lineweaver-\nBurk", "Eadie-\nHofstee", "Nonlinear"])
    axes[2].axhline(KM, color="k", ls="--", lw=1, label="true Km")
    axes[2].set_ylabel("estimated Km (mM)")
    axes[2].set_title("500 simulated experiments")
    axes[2].legend()
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
