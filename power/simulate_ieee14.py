"""Online IEEE 14-bus comparison using the offline Hankel dictionary.

Run from the repository root, after ``python power/prepare_ieee14.py``:

    python power/simulate_ieee14.py
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ddmhe import load_hankel, make_dd_mhe_pair, run_windowed  # noqa: E402
from power.ieee14_model import load_plant  # noqa: E402
from power.prepare_ieee14 import L, OFFLINE_PATH, R_WEIGHT, RHO, P_WEIGHT, TS  # noqa: E402

NOISE_STD = 0.1
ATTACK_BIAS = 40.0
U_SCALE = 0.02
T_FINAL = 6.0
T_ATTACK = 1.5


def error_metrics(x, x_hat):
    err = norm(x - x_hat, axis=1)
    return float(np.sqrt(np.mean(err**2))), float(np.max(err))


def attack_matrix(n_steps, m, support, attack_index, bias=ATTACK_BIAS):
    attack = np.zeros((n_steps, m))
    if support:
        attack[attack_index:, list(support)] = bias
    return attack


def plot_results(t, x, hats, attack_time, out_path, title):
    idx = [0, 1, 5, 6]
    labels = [r"$\delta_1$", r"$\delta_2$", r"$\omega_1$", r"$\omega_2$"]
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 4.6), sharex=True)
    colors = {"Robust DD-MHE": "#d97706", "Resilient DD-MHE": "#15803d"}
    for ax, j, lab in zip(axes.ravel(), idx, labels):
        ax.plot(t, x[:, j], color="#1d4ed8", lw=1.4, label="true")
        for name, hat in hats.items():
            ax.plot(t, hat[:, j], color=colors[name], ls="--", lw=1.2, label=name)
        ax.axvline(attack_time, color="0.35", ls=":", lw=0.8)
        ax.set_ylabel(lab)
        ax.grid(True, alpha=0.25)
    axes[0, 0].legend(fontsize=7, loc="best")
    axes[1, 0].set_xlabel("t [s]")
    axes[1, 1].set_xlabel("t [s]")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = out_path.with_name(out_path.stem + ".writing" + out_path.suffix)
    fig.savefig(temporary, dpi=160)
    os.replace(temporary, out_path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-k", type=int, default=5)
    parser.add_argument(
        "--support",
        choices=("power", "omega", "mixed"),
        default="power",
        help="Attack bus-power sensors (default), generator speeds, or a mixed set.",
    )
    args = parser.parse_args()

    _u_d, _y_d, _x_d, hu, hy, hx, meta = load_hankel(OFFLINE_PATH)
    horizon = int(meta["horizon"])
    scale = np.asarray(meta["scale"], dtype=float)
    plant = load_plant()
    A, B, C, D = plant["A"], plant["B"], plant["C"], plant["D"]
    n, p, m = plant["n"], plant["p"], plant["m"]
    n_gen = plant["n_gen"]
    if horizon != L:
        raise RuntimeError(f"Offline horizon {horizon} does not match L={L}. Rebuild the dictionary.")

    n_steps = int(round(T_FINAL / TS))
    attack_index = int(round(T_ATTACK / TS))
    if args.support == "power":
        base_support = list(range(n_gen, m))
    elif args.support == "omega":
        base_support = list(range(n_gen))
    else:
        base_support = list(range(0, m, max(m // 5, 1)))[:5]
    base_support = base_support[: max(args.max_k, 1)]

    rng = np.random.default_rng(1)
    u = rng.uniform(-U_SCALE, U_SCALE, size=(n_steps, p))
    x = np.zeros((n_steps, n))
    x[0] = rng.normal(0.0, 0.01, size=n)
    for t in range(n_steps - 1):
        x[t + 1] = A @ x[t] + B @ u[t]
    y_clean = (x @ C.T + u @ D.T) / scale
    noise = rng.normal(0.0, NOISE_STD, size=y_clean.shape)
    prior0 = np.zeros(n)

    solvers = make_dd_mhe_pair(
        hu, hy, hx, horizon=horizon, p_weight=P_WEIGHT, r_weight=R_WEIGHT, rho=RHO
    )
    estimators = {name: solver.estimate for name, solver in solvers.items()}
    print(f"attack channels ({args.support}): {base_support}", flush=True)

    rows = []
    figure_hats = None
    figure_k = min(3, args.max_k)
    t_axis = np.arange(n_steps) * TS
    for k in range(0, args.max_k + 1):
        support = tuple(base_support[:k])
        y = y_clean + noise + attack_matrix(n_steps, m, support, attack_index)
        metrics = {}
        hats = {}
        t0 = time.time()
        for name, estimate in estimators.items():
            x_hat = run_windowed(estimate, u, y, prior0, horizon)
            metrics[name] = error_metrics(x, x_hat)
            if k == figure_k:
                hats[name] = x_hat
        elapsed = time.time() - t0
        rows.append({"k": k, "support": support, "metrics": metrics, "seconds": elapsed})
        print(
            f"k={k} support={support} "
            + " ".join(f"{name[:3]}RMSE={metrics[name][0]:.3f}" for name in estimators)
            + f" ({elapsed:.1f}s)",
            flush=True,
        )
        if k == figure_k:
            figure_hats = hats

    out_dir = Path(__file__).resolve().parent
    fig_path = out_dir / "ieee14_mhe.png"
    plot_results(
        t_axis,
        x,
        figure_hats,
        attack_index * TS,
        fig_path,
        title=f"IEEE 14: {figure_k} sensor bias(es) from t={attack_index * TS:.1f}s",
    )
    summary = {
        "plant": {"n": n, "p": p, "m": m, "ts": TS, "L": horizon},
        "scale": scale.tolist(),
        "pe_rank": meta["rank"]["input_rank"],
        "attack_bias_scaled": ATTACK_BIAS,
        "rows": [
            {
                "k": row["k"],
                "support": list(row["support"]),
                "rmse": {name: row["metrics"][name][0] for name in estimators},
                "mae": {name: row["metrics"][name][1] for name in estimators},
                "seconds": row["seconds"],
            }
            for row in rows
        ],
        "figure": str(fig_path),
    }
    summary_path = out_dir / "ieee14_results.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"wrote {fig_path}", flush=True)
    print(f"wrote {summary_path}", flush=True)
    print("\nRMSE vs number of attacked sensors")
    print("k " + " ".join(f"{name:>16}" for name in estimators))
    for row in rows:
        print(f"{row['k']} " + " ".join(f"{row['metrics'][name][0]:16.3f}" for name in estimators))


if __name__ == "__main__":
    main()
