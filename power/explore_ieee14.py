"""Resilient vs robust MHE on the linearized IEEE 14-bus plant.

Uses the discrete (A, B, C, D) from power/ieee14/, taken from
ZYblend/Moving-horizon-FDIA (demo_14_bus_system). The plant has n=10,
p=19, m=19, so the sparse-recovery bound k < m/2 allows multi-sensor
attacks. Offline Hankel data come from a persistently exciting rollout of
the linear model; online runs add constant biases on an increasing number
of sensors.

Run from the repository root:

    python3 power/explore_ieee14.py
    python3 power/explore_ieee14.py --quick
    python3 power/explore_ieee14.py --paper

--paper writes image/ieee14_mhe.png and prints the table used in the text.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import cvxpy as cp
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import lstsq, norm

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from power.ieee14_model import load_plant  # noqa: E402

TS = 0.01
L = 6
P_WEIGHT = 1.0
R_WEIGHT = 10.0
RHO = 0.95
NOISE_STD = 0.1
ATTACK_BIAS = 40.0  # on the scaled measurement
U_SCALE = 0.02
N_OFFLINE = 800
HOLD = 2
T_FINAL = 6.0
T_ATTACK = 1.5
PRIOR0 = None  # filled after n is known
COMPRESS_TOL = 1e-8


def hankel(sig, depth):
    n_samples, _ = sig.shape
    cols = n_samples - depth + 1
    return np.column_stack([sig[j : j + depth].reshape(-1) for j in range(cols)])


def compress(hu, hy, hx, tol=COMPRESS_TOL):
    out_u, out_y, out_x = {}, {}, {}
    for depth in hu:
        stacked = np.vstack([hu[depth], hy[depth], hx[depth]])
        _, singular, vt = np.linalg.svd(stacked, full_matrices=False)
        keep = max(int(np.sum(singular > tol * singular[0])), 1)
        basis = vt[:keep].T
        out_u[depth] = hu[depth] @ basis
        out_y[depth] = hy[depth] @ basis
        out_x[depth] = hx[depth] @ basis
    return out_u, out_y, out_x


def offline_trajectory(A, B, C, D, n_samples, seed, u_scale=U_SCALE):
    rng = np.random.default_rng(seed)
    raw = rng.uniform(-u_scale, u_scale, size=((n_samples + HOLD - 1) // HOLD, B.shape[1]))
    u = np.repeat(raw, HOLD, axis=0)[:n_samples]
    x = np.zeros((n_samples, A.shape[0]))
    for t in range(n_samples - 1):
        x[t + 1] = A @ x[t] + B @ u[t]
    y = x @ C.T + u @ D.T
    return u, x, y


def measurement_scale(y):
    """One scale per channel from the offline trajectory."""
    scale = np.std(y, axis=0)
    scale = np.maximum(scale, 0.05 * np.median(scale[scale > 0]))
    scale = np.maximum(scale, 1e-3)
    return scale


def prepare_hankel(u_d, y_d, x_d, horizon):
    hu, hy, hx = {}, {}, {}
    for depth in range(1, horizon + 1):
        hu[depth] = hankel(u_d, depth)
        hy[depth] = hankel(y_d, depth)
        hx[depth] = hankel(x_d, depth)
    return compress(hu, hy, hx)


def _solve(problem, label):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Solution may be inaccurate")
        try:
            problem.solve(solver=cp.CLARABEL, warm_start=True, verbose=False)
        except cp.SolverError:
            problem.solve(solver=cp.CLARABEL, warm_start=False, verbose=False)
        if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
            try:
                problem.solve(solver=cp.SCS, warm_start=False, verbose=False, eps=1e-5, max_iters=8000)
            except cp.SolverError as exc:
                raise RuntimeError(f"{label} solve failed: {problem.status}") from exc
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise RuntimeError(f"{label} solve failed: {problem.status}")


class DataDrivenSolver:
    def __init__(self, hu, hy, hx, n, m, p, horizon, resilient=True):
        self.resilient = resilient
        self.hx = hx
        self.problems = {}
        for depth in range(1, horizon + 1):
            cols = hu[depth].shape[1]
            alpha = cp.Variable(cols)
            u_par = cp.Parameter(p * depth)
            y_par = cp.Parameter(m * depth)
            prior = cp.Parameter(n)
            x0 = hx[depth][:n] @ alpha
            residual = y_par - hy[depth] @ alpha
            tau = cp.Variable(nonneg=True)
            constraints = [hu[depth] @ alpha == u_par, cp.norm(x0 - prior, 2) <= tau]
            arrival = (P_WEIGHT * RHO ** (depth - 1)) * tau
            if resilient:
                slack = cp.Variable(m * depth, nonneg=True)
                weights = np.repeat(
                    [R_WEIGHT * RHO ** (depth - 1 - i) for i in range(depth)], m
                )
                constraints.append(cp.abs(residual) <= slack)
                cost = arrival + slack @ weights
            else:
                weights = np.array(
                    [R_WEIGHT * RHO ** (depth - 1 - i) for i in range(depth)]
                )
                stage = 0
                for i in range(depth):
                    sl = slice(i * m, (i + 1) * m)
                    stage += weights[i] * cp.sum_squares(residual[sl])
                cost = arrival + stage
            problem = cp.Problem(cp.Minimize(cost), constraints)
            self.problems[depth] = (problem, alpha, u_par, y_par, prior, hx[depth])

    def estimate(self, u_win, y_win, prior):
        depth = u_win.shape[0]
        problem, alpha, u_par, y_par, prior_par, hx = self.problems[depth]
        u_par.value = u_win.ravel()
        y_par.value = y_win.reshape(-1)
        prior_par.value = prior
        _solve(problem, "data-driven")
        return (hx @ alpha.value).reshape(depth, -1)


class ModelSolver:
    def __init__(self, A, B, C, D, horizon, resilient=True):
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.resilient = resilient
        self.powers = [np.eye(A.shape[0])]
        for _ in range(horizon - 1):
            self.powers.append(A @ self.powers[-1])
        self.problems = {}
        n, m = A.shape[0], C.shape[0]
        for depth in range(1, horizon + 1):
            x0 = cp.Variable(n)
            tau = cp.Variable(nonneg=True)
            g_par = cp.Parameter((depth, n))
            y_par = cp.Parameter((depth, m))
            prior = cp.Parameter(n)
            constraints = [cp.norm(x0 - prior, 2) <= tau]
            arrival = (P_WEIGHT * RHO ** (depth - 1)) * tau
            stage = 0
            if resilient:
                slack = cp.Variable((depth, m), nonneg=True)
                for i in range(depth):
                    state = self.powers[i] @ x0 + g_par[i]
                    weight = R_WEIGHT * RHO ** (depth - 1 - i)
                    constraints.append(cp.abs(y_par[i] - C @ state) <= slack[i])
                    stage += weight * cp.sum(slack[i])
            else:
                for i in range(depth):
                    state = self.powers[i] @ x0 + g_par[i]
                    weight = R_WEIGHT * RHO ** (depth - 1 - i)
                    stage += weight * cp.sum_squares(y_par[i] - C @ state)
            problem = cp.Problem(cp.Minimize(arrival + stage), constraints)
            self.problems[depth] = (problem, x0, g_par, y_par, prior)

    def estimate(self, u_win, y_win, prior):
        depth = u_win.shape[0]
        problem, x0, g_par, y_par, prior_par = self.problems[depth]
        g = np.zeros((depth, self.A.shape[0]))
        acc = np.zeros(self.A.shape[0])
        for i in range(depth):
            g[i] = acc
            acc = self.A @ acc + self.B @ u_win[i]
        g_par.value = g
        y_par.value = y_win - u_win @ self.D.T
        prior_par.value = prior
        _solve(problem, "model")
        return np.vstack([self.powers[i] @ x0.value + g[i] for i in range(depth)])


def run_estimator(estimate, u, y, prior0, horizon):
    n = prior0.shape[0]
    x_hat = np.zeros((u.shape[0], n))
    latest = {0: prior0.copy()}
    for t in range(u.shape[0]):
        if t + 1 < horizon:
            t0, prior = 0, prior0
        else:
            t0, prior = t - horizon + 1, latest[t0]
        window = estimate(u[t0 : t + 1], y[t0 : t + 1], prior)
        for k, time_index in enumerate(range(t0, t + 1)):
            latest[time_index] = window[k]
        x_hat[t] = window[-1]
    return x_hat


def error_metrics(x, x_hat):
    err = norm(x - x_hat, axis=1)
    return float(np.sqrt(np.mean(err**2))), float(np.max(err))


def attack_matrix(n_steps, m, support, attack_index, bias=ATTACK_BIAS):
    attack = np.zeros((n_steps, m))
    if support:
        attack[attack_index:, list(support)] = bias
    return attack


def plot_results(t, x, hats, attack_time, out_path, title):
    # Show two generator angles and two speeds.
    idx = [0, 1, 5, 6]
    labels = [r"$\delta_1$", r"$\delta_2$", r"$\omega_1$", r"$\omega_2$"]
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 4.6), sharex=True)
    colors = {
        "Robust DD-MHE": "#d97706",
        "Resilient DD-MHE": "#15803d",
    }
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
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--paper", action="store_true")
    parser.add_argument("--max-k", type=int, default=5)
    parser.add_argument(
        "--with-md",
        action="store_true",
        help="Also run model-driven estimators (matches DD when the Hankel is PE).",
    )
    parser.add_argument(
        "--support",
        choices=("power", "omega", "mixed"),
        default="power",
        help="Attack bus-power sensors (default), generator speeds, or a mixed set.",
    )
    args = parser.parse_args()

    plant = load_plant()
    A, B, C, D = plant["A"], plant["B"], plant["C"], plant["D"]
    n, p, m = plant["n"], plant["p"], plant["m"]
    n_gen = plant["n_gen"]
    horizon = 4 if args.quick else L
    # Willems needs N >= (p+1)(L+n)-1; keep a margin above that floor.
    pe_floor = (p + 1) * (horizon + n)
    n_offline = max(320 if args.quick else N_OFFLINE, pe_floor + 40)
    t_final = 2.5 if args.quick else T_FINAL
    n_steps = int(round(t_final / TS))
    attack_index = int(round((0.8 if args.quick else T_ATTACK) / TS))
    max_k = min(args.max_k, 3 if args.quick else args.max_k)
    if args.support == "power":
        base_support = list(range(n_gen, m))  # net bus-power channels
    elif args.support == "omega":
        base_support = list(range(n_gen))
    else:
        base_support = list(range(0, m, max(m // 5, 1)))[:5]

    print(
        f"IEEE 14 linearized plant: n={n}, p={p}, m={m}, L={horizon}, "
        f"N_offline={n_offline}, steps={n_steps}, max_k={max_k}",
        flush=True,
    )

    u_d, x_d, y_d_phys = offline_trajectory(A, B, C, D, n_offline, seed=0)
    scale = measurement_scale(y_d_phys)
    y_d = y_d_phys / scale
    # Estimators see scaled outputs; C_s = S^{-1} C, D_s = S^{-1} D.
    Sinv = np.diag(1.0 / scale)
    C_s = Sinv @ C
    D_s = Sinv @ D

    hu, hy, hx = prepare_hankel(u_d, y_d, x_d, horizon)
    pe_rank = int(np.linalg.matrix_rank(hankel(u_d, horizon + n)))
    print(
        f"input Hankel rank (order L+n={horizon + n}): {pe_rank} / {(horizon + n) * p}",
        flush=True,
    )
    print(f"compressed columns at L: {hu[horizon].shape[1]}", flush=True)

    rng = np.random.default_rng(1)
    u = rng.uniform(-U_SCALE, U_SCALE, size=(n_steps, p))
    x = np.zeros((n_steps, n))
    x[0] = rng.normal(0.0, 0.01, size=n)
    for t in range(n_steps - 1):
        x[t + 1] = A @ x[t] + B @ u[t]
    y_clean = (x @ C.T + u @ D.T) / scale
    noise = rng.normal(0.0, NOISE_STD, size=y_clean.shape)
    prior0 = np.zeros(n)

    resilient_dd = DataDrivenSolver(hu, hy, hx, n, m, p, horizon, resilient=True)
    robust_dd = DataDrivenSolver(hu, hy, hx, n, m, p, horizon, resilient=False)
    estimators = {
        "Robust DD-MHE": robust_dd.estimate,
        "Resilient DD-MHE": resilient_dd.estimate,
    }
    # Model-driven pair is optional: with a PE Hankel it matches DD at the
    # printed precision, and doubles the wall-clock cost.
    if args.with_md:
        resilient_md = ModelSolver(A, B, C_s, D_s, horizon, resilient=True)
        robust_md = ModelSolver(A, B, C_s, D_s, horizon, resilient=False)
        estimators = {
            "Robust MD-MHE": robust_md.estimate,
            "Robust DD-MHE": robust_dd.estimate,
            "Resilient MD-MHE": resilient_md.estimate,
            "Resilient DD-MHE": resilient_dd.estimate,
        }

    # Nested supports so columns of the table are comparable.
    base_support = base_support[: max(max_k, 1)]
    print(f"attack channels ({args.support}): {base_support}", flush=True)
    rows = []
    figure_hats = None
    figure_k = min(3, max_k)
    t = np.arange(n_steps) * TS

    for k in range(0, max_k + 1):
        support = tuple(base_support[:k])
        y = y_clean + noise + attack_matrix(n_steps, m, support, attack_index)
        metrics = {}
        hats = {}
        t0 = time.time()
        for name, estimate in estimators.items():
            x_hat = run_estimator(estimate, u, y, prior0, horizon)
            metrics[name] = error_metrics(x, x_hat)
            if k == figure_k and name in ("Robust DD-MHE", "Resilient DD-MHE"):
                hats[name] = x_hat
        elapsed = time.time() - t0
        rows.append({"k": k, "support": support, "metrics": metrics, "seconds": elapsed})
        print(
            f"k={k} support={support} "
            + " ".join(
                f"{name[:3]}RMSE={metrics[name][0]:.3f}" for name in estimators
            )
            + f" ({elapsed:.1f}s)",
            flush=True,
        )
        if k == figure_k:
            figure_hats = hats

    out_dir = Path(__file__).resolve().parent
    fig_path = out_dir / "ieee14_mhe.png"
    if args.paper:
        fig_path = ROOT / "image" / "ieee14_mhe.png"
    plot_results(
        t,
        x,
        figure_hats,
        attack_index * TS,
        fig_path,
        title=f"IEEE 14: {figure_k} sensor bias(es) from t={attack_index * TS:.1f}s",
    )

    summary = {
        "plant": {"n": n, "p": p, "m": m, "ts": TS, "L": horizon},
        "scale": scale.tolist(),
        "pe_rank": pe_rank,
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

    # Compact markdown-style table for the paper draft.
    print("\nRMSE vs number of attacked sensors")
    header = "k " + " ".join(f"{name:>16}" for name in estimators)
    print(header)
    for row in rows:
        print(
            f"{row['k']} "
            + " ".join(f"{row['metrics'][name][0]:16.3f}" for name in estimators)
        )


if __name__ == "__main__":
    main()
