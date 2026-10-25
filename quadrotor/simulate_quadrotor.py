"""Estimation comparison on the linearized planar quadrotor.

The plant is the safe-control-gym TWO_D Crazyflie, linearized at hover.
Outputs are horizontal position, altitude, pitch, pitch rate, and the
vertical specific force at the center of mass. The default run closes the
loop around Robust DD-MHE and Resilient DD-MHE while the vehicle tracks a
figure-8. A pitch-rate bias starts partway through that loop.

Run from the repository root:
    python3 -m quadrotor.simulate_quadrotor
    python3 -m quadrotor.simulate_quadrotor --open-loop
"""

import argparse
import json
import os
import time
import warnings
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("RAYON_NUM_THREADS", "1")

import cvxpy as cp
import numpy as np
from numpy.linalg import norm
from scipy.linalg import solve_discrete_are

try:
    from quadrotor.quadrotor_plant import (
        ARM,
        GRAVITY,
        IYY,
        MASS,
        TS,
        discrete_plant,
    )
except ImportError:
    from quadrotor_plant import ARM, GRAVITY, IYY, MASS, TS, discrete_plant


def hankel(sig, depth):
    """Block-Hankel matrix with ``depth`` block rows."""
    n_samples, _ = sig.shape
    cols = n_samples - depth + 1
    return np.column_stack([sig[j : j + depth].reshape(-1) for j in range(cols)])


A, B, C_PHYS, D_PHYS = discrete_plant()
N_STATE, N_INPUT = A.shape[0], B.shape[1]
N_MEAS = C_PHYS.shape[0]

# One normalized unit. Sensor noise in the estimator is N(0, 0.1^2).
SCALE = np.array([0.05, 0.05, 0.05, 1.0, 1.5])  # m, m, rad, rad/s, m/s^2
C = C_PHYS / SCALE[:, None]
D = D_PHYS / SCALE[:, None]

L = 20
P_WEIGHT = 1.0
R_WEIGHT = 10.0
RHO = 0.95
NOISE_STD = 0.1
N_OFFLINE = 160
DITHER = 0.015  # Newton, about each rotor
# One vertical-plane figure-8: x = AX sin(wt), z = AZ sin(2wt).
FIG8_PERIOD = 8.0
FIG8_AX = 0.45
FIG8_AZ = 0.22
T_FINAL = 12.0
T_ATTACK = 4.0
N_STEPS = int(round(T_FINAL / TS))
ATTACK_INDEX = int(round(T_ATTACK / TS))
# Physical bias added to one output. The gyro entry is rad/s.
ATTACKS = {
    "none": None,
    "x": (0, 1.0),
    "theta": (2, 0.35),
    "gyro": (3, 8.0),
    "accel": (4, 5.0),
}
PRIOR0 = np.array([0.15, 0.0, -0.08, 0.0, 0.08, 0.0])
X0 = np.zeros(N_STATE)
# Each rotor's deviation is clipped near the hover thrust, mg/2.
U_LIMIT = 0.5 * MASS * GRAVITY

NAMES = (
    "Robust MD-MHE",
    "Robust DD-MHE",
    "Resilient MD-MHE",
    "Resilient DD-MHE",
)


def lqr_gain():
    q = np.diag([6.0, 1.0, 8.0, 1.0, 0.4, 0.02])
    r = np.diag([40.0, 40.0])
    p = solve_discrete_are(A, B, q, r)
    k = np.linalg.solve(r + B.T @ p @ B, B.T @ p @ A)
    return k


def reference(step):
    """Figure-8 state reference, including the velocity and pitch that track it."""
    t = step * TS
    w = 2.0 * np.pi / FIG8_PERIOD
    s1 = np.sin(w * t)
    c1 = np.cos(w * t)
    ref = np.zeros(N_STATE)
    ref[0] = FIG8_AX * s1
    ref[1] = FIG8_AX * w * c1
    ref[2] = FIG8_AZ * np.sin(2.0 * w * t)
    ref[3] = FIG8_AZ * 2.0 * w * np.cos(2.0 * w * t)
    # x_ddot = g * theta on this linearization, and x_ddot = -AX w^2 sin(wt).
    ref[4] = -FIG8_AX * w * w * s1 / GRAVITY
    ref[5] = -FIG8_AX * (w ** 3) * c1 / GRAVITY
    return ref


def figure8_feedforward(step):
    """Thrust deviations that produce the figure-8 acceleration on the linear model."""
    t = step * TS
    w = 2.0 * np.pi / FIG8_PERIOD
    z_ddot = -FIG8_AZ * (2.0 * w) ** 2 * np.sin(2.0 * w * t)
    theta_ddot = FIG8_AX * (w ** 4) * np.sin(w * t) / GRAVITY
    moment = ARM / (IYY * np.sqrt(2.0))
    u_sum = MASS * z_ddot
    u_diff = theta_ddot / moment
    return np.array([0.5 * (u_sum - u_diff), 0.5 * (u_sum + u_diff)])


def simulate(u, x0):
    x = np.zeros((u.shape[0], N_STATE))
    x[0] = x0
    for t in range(u.shape[0] - 1):
        x[t + 1] = A @ x[t] + B @ u[t]
    return x


def output_scaled(x, u):
    return x @ C.T + u @ D.T


def closed_loop(n_steps, gain, seed, dither_scale, x0):
    rng = np.random.default_rng(seed)
    u = np.zeros((n_steps, N_INPUT))
    x = np.zeros((n_steps, N_STATE))
    x[0] = x0
    for t in range(n_steps - 1):
        dither = rng.uniform(-dither_scale, dither_scale, size=N_INPUT)
        u[t] = np.clip(
            -gain @ (x[t] - reference(t)) + figure8_feedforward(t) + dither,
            -U_LIMIT,
            U_LIMIT,
        )
        x[t + 1] = A @ x[t] + B @ u[t]
    u[-1] = np.clip(
        -gain @ (x[-1] - reference(n_steps - 1)) + figure8_feedforward(n_steps - 1),
        -U_LIMIT,
        U_LIMIT,
    )
    return u, x


def offline_trajectory(gain, seed):
    """Hover about the origin with a dithered input. The state is recorded."""
    rng = np.random.default_rng(seed)
    u = np.zeros((N_OFFLINE, N_INPUT))
    x = np.zeros((N_OFFLINE, N_STATE))
    for t in range(N_OFFLINE - 1):
        dither = rng.uniform(-DITHER, DITHER, size=N_INPUT)
        u[t] = -gain @ x[t] + dither
        x[t + 1] = A @ x[t] + B @ u[t]
    u[-1] = -gain @ x[-1]
    y = output_scaled(x, u)
    return u, x, y


def compress_hankel(hu, hy, hx, tol=1e-8):
    """Keep the numerical column space shared by the three Hankel matrices.

    The output window of this plant repeats the pitch kinematics, so the raw
    Hankel matrix has a large numerical nullspace. Restricting the coefficient
    to the joint range leaves the represented trajectories unchanged and keeps
    the cone program well scaled.
    """
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


def prepare_hankel(u_d, y_d, x_d):
    hu, hy, hx = {}, {}, {}
    for depth in range(1, L + 1):
        hu[depth] = hankel(u_d, depth)
        hy[depth] = hankel(y_d, depth)
        hx[depth] = hankel(x_d, depth)
    return hu, hy, hx


def input_excitation_rank(u):
    order = L + N_STATE
    hu = hankel(u, order)
    return int(np.linalg.matrix_rank(hu, tol=1e-8)), hu.shape[0]


def robust_dd(hu, hy, hx, u_win, y_win, prior):
    depth = u_win.shape[0]
    n = prior.shape[0]
    m = y_win.shape[1]
    alpha_dim = hu.shape[1]
    sw = np.sqrt(P_WEIGHT * RHO ** (depth - 1))
    rows = [sw * hx[:n]]
    rhs = [sw * prior]
    for i in range(depth):
        weight = np.sqrt(R_WEIGHT * RHO ** (depth - 1 - i))
        sl = slice(i * m, (i + 1) * m)
        rows.append(-weight * hy[sl])
        rhs.append(-weight * y_win.reshape(-1)[sl])
    design = np.vstack(rows)
    target = np.concatenate(rhs)
    gram = design.T @ design
    kkt = np.block(
        [
            [gram, hu.T],
            [hu, np.zeros((hu.shape[0], hu.shape[0]))],
        ]
    )
    sol = np.linalg.solve(kkt, np.concatenate([design.T @ target, u_win.ravel()]))
    alpha = sol[:alpha_dim]
    return (hx @ alpha).reshape(depth, n)


def robust_md(u_win, y_win, prior, powers):
    depth = u_win.shape[0]
    n = prior.shape[0]
    g = np.zeros((depth, n))
    acc = np.zeros(n)
    for i in range(depth):
        g[i] = acc
        acc = A @ acc + B @ u_win[i]
    y_model = y_win - u_win @ D.T
    sw = np.sqrt(P_WEIGHT * RHO ** (depth - 1))
    rows = [sw * np.eye(n)]
    rhs = [sw * prior]
    for i in range(depth):
        weight = np.sqrt(R_WEIGHT * RHO ** (depth - 1 - i))
        rows.append(-weight * (C @ powers[i]))
        rhs.append(-weight * (y_model[i] - C @ g[i]))
    design = np.vstack(rows)
    target = np.concatenate(rhs)
    x0, *_ = np.linalg.lstsq(design, target, rcond=None)
    return np.vstack([powers[i] @ x0 + g[i] for i in range(depth)])


def _solve(problem, label):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Solution may be inaccurate")
        try:
            problem.solve(solver=cp.CLARABEL, warm_start=True, verbose=False)
            if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
                problem.solve(solver=cp.CLARABEL, warm_start=False, verbose=False)
        except cp.SolverError:
            problem.solve(
                solver=cp.SCS,
                warm_start=False,
                verbose=False,
                eps=1e-5,
                max_iters=8000,
            )
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise RuntimeError(f"{label} solve failed: {problem.status}")
    return problem.status


class DataDrivenSolver:
    def __init__(self, hu, hy, hx):
        self.problems = {}
        n = N_STATE
        m = N_MEAS
        p = N_INPUT
        for depth in range(1, L + 1):
            cols = hu[depth].shape[1]
            alpha = cp.Variable(cols)
            u_par = cp.Parameter(p * depth)
            y_par = cp.Parameter(m * depth)
            prior = cp.Parameter(n)
            x0 = hx[depth][:n] @ alpha
            residual = y_par - hy[depth] @ alpha
            tau = cp.Variable(nonneg=True)
            slack = cp.Variable(m * depth, nonneg=True)
            weights = np.repeat(
                [R_WEIGHT * RHO ** (depth - 1 - i) for i in range(depth)], m
            )
            constraints = [
                hu[depth] @ alpha == u_par,
                cp.norm(x0 - prior, 2) <= tau,
                cp.abs(residual) <= slack,
            ]
            arrival = (P_WEIGHT * RHO ** (depth - 1)) * tau
            problem = cp.Problem(cp.Minimize(arrival + slack @ weights), constraints)
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
    def __init__(self, powers):
        self.powers = powers
        self.problems = {}
        n = N_STATE
        m = N_MEAS
        for depth in range(1, L + 1):
            x0 = cp.Variable(n)
            tau = cp.Variable(nonneg=True)
            g_par = cp.Parameter((depth, n))
            y_par = cp.Parameter((depth, m))
            prior = cp.Parameter(n)
            slack = cp.Variable((depth, m), nonneg=True)
            constraints = [cp.norm(x0 - prior, 2) <= tau]
            stage = 0
            for i in range(depth):
                state = powers[i] @ x0 + g_par[i]
                weight = R_WEIGHT * RHO ** (depth - 1 - i)
                constraints.append(cp.abs(y_par[i] - C @ state) <= slack[i])
                stage += weight * cp.sum(slack[i])
            arrival = (P_WEIGHT * RHO ** (depth - 1)) * tau
            problem = cp.Problem(cp.Minimize(arrival + stage), constraints)
            self.problems[depth] = (problem, x0, g_par, y_par, prior)

    def estimate(self, u_win, y_win, prior):
        depth = u_win.shape[0]
        problem, x0, g_par, y_par, prior_par = self.problems[depth]
        g = np.zeros((depth, N_STATE))
        acc = np.zeros(N_STATE)
        for i in range(depth):
            g[i] = acc
            acc = A @ acc + B @ u_win[i]
        g_par.value = g
        y_par.value = y_win - u_win @ D.T
        prior_par.value = prior
        _solve(problem, "model")
        return np.vstack([self.powers[i] @ x0.value + g[i] for i in range(depth)])


def run_estimator(estimate, u, y):
    x_hat = np.zeros((u.shape[0], N_STATE))
    latest = {0: PRIOR0.copy()}
    for t in range(u.shape[0]):
        if t + 1 < L:
            t0 = 0
            prior = PRIOR0
        else:
            t0 = t - L + 1
            prior = latest[t0]
        window = estimate(u[t0 : t + 1], y[t0 : t + 1], prior)
        for k, time_index in enumerate(range(t0, t + 1)):
            latest[time_index] = window[k]
        x_hat[t] = window[-1]
    return x_hat


def channel_rmse(x, x_hat):
    err = x - x_hat
    return {
        "x": float(np.sqrt(np.mean(err[:, 0] ** 2))),
        "z": float(np.sqrt(np.mean(err[:, 2] ** 2))),
        "pitch_deg": float(np.degrees(np.sqrt(np.mean(err[:, 4] ** 2)))),
        "x_max": float(np.max(np.abs(err[:, 0]))),
        "z_max": float(np.max(np.abs(err[:, 2]))),
        "pitch_max_deg": float(np.degrees(np.max(np.abs(err[:, 4])))),
    }


def attacked_output(y_clean, attack):
    y = y_clean.copy()
    if attack is None:
        return y
    channel, bias_phys = attack
    y[ATTACK_INDEX:, channel] += bias_phys / SCALE[channel]
    return y


def evaluate(survey=False, horizon=None, n_steps=None, cases=None):
    global L, N_OFFLINE, ATTACK_INDEX
    if horizon is not None:
        L = int(horizon)
    if n_steps is None:
        n_steps = 80 if survey else N_STEPS
    else:
        n_steps = int(n_steps)
    # Willems' length for an input of order L+n, plus a margin of columns.
    need = (N_INPUT + 1) * (L + N_STATE) + 40
    if N_OFFLINE < need:
        N_OFFLINE = need
    ctrl = np.hstack([np.linalg.matrix_power(A, i) @ B for i in range(N_STATE)])
    obsv = np.vstack([C_PHYS @ np.linalg.matrix_power(A, i) for i in range(N_STATE)])
    print(
        f"rank ctrl {np.linalg.matrix_rank(ctrl, tol=1e-8)}/{N_STATE}, "
        f"obsv {np.linalg.matrix_rank(obsv, tol=1e-8)}/{N_STATE}",
        flush=True,
    )
    gain = lqr_gain()
    u_d, x_d, y_d = offline_trajectory(gain, seed=0)
    pe_rank, pe_rows = input_excitation_rank(u_d)
    print(
        f"offline PE rank {pe_rank}/{pe_rows}, "
        f"max |x|={np.max(np.abs(x_d[:, 0])):.3f} m, "
        f"max |z|={np.max(np.abs(x_d[:, 2])):.3f} m, "
        f"max |theta|={np.degrees(np.max(np.abs(x_d[:, 4]))):.2f} deg, "
        f"max |u|={np.max(np.abs(u_d)):.4f} N",
        flush=True,
    )
    if pe_rank != pe_rows:
        raise RuntimeError("Offline input is not PE of order L+n.")
    hu, hy, hx = compress_hankel(*prepare_hankel(u_d, y_d, x_d))
    powers = [np.eye(N_STATE)]
    for _ in range(L - 1):
        powers.append(A @ powers[-1])

    u_online, x_true = closed_loop(n_steps, gain, seed=1, dither_scale=0.0, x0=X0)
    print(
        f"online max |x|={np.max(np.abs(x_true[:, 0])):.3f}, "
        f"max |z|={np.max(np.abs(x_true[:, 2])):.3f}, "
        f"max |theta|={np.degrees(np.max(np.abs(x_true[:, 4]))):.2f} deg, "
        f"max |u|={np.max(np.abs(u_online)):.4f} N",
        flush=True,
    )
    rng = np.random.default_rng(2)
    y_clean = output_scaled(x_true, u_online) + rng.normal(
        0.0, NOISE_STD, size=(n_steps, N_MEAS)
    )
    md = ModelSolver(powers)
    dd = None if survey else DataDrivenSolver(hu, hy, hx)

    def estimators_for(y):
        bank = {
            "Robust MD-MHE": lambda u_win, y_win, prior: robust_md(
                u_win, y_win, prior, powers
            ),
            "Resilient MD-MHE": md.estimate,
        }
        if not survey:
            bank["Robust DD-MHE"] = lambda u_win, y_win, prior: robust_dd(
                hu[u_win.shape[0]],
                hy[u_win.shape[0]],
                hx[u_win.shape[0]],
                u_win,
                y_win,
                prior,
            )
            bank["Resilient DD-MHE"] = dd.estimate
        return bank

    names = ("Robust MD-MHE", "Resilient MD-MHE") if survey else NAMES
    if cases is None:
        cases = ("none", "x", "theta", "gyro", "accel") if survey else ("none", "gyro")
    saved_attack = ATTACK_INDEX
    if survey:
        ATTACK_INDEX = max(n_steps // 3, L)
    results = {}
    hats = {}
    started = time.time()
    try:
        for case in cases:
            y = attacked_output(y_clean, ATTACKS[case])
            bank = estimators_for(y)
            results[case] = {}
            if case == "gyro" and not survey:
                hats[case] = {"x": x_true}
            for name in names:
                t0 = time.time()
                x_hat = run_estimator(bank[name], u_online, y)
                metrics = channel_rmse(x_true, x_hat)
                results[case][name] = metrics
                if case == "gyro" and not survey:
                    hats[case][name] = x_hat
                print(
                    f"{case:6s} {name:18s} "
                    f"x {metrics['x']:.4f}  z {metrics['z']:.4f}  "
                    f"pitch {metrics['pitch_deg']:.3f} deg  "
                    f"({time.time() - t0:.1f}s)",
                    flush=True,
                )
    finally:
        ATTACK_INDEX = saved_attack
    print(f"done in {time.time() - started:.1f}s", flush=True)

    out_dir = Path(__file__).resolve().parent
    report = {
        "mass": MASS,
        "gravity": GRAVITY,
        "arm": ARM,
        "iyy": IYY,
        "ts": TS,
        "horizon": L,
        "scale": SCALE.tolist(),
        "noise_std": NOISE_STD,
        "lambda": R_WEIGHT,
        "rho": RHO,
        "pe_rank": [pe_rank, pe_rows],
        "attack_time": T_ATTACK,
        "gyro_bias_rad_s": ATTACKS["gyro"][1],
        "results": results,
    }
    if not survey:
        figure_path = out_dir.parents[0] / "image" / "quadrotor_mhe.png"
        save_figure(figure_path, hats["gyro"]["x"], hats["gyro"])
        (out_dir / "quadrotor_results.json").write_text(json.dumps(report, indent=2))
        print(f"wrote {figure_path}", flush=True)
    return report


def save_figure(path, x, hats):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    time_axis = np.arange(x.shape[0]) * TS
    series = (
        ("Robust MD-MHE", "#d95f02", "--"),
        ("Robust DD-MHE", "#d95f02", ":"),
        ("Resilient MD-MHE", "#1b9e77", "--"),
        ("Resilient DD-MHE", "#1b9e77", ":"),
    )
    short = {
        "Robust MD-MHE": "Robust MD",
        "Robust DD-MHE": "Robust DD",
        "Resilient MD-MHE": "Resilient MD",
        "Resilient DD-MHE": "Resilient DD",
    }
    columns = ((0, "x (m)"), (2, "z (m)"), (4, "pitch (deg)"))
    fig, axes = plt.subplots(
        1, 3, figsize=(7.1, 1.68), sharex=True, constrained_layout=True
    )
    handles = []
    labels = []
    true_line = None
    for ax, (col, ylabel) in zip(axes, columns):
        scale = 180.0 / np.pi if col == 4 else 1.0
        true_line, = ax.plot(time_axis, x[:, col] * scale, color="#1f77b4", lw=1.8)
        for name, color, style in series:
            line, = ax.plot(
                time_axis,
                hats[name][:, col] * scale,
                color=color,
                lw=1.2,
                ls=style,
            )
            if ax is axes[0]:
                handles.append(line)
                labels.append(short[name])
        ax.axvline(T_ATTACK, color="0.35", ls=":", lw=0.8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(labelsize=8)
        ax.set_xlabel("Time (s)", fontsize=8)
        ax.grid(True, color="0.9", lw=0.6)
    fig.legend(
        (true_line, *handles),
        ("True", *labels),
        frameon=False,
        fontsize=8,
        loc="outside upper center",
        ncol=5,
    )
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _window_estimate(estimate, u_hist, y_hist, latest, prior0):
    t = len(y_hist) - 1
    if t + 1 < L:
        t0 = 0
        prior = prior0
    else:
        t0 = t - L + 1
        prior = latest[t0]
    window = estimate(
        np.asarray(u_hist[t0 : t + 1]),
        np.asarray(y_hist[t0 : t + 1]),
        prior,
    )
    for k, time_index in enumerate(range(t0, t + 1)):
        latest[time_index] = window[k]
    return window[-1]


def fly_with_estimator(n_steps, gain, estimate, seed):
    """Close the linear plant around one estimator on the figure-8."""
    rng = np.random.default_rng(seed)
    x = np.zeros((n_steps, N_STATE))
    x_hat_log = np.zeros((n_steps, N_STATE))
    u_hist, y_hist = [], []
    latest = {}
    x_hat = np.zeros(N_STATE)
    prior0 = np.zeros(N_STATE)
    for t in range(n_steps):
        ref = reference(t)
        u_dev = np.clip(
            -gain @ (x_hat - ref) + figure8_feedforward(t),
            -U_LIMIT,
            U_LIMIT,
        )
        y = output_scaled(x[t], u_dev)
        y = y + rng.normal(0.0, NOISE_STD, size=N_MEAS)
        if t >= ATTACK_INDEX:
            y[3] += ATTACKS["gyro"][1] / SCALE[3]
        u_hist.append(u_dev.copy())
        y_hist.append(y)
        x_hat = _window_estimate(estimate, u_hist, y_hist, latest, prior0)
        x_hat_log[t] = x_hat
        if t + 1 < n_steps:
            x[t + 1] = A @ x[t] + B @ u_dev
        if t % 100 == 0:
            print(
                f"  t={t * TS:.1f}s x={x[t, 0]:+.3f} z={x[t, 2]:+.3f} "
                f"pitch={np.degrees(x[t, 4]):+.1f}",
                flush=True,
            )
    return x, x_hat_log


def tracking_metrics(x):
    ref = np.vstack([reference(t) for t in range(x.shape[0])])
    err = x - ref

    def pack(sl):
        return {
            "x": float(np.sqrt(np.mean(err[sl, 0] ** 2))),
            "z": float(np.sqrt(np.mean(err[sl, 2] ** 2))),
            "pitch_deg": float(np.degrees(np.sqrt(np.mean(err[sl, 4] ** 2)))),
        }

    return {
        "before": pack(slice(0, ATTACK_INDEX)),
        "after": pack(slice(ATTACK_INDEX, None)),
        "reference": ref,
    }


def save_figure8_figure(path, flights):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ref = next(iter(flights.values()))["reference"]
    time_axis = np.arange(ref.shape[0]) * TS
    colors = {"Robust DD-MHE": "#d97706", "Resilient DD-MHE": "#15803d"}
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 5.2), constrained_layout=True)
    path_ax, err_ax, x_ax, z_ax = axes.ravel()
    path_ax.plot(ref[:, 0], ref[:, 2], color="#1d4ed8", lw=1.6, label="reference")
    path_ax.scatter(
        [ref[ATTACK_INDEX, 0]],
        [ref[ATTACK_INDEX, 2]],
        s=28,
        c="0.1",
        zorder=5,
        label="attack starts",
    )
    for name, flight in flights.items():
        path_ax.plot(
            flight["x"][:, 0],
            flight["x"][:, 2],
            color=colors[name],
            lw=1.3,
            label=name,
        )
        err_ax.plot(
            time_axis,
            flight["x"][:, 0] - ref[:, 0],
            color=colors[name],
            lw=1.2,
        )
        x_ax.plot(time_axis, flight["x"][:, 0], color=colors[name], lw=1.2)
        z_ax.plot(time_axis, flight["x"][:, 2], color=colors[name], lw=1.2)
    err_ax.axhline(0.0, color="#1d4ed8", lw=1.0)
    x_ax.plot(time_axis, ref[:, 0], color="#1d4ed8", lw=1.4)
    z_ax.plot(time_axis, ref[:, 2], color="#1d4ed8", lw=1.4)
    for ax in (err_ax, x_ax, z_ax):
        ax.axvline(T_ATTACK, color="0.2", ls=":", lw=0.9)
        ax.set_xlabel("Time (s)")
        ax.grid(True, alpha=0.3)
    span = 1.15 * max(FIG8_AX, FIG8_AZ)
    path_ax.set_xlim(-2.2 * span, 2.2 * span)
    path_ax.set_ylim(-1.6 * span, 1.6 * span)
    path_ax.set_aspect("equal", adjustable="box")
    path_ax.set_xlabel("x (m)")
    path_ax.set_ylabel("z (m)")
    path_ax.grid(True, alpha=0.3)
    path_ax.legend(fontsize=7, loc="best")
    x_ax.set_ylabel("x (m)")
    z_ax.set_ylabel("z (m)")
    err_ax.set_ylabel("x error (m)")
    fig.suptitle(
        f"Figure-8, {ATTACKS['gyro'][1]:.0f} rad/s gyro bias from t={T_ATTACK:.0f} s"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def save_figure8_gif(path, flights, stride=4):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image

    names = ("Robust DD-MHE", "Resilient DD-MHE")
    colors = {"Robust DD-MHE": "#d97706", "Resilient DD-MHE": "#15803d"}
    ref = flights[names[0]]["reference"]
    n_steps = ref.shape[0]
    frames = []
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.6), constrained_layout=True)
    span = 1.15 * max(FIG8_AX, FIG8_AZ)
    for ax, name in zip(axes, names):
        ax.plot(ref[:, 0], ref[:, 2], color="#cbd5e1", lw=1.2)
        ax.set_xlim(-2.2 * span, 2.2 * span)
        ax.set_ylim(-1.6 * span, 1.6 * span)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x (m)")
        ax.set_title(name)
        ax.grid(True, alpha=0.25)
    axes[0].set_ylabel("z (m)")
    trail = {name: ax.plot([], [], color=colors[name], lw=1.4)[0] for name, ax in zip(names, axes)}
    body = {name: ax.plot([], [], color="#111827", lw=2.4)[0] for name, ax in zip(names, axes)}
    hubs = {
        name: ax.scatter([], [], s=18, c=colors[name], zorder=4)
        for name, ax in zip(names, axes)
    }
    title = fig.suptitle("")
    indices = list(range(0, n_steps, stride))
    for t in indices:
        attacked = t >= ATTACK_INDEX
        title.set_text(
            f"t={t * TS:.1f}s"
            + ("   gyro bias on" if attacked else "   no attack")
        )
        for name in names:
            state = flights[name]["x"][t]
            trail[name].set_data(
                flights[name]["x"][: t + 1, 0],
                flights[name]["x"][: t + 1, 2],
            )
            arm = 0.16
            dx = arm * np.cos(state[4])
            dz = arm * np.sin(state[4])
            body[name].set_data(
                [state[0] - dx, state[0] + dx],
                [state[2] - dz, state[2] + dz],
            )
            hubs[name].set_offsets(
                np.array(
                    [
                        [state[0] - dx, state[2] - dz],
                        [state[0] + dx, state[2] + dz],
                    ]
                )
            )
        fig.canvas.draw()
        rgba = np.asarray(fig.canvas.buffer_rgba())
        frames.append(Image.fromarray(rgba[..., :3]))
    plt.close(fig)
    duration = int(round(1000 * stride * TS))
    path.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(
        path,
        save_all=True,
        append_images=frames[1:],
        duration=duration,
        loop=0,
        optimize=True,
    )


def figure8_flights():
    """Fly the figure-8 twice, once with each data-driven estimator in the loop."""
    gain = lqr_gain()
    need = (N_INPUT + 1) * (L + N_STATE) + 40
    global N_OFFLINE
    if N_OFFLINE < need:
        N_OFFLINE = need
    u_d, x_d, y_d = offline_trajectory(gain, seed=0)
    pe_rank, pe_rows = input_excitation_rank(u_d)
    print(f"offline PE rank {pe_rank}/{pe_rows}", flush=True)
    if pe_rank != pe_rows:
        raise RuntimeError("Offline input is not PE of order L+n.")
    hu, hy, hx = compress_hankel(*prepare_hankel(u_d, y_d, x_d))
    dd = DataDrivenSolver(hu, hy, hx)

    def robust(u_win, y_win, prior):
        depth = u_win.shape[0]
        return robust_dd(hu[depth], hy[depth], hx[depth], u_win, y_win, prior)

    bank = {
        "Robust DD-MHE": robust,
        "Resilient DD-MHE": dd.estimate,
    }
    flights = {}
    report = {
        "trajectory": "figure8",
        "period_s": FIG8_PERIOD,
        "ax_m": FIG8_AX,
        "az_m": FIG8_AZ,
        "ts": TS,
        "horizon": L,
        "attack_time": T_ATTACK,
        "gyro_bias_rad_s": ATTACKS["gyro"][1],
        "results": {},
    }
    for name, estimate in bank.items():
        print(f"--- {name} ---", flush=True)
        started = time.time()
        x, x_hat = fly_with_estimator(N_STEPS, gain, estimate, seed=2)
        metrics = tracking_metrics(x)
        flights[name] = {"x": x, "x_hat": x_hat, "reference": metrics.pop("reference")}
        report["results"][name] = metrics
        print(
            f"{name}: before x {metrics['before']['x']:.3f} m, "
            f"after x {metrics['after']['x']:.3f} m, "
            f"after z {metrics['after']['z']:.3f} m, "
            f"after pitch {metrics['after']['pitch_deg']:.2f} deg "
            f"({time.time() - started:.1f}s)",
            flush=True,
        )
    out_dir = Path(__file__).resolve().parent
    figure_path = out_dir / "quadrotor_figure8.png"
    gif_path = out_dir / "quadrotor_figure8.gif"
    (out_dir / "quadrotor_results.json").write_text(json.dumps(report, indent=2))
    save_figure8_figure(figure_path, flights)
    print(f"wrote {figure_path}", flush=True)
    save_figure8_gif(gif_path, flights)
    print(f"wrote {gif_path}", flush=True)
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--survey",
        action="store_true",
        help="Model-driven only, for comparing attack channels.",
    )
    parser.add_argument(
        "--open-loop",
        action="store_true",
        help="Estimate a recorded figure-8 flown with the true state.",
    )
    parser.add_argument("--horizon", type=int, default=None)
    parser.add_argument("--steps", type=int, default=None)
    parser.add_argument("--cases", default=None, help="Comma-separated attack names.")
    args = parser.parse_args()
    if args.survey or args.open_loop:
        cases = tuple(args.cases.split(",")) if args.cases else None
        evaluate(
            survey=args.survey,
            horizon=args.horizon,
            n_steps=args.steps,
            cases=cases,
        )
    else:
        figure8_flights()


if __name__ == "__main__":
    main()
