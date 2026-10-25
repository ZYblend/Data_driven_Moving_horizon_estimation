"""Online figure-8 comparison using the offline Hankel dictionary.

Closes the linear quadrotor around Robust DD-MHE and Resilient DD-MHE.
A pitch-rate bias starts partway through the loop. Run from the repository
root, after ``python -m quadrotor.prepare_offline``:

    python -m quadrotor.simulate_quadrotor
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ddmhe import load_hankel, make_dd_mhe_pair  # noqa: E402
from quadrotor import trajectory as traj  # noqa: E402


def fly_with_estimator(n_steps, gain, estimate, seed):
    """Close the linear plant around one estimator on the figure-8."""
    rng = np.random.default_rng(seed)
    x = np.zeros((n_steps, traj.N_STATE))
    x_hat_log = np.zeros((n_steps, traj.N_STATE))
    u_hist, y_hist = [], []
    latest = {}
    x_hat = np.zeros(traj.N_STATE)
    prior0 = np.zeros(traj.N_STATE)
    for t in range(n_steps):
        ref = traj.reference(t)
        u_dev = np.clip(
            -gain @ (x_hat - ref) + traj.figure8_feedforward(t),
            -traj.U_LIMIT,
            traj.U_LIMIT,
        )
        y = traj.output_scaled(x[t], u_dev)
        y = y + rng.normal(0.0, traj.NOISE_STD, size=traj.N_MEAS)
        if t >= traj.ATTACK_INDEX:
            y[3] += traj.ATTACKS["gyro"][1] / traj.SCALE[3]
        u_hist.append(u_dev.copy())
        y_hist.append(y)
        x_hat = _window_estimate(estimate, u_hist, y_hist, latest, prior0)
        x_hat_log[t] = x_hat
        if t + 1 < n_steps:
            x[t + 1] = traj.A @ x[t] + traj.B @ u_dev
        if t % 100 == 0:
            print(
                f"  t={t * traj.TS:.1f}s x={x[t, 0]:+.3f} z={x[t, 2]:+.3f} "
                f"pitch={np.degrees(x[t, 4]):+.1f}",
                flush=True,
            )
    return x, x_hat_log


def _window_estimate(estimate, u_hist, y_hist, latest, prior0):
    t = len(y_hist) - 1
    if t + 1 < traj.L:
        t0 = 0
        prior = prior0
    else:
        t0 = t - traj.L + 1
        prior = latest[t0]
    window = estimate(np.asarray(u_hist[t0 : t + 1]), np.asarray(y_hist[t0 : t + 1]), prior)
    for k, time_index in enumerate(range(t0, t + 1)):
        latest[time_index] = window[k]
    return window[-1]


def tracking_metrics(x):
    ref = np.vstack([traj.reference(t) for t in range(x.shape[0])])
    err = x - ref

    def pack(sl):
        return {
            "x": float(np.sqrt(np.mean(err[sl, 0] ** 2))),
            "z": float(np.sqrt(np.mean(err[sl, 2] ** 2))),
            "pitch_deg": float(np.degrees(np.sqrt(np.mean(err[sl, 4] ** 2)))),
        }

    return {
        "before": pack(slice(0, traj.ATTACK_INDEX)),
        "after": pack(slice(traj.ATTACK_INDEX, None)),
        "reference": ref,
    }


def estimate_metrics(x, x_hat):
    err = x - x_hat

    def pack(sl):
        return {
            "x": float(np.sqrt(np.mean(err[sl, 0] ** 2))),
            "pitch_deg": float(np.degrees(np.sqrt(np.mean(err[sl, 4] ** 2)))),
            "pitch_rate": float(np.sqrt(np.mean(err[sl, 5] ** 2))),
        }

    return {
        "before": pack(slice(0, traj.ATTACK_INDEX)),
        "after": pack(slice(traj.ATTACK_INDEX, None)),
    }


def _write_figure(fig, path):
    """Save beside the destination, then replace it.

    OneDrive can reject a truncating write of a file it is syncing.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.stem + ".writing" + path.suffix)
    fig.savefig(temporary, dpi=160)
    os.replace(temporary, path)


def save_figure8_figure(path, flights):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ref = next(iter(flights.values()))["reference"]
    time_axis = np.arange(ref.shape[0]) * traj.TS
    colors = {"Robust DD-MHE": "#d97706", "Resilient DD-MHE": "#15803d"}
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 5.2), constrained_layout=True)
    path_ax, err_ax, x_ax, z_ax = axes.ravel()
    path_ax.plot(ref[:, 0], ref[:, 2], color="#1d4ed8", lw=1.6, label="reference")
    path_ax.scatter(
        [ref[traj.ATTACK_INDEX, 0]],
        [ref[traj.ATTACK_INDEX, 2]],
        s=28,
        c="0.1",
        zorder=5,
        label="attack starts",
    )
    for name, flight in flights.items():
        path_ax.plot(flight["x"][:, 0], flight["x"][:, 2], color=colors[name], lw=1.3, label=name)
        err_ax.plot(time_axis, flight["x"][:, 0] - ref[:, 0], color=colors[name], lw=1.2)
        x_ax.plot(time_axis, flight["x"][:, 0], color=colors[name], lw=1.2)
        z_ax.plot(time_axis, flight["x"][:, 2], color=colors[name], lw=1.2)
    err_ax.axhline(0.0, color="#1d4ed8", lw=1.0)
    x_ax.plot(time_axis, ref[:, 0], color="#1d4ed8", lw=1.4)
    z_ax.plot(time_axis, ref[:, 2], color="#1d4ed8", lw=1.4)
    for ax in (err_ax, x_ax, z_ax):
        ax.axvline(traj.T_ATTACK, color="0.2", ls=":", lw=0.9)
        ax.set_xlabel("Time (s)")
        ax.grid(True, alpha=0.3)
    span = 1.15 * max(traj.FIG8_AX, traj.FIG8_AZ)
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
    fig.suptitle(f"Figure-8, {traj.ATTACKS['gyro'][1]:.0f} rad/s gyro bias from t={traj.T_ATTACK:.0f} s")
    _write_figure(fig, path)
    plt.close(fig)


def save_estimate_figure(path, flights):
    """True state against each estimator. The pitch-rate panel is the biased channel."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    time_axis = np.arange(next(iter(flights.values()))["x"].shape[0]) * traj.TS
    colors = {"Robust DD-MHE": "#d97706", "Resilient DD-MHE": "#15803d"}
    columns = (
        (0, 1.0, "x (m)"),
        (4, 180.0 / np.pi, "pitch (deg)"),
        (5, 1.0, "pitch rate (rad/s)"),
    )
    fig, axes = plt.subplots(
        2, 3, figsize=(8.6, 4.8), sharex=True, sharey="col", constrained_layout=True
    )
    for row, name in enumerate(("Robust DD-MHE", "Resilient DD-MHE")):
        flight = flights[name]
        for col, (index, scale, ylabel) in enumerate(columns):
            ax = axes[row, col]
            ax.plot(time_axis, flight["x"][:, index] * scale, color="#1d4ed8", lw=1.5, label="true")
            ax.plot(
                time_axis,
                flight["x_hat"][:, index] * scale,
                color=colors[name],
                lw=1.2,
                ls="--",
                label="estimate",
            )
            ax.axvline(traj.T_ATTACK, color="0.2", ls=":", lw=0.9)
            ax.grid(True, alpha=0.3)
            if row == 1:
                ax.set_xlabel("Time (s)")
            if col == 0:
                ax.set_ylabel(f"{name}\n{ylabel}")
            else:
                ax.set_ylabel(ylabel)
            if row == 0 and col == 0:
                ax.legend(fontsize=7, loc="best")
    fig.suptitle("Estimate versus true state. The gyro bias enters the pitch-rate measurement.")
    _write_figure(fig, path)
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
    span = 1.15 * max(traj.FIG8_AX, traj.FIG8_AZ)
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
        name: ax.scatter([], [], s=18, c=colors[name], zorder=4) for name, ax in zip(names, axes)
    }
    title = fig.suptitle("")
    for t in range(0, n_steps, stride):
        title.set_text(
            f"t={t * traj.TS:.1f}s" + ("   gyro bias on" if t >= traj.ATTACK_INDEX else "   no attack")
        )
        for name in names:
            state = flights[name]["x"][t]
            trail[name].set_data(flights[name]["x"][: t + 1, 0], flights[name]["x"][: t + 1, 2])
            arm = 0.16
            dx = arm * np.cos(state[4])
            dz = arm * np.sin(state[4])
            body[name].set_data([state[0] - dx, state[0] + dx], [state[2] - dz, state[2] + dz])
            hubs[name].set_offsets(
                np.array([[state[0] - dx, state[2] - dz], [state[0] + dx, state[2] + dz]])
            )
        fig.canvas.draw()
        rgba = np.asarray(fig.canvas.buffer_rgba())
        frames.append(Image.fromarray(rgba[..., :3]))
    plt.close(fig)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.stem + ".writing" + path.suffix)
    frames[0].save(
        temporary,
        save_all=True,
        append_images=frames[1:],
        duration=int(round(1000 * stride * traj.TS)),
        loop=0,
        optimize=True,
    )
    os.replace(temporary, path)


def figure8_flights():
    _u, _y, _x, hu, hy, hx, meta = load_hankel(traj.OFFLINE_PATH)
    horizon = int(meta["horizon"])
    if horizon != traj.L:
        raise RuntimeError(
            f"Offline horizon {horizon} does not match L={traj.L}. Rebuild with prepare_offline."
        )
    solvers = make_dd_mhe_pair(
        hu,
        hy,
        hx,
        horizon=horizon,
        p_weight=traj.P_WEIGHT,
        r_weight=traj.R_WEIGHT,
        rho=traj.RHO,
    )
    gain = traj.lqr_gain()
    flights = {}
    report = {
        "trajectory": "figure8",
        "period_s": traj.FIG8_PERIOD,
        "ax_m": traj.FIG8_AX,
        "az_m": traj.FIG8_AZ,
        "ts": traj.TS,
        "horizon": horizon,
        "attack_time": traj.T_ATTACK,
        "gyro_bias_rad_s": traj.ATTACKS["gyro"][1],
        "results": {},
        "estimates": {},
    }
    for name, solver in solvers.items():
        print(f"--- {name} ---", flush=True)
        started = time.time()
        x, x_hat = fly_with_estimator(traj.N_STEPS, gain, solver.estimate, seed=2)
        metrics = tracking_metrics(x)
        estimates = estimate_metrics(x, x_hat)
        flights[name] = {"x": x, "x_hat": x_hat, "reference": metrics.pop("reference")}
        report["results"][name] = metrics
        report["estimates"][name] = estimates
        print(
            f"{name}: vehicle after x {metrics['after']['x']:.3f} m, "
            f"pitch-rate estimate error after {estimates['after']['pitch_rate']:.3f} rad/s "
            f"({time.time() - started:.1f}s)",
            flush=True,
        )
    out_dir = Path(__file__).resolve().parent
    figure_path = out_dir / "quadrotor_figure8.png"
    estimate_path = out_dir / "quadrotor_estimate.png"
    gif_path = out_dir / "quadrotor_figure8.gif"
    (out_dir / "quadrotor_results.json").write_text(json.dumps(report, indent=2))
    save_figure8_figure(figure_path, flights)
    save_estimate_figure(estimate_path, flights)
    print(f"wrote {figure_path}", flush=True)
    print(f"wrote {estimate_path}", flush=True)
    save_figure8_gif(gif_path, flights)
    print(f"wrote {gif_path}", flush=True)
    return report


def main():
    figure8_flights()


if __name__ == "__main__":
    main()
