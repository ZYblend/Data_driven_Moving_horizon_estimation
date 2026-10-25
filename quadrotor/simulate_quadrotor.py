"""Online figure-8 comparison using the offline Hankel dictionary.

Closes the linear quadrotor around Robust DD-MHE and Resilient DD-MHE.
A pitch-rate bias starts partway through the loop. Run from the repository
root, after ``python -m quadrotor.prepare_offline``:

    python -m quadrotor.simulate_quadrotor
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import threading
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ddmhe import load_hankel, make_dd_mhe_pair  # noqa: E402
from quadrotor import trajectory as traj  # noqa: E402


class _Flight:
    """One closed-loop figure-8. ``step`` advances a single sample."""

    def __init__(self, n_steps, gain, estimate, seed):
        self.n_steps = n_steps
        self.gain = gain
        self.estimate = estimate
        self.rng = np.random.default_rng(seed)
        self.x = np.zeros((n_steps, traj.N_STATE))
        self.x_hat_log = np.zeros((n_steps, traj.N_STATE))
        self.u_hist, self.y_hist = [], []
        self.latest = {}
        self.x_hat = np.zeros(traj.N_STATE)
        self.prior0 = np.zeros(traj.N_STATE)
        self.t = 0

    def step(self):
        t = self.t
        ref = traj.reference(t)
        u_dev = np.clip(
            -self.gain @ (self.x_hat - ref) + traj.figure8_feedforward(t),
            -traj.U_LIMIT,
            traj.U_LIMIT,
        )
        y = traj.output_scaled(self.x[t], u_dev)
        y = y + self.rng.normal(0.0, traj.NOISE_STD, size=traj.N_MEAS)
        if t >= traj.ATTACK_INDEX:
            y[3] += traj.ATTACKS["gyro"][1] / traj.SCALE[3]
        self.u_hist.append(u_dev.copy())
        self.y_hist.append(y)
        self.x_hat = _window_estimate(self.estimate, self.u_hist, self.y_hist, self.latest, self.prior0)
        self.x_hat_log[t] = self.x_hat
        if t + 1 < self.n_steps:
            self.x[t + 1] = traj.A @ self.x[t] + traj.B @ u_dev
        self.t = t + 1


def fly_with_estimator(n_steps, gain, estimate, seed):
    """Close the linear plant around one estimator on the figure-8."""
    flight = _Flight(n_steps, gain, estimate, seed)
    for t in range(n_steps):
        flight.step()
        if t % 100 == 0:
            print(
                f"  t={t * traj.TS:.1f}s x={flight.x[t, 0]:+.3f} z={flight.x[t, 2]:+.3f} "
                f"pitch={np.degrees(flight.x[t, 4]):+.1f}",
                flush=True,
            )
    return flight.x, flight.x_hat_log


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


def figure8_flights(show_gui=True):
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
    runners = [(name, _Flight(traj.N_STEPS, gain, solver.estimate, seed=2)) for name, solver in solvers.items()]
    view = None
    if show_gui:
        from quadrotor.gym_view import RESILIENT, ROBUST, LiveView

        reference = np.vstack([traj.reference(t) for t in range(traj.N_STEPS)])
        view = LiveView(traj.TS, reference)

    def run_estimators():
        started = time.time()
        try:
            for t in range(traj.N_STEPS):
                states = {}
                for name, flight in runners:
                    flight.step()
                    states[name] = flight.x[t]
                if view is not None:
                    view.push(states[ROBUST], states[RESILIENT], traj.reference(t))
                if t % 100 == 0:
                    print(
                        f"  t={t * traj.TS:.1f}s  "
                        + "  ".join(
                            f"{name.split()[0]} x={states[name][0]:+.3f} z={states[name][2]:+.3f}"
                            for name, _flight in runners
                        ),
                        flush=True,
                    )
            print(f"estimators finished in {time.time() - started:.1f}s", flush=True)
        except Exception as exc:
            errors.append(exc)

    errors = []

    if view is None:
        run_estimators()
        if errors:
            raise errors[0]
    else:
        worker = threading.Thread(target=run_estimators, name="dd-mhe")
        worker.start()
        window_open = True
        try:
            while worker.is_alive():
                if not view.pump():
                    window_open = False
                    break
                time.sleep(0.016)
            worker.join()
        finally:
            if not window_open:
                view.close()
        if errors:
            if window_open:
                view.close()
            raise errors[0]

    for name, flight in runners:
        metrics = tracking_metrics(flight.x)
        estimates = estimate_metrics(flight.x, flight.x_hat_log)
        flights[name] = {"x": flight.x, "x_hat": flight.x_hat_log, "reference": metrics.pop("reference")}
        report["results"][name] = metrics
        report["estimates"][name] = estimates
        print(
            f"{name}: vehicle after x {metrics['after']['x']:.3f} m, "
            f"pitch-rate estimate error after {estimates['after']['pitch_rate']:.3f} rad/s",
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
    if show_gui and view is not None and view.pump():
        print("Plots are saved. Close the PyBullet window to exit.", flush=True)
        try:
            while view.pump():
                time.sleep(0.016)
        finally:
            view.close()
    return report


def main():
    parser = argparse.ArgumentParser(description="Figure-8 comparison of the two DD-MHE estimators.")
    parser.add_argument(
        "--no-gui",
        action="store_true",
        help="Skip the live safe-control-gym window.",
    )
    args = parser.parse_args()
    figure8_flights(show_gui=not args.no_gui)


if __name__ == "__main__":
    main()
