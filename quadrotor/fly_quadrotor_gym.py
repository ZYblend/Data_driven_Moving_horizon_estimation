"""Fly the safe-control-gym planar Crazyflie with the paper's estimators in the loop.

The vehicle is the TWO_D quadrotor, integrated with that suite's explicit
dynamics and rendered by its PyBullet camera. Each flight is closed around
one estimator: the thrust is the same linear-quadratic regulator used in the
paper, evaluated at the estimated deviation from hover. A bias is added to
the pitch-rate measurement that the estimator sees. The two flights are the
robust and resilient data-driven estimators.

safe-control-gym is the submodule in this folder. Its explicit dynamics
read an attribute named TIMESTEP that the class never sets; this script sets
it to the PyBullet period.

Run from the repository root, after `git submodule update --init quadrotor/safe-control-gym`:
    python3 -m quadrotor.fly_quadrotor_gym --paper
"""

import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

# Their control rate is 60 Hz, and the physics rate must divide it.
import quadrotor.quadrotor_plant as plant

plant.TS = 1.0 / 60.0

import numpy as np
import pybullet as p
import quadrotor.simulate_quadrotor as sim

_HERE = os.path.dirname(os.path.abspath(__file__))
GYM_ROOT = os.environ.get("SAFE_CONTROL_GYM", os.path.join(_HERE, "safe-control-gym"))
if GYM_ROOT not in sys.path:
    sys.path.insert(0, GYM_ROOT)

# normalization.py imports torch at module scope. The planar flight never
# calls the torch path, so a stand-in is enough when torch is absent.
try:
    import torch  # noqa: F401
except ImportError:
    import types

    _torch = types.ModuleType("torch")

    class _Tensor:
        pass

    _torch.Tensor = _Tensor
    sys.modules["torch"] = _torch

from safe_control_gym.envs.gym_pybullet_drones.quadrotor import Quadrotor

# Window length in seconds stays near the paper's 0.4 s.
sim.L = 24
HOVER_Z = 1.0
DURATION = 6.0
ATTACK_TIME = 3.0
GYRO_BIAS = 8.0  # rad/s, same bias as the paper comparison
PRIOR = np.array([0.05, 0.0, 0.0, 0.0, 0.02, 0.0])
FRAME_STRIDE = 3


def make_env():
    env = Quadrotor(
        quad_type=2,
        physics="dyn",
        gui=False,
        ctrl_freq=60,
        pyb_freq=240,
        episode_len_sec=int(DURATION) + 2,
        randomized_init=False,
        done_on_out_of_bound=False,
        init_state=np.array([0.0, 0.0, HOVER_Z, 0.0, 0.0, 0.0]),
    )
    env.reset()
    # base_aviary._dynamics integrates with self.TIMESTEP, which is unset.
    env.TIMESTEP = env.PYB_TIMESTEP
    env.RENDER_WIDTH = 480
    env.RENDER_HEIGHT = 360
    env.CAM_PRO = p.computeProjectionMatrixFOV(
        fov=42.0,
        aspect=env.RENDER_WIDTH / env.RENDER_HEIGHT,
        nearVal=0.02,
        farVal=20.0,
        physicsClientId=env.PYB_CLIENT,
    )
    # The waypoint the regulator is flying toward, drawn in both panels.
    marker = p.createVisualShape(
        p.GEOM_SPHERE,
        radius=0.018,
        rgbaColor=[0.85, 0.15, 0.12, 1.0],
        physicsClientId=env.PYB_CLIENT,
    )
    p.createMultiBody(
        baseMass=0,
        baseVisualShapeIndex=marker,
        basePosition=[sim.X_GOAL, 0.0, HOVER_Z + sim.Z_GOAL],
        physicsClientId=env.PYB_CLIENT,
    )
    aim_camera(env)
    return env


def aim_camera(env):
    """Close side view of the short move, so a drift of a few decimetres reads clearly."""
    env.CAM_VIEW = p.computeViewMatrix(
        cameraEyePosition=[0.20, -0.52, 1.10],
        cameraTargetPosition=[0.20, 0.0, 1.10],
        cameraUpVector=[0.0, 0.0, 1.0],
        physicsClientId=env.PYB_CLIENT,
    )


def deviation(state):
    dev = np.array(state, dtype=float).copy()
    dev[2] -= HOVER_Z
    return dev


def measurement(dev, u_dev, rng, attacked):
    y = sim.output_scaled(dev[None, :], u_dev[None, :])[0]
    y = y + rng.normal(0.0, sim.NOISE_STD, size=y.shape)
    if attacked:
        y[3] += GYRO_BIAS / sim.SCALE[3]
    return y


def estimate_now(estimate, u_hist, y_hist, latest):
    t = len(y_hist) - 1
    if t + 1 < sim.L:
        t0 = 0
        prior = PRIOR
    else:
        t0 = t - sim.L + 1
        prior = latest[t0]
    window = estimate(
        np.asarray(u_hist[t0 : t + 1]),
        np.asarray(y_hist[t0 : t + 1]),
        prior,
    )
    for k, time_index in enumerate(range(t0, t + 1)):
        latest[time_index] = window[k]
    return window[-1]


def fly(name, estimate, gain, hover, attack=True):
    env = make_env()
    n_steps = int(round(DURATION / sim.TS))
    attack_index = int(round(ATTACK_TIME / sim.TS))
    rng = np.random.default_rng(2)
    u_hist, y_hist = [], []
    states = []
    frames = []
    latest = {}
    # The command at this step is computed from the previous state estimate,
    # then the current state is measured and the vehicle is advanced. Using the
    # just-computed terminal estimate on the state it describes excites a pitch
    # oscillation in the resilient loop.
    x_hat = PRIOR.copy()
    low = np.array(env.action_space.low, dtype=float)
    high = np.array(env.action_space.high, dtype=float)
    for t in range(n_steps):
        dev = deviation(env.state)
        ref = sim.reference(t)
        u_dev = -gain @ (x_hat - ref)
        thrust = np.clip(hover + u_dev, low, high)
        u_dev = thrust - hover
        y = measurement(dev, u_dev, rng, attack and t >= attack_index)
        u_hist.append(u_dev.copy())
        y_hist.append(y)
        x_hat = estimate_now(estimate, u_hist, y_hist, latest)
        states.append(dev.copy())
        if t % FRAME_STRIDE == 0:
            frames.append(np.ascontiguousarray(env.render()[..., :3], dtype=np.uint8))
        env.step(thrust)
        if t % 60 == 0:
            print(
                f"  {name} t={t * sim.TS:.1f}s  x={dev[0]:+.3f}  z={dev[2] + HOVER_Z:.3f}  "
                f"pitch={np.degrees(dev[4]):+.1f}  rate={dev[5]:+.2f}  "
                f"hat_pitch={np.degrees(x_hat[4]):+.1f}  hat_rate={x_hat[5]:+.2f}",
                flush=True,
            )
    env.close()
    states = np.asarray(states)
    print(
        f"{name}: final x {states[-1, 0]:.3f} m, z {states[-1, 2] + HOVER_Z:.3f} m, "
        f"pitch {np.degrees(states[-1, 4]):.1f} deg, "
        f"max |pitch| {np.degrees(np.max(np.abs(states[:, 4]))):.1f} deg",
        flush=True,
    )
    return states, frames


def annotate(frame, title, dev, t):
    from PIL import Image, ImageDraw

    image = Image.fromarray(frame)
    draw = ImageDraw.Draw(image)
    draw.rectangle((0, 0, image.width, 36), fill=(255, 255, 255))
    text = (
        f"{title}   t={t:0.1f}s   "
        f"x={dev[0]:+0.2f} m   z={dev[2] + HOVER_Z:0.2f} m   "
        f"pitch={np.degrees(dev[4]):+0.1f} deg"
    )
    draw.text((8, 10), text, fill=(20, 20, 20))
    return np.asarray(image)


def fly_truth(gain, hover):
    """Same regulator on the measured state, so the controller itself can be judged."""
    env = make_env()
    n_steps = int(round(DURATION / sim.TS))
    low = np.array(env.action_space.low, dtype=float)
    high = np.array(env.action_space.high, dtype=float)
    states = []
    for t in range(n_steps):
        dev = deviation(env.state)
        thrust = np.clip(hover - gain @ (dev - sim.reference(t)), low, high)
        env.step(thrust)
        states.append(deviation(env.state))
    env.close()
    states = np.asarray(states)
    print(
        f"True state: final x {states[-1, 0]:.3f} m, z {states[-1, 2] + HOVER_Z:.3f} m, "
        f"pitch {np.degrees(states[-1, 4]):.1f} deg",
        flush=True,
    )
    return states


def write_video(path, robust_frames, resilient_frames, robust_states, resilient_states):
    import subprocess

    frames = []
    n = min(len(robust_frames), len(resilient_frames))
    for i in range(n):
        t = i * FRAME_STRIDE * sim.TS
        left = annotate(robust_frames[i], "Robust DD-MHE", robust_states[i * FRAME_STRIDE], t)
        right = annotate(
            resilient_frames[i],
            "Resilient DD-MHE",
            resilient_states[i * FRAME_STRIDE],
            t,
        )
        gap = np.full((left.shape[0], 8, 3), 255, dtype=np.uint8)
        frames.append(np.concatenate([left, gap, right], axis=1))
    height, width = frames[0].shape[:2]
    proc = subprocess.Popen(
        [
            "ffmpeg", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24",
            "-s", f"{width}x{height}", "-r", "20", "-i", "-",
            "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
            "-movflags", "+faststart",
            path,
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    for frame in frames:
        proc.stdin.write(frame.tobytes())
    proc.stdin.close()
    if proc.wait() != 0:
        raise RuntimeError("ffmpeg failed while writing the flight video")
    print(f"wrote {path}", flush=True)


def tracking_report(name, states):
    n = states.shape[0]
    attack_index = int(round(ATTACK_TIME / sim.TS))
    ref = np.vstack([sim.reference(t) for t in range(n)])
    err = states - ref
    post = slice(attack_index, None)
    stats = {
        "x_rmse": float(np.sqrt(np.mean(err[:, 0] ** 2))),
        "z_rmse": float(np.sqrt(np.mean(err[:, 2] ** 2))),
        "x_rmse_post": float(np.sqrt(np.mean(err[post, 0] ** 2))),
        "z_rmse_post": float(np.sqrt(np.mean(err[post, 2] ** 2))),
        "final_x": float(states[-1, 0]),
        "final_z": float(states[-1, 2]),
        "max_pitch_deg": float(np.degrees(np.max(np.abs(states[:, 4])))),
    }
    print(
        f"{name:18s} x {stats['x_rmse']:.4f} ({stats['x_rmse_post']:.4f} after)  "
        f"z {stats['z_rmse']:.4f} ({stats['z_rmse_post']:.4f} after)  "
        f"final x {stats['final_x']:.3f} z {stats['final_z']:.3f}  "
        f"max |pitch| {stats['max_pitch_deg']:.2f} deg",
        flush=True,
    )
    return stats


def save_flight_figure(path, flights):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_steps = next(iter(flights.values())).shape[0]
    time_axis = np.arange(n_steps) * sim.TS
    ref = np.vstack([sim.reference(t) for t in range(n_steps)])
    attack_index = int(round(ATTACK_TIME / sim.TS))
    # Model-driven paths coincide with these data-driven paths, so the figure
    # draws the data-driven flights only.
    series = (
        ("Robust DD-MHE", "#d95f02", "-"),
        ("Resilient DD-MHE", "#1b9e77", "-"),
    )
    short = {
        "Robust DD-MHE": "Robust DD-MHE",
        "Resilient DD-MHE": "Resilient DD-MHE",
    }
    fig = plt.figure(figsize=(7.1, 2.15), constrained_layout=True)
    path_ax, pitch_ax = fig.subplots(1, 2, width_ratios=(1.15, 1.0))
    handles = []
    labels = []
    ref_line, = path_ax.plot(ref[:, 0], ref[:, 2], color="#1f77b4", lw=1.8)
    for name, color, style in series:
        line, = path_ax.plot(
            flights[name][:, 0],
            flights[name][:, 2],
            color=color,
            lw=1.2,
            ls=style,
        )
        handles.append(line)
        labels.append(short[name])
    path_ax.scatter(
        [ref[attack_index, 0]],
        [ref[attack_index, 2]],
        s=16,
        c="0.15",
        zorder=5,
    )
    path_ax.set_aspect("equal", adjustable="box")
    path_ax.set_xlabel("x (m)", fontsize=8)
    path_ax.set_ylabel("z (m)", fontsize=8)
    path_ax.tick_params(labelsize=8)
    path_ax.grid(True, color="0.9", lw=0.6)

    pitch_ax.plot(time_axis, np.degrees(ref[:, 4]), color="#1f77b4", lw=1.8)
    for name, color, style in series:
        pitch_ax.plot(
            time_axis,
            np.degrees(flights[name][:, 4]),
            color=color,
            lw=1.2,
            ls=style,
        )
    pitch_ax.axvline(ATTACK_TIME, color="0.35", ls=":", lw=0.8)
    pitch_ax.set_xlabel("Time (s)", fontsize=8)
    pitch_ax.set_ylabel("pitch (deg)", fontsize=8)
    pitch_ax.tick_params(labelsize=8)
    pitch_ax.grid(True, color="0.9", lw=0.6)
    fig.legend(
        (ref_line, *handles),
        ("Reference", *labels),
        frameon=False,
        fontsize=8,
        loc="outside upper center",
        ncol=3,
    )
    fig.savefig(path, dpi=200)
    plt.close(fig)


def paper_flights():
    """Closed-loop flights for the quadrotor subsection, with and without the bias."""
    global FRAME_STRIDE
    FRAME_STRIDE = 10**9
    gain = sim.lqr_gain()
    sim.N_OFFLINE = max(sim.N_OFFLINE, (sim.N_INPUT + 1) * (sim.L + sim.N_STATE) + 40)
    u_d, x_d, y_d = sim.offline_trajectory(gain, seed=0)
    pe_rank, pe_rows = sim.input_excitation_rank(u_d)
    print(f"sample {sim.TS:.6f}s, L={sim.L}, offline rank {pe_rank}/{pe_rows}", flush=True)
    hu, hy, hx = sim.compress_hankel(*sim.prepare_hankel(u_d, y_d, x_d))
    powers = [np.eye(sim.N_STATE)]
    for _ in range(sim.L - 1):
        powers.append(sim.A @ powers[-1])
    dd = sim.DataDrivenSolver(hu, hy, hx)
    md = sim.ModelSolver(powers)

    def robust_dd(u_win, y_win, prior):
        depth = u_win.shape[0]
        return sim.robust_dd(hu[depth], hy[depth], hx[depth], u_win, y_win, prior)

    def robust_md(u_win, y_win, prior):
        return sim.robust_md(u_win, y_win, prior, powers)

    bank = {
        "Robust MD-MHE": robust_md,
        "Robust DD-MHE": robust_dd,
        "Resilient MD-MHE": md.estimate,
        "Resilient DD-MHE": dd.estimate,
    }
    env = make_env()
    hover = np.full(2, env.MASS * env.GRAVITY_ACC / 2.0)
    env.close()
    report = {"pe_rank": [pe_rank, pe_rows], "cases": {}}
    attacked = {}
    for attack, case in (False, "none"), (True, "gyro"):
        print(f"--- {case} ---", flush=True)
        report["cases"][case] = {}
        for name, estimate in bank.items():
            states, _ = fly(name, estimate, gain, hover, attack=attack)
            report["cases"][case][name] = tracking_report(name, states)
            if attack:
                attacked[name] = states
    figure_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "image",
        "quadrotor_mhe.png",
    )
    save_flight_figure(figure_path, attacked)
    print(f"wrote {figure_path}", flush=True)
    return report


def main():
    gain = sim.lqr_gain()
    # Offline excitation of the linearization at this sample time.
    sim.N_OFFLINE = max(sim.N_OFFLINE, (sim.N_INPUT + 1) * (sim.L + sim.N_STATE) + 40)
    u_d, x_d, y_d = sim.offline_trajectory(gain, seed=0)
    pe_rank, pe_rows = sim.input_excitation_rank(u_d)
    print(f"sample {sim.TS:.4f}s, L={sim.L}, offline rank {pe_rank}/{pe_rows}", flush=True)
    if pe_rank != pe_rows:
        raise RuntimeError("Offline input is not persistently exciting.")
    hu, hy, hx = sim.compress_hankel(*sim.prepare_hankel(u_d, y_d, x_d))
    powers = [np.eye(sim.N_STATE)]
    for _ in range(sim.L - 1):
        powers.append(sim.A @ powers[-1])
    dd = sim.DataDrivenSolver(hu, hy, hx)

    def robust(u_win, y_win, prior):
        depth = u_win.shape[0]
        return sim.robust_dd(hu[depth], hy[depth], hx[depth], u_win, y_win, prior)

    # One window, gyro bias only, to confirm the resilient cost still rejects it.
    depth = sim.L
    u0 = np.zeros((depth, sim.N_INPUT))
    y0 = np.zeros((depth, sim.N_MEAS))
    y0[:, 3] = GYRO_BIAS / sim.SCALE[3]
    md = sim.ModelSolver(powers)
    hat_s = md.estimate(u0, y0, np.zeros(sim.N_STATE))
    hat_r = sim.robust_md(u0, y0, np.zeros(sim.N_STATE), powers)
    print(
        f"one-window pitch-rate: resilient {hat_s[:, 5].mean():.3f}, "
        f"robust {hat_r[:, 5].mean():.3f} rad/s",
        flush=True,
    )

    env = make_env()
    hover = np.full(2, env.MASS * env.GRAVITY_ACC / 2.0)
    env.close()
    print("hover thrust", hover, "gain theta column", gain[:, 4], flush=True)
    fly_truth(gain, hover)

    robust_states, robust_frames = fly("Robust DD-MHE", robust, gain, hover)
    resilient_states, resilient_frames = fly(
        "Resilient DD-MHE", dd.estimate, gain, hover
    )
    out = os.environ.get(
        "FLIGHT_VIDEO",
        "/opt/cursor/artifacts/quadrotor_feedback_corrected.mp4",
    )
    os.makedirs(os.path.dirname(out), exist_ok=True)
    write_video(out, robust_frames, resilient_frames, robust_states, resilient_states)
    # Position error against the reference after the move, before and after the attack.
    n = robust_states.shape[0]
    attack_index = int(round(ATTACK_TIME / sim.TS))
    ref = np.vstack([sim.reference(t) for t in range(n)])
    for label, states in ("robust", robust_states), ("resilient", resilient_states):
        err = states[:, 0] - ref[:, 0]
        post = slice(attack_index, None)
        print(
            f"{label} x RMSE after the bias {np.sqrt(np.mean(err[post] ** 2)):.3f} m, "
            f"final pitch {np.degrees(states[-1, 4]):.1f} deg",
            flush=True,
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--paper",
        action="store_true",
        help="Run the four closed-loop flights and write the paper figure.",
    )
    args = parser.parse_args()
    if args.paper:
        paper_flights()
    else:
        main()
