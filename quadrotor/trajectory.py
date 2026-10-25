"""Planar quadrotor linearization, figure-8 reference, and offline rollout."""

from __future__ import annotations

from pathlib import Path

import numpy as np
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

A, B, C_PHYS, D_PHYS = discrete_plant()
N_STATE, N_INPUT = A.shape[0], B.shape[1]
N_MEAS = C_PHYS.shape[0]

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
FIG8_PERIOD = 8.0
FIG8_AX = 0.45
FIG8_AZ = 0.22
T_FINAL = 12.0
T_ATTACK = 4.0
N_STEPS = int(round(T_FINAL / TS))
ATTACK_INDEX = int(round(T_ATTACK / TS))
ATTACKS = {
    "none": None,
    "x": (0, 1.0),
    "theta": (2, 0.35),
    "gyro": (3, 8.0),
    "accel": (4, 5.0),
}
X0 = np.zeros(N_STATE)
U_LIMIT = 0.5 * MASS * GRAVITY
OFFLINE_PATH = Path(__file__).resolve().parent / "offline" / "hankel.npz"


def lqr_gain():
    q = np.diag([6.0, 1.0, 8.0, 1.0, 0.4, 0.02])
    r = np.diag([40.0, 40.0])
    p = solve_discrete_are(A, B, q, r)
    return np.linalg.solve(r + B.T @ p @ B, B.T @ p @ A)


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


def output_scaled(x, u):
    return x @ C.T + u @ D.T


def offline_trajectory(gain, seed, n_samples=None):
    """Hover about the origin with a dithered input. The state is recorded."""
    n_samples = N_OFFLINE if n_samples is None else n_samples
    rng = np.random.default_rng(seed)
    u = np.zeros((n_samples, N_INPUT))
    x = np.zeros((n_samples, N_STATE))
    for t in range(n_samples - 1):
        dither = rng.uniform(-DITHER, DITHER, size=N_INPUT)
        u[t] = -gain @ x[t] + dither
        x[t + 1] = A @ x[t] + B @ u[t]
    u[-1] = -gain @ x[-1]
    y = output_scaled(x, u)
    return u, x, y
