"""Build the IEEE 14-bus offline Hankel dictionary.

Run from the repository root, before the online attack sweep:

    python power/prepare_ieee14.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ddmhe import (  # noqa: E402
    build_compressed_hankel,
    print_hankel_analysis,
    rank_analysis,
    save_hankel,
)
from power.ieee14_model import load_plant  # noqa: E402

TS = 0.01
L = 6
P_WEIGHT = 1.0
R_WEIGHT = 10.0
RHO = 0.95
U_SCALE = 0.02
N_OFFLINE = 800
HOLD = 2
OFFLINE_PATH = Path(__file__).resolve().parent / "offline" / "ieee14_hankel.npz"


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
    scale = np.std(y, axis=0)
    scale = np.maximum(scale, 0.05 * np.median(scale[scale > 0]))
    scale = np.maximum(scale, 1e-3)
    return scale


def main():
    plant = load_plant()
    A, B, C, D = plant["A"], plant["B"], plant["C"], plant["D"]
    n, p = plant["n"], plant["p"]
    n_offline = max(N_OFFLINE, (p + 1) * (L + n) + 40)
    print(
        f"IEEE 14 offline rollout: n={n}, p={p}, m={plant['m']}, L={L}, N={n_offline}",
        flush=True,
    )
    u, x, y_phys = offline_trajectory(A, B, C, D, n_offline, seed=0)
    scale = measurement_scale(y_phys)
    y = y_phys / scale
    report = rank_analysis(u, y, x, L)
    hu, hy, hx = build_compressed_hankel(u, y, x, L)
    print_hankel_analysis(report, compressed_columns=hu[L].shape[1])
    if not report["persistently_exciting"]:
        raise SystemExit("Offline input is not persistently exciting of order L+n.")
    meta = {
        "horizon": L,
        "ts": TS,
        "n": n,
        "p": p,
        "m": plant["m"],
        "p_weight": P_WEIGHT,
        "r_weight": R_WEIGHT,
        "rho": RHO,
        "scale": scale.tolist(),
        "rank": report,
    }
    save_hankel(OFFLINE_PATH, u, y, x, hu, hy, hx, meta)


if __name__ == "__main__":
    main()
