"""Discrete linearized IEEE 14-bus plant from Zheng/Anubi FDIA demos.

Matrices in ``power/ieee14/*.csv`` are the discrete plant from
https://github.com/ZYblend/Moving-horizon-FDIA (branch ``demo_14_bus_system``):

    A_bar_d, B_bar_d, C_obsv_d, D_obsv_d

with sample time ``T = 0.01`` s. The continuous-time construction and the
discretizer ``B_d = 0.5 T (expm(A T) + I) B`` are documented in
``SIM_Cyber-Physical-Systems/IEEE_14_bus_System/run_model.m``.

States (n = 10): generator angle and speed deviations.
Inputs (p = 19): generator mechanical power and bus load deviations.
Outputs (m = 19): generator speeds and net bus power injections.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

TS = 0.01
N_GEN = 5
N_BUS = 14
N_STATE = 2 * N_GEN
N_INPUT = N_GEN + N_BUS
N_MEAS = N_GEN + N_BUS

DATA_DIR = Path(__file__).resolve().parent / "ieee14"


def load_plant(data_dir: Path | None = None):
    data_dir = Path(data_dir or DATA_DIR)
    plant = {
        "A": np.loadtxt(data_dir / "A.csv", delimiter=","),
        "B": np.loadtxt(data_dir / "B.csv", delimiter=","),
        "C": np.loadtxt(data_dir / "C.csv", delimiter=","),
        "D": np.loadtxt(data_dir / "D.csv", delimiter=","),
        "ts": TS,
        "n": N_STATE,
        "p": N_INPUT,
        "m": N_MEAS,
        "n_gen": N_GEN,
        "n_bus": N_BUS,
    }
    assert plant["A"].shape == (N_STATE, N_STATE)
    assert plant["B"].shape == (N_STATE, N_INPUT)
    assert plant["C"].shape == (N_MEAS, N_STATE)
    assert plant["D"].shape == (N_MEAS, N_INPUT)
    return plant


def plant_summary(plant=None):
    plant = plant or load_plant()
    a, b, c = plant["A"], plant["B"], plant["C"]
    ctrb = np.linalg.matrix_rank(
        np.hstack([np.linalg.matrix_power(a, i) @ b for i in range(a.shape[0])])
    )
    obsv = np.linalg.matrix_rank(
        np.vstack([c @ np.linalg.matrix_power(a, i) for i in range(a.shape[0])])
    )
    return {
        "n": plant["n"],
        "p": plant["p"],
        "m": plant["m"],
        "rank_C": int(np.linalg.matrix_rank(c)),
        "ctrb": int(ctrb),
        "obsv": int(obsv),
        "rho_A": float(np.max(np.abs(np.linalg.eigvals(a)))),
    }


if __name__ == "__main__":
    info = plant_summary()
    for key, value in info.items():
        print(f"{key}: {value}")
