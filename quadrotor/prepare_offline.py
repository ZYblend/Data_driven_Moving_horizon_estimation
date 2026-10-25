"""Build the quadrotor offline Hankel dictionary.

Run from the repository root, before the online figure-8 simulation:

    python -m quadrotor.prepare_offline
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ddmhe import (  # noqa: E402
    build_compressed_hankel,
    print_hankel_analysis,
    rank_analysis,
    save_hankel,
)
from quadrotor.trajectory import (  # noqa: E402
    L,
    N_INPUT,
    N_OFFLINE,
    N_STATE,
    OFFLINE_PATH,
    P_WEIGHT,
    R_WEIGHT,
    RHO,
    TS,
    lqr_gain,
    offline_trajectory,
)


def main():
    need = (N_INPUT + 1) * (L + N_STATE) + 40
    n_samples = max(N_OFFLINE, need)
    gain = lqr_gain()
    u, x, y = offline_trajectory(gain, seed=0, n_samples=n_samples)
    report = rank_analysis(u, y, x, L)
    hu, hy, hx = build_compressed_hankel(u, y, x, L)
    print_hankel_analysis(report, compressed_columns=hu[L].shape[1])
    if not report["persistently_exciting"]:
        raise SystemExit("Offline input is not persistently exciting of order L+n.")
    meta = {
        "horizon": L,
        "ts": TS,
        "n_state": N_STATE,
        "n_input": N_INPUT,
        "p_weight": P_WEIGHT,
        "r_weight": R_WEIGHT,
        "rho": RHO,
        "rank": report,
    }
    save_hankel(OFFLINE_PATH, u, y, x, hu, hy, hx, meta)


if __name__ == "__main__":
    main()
