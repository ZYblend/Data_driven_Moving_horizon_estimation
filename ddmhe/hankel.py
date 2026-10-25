"""Offline Hankel matrices and the persistence-of-excitation rank check."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def block_hankel(sig, depth):
    """Block-Hankel matrix with ``depth`` block rows."""
    n_samples, _ = sig.shape
    cols = n_samples - depth + 1
    if cols < 1:
        raise ValueError(f"Need at least {depth} samples to build a depth-{depth} Hankel matrix.")
    return np.column_stack([sig[j : j + depth].reshape(-1) for j in range(cols)])


def compress_hankel(hu, hy, hx, tol=1e-8):
    """Keep the joint numerical column space of the input, output, and state Hankel matrices."""
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


def build_compressed_hankel(u, y, x, horizon, tol=1e-8):
    hu, hy, hx = {}, {}, {}
    for depth in range(1, horizon + 1):
        hu[depth] = block_hankel(u, depth)
        hy[depth] = block_hankel(y, depth)
        hx[depth] = block_hankel(x, depth)
    return compress_hankel(hu, hy, hx, tol=tol)


def rank_analysis(u, y, x, horizon, tol=1e-8):
    """Rank of the input Hankel matrix of order L+n, and of the joint data matrix at depth L."""
    n_state = x.shape[1]
    n_input = u.shape[1]
    order = horizon + n_state
    hu = block_hankel(u, order)
    rank = int(np.linalg.matrix_rank(hu, tol=tol))
    rows = int(hu.shape[0])
    joint = np.vstack(
        [block_hankel(u, horizon), block_hankel(y, horizon), block_hankel(x, horizon)]
    )
    joint_rank = int(np.linalg.matrix_rank(joint, tol=tol))
    return {
        "samples": int(u.shape[0]),
        "horizon": int(horizon),
        "n_state": int(n_state),
        "n_input": int(n_input),
        "order": int(order),
        "input_rows": rows,
        "input_columns": int(hu.shape[1]),
        "input_rank": rank,
        "persistently_exciting": rank == rows,
        "joint_rows": int(joint.shape[0]),
        "joint_columns": int(joint.shape[1]),
        "joint_rank": joint_rank,
    }


def print_hankel_analysis(report, compressed_columns=None):
    """Print the rank check that decides whether the offline input is persistently exciting."""
    print("Hankel rank analysis", flush=True)
    print(
        f"  samples {report['samples']}, horizon L={report['horizon']}, "
        f"n={report['n_state']}, p={report['n_input']}",
        flush=True,
    )
    print(
        f"  input Hankel of order L+n={report['order']}: "
        f"rank {report['input_rank']} / {report['input_rows']} rows, "
        f"{report['input_columns']} columns",
        flush=True,
    )
    if report["persistently_exciting"]:
        print("  input is persistently exciting of order L+n", flush=True)
    else:
        deficit = report["input_rows"] - report["input_rank"]
        print(
            f"  input is rank deficient by {deficit}. "
            "A data-driven MHE dictionary built from it does not span every length-L trajectory.",
            flush=True,
        )
    print(
        f"  joint [Hu; Hy; Hx] at depth L: "
        f"rank {report['joint_rank']} / {report['joint_rows']} rows, "
        f"{report['joint_columns']} columns",
        flush=True,
    )
    if compressed_columns is not None:
        print(f"  compressed dictionary columns at depth L: {compressed_columns}", flush=True)


def save_hankel(path, u, y, x, hu, hy, hx, meta):
    """Save the offline trajectory, the compressed Hankel dictionary, and a JSON sidecar."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    horizon = int(meta["horizon"])
    payload = {"u": u, "y": y, "x": x, "horizon": np.array([horizon])}
    for key in ("scale",):
        if key in meta:
            payload[key] = np.asarray(meta[key], dtype=float)
    for depth in range(1, horizon + 1):
        payload[f"hu_{depth}"] = hu[depth]
        payload[f"hy_{depth}"] = hy[depth]
        payload[f"hx_{depth}"] = hx[depth]
    np.savez(path, **payload)
    sidecar = {k: v for k, v in meta.items() if k != "scale"}
    path.with_suffix(".json").write_text(json.dumps(sidecar, indent=2))
    print(f"wrote {path}", flush=True)


def load_hankel(path):
    """Load a dictionary written by :func:`save_hankel`."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Missing offline Hankel file {path}. Run the example's offline preparation script first."
        )
    data = np.load(path)
    horizon = int(data["horizon"][0])
    hu, hy, hx = {}, {}, {}
    for depth in range(1, horizon + 1):
        hu[depth] = data[f"hu_{depth}"]
        hy[depth] = data[f"hy_{depth}"]
        hx[depth] = data[f"hx_{depth}"]
    meta = json.loads(path.with_suffix(".json").read_text())
    if "scale" in data:
        meta["scale"] = data["scale"]
    return data["u"], data["y"], data["x"], hu, hy, hx, meta
