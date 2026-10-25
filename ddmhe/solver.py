"""Unified data-driven MHE. One factory builds the robust and resilient solvers."""

from __future__ import annotations

import warnings

import cvxpy as cp
import numpy as np


def _solve(problem, label):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Solution may be inaccurate")
        try:
            problem.solve(solver=cp.CLARABEL, warm_start=True, verbose=False)
        except cp.SolverError:
            problem.solve(solver=cp.CLARABEL, warm_start=False, verbose=False)
        if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
            try:
                problem.solve(
                    solver=cp.SCS,
                    warm_start=False,
                    verbose=False,
                    eps=1e-5,
                    max_iters=8000,
                )
            except cp.SolverError as exc:
                raise RuntimeError(f"{label} solve failed: {problem.status}") from exc
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise RuntimeError(f"{label} solve failed: {problem.status}")


class DataDrivenMHE:
    """Moving-horizon fit of a window to a Hankel dictionary.

    ``resilient=False`` penalizes the squared measurement residual (robust DD-MHE).
    ``resilient=True`` penalizes the absolute residual (resilient DD-MHE). Both use
    the same forgetting factor and the same arrival cost on the initial state.
    """

    def __init__(self, hu, hy, hx, horizon, resilient, p_weight=1.0, r_weight=10.0, rho=0.95):
        self.resilient = resilient
        self.horizon = horizon
        n = hx[1].shape[0]
        m = hy[1].shape[0]
        p = hu[1].shape[0]
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
            arrival = (p_weight * rho ** (depth - 1)) * tau
            if resilient:
                slack = cp.Variable(m * depth, nonneg=True)
                weights = np.repeat([r_weight * rho ** (depth - 1 - i) for i in range(depth)], m)
                constraints.append(cp.abs(residual) <= slack)
                cost = arrival + slack @ weights
            else:
                weights = np.array([r_weight * rho ** (depth - 1 - i) for i in range(depth)])
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
        u_par.value = np.asarray(u_win, dtype=float).ravel()
        y_par.value = np.asarray(y_win, dtype=float).reshape(-1)
        prior_par.value = np.asarray(prior, dtype=float)
        label = "resilient-dd" if self.resilient else "robust-dd"
        _solve(problem, label)
        return (hx @ alpha.value).reshape(depth, -1)


def make_dd_mhe(hu, hy, hx, *, horizon, resilient, p_weight=1.0, r_weight=10.0, rho=0.95):
    """Build one data-driven MHE solver.

    ``resilient=False`` is Robust DD-MHE. ``resilient=True`` is Resilient DD-MHE.
    """
    return DataDrivenMHE(
        hu,
        hy,
        hx,
        horizon=horizon,
        resilient=resilient,
        p_weight=p_weight,
        r_weight=r_weight,
        rho=rho,
    )


def make_dd_mhe_pair(hu, hy, hx, *, horizon, p_weight=1.0, r_weight=10.0, rho=0.95):
    """Robust and resilient solvers for the same offline dictionary."""
    shared = dict(hu=hu, hy=hy, hx=hx, horizon=horizon, p_weight=p_weight, r_weight=r_weight, rho=rho)
    return {
        "Robust DD-MHE": make_dd_mhe(resilient=False, **shared),
        "Resilient DD-MHE": make_dd_mhe(resilient=True, **shared),
    }


def run_windowed(estimate, u, y, prior0, horizon):
    """Replay ``estimate`` on a recorded input and output sequence."""
    n = prior0.shape[0]
    x_hat = np.zeros((u.shape[0], n))
    latest = {0: np.asarray(prior0, dtype=float).copy()}
    for t in range(u.shape[0]):
        if t + 1 < horizon:
            t0 = 0
            prior = prior0
        else:
            t0 = t - horizon + 1
            prior = latest[t0]
        window = estimate(u[t0 : t + 1], y[t0 : t + 1], prior)
        for k, time_index in enumerate(range(t0, t + 1)):
            latest[time_index] = window[k]
        x_hat[t] = window[-1]
    return x_hat
