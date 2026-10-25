"""Open-loop MHE on the ANDES single-machine infinite-bus case.

The plant is the textbook SMIB shipped with ANDES (classical machine against
an infinite bus). The input is a known deviation of mechanical torque. Outputs
are rotor angle, rotor speed, and electrical power. A constant bias is added
to the speed measurement partway through the run.

ANDES is imported from the installed package. Its source is GPL-3 and is not
vendored here. Run from the repository root:

    python3 power/explore_smib.py

The script writes two figures. power/smib_mhe.png is the open-loop estimate.
power/smib_closed_loop.png, also copied to image/smib_mhe.png, is the machine
trajectory when that estimate drives the torque command: a feedforward torque
for a slow angle move, minus a linear gain on the estimated tracking error,
applied one step after the estimate. Model-driven and data-driven results
agree at the printed precision. The closed-loop figure shows the data-driven
pair. The data-driven dictionary is the
linearized SMIB, because a Hankel matrix built from the nonlinear samples does
not have a sharp rank gap and will absorb a speed bias as if it were real.
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import json
import sys
import types
import warnings
from pathlib import Path

import cvxpy as cp
import matplotlib.pyplot as plt
import numpy as np
import andes
from numpy.linalg import lstsq
from scipy.linalg import expm

andes.config_logger(stream_level=40)

if "smib_cl_hook" not in sys.modules:
    sys.modules["smib_cl_hook"] = types.ModuleType("smib_cl_hook")
sys.modules["smib_cl_hook"].HOOK = None

TS = 0.02
L = 25
P_WEIGHT = 1.0
R_WEIGHT = 10.0
RHO = 0.95
NOISE_STD = 0.1
# Physical units of one normalized measurement: rad, pu speed, pu power.
SCALE = np.array([0.01, 2.0e-4, 0.02])
SPEED_BIAS = 1.0e-3  # pu, about 0.06 Hz
ATTACK_TIME = 2.0
STEP_TIME = 0.5
STEP = 0.02  # pu mechanical torque
T_FINAL = 6.0
# Closed loop: a slow angle move, then a speed bias while torque is in feedback.
CL_T_FINAL = 8.0
CL_ATTACK = 3.0
MOVE_START = 0.5
MOVE_END = 2.2
DELTA_GOAL = 0.012  # rad, about a 0.02 pu torque step
U_LIMIT = 0.08
N_OFFLINE = 400
HOLD = 2  # samples per PRBS level
PRIOR = np.zeros(2)


def hankel(sig, depth):
    n_samples, _ = sig.shape
    cols = n_samples - depth + 1
    return np.column_stack([sig[j : j + depth].reshape(-1) for j in range(cols)])


def load_case():
    data = json.loads(Path(andes.get_case("smib/SMIB.json")).read_text())
    data["Fault"][0]["u"] = 0.0
    path = Path("/tmp/smib_nofault.json")
    path.write_text(json.dumps(data))
    ss = andes.load(str(path), setup=True, no_output=True, default_config=True)
    ss.PFlow.run()
    if not ss.PFlow.converged:
        raise RuntimeError("SMIB power flow did not converge")
    return ss


def linear_model(ss):
    ss.EIG.run()
    full = np.array(ss.EIG.As, dtype=float)
    # GENCLS 1 is the finite machine. GENCLS 2 is the infinite bus.
    index = [0, 2]
    A = full[np.ix_(index, index)]
    inertia = float(ss.GENCLS.M.v[0])
    B = np.array([[0.0], [1.0 / inertia]])
    return A, B


def write_schedule(u, tm_eq):
    path = "/tmp/smib_schedule.npz"
    np.savez(path, u=np.asarray(u, dtype=float), ts=TS, tm_eq=tm_eq)
    os.environ["SMIB_U"] = path
    pert = Path("/tmp/smib_pert.py")
    pert.write_text(
        "import os\n"
        "import numpy as np\n"
        "def pert(t, system):\n"
        "    z = np.load(os.environ['SMIB_U'])\n"
        "    k = int(round(float(t) / float(z['ts'])))\n"
        "    k = max(0, min(k, len(z['u']) - 1))\n"
        "    system.GENCLS.tm0.v[0] = float(z['tm_eq'] + z['u'][k])\n"
    )
    return pert


def equilibrium_torque():
    ss = load_case()
    ss.TDS.init()
    return float(ss.GENCLS.tm0.v[0])


def simulate(u, tm_eq):
    """Run ANDES and resample states and power onto the uniform input grid."""
    ss = load_case()
    pert = write_schedule(u, tm_eq)
    ss.files.pert = str(pert)
    ss.TDS.config.tf = TS * (len(u) - 1)
    ss.TDS.config.tstep = TS
    ss.TDS.config.no_tqdm = 1
    ss.TDS.config.criteria = 0
    ss.TDS.run()
    if ss.TDS.busted:
        raise RuntimeError(ss.TDS.err_msg)
    ss.dae.ts.unpack()
    t = ss.dae.ts.t
    grid = np.arange(len(u)) * TS
    delta = np.interp(grid, t, ss.dae.ts.x[:, 0])
    omega = np.interp(grid, t, ss.dae.ts.x[:, 2])
    pe_index = list(ss.dae.y_name).index("Pe GENCLS 1")
    pe = np.interp(grid, t, ss.dae.ts.y[:, pe_index])
    x = np.column_stack([delta - delta[0], omega - omega[0]])
    y_phys = np.column_stack([x[:, 0], x[:, 1], pe - pe[0]])
    dt = np.diff(t)
    return x, y_phys, tm_eq, (float(dt.min()), float(dt.max()), len(t))


def output_matrices(x, u, y_phys):
    """Electrical power is an algebraic output. Fit its row on the offline run."""
    phi = np.column_stack([x, u])
    coef, residual, _, _ = lstsq(phi, y_phys[:, 2], rcond=None)
    fit = phi @ coef
    rmse = float(np.sqrt(np.mean((y_phys[:, 2] - fit) ** 2)))
    signal = float(np.max(np.abs(y_phys[:, 2])))
    C = np.array([[1.0, 0.0], [0.0, 1.0], coef[:2]], dtype=float)
    D = np.array([[0.0], [0.0], [coef[2]]], dtype=float)
    return C, D, rmse, signal


def discretize(A, B):
    n = A.shape[0]
    block = np.zeros((n + 1, n + 1))
    block[:n, :n] = A * TS
    block[:n, n:] = B * TS
    discrete = expm(block)
    return discrete[:n, :n], discrete[:n, n:]


def compress(hu, hy, hx, tol=1e-8):
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


def solve(problem, label):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Solution may be inaccurate")
        try:
            problem.solve(solver=cp.CLARABEL, warm_start=True, verbose=False)
            if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
                problem.solve(solver=cp.CLARABEL, warm_start=False, verbose=False)
        except cp.SolverError:
            problem.solve(solver=cp.SCS, warm_start=False, verbose=False, eps=1e-5, max_iters=8000)
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise RuntimeError(f"{label} solve failed: {problem.status}")


class DataDriven:
    def __init__(self, hu, hy, hx, n, m):
        self.problems = {}
        for depth in range(1, L + 1):
            alpha = cp.Variable(hu[depth].shape[1])
            u_par = cp.Parameter(depth)  # one input
            y_par = cp.Parameter(m * depth)
            prior = cp.Parameter(n)
            x0 = hx[depth][:n] @ alpha
            residual = y_par - hy[depth] @ alpha
            tau = cp.Variable(nonneg=True)
            slack = cp.Variable(m * depth, nonneg=True)
            weights = np.repeat([R_WEIGHT * RHO ** (depth - 1 - i) for i in range(depth)], m)
            constraints = [
                hu[depth] @ alpha == u_par,
                cp.norm(x0 - prior, 2) <= tau,
                cp.abs(residual) <= slack,
            ]
            problem = cp.Problem(
                cp.Minimize((P_WEIGHT * RHO ** (depth - 1)) * tau + slack @ weights),
                constraints,
            )
            self.problems[depth] = (problem, alpha, u_par, y_par, prior, hx[depth])

    def estimate(self, u_win, y_win, prior):
        depth = u_win.shape[0]
        problem, alpha, u_par, y_par, prior_par, hx = self.problems[depth]
        u_par.value = u_win.ravel()
        y_par.value = y_win.reshape(-1)
        prior_par.value = prior
        solve(problem, "data-driven")
        return (hx @ alpha.value).reshape(depth, -1)


class ModelDriven:
    def __init__(self, Ad, Bd, C, D, powers):
        self.Ad = Ad
        self.Bd = Bd
        self.C = C
        self.D = D
        self.powers = powers
        self.problems = {}
        n = Ad.shape[0]
        m = C.shape[0]
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
            problem = cp.Problem(
                cp.Minimize((P_WEIGHT * RHO ** (depth - 1)) * tau + stage),
                constraints,
            )
            self.problems[depth] = (problem, x0, g_par, y_par, prior)

    def estimate(self, u_win, y_win, prior):
        depth = u_win.shape[0]
        problem, x0, g_par, y_par, prior_par = self.problems[depth]
        g = np.zeros((depth, self.Ad.shape[0]))
        acc = np.zeros(self.Ad.shape[0])
        for i in range(depth):
            g[i] = acc
            acc = self.Ad @ acc + self.Bd.ravel() * u_win[i, 0]
        g_par.value = g
        y_par.value = y_win - u_win @ self.D.T
        prior_par.value = prior
        solve(problem, "model")
        return np.vstack([self.powers[i] @ x0.value + g[i] for i in range(depth)])


def robust_dd(hu, hy, hx, u_win, y_win, prior):
    depth = u_win.shape[0]
    n = prior.shape[0]
    m = y_win.shape[1]
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
    kkt = np.block([[gram, hu.T], [hu, np.zeros((hu.shape[0], hu.shape[0]))]])
    sol = np.linalg.solve(kkt, np.concatenate([design.T @ target, u_win.ravel()]))
    return (hx @ sol[: hu.shape[1]]).reshape(depth, n)


def robust_md(Ad, Bd, C, D, powers, u_win, y_win, prior):
    depth = u_win.shape[0]
    n = prior.shape[0]
    g = np.zeros((depth, n))
    acc = np.zeros(n)
    for i in range(depth):
        g[i] = acc
        acc = Ad @ acc + Bd.ravel() * u_win[i, 0]
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
    x0, *_ = lstsq(design, target, rcond=None)
    return np.vstack([powers[i] @ x0 + g[i] for i in range(depth)])


def run_estimator(estimate, u, y):
    """u is physical torque deviation. y is scaled. The estimate is physical."""
    x_hat = np.zeros((u.shape[0], 2))
    latest = {0: PRIOR.copy()}
    for t in range(u.shape[0]):
        if t + 1 < L:
            t0, prior = 0, PRIOR
        else:
            t0, prior = t - L + 1, latest[t0]
        window = estimate(u[t0 : t + 1], y[t0 : t + 1], prior)
        for k, time_index in enumerate(range(t0, t + 1)):
            latest[time_index] = window[k]
        x_hat[t] = window[-1]
    return x_hat


def rmse(err):
    return float(np.sqrt(np.mean(err ** 2)))


def lqr_gain(Ad, Bd):
    """Discrete LQR. Weights are Bryson scales: 0.01 rad, 5e-4 pu, 0.02 pu."""
    from scipy.linalg import solve_discrete_are

    weight = np.diag([(1.0 / 0.01) ** 2, (1.0 / 5.0e-4) ** 2])
    effort = np.array([[(1.0 / 0.02) ** 2]])
    gram = solve_discrete_are(Ad, Bd, weight, effort)
    return np.linalg.solve(effort + Bd.T @ gram @ Bd, Bd.T @ gram @ Ad).ravel()


def reference_trajectory(A, B, n_steps):
    """Slow angle move that the linear model can track with a torque feedforward."""
    t = np.arange(n_steps) * TS
    span = MOVE_END - MOVE_START
    z = np.clip((t - MOVE_START) / span, 0.0, 1.0)
    shape = z * z * (3.0 - 2.0 * z)
    dshape = np.where((t > MOVE_START) & (t < MOVE_END), 6.0 * z * (1.0 - z) / span, 0.0)
    ddshape = np.where((t > MOVE_START) & (t < MOVE_END), (6.0 - 12.0 * z) / span ** 2, 0.0)
    delta = DELTA_GOAL * shape
    omega = DELTA_GOAL * dshape / A[0, 1]
    domega = DELTA_GOAL * ddshape / A[0, 1]
    torque = (domega - A[1, 0] * delta - A[1, 1] * omega) / B[1, 0]
    return np.column_stack([delta, omega]), torque


def _closed_loop_pert_path():
    path = Path("/tmp/smib_cl_pert.py")
    path.write_text(
        "def pert(t, system):\n"
        "    import smib_cl_hook\n"
        "    hook = smib_cl_hook.HOOK\n"
        "    if hook is not None:\n"
        "        hook(t, system)\n"
    )
    return path


def closed_loop_andes(estimate, gain, x_ref, u_ff, C_phys, D_phys, attack, noise):
    """Torque command from the previous estimate, one control step of delay.

    The command is the feedforward torque minus a linear gain on the estimated
    tracking error. Measurements come from ANDES. The estimate is updated after
    the command is chosen, then the integrator advances one step.
    """
    n_steps = x_ref.shape[0]
    records = {
        "t": [],
        "x": [],
        "u": [],
        "x_hat": [],
    }
    hold = {"torque": 0.0, "k": 0, "next_t": 0.0, "origin": None, "hat": np.zeros(2), "latest": {0: np.zeros(2)}}
    applied = np.zeros((n_steps, 1))
    measured = np.zeros((n_steps, 3))

    def hook(t, system):
        if hold["origin"] is None:
            names_x = list(system.dae.x_name)
            names_y = list(system.dae.y_name)
            hold["origin"] = {
                "uid": system.GENCLS.idx2uid("GENCLS_1"),
                "tm": float(system.GENCLS.tm0.v[system.GENCLS.idx2uid("GENCLS_1")]),
                "i_delta": names_x.index("delta GENCLS 1"),
                "i_omega": names_x.index("omega GENCLS 1"),
                "i_pe": names_y.index("Pe GENCLS 1"),
                "delta": float(system.dae.x[names_x.index("delta GENCLS 1")]),
                "omega": float(system.dae.x[names_x.index("omega GENCLS 1")]),
                "pe": float(system.dae.y[names_y.index("Pe GENCLS 1")]),
            }
        origin = hold["origin"]
        system.GENCLS.tm0.v[origin["uid"]] = origin["tm"] + hold["torque"]
        if hold["k"] >= n_steps or t < hold["next_t"] - 1e-4:
            return
        k = hold["k"]
        delta = float(system.dae.x[origin["i_delta"]]) - origin["delta"]
        omega = float(system.dae.x[origin["i_omega"]]) - origin["omega"]
        power = float(system.dae.y[origin["i_pe"]]) - origin["pe"]
        state = np.array([delta, omega])
        # Command from the estimate of the previous sample.
        torque = float(u_ff[k] - gain @ (hold["hat"] - x_ref[k]))
        torque = float(np.clip(torque, -U_LIMIT, U_LIMIT))
        hold["torque"] = torque
        system.GENCLS.tm0.v[origin["uid"]] = origin["tm"] + torque
        y = np.array([delta, omega, power]) / SCALE + noise[k]
        if attack and t >= CL_ATTACK - 1e-9:
            y = y.copy()
            y[1] += SPEED_BIAS / SCALE[1]
        applied[k, 0] = torque
        measured[k] = y
        if k + 1 < L:
            t0, prior = 0, np.zeros(2)
        else:
            t0 = k - L + 1
            prior = hold["latest"][t0]
        window = estimate(applied[t0 : k + 1], measured[t0 : k + 1], prior)
        for j, time_index in enumerate(range(t0, k + 1)):
            hold["latest"][time_index] = window[j]
        hold["hat"] = window[-1].copy()
        records["t"].append(float(t))
        records["x"].append(state.copy())
        records["u"].append(torque)
        records["x_hat"].append(hold["hat"].copy())
        hold["k"] = k + 1
        hold["next_t"] += TS

    sys.modules["smib_cl_hook"].HOOK = hook
    ss = load_case()
    ss.files.pert = str(_closed_loop_pert_path())
    ss.TDS.config.tf = CL_T_FINAL
    ss.TDS.config.tstep = TS
    ss.TDS.config.no_tqdm = 1
    ss.TDS.config.criteria = 0
    ss.TDS.run()
    sys.modules["smib_cl_hook"].HOOK = None
    if ss.TDS.busted:
        raise RuntimeError(ss.TDS.err_msg or "ANDES time-domain simulation failed")
    if hold["k"] < n_steps:
        raise RuntimeError(f"closed loop recorded {hold['k']} of {n_steps} samples")
    return {
        "t": np.asarray(records["t"]),
        "x": np.vstack(records["x"]),
        "u": np.asarray(records["u"]),
        "x_hat": np.vstack(records["x_hat"]),
        "x_ref": x_ref,
    }


def save_closed_loop_figure(t, x_ref, flights, path):
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.4), constrained_layout=True)
    series = [("Reference", x_ref, "#1f77b4", 1.6)]
    colors = {
        "Robust MD-MHE": "#d95f02",
        "Resilient MD-MHE": "#1b9e77",
        "Robust DD-MHE": "#d95f02",
        "Resilient DD-MHE": "#1b9e77",
    }
    for name, traj in flights:
        series.append((name, traj, colors[name], 1.2))
    for name, traj, color, width in series:
        axes[0].plot(t, 1e3 * traj[:, 0], color=color, lw=width, label=name)
        axes[1].plot(t, 60e3 * traj[:, 1], color=color, lw=width, label=name)
    for ax in axes:
        ax.axvline(CL_ATTACK, color="0.4", ls="--", lw=0.8)
        ax.set_xlabel("Time (s)")
        ax.set_xlim(0.0, t[-1])
    axes[0].set_ylabel("Rotor angle (mrad)")
    axes[1].set_ylabel("Speed (mHz)")
    axes[0].legend(loc="upper center", bbox_to_anchor=(1.15, 1.38), ncol=3, frameon=False)
    fig.savefig(path, dpi=160)
    print(f"wrote {path}")


def main():
    probe = load_case()
    A, B = linear_model(probe)
    Ad, Bd = discretize(A, B)
    print("A", A)
    print("B", B.ravel())
    print("Ad eigs", np.linalg.eigvals(Ad))

    rng = np.random.default_rng(0)
    levels = rng.uniform(-STEP, STEP, size=int(np.ceil(N_OFFLINE / HOLD)))
    u_off = np.repeat(levels, HOLD)[:N_OFFLINE]
    tm_eq = equilibrium_torque()
    x_off, y_off, _, dt_info = simulate(u_off, tm_eq)
    print(f"offline tm_eq {tm_eq:.4f} raw dt range {dt_info}")
    print(
        f"offline max |delta| {np.max(np.abs(x_off[:, 0])):.4e} rad, "
        f"|omega| {np.max(np.abs(x_off[:, 1])):.4e} pu"
    )
    C_phys, D_phys, pe_rmse, pe_span = output_matrices(x_off, u_off, y_off)
    print(f"power fit rmse {pe_rmse:.4e} on a span of {pe_span:.4e}")
    print("C", C_phys)
    print("D", D_phys.ravel())

    # One-step linear prediction of the offline trajectory.
    x_lin = np.zeros_like(x_off)
    for k in range(len(u_off) - 1):
        x_lin[k + 1] = Ad @ x_lin[k] + Bd.ravel() * u_off[k]
    print(
        "linear one-step is free-run from 0; "
        f"rmse delta {rmse(x_off[:, 0] - x_lin[:, 0]):.4e}, "
        f"omega {rmse(x_off[:, 1] - x_lin[:, 1]):.4e}"
    )

    order = L + 2
    rank = np.linalg.matrix_rank(hankel(u_off.reshape(-1, 1), order), tol=1e-8)
    print(f"input Hankel rank {rank}/{order}")

    C = C_phys / SCALE[:, None]
    D = D_phys / SCALE[:, None]
    # The nonlinear samples do not have a sharp Hankel gap, so the data-driven
    # dictionary is this linear rollout, not the ANDES samples above.
    y_lin = x_lin @ C_phys.T + u_off.reshape(-1, 1) @ D_phys.T
    y_off_s = y_lin / SCALE
    x_off_s = x_lin
    u_off_s = u_off.reshape(-1, 1)
    hu, hy, hx = {}, {}, {}
    for depth in range(1, L + 1):
        hu[depth] = hankel(u_off_s, depth)
        hy[depth] = hankel(y_off_s, depth)
        hx[depth] = hankel(x_off_s, depth)
    hu, hy, hx = compress(hu, hy, hx)
    print(f"compressed columns at L: {hu[L].shape[1]}")

    powers = [np.eye(2)]
    for _ in range(L - 1):
        powers.append(Ad @ powers[-1])
    md = ModelDriven(Ad, Bd, C, D, powers)
    dd = DataDriven(hu, hy, hx, n=2, m=3)

    n_on = int(round(T_FINAL / TS)) + 1
    u_on = np.zeros(n_on)
    u_on[int(round(STEP_TIME / TS)) :] = STEP
    x_on, y_on, _, dt_on = simulate(u_on, tm_eq)
    print(f"online raw dt range {dt_on}")
    y_clean = y_on / SCALE
    attack_index = int(round(ATTACK_TIME / TS))
    rng = np.random.default_rng(1)
    y_noise = y_clean + rng.normal(0.0, NOISE_STD, size=y_clean.shape)

    def attacked(y):
        out = y.copy()
        out[attack_index:, 1] += SPEED_BIAS / SCALE[1]
        return out

    bank_builder = {
        "Robust MD-MHE": lambda: (
            lambda u_win, y_win, prior: robust_md(
                Ad, Bd, C, D, powers, u_win, y_win, prior
            )
        ),
        "Resilient MD-MHE": lambda: md.estimate,
        "Robust DD-MHE": lambda: (
            lambda u_win, y_win, prior: robust_dd(
                hu[u_win.shape[0]], hy[u_win.shape[0]], hx[u_win.shape[0]], u_win, y_win, prior
            )
        ),
        "Resilient DD-MHE": lambda: dd.estimate,
    }

    results = {}
    hats = {}
    for case, y in (("none", y_noise), ("speed", attacked(y_noise))):
        results[case] = {}
        hats[case] = {}
        for name, build in bank_builder.items():
            x_hat = run_estimator(build(), u_on.reshape(-1, 1), y)
            sl = slice(attack_index, None)
            metrics = {
                "delta_mrad": 1e3 * rmse(x_on[sl, 0] - x_hat[sl, 0]),
                "omega_mHz": 60e3 * rmse(x_on[sl, 1] - x_hat[sl, 1]),
            }
            results[case][name] = metrics
            hats[case][name] = x_hat
            print(
                f"{case:6s} {name:18s} "
                f"angle {metrics['delta_mrad']:.2f} mrad  "
                f"speed {metrics['omega_mHz']:.3f} mHz",
                flush=True,
            )
    t = np.arange(n_on) * TS
    figure = Path(__file__).resolve().parent / "smib_mhe.png"
    save_figure(t, x_on, hats["speed"], figure)

    gain = lqr_gain(Ad, Bd)
    n_cl = int(round(CL_T_FINAL / TS)) + 1
    x_ref, u_ff = reference_trajectory(A, B, n_cl)
    print(
        f"LQR gain {gain}, feedforward peak {np.max(np.abs(u_ff)):.4f} pu",
        flush=True,
    )
    attack_cl = int(round(CL_ATTACK / TS))
    drawn = []
    for case, attack in (("none", False), ("speed", True)):
        noise = np.random.default_rng(1).normal(0.0, NOISE_STD, size=(n_cl, 3))
        md = ModelDriven(Ad, Bd, C, D, powers)
        dd = DataDriven(hu, hy, hx, n=2, m=3)
        runs = (
            ("Robust MD-MHE", lambda u, y, p: robust_md(Ad, Bd, C, D, powers, u, y, p)),
            ("Resilient MD-MHE", md.estimate),
            (
                "Robust DD-MHE",
                lambda u, y, p: robust_dd(
                    hu[u.shape[0]], hy[u.shape[0]], hx[u.shape[0]], u, y, p
                ),
            ),
            ("Resilient DD-MHE", dd.estimate),
        )
        for name, estimate in runs:
            out = closed_loop_andes(
                estimate, gain, x_ref, u_ff, C_phys, D_phys, attack, noise
            )
            sl = slice(attack_cl, None)
            track = 1e3 * rmse(out["x"][sl, 0] - x_ref[sl, 0])
            final = 1e3 * out["x"][-1, 0]
            print(
                f"closed {case:6s} {name:18s} "
                f"angle {track:.2f} mrad  final {final:.2f} mrad  "
                f"samples {len(out['t'])}  dt {np.diff(out['t']).max():.4f}",
                flush=True,
            )
            if case == "speed" and name in ("Robust DD-MHE", "Resilient DD-MHE"):
                drawn.append((name, out["x"]))
    closed_figure = Path(__file__).resolve().parent / "smib_closed_loop.png"
    paper_figure = Path(__file__).resolve().parents[1] / "image" / "smib_mhe.png"
    save_closed_loop_figure(np.arange(n_cl) * TS, x_ref, drawn, closed_figure)
    save_closed_loop_figure(np.arange(n_cl) * TS, x_ref, drawn, paper_figure)
    return results


def save_figure(t, x, hats, path):
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.4), constrained_layout=True)
    series = (
        ("True", x, "#1f77b4", 1.6),
        ("Robust MD-MHE", hats["Robust MD-MHE"], "#d95f02", 1.2),
        ("Resilient MD-MHE", hats["Resilient MD-MHE"], "#1b9e77", 1.2),
    )
    for name, traj, color, width in series:
        axes[0].plot(t, 1e3 * traj[:, 0], color=color, lw=width, label=name)
        axes[1].plot(t, 60e3 * traj[:, 1], color=color, lw=width, label=name)
    for ax in axes:
        ax.axvline(ATTACK_TIME, color="0.4", ls="--", lw=0.8)
        ax.set_xlabel("Time (s)")
        ax.set_xlim(0, t[-1])
    axes[0].set_ylabel("Rotor angle (mrad)")
    axes[1].set_ylabel("Speed (mHz)")
    axes[0].legend(loc="upper center", bbox_to_anchor=(1.15, 1.38), ncol=3, frameon=False)
    fig.savefig(path, dpi=160)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
