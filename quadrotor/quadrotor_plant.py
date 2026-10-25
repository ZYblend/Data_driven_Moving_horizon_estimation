"""Planar quadrotor linearized at hover, as in safe-control-gym.

The continuous model is the TWO_D quadrotor in
safe_control_gym.envs.gym_pybullet_drones.quadrotor, with the Crazyflie
parameters that library loads from its URDF. States are deviations from hover,

    [x, x_dot, z, z_dot, theta, theta_dot],

and the inputs are deviations of the two rotor thrusts from mg/2. The output
is horizontal position, vertical position, pitch, pitch rate, and the vertical
specific force of an accelerometer at the center of mass. That specific force
depends directly on the two thrusts.
"""

import numpy as np
from scipy.linalg import expm

# Crazyflie parameters used by safe-control-gym / gym-pybullet-drones.
MASS = 0.027
GRAVITY = 9.8  # safe-control-gym BaseAviary.GRAVITY_ACC
ARM = 0.0397
IYY = 1.4e-5
TS = 0.02  # 50 Hz control interval

# y = [x, z, theta, theta_dot, f_z]
L_HORIZON = 2


def continuous_matrices():
    g, m, length, inertia = GRAVITY, MASS, ARM, IYY
    moment = length / (inertia * np.sqrt(2.0))
    ac = np.array(
        [
            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, g, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    )
    bc = np.array(
        [
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [1.0 / m, 1.0 / m],
            [0.0, 0.0],
            [-moment, moment],
        ]
    )
    # Accelerometer row is the only direct term: f_z = (T1 + T2) / m.
    cc = np.zeros((5, 6))
    cc[0, 0] = 1.0
    cc[1, 2] = 1.0
    cc[2, 4] = 1.0
    cc[3, 5] = 1.0
    dc = np.zeros((5, 2))
    dc[4, :] = 1.0 / m
    return ac, bc, cc, dc


def discrete_plant():
    ac, bc, cc, dc = continuous_matrices()
    n, p = bc.shape
    block = np.zeros((n + p, n + p))
    block[:n, :n] = ac
    block[:n, n:] = bc
    step = expm(block * TS)
    return step[:n, :n], step[:n, n:], cc, dc
