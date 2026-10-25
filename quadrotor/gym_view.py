"""Live safe-control-gym view of the figure-8 while the estimators run.

States are deviations from hover, ``[x, x_dot, z, z_dot, theta, theta_dot]``.
The camera starts on the negative y-axis so the vertical figure-8 faces the
window. It is not moved again, so the mouse can orbit, pan, and zoom.
"""

from __future__ import annotations

import os
import sys
import threading
import types
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore", category=UserWarning, module=r"gymnasium\.spaces\.box")


def _quiet_native():
    """Hide PyBullet's C++ startup lines. They bypass ``sys.stdout``."""
    class _Quiet:
        def __enter__(self):
            self._devnull = os.open(os.devnull, os.O_WRONLY)
            self._saved = (os.dup(1), os.dup(2))
            os.dup2(self._devnull, 1)
            os.dup2(self._devnull, 2)
            return self

        def __exit__(self, *_exc):
            os.dup2(self._saved[0], 1)
            os.dup2(self._saved[1], 2)
            os.close(self._saved[0])
            os.close(self._saved[1])
            os.close(self._devnull)
            return False

    return _Quiet()

HOVER_Z = 1.0
ROBUST = "Robust DD-MHE"
RESILIENT = "Resilient DD-MHE"
ROBUST_RGBA = (0.85, 0.47, 0.02, 1.0)
RESILIENT_RGBA = (0.08, 0.50, 0.16, 1.0)


def _import_quadrotor():
    gym_root = Path(__file__).resolve().parent / "safe-control-gym"
    if not gym_root.is_dir():
        raise FileNotFoundError(
            "safe-control-gym is missing. From the repository root, run "
            "git submodule update --init quadrotor/safe-control-gym"
        )
    root = str(gym_root)
    if root not in sys.path:
        sys.path.insert(0, root)
    # normalization.py imports torch at module scope. This playback never calls it.
    if "torch" not in sys.modules:
        try:
            import torch  # noqa: F401
        except ImportError:
            stand_in = types.ModuleType("torch")

            class _Tensor:
                pass

            stand_in.Tensor = _Tensor
            sys.modules["torch"] = stand_in
    from safe_control_gym.envs.gym_pybullet_drones.quadrotor import Quadrotor

    return Quadrotor


def _paint(pybullet, body, rgba, client):
    link_count = pybullet.getNumJoints(body, physicsClientId=client)
    for link in range(-1, link_count):
        pybullet.changeVisualShape(body, link, rgbaColor=list(rgba), physicsClientId=client)


def _pose(pybullet, body, state, y, client):
    position = [float(state[0]), float(y), float(HOVER_Z + state[2])]
    orientation = pybullet.getQuaternionFromEuler([0.0, float(state[4]), 0.0])
    pybullet.resetBasePositionAndOrientation(body, position, orientation, physicsClientId=client)
    pybullet.resetBaseVelocity(body, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], physicsClientId=client)


def _draw_course(pybullet, client, reference):
    bead = pybullet.createVisualShape(
        pybullet.GEOM_SPHERE,
        radius=0.012,
        rgbaColor=[0.25, 0.45, 0.85, 0.85],
        physicsClientId=client,
    )
    count = reference.shape[0]
    for index in range(0, count, max(count // 64, 1)):
        state = reference[index]
        pybullet.createMultiBody(
            baseMass=0,
            baseVisualShapeIndex=bead,
            basePosition=[float(state[0]), 0.0, float(HOVER_Z + state[2])],
            physicsClientId=client,
        )
    marker_shape = pybullet.createVisualShape(
        pybullet.GEOM_SPHERE,
        radius=0.02,
        rgbaColor=[0.85, 0.15, 0.12, 1.0],
        physicsClientId=client,
    )
    marker = pybullet.createMultiBody(
        baseMass=0,
        baseVisualShapeIndex=marker_shape,
        basePosition=[0.0, 0.0, HOVER_Z],
        physicsClientId=client,
    )
    return marker


class LiveView:
    """PyBullet window updated from the estimator thread.

    Call :meth:`push` from the solver thread and :meth:`pump` from the thread
    that opened the window. ``pump`` steps the simulator so mouse orbit, pan,
    and zoom are applied while a solve is in progress.
    """

    def __init__(self, sample_time, reference):
        import pybullet as p

        self._p = p
        reference = np.asarray(reference, dtype=float)
        ctrl = max(int(round(1.0 / sample_time)), 1)
        with _quiet_native():
            Quadrotor = _import_quadrotor()
            self._env = Quadrotor(
                quad_type=2,
                physics="dyn",
                gui=True,
                ctrl_freq=ctrl,
                pyb_freq=ctrl * 4,
                episode_len_sec=int(reference.shape[0] * sample_time) + 2,
                randomized_init=False,
                done_on_out_of_bound=False,
                init_state=np.array([0.0, 0.0, HOVER_Z, 0.0, 0.0, 0.0], dtype=np.float32),
            )
            self._env.reset()
            self._client = self._env.PYB_CLIENT
            p.setGravity(0.0, 0.0, 0.0, physicsClientId=self._client)
            # Left-drag otherwise grabs the floor, so the camera never orbits.
            p.configureDebugVisualizer(p.COV_ENABLE_MOUSE_PICKING, 0, physicsClientId=self._client)
            p.configureDebugVisualizer(p.COV_ENABLE_SHADOWS, 0, physicsClientId=self._client)
            # yaw=0, pitch=0 looks along +y: the x-z figure-8 faces the window.
            p.resetDebugVisualizerCamera(
                cameraDistance=2.6,
                cameraYaw=0,
                cameraPitch=0,
                cameraTargetPosition=[0.0, 0.0, HOVER_Z],
                physicsClientId=self._client,
            )
            p.stepSimulation(physicsClientId=self._client)
        _paint(p, self._env.DRONE_ID, ROBUST_RGBA, self._client)
        self._other = p.loadURDF(
            self._env.URDF_PATH,
            [0.0, 0.22, HOVER_Z],
            p.getQuaternionFromEuler([0.0, 0.0, 0.0]),
            flags=p.URDF_USE_INERTIA_FROM_FILE,
            physicsClientId=self._client,
        )
        _paint(p, self._other, RESILIENT_RGBA, self._client)
        p.addUserDebugText(
            "Robust",
            [0.0, 0.0, 0.08],
            textColorRGB=ROBUST_RGBA[:3],
            textSize=1.2,
            parentObjectUniqueId=self._env.DRONE_ID,
            physicsClientId=self._client,
        )
        p.addUserDebugText(
            "Resilient",
            [0.0, 0.0, 0.08],
            textColorRGB=RESILIENT_RGBA[:3],
            textSize=1.2,
            parentObjectUniqueId=self._other,
            physicsClientId=self._client,
        )
        self._marker = _draw_course(p, self._client, reference)
        self._lock = threading.Lock()
        self._frame = None
        print(
            "PyBullet window is open. The figure-8 faces you. "
            "Drag to orbit, scroll to zoom, Ctrl+drag to pan. "
            "The vehicles move as the estimators run.",
            flush=True,
        )

    def push(self, robust, resilient, reference):
        """Store the newest states. Safe to call from the solver thread."""
        frame = (
            np.asarray(robust, dtype=float).copy(),
            np.asarray(resilient, dtype=float).copy(),
            np.asarray(reference, dtype=float).copy(),
        )
        with self._lock:
            self._frame = frame

    def pump(self):
        """Draw the newest states and let the GUI handle the mouse.

        Returns false once the window has been closed.
        """
        client = self._client
        if client < 0 or not self._p.isConnected(physicsClientId=client):
            return False
        try:
            self._p.stepSimulation(physicsClientId=client)
            with self._lock:
                frame = self._frame
            if frame is None:
                return True
            robust, resilient, reference = frame
            _pose(self._p, self._env.DRONE_ID, robust, y=-0.18, client=client)
            _pose(self._p, self._other, resilient, y=0.18, client=client)
            self._p.resetBasePositionAndOrientation(
                self._marker,
                [float(reference[0]), 0.0, float(HOVER_Z + reference[2])],
                [0.0, 0.0, 0.0, 1.0],
                physicsClientId=client,
            )
        except self._p.error:
            return False
        return True

    def close(self):
        env = self._env
        client = self._client
        self._env = None
        self._client = -1
        if env is None or client < 0:
            return
        try:
            if self._p.isConnected(physicsClientId=client):
                env.close()
        except self._p.error:
            return
