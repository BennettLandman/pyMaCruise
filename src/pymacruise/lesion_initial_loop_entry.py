"""Carry initial LesionTOADS memberships into the active loop entry.

This captures the selected wrapper seam after its initial membership update:
set GMflag, seed the previous energy and ratio, initialize iteration/stop
state, and evaluate the main-loop gate. The next selected call is
``propagateThinning``; this adapter does not run it or model the later
13-class membership update, whose post-initial GMflag behavior needs its own
source-matched kernel.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .lesion_standard13_volume import Standard13VolumeResult


@dataclass(frozen=True)
class LesionToadsInitialLoopEntry:
    """Selected loop-control values handed forward from initial memberships."""

    initial_membership: Standard13VolumeResult
    previous_energy: np.float32
    energy_ratio: np.float32
    gm_flag: bool
    iteration: int
    stop: bool
    enter_loop: bool
    next_stage: str | None
    stage_order: tuple[str, ...]


def seed_lesion_toads_initial_loop(
    initial_membership: Standard13VolumeResult,
    *,
    iterations: int,
    thread_stopped: bool = False,
    is_working: bool = True,
) -> LesionToadsInitialLoopEntry:
    """Apply the wrapper's post-initial-membership loop-entry transition.

    The initial standard-13 volume result is carried by reference. The source
    sets GMflag only after this first ``computeMemberships()`` call, stores
    that call's energy in ``Eprev``, seeds ``Eratio`` to exactly ``1.0f``,
    and enters the main loop unless the iteration count is zero, the thread
    is stopped, or the segmentation reports that it is no longer working.
    """
    if not isinstance(initial_membership, Standard13VolumeResult):
        raise TypeError("initial_membership must be a Standard13VolumeResult")
    if isinstance(iterations, (bool, np.bool_)) or not isinstance(
        iterations, (int, np.integer)
    ) or iterations < 0:
        raise ValueError("iterations must be a nonnegative integer")
    if not isinstance(thread_stopped, (bool, np.bool_)) or not isinstance(
        is_working, (bool, np.bool_)
    ):
        raise ValueError("thread_stopped and is_working must be booleans")

    energy = np.float32(initial_membership.objective)
    if not np.isfinite(energy):
        raise ValueError("initial membership objective must be finite")
    stop = int(iterations) == 0
    enter_loop = bool((not stop) and (not thread_stopped) and is_working)
    stages = [
        "computeMemberships_initial",
        "setGMflag",
        "Eprev=initial_energy",
        "Eratio=1.0f",
        "initialize_stop_and_iteration",
        "main_loop_gate",
    ]
    if enter_loop:
        stages.extend(("iteration=1", "propagateThinning_pending"))
    return LesionToadsInitialLoopEntry(
        initial_membership=initial_membership,
        previous_energy=energy,
        energy_ratio=np.float32(1.0),
        gm_flag=True,
        iteration=1 if enter_loop else 0,
        stop=stop,
        enter_loop=enter_loop,
        next_stage="propagateThinning" if enter_loop else None,
        stage_order=tuple(stages),
    )
