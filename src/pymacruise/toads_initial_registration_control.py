"""P18cyp continuation: fixed-membership rigid registration and loop control.

This selected wrapper slice consumes the exact membership map and objective
already produced by P18cyo. It does not call the membership kernel again.
The post-thinning correction-field refit in ``AlgorithmMultipleToads`` is
intentionally omitted here; channel correction fields remain fixed at their
P18cyo values and this is not a complete TOADS iteration.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .toads_control import ToadsIterationControl, toads_iteration_control
from .toads_default_registration_handoff import (
    ToadsDefaultRegistrationHandoff,
    _bounds_for_transform,
)
from .toads_initial_thinning_statistics_membership import (
    ToadsInitialStatisticsMembershipContinuation,
)
from .toads_registration_sampling import register_toads_atlas_rigid
from .toads_shapes import _rotation_matrix


@dataclass(frozen=True)
class ToadsInitialRegistrationControl:
    """P18cyo memberships passed through registration and wrapper control."""

    handoff: ToadsDefaultRegistrationHandoff
    control: ToadsIterationControl
    stage_order: tuple[str, ...]


def register_initial_memberships_and_control(
    continuation: ToadsInitialStatisticsMembershipContinuation,
    atlas_shapes: np.ndarray,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    initial_transform: tuple[float, float, float, float, float, float],
    previous_energy: float,
    iteration: int,
    maximum_iterations: int,
    convergence_threshold: float,
    do_growing: bool,
    offset: int = 0,
    subsample: int = 2,
    max_registration_iterations: int = 5,
    registration_energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
    align_shapes: bool = True,
) -> ToadsInitialRegistrationControl:
    """Register the existing 11-class membership state, then decide control.

    Defaults follow the selected MATLAB wrapper settings: rigid alignment,
    offset 0, subsample 2, five registration updates, ``minEdiff=1e-6``, and
    ``minLambda=1e-3``. Atlas class 0 is a mask and is excluded from shape
    registration. Registration uses the P18cyo exact maps directly; this
    composition never recomputes memberships or their objective.
    """
    if not isinstance(continuation, ToadsInitialStatisticsMembershipContinuation):
        raise TypeError("continuation must be a selected P18cyo continuation")
    memberships = continuation.membership_pass.memberships
    objective = np.float32(continuation.membership_pass.objective)
    shapes = np.asarray(atlas_shapes, dtype=np.float32)
    transform = np.asarray(initial_transform, dtype=np.float32)
    if memberships.ndim != 4 or memberships.shape[-1] != 11:
        raise ValueError("P18cyo memberships must have the selected 11-class layout")
    if shapes.ndim != 4 or shapes.shape[-1] != 11 or min(shapes.shape[:3]) < 2:
        raise ValueError("atlas_shapes must be a template grid with eleven classes")
    if not np.isfinite(shapes).all():
        raise ValueError("atlas_shapes must be finite")
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("initial_transform must contain six finite rigid parameters")

    current_bounds = _bounds_for_transform(
        shapes, memberships.shape[:3], image_resolutions_mm,
        template_resolutions_mm, transform,
    )
    registration = None
    final_transform = transform.copy()
    final_bounds = current_bounds
    stages = [*continuation.stage_order]
    if align_shapes:
        low = np.asarray([bound.lower_xyz for bound in current_bounds], dtype=np.int32)
        high = np.asarray([bound.upper_xyz for bound in current_bounds], dtype=np.int32)
        registered = np.ones(11, dtype=bool)
        registered[0] = False
        registration = register_toads_atlas_rigid(
            memberships,
            shapes,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            initial_transform=tuple(float(v) for v in transform),
            class_bounds_low_xyz=low,
            class_bounds_high_xyz=high,
            registered_classes=registered,
            offset=offset,
            subsample=subsample,
            max_iterations=max_registration_iterations,
            energy_tolerance=registration_energy_tolerance,
            minimum_lambda=minimum_lambda,
        )
        final_transform = registration.transform.copy()
        final_bounds = _bounds_for_transform(
            shapes, memberships.shape[:3], image_resolutions_mm,
            template_resolutions_mm, final_transform,
        )
        stages.append("atlas.registerShapes() -> registerAllShapes()")

    handoff = ToadsDefaultRegistrationHandoff(
        membership_pass=continuation.membership_pass,
        registration=registration,
        transform=final_transform,
        rotation=_rotation_matrix(final_transform[:3]),
        shape_bounds=final_bounds,
    )
    stages.append("energyRatio")
    control = toads_iteration_control(
        previous_energy,
        objective,
        iteration,
        maximum_iterations,
        convergence_threshold,
        do_growing=do_growing,
    )
    stages.extend(("convergenceCheck", "iterationLimitCheck"))
    return ToadsInitialRegistrationControl(handoff, control, tuple(stages))
