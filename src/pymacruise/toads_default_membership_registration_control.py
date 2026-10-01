"""Continue one exact default-atlas membership pass through wrapper control.

The selected ``AlgorithmMultipleToads`` wrapper calculates its float32 energy
ratio immediately after exact memberships, optionally calls rigid
``atlas.registerShapes()``, and only then checks convergence and the iteration
limit.  This adapter preserves that order for a completed membership result.
It does not compute the membership/statistics state itself or continue into
growth/classification.  Paper context: Bazin and Pham (2007, 2008),
``docs/papers/toads.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .toads_control import (
    ToadsIterationControl,
    toads_energy_ratio,
    toads_iteration_control_from_ratio,
)
from .toads_default_atlas_membership import ToadsDefaultAtlasMembershipPass
from .toads_default_registration_handoff import (
    ToadsDefaultRegistrationHandoff,
    _bounds_for_transform,
)
from .toads_registration_sampling import register_toads_atlas_rigid
from .toads_shapes import _rotation_matrix


@dataclass(frozen=True)
class ToadsDefaultMembershipRegistrationControl:
    """Membership result, optional rigid registration, and later loop control."""

    handoff: ToadsDefaultRegistrationHandoff
    control: ToadsIterationControl
    stage_order: tuple[str, ...]


def continue_default_membership_through_registration_control(
    membership_pass: ToadsDefaultAtlasMembershipPass,
    atlas_shapes: Any,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    initial_transform: tuple[float, float, float, float, float, float],
    previous_energy: float,
    iteration: int,
    maximum_iterations: int,
    convergence_threshold: float,
    do_growing: bool,
    align_shapes: bool = True,
    offset: int = 0,
    subsample: int = 2,
    max_registration_iterations: int = 5,
    registration_energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
) -> ToadsDefaultMembershipRegistrationControl:
    """Apply energy ratio, optional rigid registration, then wrapper stop checks.

    ``membership_pass`` is retained by identity in the returned handoff, so
    the registration stage receives the exact membership array whose
    objective produced ``energyRatio``.  ``toads_energy_ratio`` is shared with
    the corrected P18cyq source-order helper: compute it before registration,
    then pass that same value to the later convergence/iteration checks.
    """
    if not isinstance(membership_pass, ToadsDefaultAtlasMembershipPass):
        raise TypeError("membership_pass must be a selected default-atlas membership pass")
    memberships = np.asarray(membership_pass.memberships, dtype=np.float32)
    if memberships.ndim != 4 or memberships.shape[-1] != 11 or min(memberships.shape[:3]) < 1:
        raise ValueError("membership_pass.memberships must have a nonempty 11-class layout")

    shapes = np.asarray(atlas_shapes, dtype=np.float32)
    if shapes.ndim != 4 or shapes.shape[-1] != 11 or min(shapes.shape[:3]) < 2:
        raise ValueError("atlas_shapes must be a 4D template grid with eleven classes")
    if not np.isfinite(shapes).all():
        raise ValueError("atlas_shapes must be finite")
    transform = np.asarray(initial_transform, dtype=np.float32)
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("initial_transform must contain six finite rigid parameters")

    ratio = toads_energy_ratio(previous_energy, membership_pass.objective)
    stages = ["computeExactMemberships", "energyRatio"]
    current_bounds = _bounds_for_transform(
        shapes,
        memberships.shape[:3],
        image_resolutions_mm,
        template_resolutions_mm,
        transform,
    )
    registration = None
    final_transform = transform.copy()
    final_bounds = current_bounds
    if align_shapes:
        # The wrapper's atlas registration samples the full subject lattice;
        # bounds are refreshed outputs and do not clip this optimizer domain.
        low = np.zeros((11, 3), dtype=np.int32)
        high = np.broadcast_to(
            np.asarray(memberships.shape[:3], dtype=np.int32), (11, 3)
        ).copy()
        registered = np.ones(11, dtype=bool)
        registered[0] = False
        registration = register_toads_atlas_rigid(
            memberships,
            shapes,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            initial_transform=tuple(float(value) for value in transform),
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
            shapes,
            memberships.shape[:3],
            image_resolutions_mm,
            template_resolutions_mm,
            final_transform,
        )
        stages.append("atlas.registerShapes() -> registerAllShapes()")

    handoff = ToadsDefaultRegistrationHandoff(
        membership_pass=membership_pass,
        registration=registration,
        transform=final_transform,
        rotation=_rotation_matrix(final_transform[:3]),
        shape_bounds=final_bounds,
    )
    control = toads_iteration_control_from_ratio(
        ratio,
        iteration,
        maximum_iterations,
        convergence_threshold,
        do_growing=do_growing,
    )
    stages.extend(("convergenceCheck", "iterationLimitCheck"))
    return ToadsDefaultMembershipRegistrationControl(
        handoff=handoff,
        control=control,
        stage_order=tuple(stages),
    )
