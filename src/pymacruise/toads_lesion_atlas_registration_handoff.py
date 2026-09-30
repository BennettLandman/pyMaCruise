"""Compose a lesion-bearing 13-class TOADS pass with rigid atlas registration.

This is the ``computeExactMemberships`` -> ``registerShapes`` handoff for the
selected base-TOADS wrapper when it is configured with the lesion-bearing
atlas descriptor.  It shares the validated exact-membership registration
kernel with the default-atlas handoff, but retains the descriptor class order
and its lesion map as part of the membership result.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_lesion_atlas_type_membership import (
    LESION_TOADS_CLASS_TYPES,
    ToadsLesionAtlasTypeMembershipPass,
)
from .toads_registration_sampling import (
    ToadsRigidRegistrationResult,
    register_toads_atlas_rigid,
)
from .toads_shapes import ToadShapeBounds, _rotation_matrix, compute_toad_shape_bounds


@dataclass(frozen=True)
class ToadsLesionAtlasRegistrationHandoff:
    """Membership result and registration state after the wrapper handoff."""

    membership_pass: ToadsLesionAtlasTypeMembershipPass
    registration: ToadsRigidRegistrationResult | None
    transform: NDArray[np.float32]
    rotation: NDArray[np.float32]
    shape_bounds: tuple[ToadShapeBounds, ...]
    registered_class_indices: tuple[int, ...]


def _shape_bounds(
    atlas_shapes: NDArray[np.float32],
    image_dimensions_xyz: tuple[int, int, int],
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    transform: NDArray[np.float32],
) -> tuple[ToadShapeBounds, ...]:
    return tuple(
        compute_toad_shape_bounds(
            atlas_shapes[..., k],
            image_dimensions_xyz=image_dimensions_xyz,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            transform=tuple(float(value) for value in transform),
        )
        for k in range(atlas_shapes.shape[-1])
    )


def toads_lesion_atlas_registration_handoff(
    membership_pass: ToadsLesionAtlasTypeMembershipPass,
    atlas_shapes: Any,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    initial_transform: tuple[float, float, float, float, float, float],
    align_shapes: bool = True,
    offset: int = 0,
    subsample: int = 1,
    max_iterations: int = 50,
    energy_tolerance: float = 1e-4,
    minimum_lambda: float = 1e-3,
    source_faithful: bool = True,
) -> ToadsLesionAtlasRegistrationHandoff:
    """Register descriptor-ordered atlas maps using the fresh pass memberships.

    Atlas shapes stay on their template grid.  The membership and lesion maps
    remain on the RAS subject grid; their geometry is not changed.  The
    selected registration includes every declared class except the MASK
    background class, including LESION at class index 10.  This function stops
    after the optional rigid registration and bounds/rotation refresh.
    """
    if not isinstance(membership_pass, ToadsLesionAtlasTypeMembershipPass):
        raise TypeError("membership_pass must be a ToadsLesionAtlasTypeMembershipPass")
    if tuple(membership_pass.class_types) != LESION_TOADS_CLASS_TYPES:
        raise ValueError("membership_pass must use the selected lesion 13-class order")
    memberships = np.asarray(membership_pass.memberships, dtype=np.float32)
    lesion = np.asarray(membership_pass.lesion_membership, dtype=np.float32)
    if memberships.ndim != 4 or memberships.shape[-1] != 13:
        raise ValueError("membership_pass.memberships must have shape (x,y,z,13)")
    if lesion.shape != memberships.shape[:3] or not np.isfinite(lesion).all():
        raise ValueError("membership_pass.lesion_membership must be finite on the subject grid")
    if not np.isfinite(memberships).all():
        raise ValueError("membership_pass.memberships must be finite")

    shapes = np.asarray(atlas_shapes, dtype=np.float32)
    transform = np.asarray(initial_transform, dtype=np.float32)
    if shapes.ndim != 4 or shapes.shape[-1] != 13 or min(shapes.shape[:3]) < 2:
        raise ValueError("atlas_shapes must be a 4D template grid with thirteen class maps")
    if not np.isfinite(shapes).all():
        raise ValueError("atlas_shapes must be finite")
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("initial_transform must contain six finite rigid parameters")
    if len(image_resolutions_mm) != 3 or len(template_resolutions_mm) != 3:
        raise ValueError("image and template resolutions must contain three axes")
    if any(not np.isfinite(value) or value <= 0 for value in (*image_resolutions_mm, *template_resolutions_mm)):
        raise ValueError("image and template resolutions must be positive finite millimetres")
    if isinstance(offset, bool) or not isinstance(offset, (int, np.integer)) or offset < 0:
        raise ValueError("offset must be a nonnegative integer")
    if isinstance(subsample, bool) or not isinstance(subsample, (int, np.integer)) or subsample < 1:
        raise ValueError("subsample must be a positive integer")

    registered_indices = tuple(k for k, kind in enumerate(membership_pass.class_types) if kind != "MASK")
    registered = np.zeros(13, dtype=bool)
    registered[list(registered_indices)] = True
    current_bounds = _shape_bounds(
        shapes, memberships.shape[:3], image_resolutions_mm,
        template_resolutions_mm, transform,
    )
    if not align_shapes:
        return ToadsLesionAtlasRegistrationHandoff(
            membership_pass, None, transform.copy(), _rotation_matrix(transform[:3]),
            current_bounds, registered_indices,
        )

    low = np.asarray([item.lower_xyz for item in current_bounds], dtype=np.int32)
    high = np.asarray([item.upper_xyz for item in current_bounds], dtype=np.int32)
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
        max_iterations=max_iterations,
        energy_tolerance=energy_tolerance,
        minimum_lambda=minimum_lambda,
        source_faithful=source_faithful,
    )
    final_bounds = _shape_bounds(
        shapes, memberships.shape[:3], image_resolutions_mm,
        template_resolutions_mm, registration.transform,
    )
    return ToadsLesionAtlasRegistrationHandoff(
        membership_pass, registration, registration.transform.copy(),
        _rotation_matrix(registration.transform[:3]), final_bounds,
        registered_indices,
    )
