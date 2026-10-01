"""Coarse-to-fine rigid registration for default TOADS top-k memberships.

This composes the default 11-class initial membership state with the source
approximate-registration pyramid. It reuses the selected registration
sampling, coefficient, line-search, and block-mean scale kernels validated in
``toads_approximate_registration``; it does not run a segmentation sweep.

Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``, place
image registration and atlas priors within TOADS. This adapter connects an
11-class top-k initialization to that registration role. Pyramid levels,
sampling, line search, and parameter handoff follow selected VandyCRUISE
source behavior; matching this synthetic reference does not establish the
papers' full atlas-registration or segmentation method.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_approximate_registration import (
    _approx_registration_level,
    _expand_slots,
    subsample_toads_approximate_memberships,
)
from .toads_default_approximate_initial_membership import (
    ToadsDefaultApproximateInitialMemberships,
    initialize_toads_default_approximate_initial_memberships,
)
from .toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES
from .toads_registration_sampling import ToadsRigidRegistrationResult
from .toads_shapes import ToadShapeBounds, _rotation_matrix, compute_toad_shape_bounds


@dataclass(frozen=True)
class ToadsDefaultApproximateRegistrationLevel:
    """One default-atlas top-k pyramid scale and its optimizer trace."""

    scale: int
    initial_transform: NDArray[np.float32]
    memberships: NDArray[np.float32]
    best_class_indices: NDArray[np.uint8]
    expanded_memberships: NDArray[np.float32]
    result: ToadsRigidRegistrationResult
    iteration_coefficients: tuple[tuple[NDArray[np.float32], NDArray[np.float32]], ...]
    iteration_trials: tuple[tuple[tuple[NDArray[np.float32], np.float32], ...], ...]


@dataclass(frozen=True)
class ToadsDefaultApproximateAtlasRegistration:
    """Default 11-class approximate registration result and transformed bounds."""

    initialization: ToadsDefaultApproximateInitialMemberships
    levels: tuple[ToadsDefaultApproximateRegistrationLevel, ...]
    transform: NDArray[np.float32]
    rotation: NDArray[np.float32]
    shape_bounds: tuple[ToadShapeBounds, ...]


def register_toads_default_approximate_atlas_pyramid(
    initialization: ToadsDefaultApproximateInitialMemberships,
    atlas_shapes: ArrayLike,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    initial_transform: ArrayLike,
    levels: int,
    max_iterations: int = 50,
    energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
) -> ToadsDefaultApproximateAtlasRegistration:
    """Register default-atlas shape maps using the initial top-k state.

    The source pyramid visits scales ``2**(levels-1)`` down through two,
    floor-truncates incomplete blocks, averages candidate membership per
    coarse block, and resets the line-search lambda on each level. Input and
    output grids use RAS voxel ordering; resolutions are supplied in RAS-axis
    order in millimetres.
    """

    if tuple(initialization.class_types) != DEFAULT_ATLAS_CLASS_TYPES:
        raise ValueError("initialization class order must match the default 11-class atlas")
    memberships = np.asarray(initialization.memberships, dtype=np.float32)
    best = np.asarray(initialization.best_class_indices)
    shapes = np.asarray(atlas_shapes, dtype=np.float32)
    transform = np.asarray(initial_transform, dtype=np.float32)
    classes = len(DEFAULT_ATLAS_CLASS_TYPES)
    if memberships.ndim != 4 or best.shape != memberships.shape or not np.issubdtype(best.dtype, np.integer):
        raise ValueError("initialization must contain matching 4D membership and candidate maps")
    if shapes.ndim != 4 or shapes.shape[-1] != classes or not np.isfinite(shapes).all():
        raise ValueError("atlas_shapes must contain eleven finite template-space maps")
    if not np.isfinite(memberships).all() or np.any((best < 0) | (best >= classes)):
        raise ValueError("top-k memberships or candidate indices are invalid")
    if isinstance(levels, bool) or not isinstance(levels, (int, np.integer)) or levels < 2:
        raise ValueError("levels must be an integer >= 2 for pyramid registration")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, (int, np.integer)) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer")
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("initial_transform must contain six finite rigid parameters")
    for name, resolution in (
        ("image_resolutions_mm", image_resolutions_mm),
        ("template_resolutions_mm", template_resolutions_mm),
    ):
        values = np.asarray(resolution, dtype=np.float32)
        if values.shape != (3,) or not np.isfinite(values).all() or np.any(values <= 0):
            raise ValueError(f"{name} must contain three positive finite RAS voxel sizes")
    if not np.isfinite(energy_tolerance) or energy_tolerance < 0:
        raise ValueError("energy_tolerance must be finite and nonnegative")
    if not np.isfinite(minimum_lambda) or minimum_lambda < 0:
        raise ValueError("minimum_lambda must be finite and nonnegative")

    full_image_dims = tuple(int(size) for size in memberships.shape[:3])
    full_template_dims = tuple(int(size) for size in shapes.shape[:3])
    registered = np.asarray([kind != "MASK" for kind in DEFAULT_ATLAS_CLASS_TYPES], dtype=bool)
    level_results: list[ToadsDefaultApproximateRegistrationLevel] = []
    for scale in (2**power for power in range(int(levels) - 1, 0, -1)):
        pooled = subsample_toads_approximate_memberships(
            memberships, best, classes=classes, factor=scale,
        )
        expanded = _expand_slots(pooled, classes)
        level_initial = transform.copy()
        result, coefficient_history, trial_history = _approx_registration_level(
            pooled, shapes, transform,
            full_image_dims=full_image_dims,
            full_template_dims=full_template_dims,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            scale=scale,
            registered=registered,
            max_iterations=int(max_iterations),
            energy_tolerance=energy_tolerance,
            minimum_lambda=minimum_lambda,
        )
        level_results.append(ToadsDefaultApproximateRegistrationLevel(
            scale, level_initial, pooled.memberships, pooled.best_class_indices,
            expanded, result, coefficient_history, trial_history,
        ))
        transform = result.transform.copy()

    bounds = tuple(
        compute_toad_shape_bounds(
            shapes[..., klass],
            image_dimensions_xyz=full_image_dims,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            transform=tuple(float(value) for value in transform),
        )
        for klass in range(classes)
    )
    return ToadsDefaultApproximateAtlasRegistration(
        initialization,
        tuple(level_results),
        transform.copy(),
        _rotation_matrix(transform[:3]),
        bounds,
    )


def initialize_and_register_toads_default_approximate_atlas(
    image_values: ArrayLike,
    centroids: ArrayLike,
    atlas_shapes: ArrayLike,
    *,
    approx: int,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    initial_transform: ArrayLike,
    levels: int,
    max_iterations: int = 50,
    energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
) -> ToadsDefaultApproximateAtlasRegistration:
    """Compose the selected default-atlas initialization and pyramid call."""

    initialization = initialize_toads_default_approximate_initial_memberships(
        image_values, centroids, approx=approx,
    )
    return register_toads_default_approximate_atlas_pyramid(
        initialization,
        atlas_shapes,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        initial_transform=initial_transform,
        levels=levels,
        max_iterations=max_iterations,
        energy_tolerance=energy_tolerance,
        minimum_lambda=minimum_lambda,
    )
