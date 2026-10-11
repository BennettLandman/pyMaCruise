"""Connect updated LesionTOADS memberships to its active atlas registration.

Paper context: Shiee et al. (2010), ``docs/papers/lesion-toads.md`` describes
statistical/topological atlas priors. This is the selected active source
handoff with supplied memberships; it does not establish full workflow
registration equivalence.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_shapes import (
    LesionRegistrationOptimizationResult,
    _linear_interpolation_default_grid,
    optimize_lesion_shapes_scaled,
)
from .toads_shapes import ToadShapeBounds, _rotation_matrix


@dataclass(frozen=True)
class LesionRegistrationHandoff:
    """Registration state after the wrapper's membership-to-atlas handoff."""

    memberships: NDArray[np.float32]
    optimization: LesionRegistrationOptimizationResult | None
    transform: NDArray[np.float32]
    rotation: NDArray[np.float32]
    shape_bounds: tuple[ToadShapeBounds, ...]


def _bounds_for_transform(
    shape: ArrayLike,
    image_dimensions_xyz: tuple[int, int, int],
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    transform: NDArray[np.float32],
) -> ToadShapeBounds:
    values = np.asarray(shape, dtype=np.float32)
    image_dims = tuple(int(value) for value in image_dimensions_xyz)
    template_res = np.asarray(template_resolutions_mm, dtype=np.float32)
    image_res = np.asarray(image_resolutions_mm, dtype=np.float32)
    rotation = _rotation_matrix(transform[:3], java_float32_norm=True)
    centers = np.asarray(image_dims, dtype=np.float32) / np.float32(2.0)
    template_centers = np.asarray(values.shape, dtype=np.float32) / np.float32(2.0)
    lower = list(image_dims)
    upper = [0, 0, 0]
    support_count = 0
    y_indices = np.arange(image_dims[1], dtype=np.float32)
    z_indices = np.arange(image_dims[2], dtype=np.float32)
    y_grid, z_grid = np.meshgrid(y_indices, z_indices, indexing="ij")
    coordinates = np.empty((image_dims[1], image_dims[2], 3), dtype=np.float32)
    for x in range(image_dims[0]):
        image_x = np.float32(np.float32(x - centers[0]) * image_res[0])
        image_y = np.float32((y_grid - centers[1]) * image_res[1])
        image_z = np.float32((z_grid - centers[2]) * image_res[2])
        for axis in range(3):
            # Match ParametricTransform.imageToTemplate's float expression:
            # rotate physical offsets, add translation, divide by template
            # resolution, then add the template center. Precomputed affine
            # coefficients alter float32 rounding at positive-support edges.
            coordinate = np.float32(rotation[axis, 0] * image_x)
            coordinate = np.float32(
                coordinate + np.float32(rotation[axis, 1] * image_y)
            )
            coordinate = np.float32(
                coordinate + np.float32(rotation[axis, 2] * image_z)
            )
            coordinate = np.float32(coordinate + transform[axis + 3])
            coordinate = np.float32(coordinate / template_res[axis])
            coordinates[..., axis] = np.float32(coordinate + template_centers[axis])
        sampled = _linear_interpolation_default_grid(
            values, coordinates, np.float32(0.0)
        )
        selected = sampled > np.float32(0.0)
        count = int(np.count_nonzero(selected))
        if count:
            selected_y = y_grid[selected].astype(np.int64)
            selected_z = z_grid[selected].astype(np.int64)
            support_count += count
            lower[0] = min(lower[0], x)
            lower[1] = min(lower[1], int(selected_y.min()))
            lower[2] = min(lower[2], int(selected_z.min()))
            upper[0] = max(upper[0], x)
            upper[1] = max(upper[1], int(selected_y.max()))
            upper[2] = max(upper[2], int(selected_z.max()))
    return ToadShapeBounds(tuple(lower), tuple(upper), support_count, None)


def register_lesion_membership_handoff(
    membership_result: Any,
    shape_priors: Sequence[ArrayLike],
    initial_transform: ArrayLike,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    align_shapes: bool = True,
    init_align_iter: int = 10,
    progress_callback: Callable[[str], None] | None = None,
) -> LesionRegistrationHandoff:
    """Apply the active single-rigid LesionTOADS wrapper alignment branch.

    ``membership_result`` is the membership-pass result returned by the
    P39zzzzs statistics bridge. The selected wrapper uses the updated class
    membership array, scale 1, offset 0, and subsample 2. The optimizer's
    source defaults are made explicit here: lambda 1, multiplier 1.5,
    ``minEdiff=1e-6``, ``minLambda=0.001``, and at most 10 rejected trials per
    line search. Rotation and transformed positive-support bounds are
    refreshed after the optimizer, as in ``LesionToadAtlas.registerAllShapes``.
    """
    membership_source = getattr(membership_result, "membership", membership_result)
    membership_source = getattr(membership_source, "class_memberships", membership_source)
    memberships = np.asarray(membership_source, dtype=np.float32)
    transform = np.asarray(initial_transform, dtype=np.float32)
    if memberships.ndim != 4 or memberships.shape[-1] != len(shape_priors):
        raise ValueError("updated memberships and shape-prior class counts must match")
    if not np.isfinite(memberships).all():
        raise ValueError("updated memberships must be finite")
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("initial_transform must contain six finite rigid parameters")
    if (
        isinstance(init_align_iter, bool)
        or not isinstance(init_align_iter, (int, np.integer))
        or init_align_iter < 0
    ):
        raise ValueError("init_align_iter must be a nonnegative integer")
    if progress_callback is not None and not callable(progress_callback):
        raise TypeError("progress_callback must be callable or None")

    optimization = None
    if align_shapes:
        if progress_callback is not None:
            progress_callback("shape_optimizer:running")
        optimization = optimize_lesion_shapes_scaled(
            shape_priors,
            memberships,
            transform,
            image_dimensions_xyz=tuple(int(size) for size in memberships.shape[:3]),
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            scale=1.0,
            offset=0,
            subsample=2,
            max_iterations=int(init_align_iter),
            min_energy_change=1.0e-6,
            min_lambda=0.001,
            initial_lambda=1.0,
            lambda_factor=1.5,
            max_rejections=10,
        )
        transform = optimization.transform.copy()
        if progress_callback is not None:
            progress_callback("shape_optimizer:complete")

    rotation = _rotation_matrix(
        transform[:3], java_float32_norm=True
    ).astype(np.float32, copy=True)
    bounds_list = []
    for index, shape in enumerate(shape_priors):
        if progress_callback is not None:
            progress_callback(f"shape_bounds:{index}:running")
        bounds_list.append(_bounds_for_transform(
            shape,
            image_dimensions_xyz=tuple(int(size) for size in memberships.shape[:3]),
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            transform=tuple(float(value) for value in transform),
        ))
        if progress_callback is not None:
            progress_callback(f"shape_bounds:{index}:complete")
    bounds = tuple(bounds_list)
    return LesionRegistrationHandoff(
        memberships.copy(), optimization, transform.copy(), rotation, bounds
    )
