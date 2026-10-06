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
    _linear_interpolation_with_default,
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
    matrix = np.empty((3, 4), dtype=np.float32)
    centers = np.asarray(image_dims, dtype=np.float32) / np.float32(2.0)
    template_centers = np.asarray(values.shape, dtype=np.float32) / np.float32(2.0)
    for axis in range(3):
        for source_axis in range(3):
            matrix[axis, source_axis] = np.float32(
                np.float32(rotation[axis, source_axis] * image_res[source_axis])
                / template_res[axis]
            )
        translated = np.float32(0.0)
        for source_axis in range(3):
            term = np.float32(
                np.float32(
                    np.float32(rotation[axis, source_axis] * np.float32(-centers[source_axis]))
                    * image_res[source_axis]
                )
            )
            translated = np.float32(translated + term)
        translated = np.float32(translated + transform[axis + 3])
        translated = np.float32(translated / template_res[axis])
        matrix[axis, 3] = np.float32(translated + template_centers[axis])

    lower = list(image_dims)
    upper = [0, 0, 0]
    support_count = 0
    for x in range(image_dims[0]):
        for y in range(image_dims[1]):
            for z in range(image_dims[2]):
                point = (x, y, z)
                mapped = []
                for axis in range(3):
                    coordinate = np.float32(matrix[axis, 0] * np.float32(x))
                    coordinate = np.float32(
                        coordinate + np.float32(matrix[axis, 1] * np.float32(y))
                    )
                    coordinate = np.float32(
                        coordinate + np.float32(matrix[axis, 2] * np.float32(z))
                    )
                    coordinate = np.float32(coordinate + matrix[axis, 3])
                    mapped.append(coordinate)
                if _linear_interpolation_with_default(
                    values, tuple(mapped), np.float32(0.0)
                ) > np.float32(0.0):
                    support_count += 1
                    for axis in range(3):
                        lower[axis] = min(lower[axis], point[axis])
                        upper[axis] = max(upper[axis], point[axis])
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
