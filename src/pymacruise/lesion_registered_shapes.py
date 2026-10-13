"""Sample atlas priors after LesionTOADS' selected shape registration.

This is the precomputed rigid ``getFastTransformedShape`` path used when the
wrapper reaches its next thinning pass. It consumes the bounds and transform
refreshed by ``registerShapes``; it does not model other transform modes or
claim equivalence of the full LesionTOADS workflow.
"""

from __future__ import annotations

from typing import Literal, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_shapes import _linear_interpolation_with_default
from .toads_shapes import _rotation_matrix
from .lesion_registered_shapes_numba import registered_shapes_numba


def transform_registered_lesion_shape_priors(
    registration_result,
    atlas_shape_priors: Sequence[ArrayLike],
    *,
    image_dimensions_xyz: tuple[int, int, int] | None = None,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    backend: Literal["numba", "python"] = "numba",
) -> NDArray[np.float32]:
    """Return selected registered atlas class maps on the image grid.

    The selected source defaults to its precomputed path for linear transforms.
    It builds a rigid image-to-template matrix, uses the post-registration
    support bounds as lower-inclusive / upper-exclusive tests, and samples each
    atlas prior with ``ImageFunctions.linearInterpolation``. Class zero uses a
    default value of one outside its support; all other classes use zero.
    """
    handoff = getattr(registration_result, "registration", registration_result)
    transform = np.asarray(getattr(handoff, "transform", None), dtype=np.float32)
    bounds = getattr(handoff, "shape_bounds", None)
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("registration must provide six finite rigid parameters")
    if bounds is None or len(bounds) != len(atlas_shape_priors):
        raise ValueError("registration bounds and atlas class count must match")

    shapes = [np.asarray(prior, dtype=np.float32) for prior in atlas_shape_priors]
    if not shapes or any(shape.ndim != 3 for shape in shapes):
        raise ValueError("atlas priors must be a nonempty sequence of 3D arrays")
    membership_grid = getattr(handoff, "memberships", None)
    if membership_grid is None:
        if image_dimensions_xyz is None:
            raise ValueError(
                "image_dimensions_xyz is required when registration has no membership grid"
            )
        image_dims = tuple(int(size) for size in image_dimensions_xyz)
    else:
        image_dims = tuple(int(size) for size in membership_grid.shape[:3])
        if image_dimensions_xyz is not None and tuple(image_dimensions_xyz) != image_dims:
            raise ValueError("image_dimensions_xyz must match the registration membership grid")
    if len(image_dims) != 3 or any(size < 2 for size in image_dims):
        raise ValueError("image_dimensions_xyz must contain three dimensions >= 2")
    if any(shape.shape != shapes[0].shape for shape in shapes):
        raise ValueError("atlas priors must share a 3D template grid")
    if any(not np.isfinite(shape).all() for shape in shapes):
        raise ValueError("atlas priors must be finite")
    if backend not in ("numba", "python"):
        raise ValueError("backend must be 'numba' or 'python'")

    image_res = np.asarray(image_resolutions_mm, dtype=np.float32)
    template_res = np.asarray(template_resolutions_mm, dtype=np.float32)
    if image_res.shape != (3,) or template_res.shape != (3,):
        raise ValueError("image and template resolutions must have three values")
    if np.any(image_res == 0) or np.any(template_res <= 0) or not np.isfinite(image_res).all() or not np.isfinite(template_res).all():
        raise ValueError("resolutions must be finite and template values positive")

    rotation = _rotation_matrix(transform[:3], java_float32_norm=True)
    image_center = np.asarray(image_dims, dtype=np.float32) / np.float32(2.0)
    template_center = np.asarray(shapes[0].shape, dtype=np.float32) / np.float32(2.0)
    matrix = np.empty((3, 4), dtype=np.float32)
    # This operation order mirrors ParametricTransform's Java float arithmetic.
    for axis in range(3):
        for source_axis in range(3):
            matrix[axis, source_axis] = np.float32(
                np.float32(rotation[axis, source_axis] * image_res[source_axis])
                / template_res[axis]
            )
        centered = np.float32(0.0)
        for source_axis in range(3):
            centered = np.float32(
                centered
                + np.float32(
                    np.float32(rotation[axis, source_axis] * np.float32(-image_center[source_axis]))
                    * image_res[source_axis]
                )
            )
        centered = np.float32(centered + transform[axis + 3])
        centered = np.float32(centered / template_res[axis])
        matrix[axis, 3] = np.float32(centered + template_center[axis])

    if backend == "numba":
        source_maps = np.ascontiguousarray(np.stack(shapes, axis=0), dtype=np.float32)
        source_bounds = np.asarray(
            [(*box.lower_xyz, *box.upper_xyz) for box in bounds], dtype=np.int64
        )
        return registered_shapes_numba(source_maps, matrix, source_bounds, image_dims)

    # Retained scalar Python oracle for source-order and edge-case comparisons.
    output = np.empty((*image_dims, len(shapes)), dtype=np.float32)
    for x in range(image_dims[0]):
        for y in range(image_dims[1]):
            for z in range(image_dims[2]):
                point = (x, y, z)
                mapped = None
                for k, (shape, box) in enumerate(zip(shapes, bounds, strict=True)):
                    lower, upper = box.lower_xyz, box.upper_xyz
                    if any(point[a] < lower[a] or point[a] >= upper[a] for a in range(3)):
                        output[x, y, z, k] = np.float32(1.0 if k == 0 else 0.0)
                        continue
                    if mapped is None:
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
                    output[x, y, z, k] = _linear_interpolation_with_default(
                        shape, tuple(mapped), np.float32(1.0 if k == 0 else 0.0)
                    )
    return output
