"""Resample TOADS shape priors over the current RAS subject grid."""

from __future__ import annotations

from typing import Any, Sequence

import numpy as np
from numpy.typing import NDArray

from .toads_shapes import ToadShapeBounds, _rotation_matrix


def _source_shape_matrix(
    image_dimensions_xyz: tuple[int, int, int],
    template_dimensions_xyz: tuple[int, int, int],
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    transform: NDArray[np.float32],
) -> NDArray[np.float32]:
    """Build the selected precomputed rigid image-to-template matrix."""
    image_center = tuple(np.float32(size / 2.0) for size in image_dimensions_xyz)
    template_center = tuple(np.float32(size / 2.0) for size in template_dimensions_xyz)
    image_res = np.asarray(image_resolutions_mm, dtype=np.float32)
    template_res = np.asarray(template_resolutions_mm, dtype=np.float32)
    rotation = _rotation_matrix(transform[:3])
    matrix = np.empty((3, 4), dtype=np.float32)
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
                    np.float32(rotation[axis, source_axis] * np.float32(-image_center[source_axis]))
                    * image_res[source_axis]
                )
            )
            translated = np.float32(translated + term)
        translated = np.float32(translated + transform[axis + 3])
        translated = np.float32(translated / template_res[axis])
        matrix[axis, 3] = np.float32(translated + template_center[axis])
    return matrix


def _linear_closest_grid(
    values: NDArray[np.float32],
    x: NDArray[np.float32],
    y: NDArray[np.float32],
    z: NDArray[np.float32],
) -> NDArray[np.float32]:
    """Vectorize ImageFunctions.linearClosestInterpolation in float32 order."""
    coordinates = (x, y, z)
    indexes = tuple(
        np.clip(np.floor(value).astype(np.int64), 0, values.shape[axis] - 2)
        for axis, value in enumerate(coordinates)
    )
    alpha = tuple(
        np.clip(
            np.subtract(value, index.astype(np.float32), dtype=np.float32),
            np.float32(0.0),
            np.float32(1.0),
        ).astype(np.float32, copy=False)
        for value, index in zip(coordinates, indexes, strict=True)
    )
    complement = tuple(
        np.subtract(np.float32(1.0), value, dtype=np.float32) for value in alpha
    )
    ix, iy, iz = indexes
    ax, ay, az = alpha
    bx, by, bz = complement
    terms = (
        np.multiply(np.multiply(np.multiply(bx, by, dtype=np.float32), bz, dtype=np.float32), values[ix, iy, iz], dtype=np.float32),
        np.multiply(np.multiply(np.multiply(ax, by, dtype=np.float32), bz, dtype=np.float32), values[ix + 1, iy, iz], dtype=np.float32),
        np.multiply(np.multiply(np.multiply(bx, ay, dtype=np.float32), bz, dtype=np.float32), values[ix, iy + 1, iz], dtype=np.float32),
        np.multiply(np.multiply(np.multiply(bx, by, dtype=np.float32), az, dtype=np.float32), values[ix, iy, iz + 1], dtype=np.float32),
        np.multiply(np.multiply(np.multiply(ax, ay, dtype=np.float32), bz, dtype=np.float32), values[ix + 1, iy + 1, iz], dtype=np.float32),
        np.multiply(np.multiply(np.multiply(bx, ay, dtype=np.float32), az, dtype=np.float32), values[ix, iy + 1, iz + 1], dtype=np.float32),
        np.multiply(np.multiply(np.multiply(ax, by, dtype=np.float32), az, dtype=np.float32), values[ix + 1, iy, iz + 1], dtype=np.float32),
        np.multiply(np.multiply(np.multiply(ax, ay, dtype=np.float32), az, dtype=np.float32), values[ix + 1, iy + 1, iz + 1], dtype=np.float32),
    )
    total = np.add(terms[0], terms[1], dtype=np.float32)
    for term in terms[2:]:
        total = np.add(total, term, dtype=np.float32)
    return total


def transform_toads_shape_priors(
    atlas_shapes: Any,
    *,
    image_dimensions_xyz: tuple[int, int, int],
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    transform: Any,
    class_bounds: Sequence[ToadShapeBounds],
    shape_slope: float = 2.0,
) -> NDArray[np.float32]:
    """Sample descriptor-ordered shape maps onto an image grid.

    Class bounds are the selected source's current-transform bounds, whose
    upper coordinates are intentionally exclusive. Outside each bound, class
    zero defaults to one and all other classes default to zero. Inside, the
    source linear-closest interpolation and slope restoration are applied.
    The voxel arithmetic is vectorized but keeps the selected float32
    interpolation-term and accumulation order. All voxel arrays are RAS
    `(x,y,z,class)`; spacing is signed as required by the source RAS mapping.
    """
    shapes = np.asarray(atlas_shapes, dtype=np.float32)
    dimensions = tuple(int(value) for value in image_dimensions_xyz)
    parameters = np.asarray(transform, dtype=np.float32)
    slope = np.float32(shape_slope)
    image_res = np.asarray(image_resolutions_mm, dtype=np.float32)
    template_res = np.asarray(template_resolutions_mm, dtype=np.float32)
    if shapes.ndim != 4 or any(size < 2 for size in shapes.shape[:3]) or shapes.shape[-1] < 1:
        raise ValueError("atlas_shapes must have shape (template_x,template_y,template_z,classes)")
    if not np.isfinite(shapes).all():
        raise ValueError("atlas_shapes must be finite")
    if len(dimensions) != 3 or min(dimensions) < 2:
        raise ValueError("image_dimensions_xyz must contain three dimensions >= 2")
    if len(class_bounds) != shapes.shape[-1]:
        raise ValueError("class_bounds must contain one bound for every descriptor class")
    if parameters.shape != (6,) or not np.isfinite(parameters).all():
        raise ValueError("transform must contain six finite rigid parameters")
    if image_res.shape != (3,) or not np.isfinite(image_res).all() or np.any(image_res == 0):
        raise ValueError("image_resolutions_mm must contain three finite nonzero signed values")
    if template_res.shape != (3,) or not np.isfinite(template_res).all() or np.any(template_res <= 0):
        raise ValueError("template_resolutions_mm must contain three finite positive values")
    if not np.isfinite(slope):
        raise ValueError("shape_slope must be finite")

    output = np.zeros(dimensions + (shapes.shape[-1],), dtype=np.float32)
    output[..., 0] = np.float32(1.0)
    matrix = _source_shape_matrix(
        dimensions,
        tuple(int(value) for value in shapes.shape[:3]),
        tuple(float(value) for value in image_res),
        tuple(float(value) for value in template_res),
        parameters,
    )
    y_grid, z_grid = np.meshgrid(
        np.arange(dimensions[1], dtype=np.float32),
        np.arange(dimensions[2], dtype=np.float32),
        indexing="ij",
    )
    for class_index, bounds in enumerate(class_bounds):
        low = bounds.lower_xyz
        high = bounds.upper_xyz
        if len(low) != 3 or len(high) != 3:
            raise ValueError("each class bound must contain lower and upper xyz triples")
        lower = tuple(max(0, int(value)) for value in low)
        upper = tuple(min(dimensions[axis], int(high[axis])) for axis in range(3))
        if any(lower[axis] >= upper[axis] for axis in range(3)):
            continue
        shape = shapes[..., class_index]
        for x_index in range(lower[0], upper[0]):
            ys = y_grid[lower[1]:upper[1], lower[2]:upper[2]]
            zs = z_grid[lower[1]:upper[1], lower[2]:upper[2]]
            xs = np.full(ys.shape, np.float32(x_index), dtype=np.float32)
            coordinates: list[NDArray[np.float32]] = []
            for axis in range(3):
                coordinate = np.multiply(matrix[axis, 0], xs, dtype=np.float32)
                coordinate = np.add(
                    coordinate,
                    np.multiply(matrix[axis, 1], ys, dtype=np.float32),
                    dtype=np.float32,
                )
                coordinate = np.add(
                    coordinate,
                    np.multiply(matrix[axis, 2], zs, dtype=np.float32),
                    dtype=np.float32,
                )
                coordinate = np.add(coordinate, matrix[axis, 3], dtype=np.float32)
                coordinates.append(coordinate)
            sampled = _linear_closest_grid(shape, *coordinates)
            recovered = np.add(
                np.multiply(slope, np.subtract(sampled, np.float32(0.5), dtype=np.float32), dtype=np.float32),
                np.float32(0.5),
                dtype=np.float32,
            )
            output[
                x_index,
                lower[1]:upper[1],
                lower[2]:upper[2],
                class_index,
            ] = np.clip(recovered, np.float32(0.0), np.float32(1.0))
    return output
