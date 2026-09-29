"""Source-shaped sampling of TOADS atlas membership priors.

This supplies the statistical-atlas shape-prior term used alongside tissue
likelihoods in the framework of Bazin & Pham (2008). The transform, bounds,
interpolation, and clipping details are ported from the selected VandyCRUISE
source, so this sampling kernel is not evidence of the paper's full method or
of segmentation equivalence. See ``docs/papers/toads.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ToadShapeSample:
    """One transformed atlas membership sample and its template coordinate."""

    value: np.float32
    template_xyz: tuple[np.float32, np.float32, np.float32] | None


@dataclass(frozen=True)
class ToadShapeBounds:
    """Per-class support bounds produced by the selected atlas scan.

    ``upper_xyz`` intentionally stores the largest positive-support voxel
    index, not that index plus one. The selected source later treats it as an
    exclusive bound, so the maximum-support plane is considered outside.
    """

    lower_xyz: tuple[int, int, int]
    upper_xyz: tuple[int, int, int]
    support_count: int
    positive_support_xyz: tuple[tuple[int, int, int], ...] | None


def compute_toad_shape_bounds(
    shape: np.ndarray,
    *,
    image_dimensions_xyz: tuple[int, int, int],
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    transform: tuple[float, float, float, float, float, float],
    image_origin_xyz: tuple[float, float, float] | None = None,
    template_origin_xyz: tuple[float, float, float] | None = None,
    record_support: bool = False,
) -> ToadShapeBounds:
    """Scan positive atlas-shape support using TOADS' precomputed rigid map.

    This reproduces the selected `computeTransformedShapeBoundingBox` scan
    and its cached `ParametricTransform` matrix path for one rigid transform.
    The source initializes minima to the input dimensions and maxima to zero,
    then stores the greatest positive-support voxel index directly; it does
    not add one to make the upper bound half-open.

    For production use, only the bounds and support count are retained. Set
    ``record_support`` when a test or diagnostic needs every support coordinate.
    """
    values = np.asarray(shape, dtype=np.float32)
    if values.ndim != 3 or any(size < 2 for size in values.shape):
        raise ValueError("TOADS shape prior must be 3D with dimensions >= 2")
    if not np.isfinite(values).all():
        raise ValueError("TOADS shape prior must contain only finite values")
    image_dims = _triple_int(image_dimensions_xyz, "image_dimensions_xyz", positive=True)
    image_res = _triple_float(image_resolutions_mm, "image_resolutions_mm", positive=True)
    template_res = _triple_float(
        template_resolutions_mm, "template_resolutions_mm", positive=True
    )
    parameters = np.asarray(transform, dtype=np.float32)
    if parameters.shape != (6,) or not np.isfinite(parameters).all():
        raise ValueError("transform must contain six finite rigid parameters")

    image_origin = _triple_float(
        image_origin_xyz
        if image_origin_xyz is not None
        else tuple(np.float32(size / 2.0) for size in image_dims),
        "image_origin_xyz",
    )
    template_origin = _triple_float(
        template_origin_xyz
        if template_origin_xyz is not None
        else tuple(np.float32(size / 2.0) for size in values.shape),
        "template_origin_xyz",
    )
    rotation = _rotation_matrix(parameters[:3])
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
                    np.float32(rotation[axis, source_axis] * np.float32(-image_origin[source_axis]))
                    * image_res[source_axis]
                )
            )
            translated = np.float32(translated + term)
        translated = np.float32(translated + parameters[axis + 3])
        translated = np.float32(translated / template_res[axis])
        matrix[axis, 3] = np.float32(translated + template_origin[axis])

    lower = list(image_dims)
    upper = [0, 0, 0]
    support: list[tuple[int, int, int]] | None = [] if record_support else None
    support_count = 0
    for x in range(image_dims[0]):
        for y in range(image_dims[1]):
            for z in range(image_dims[2]):
                point = (x, y, z)
                mapped: list[np.float32] = []
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
                if _linear_closest(values, mapped) > np.float32(0.0):
                    support_count += 1
                    if support is not None:
                        support.append(point)
                    for axis in range(3):
                        lower[axis] = min(lower[axis], point[axis])
                        upper[axis] = max(upper[axis], point[axis])
    return ToadShapeBounds(
        tuple(lower),
        tuple(upper),
        support_count,
        None if support is None else tuple(support),
    )


def sample_toad_shape_prior(
    shape: np.ndarray,
    image_point_xyz: tuple[int, int, int],
    *,
    image_dimensions_xyz: tuple[int, int, int],
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    transform: tuple[float, float, float, float, float, float],
    class_index: int,
    bbox_low_xyz: tuple[int, int, int],
    bbox_high_xyz: tuple[int, int, int],
    shape_slope: float = 2.0,
    scale: float = 1.0,
    image_center_xyz: tuple[float, float, float] | None = None,
    template_center_xyz: tuple[float, float, float] | None = None,
) -> ToadShapeSample:
    """Sample a source TOADS shape prior at one input-grid voxel.

    The transform is VandyCRUISE's six-float rigid parameter vector: the first
    three values define its inverse-Rodrigues rotation and the last three are
    translations in millimetres. Input and template arrays use the source
    x/y/z order; callers must pass RAS-oriented image arrays. The atlas map is
    linearly interpolated with the source closest-edge rule, then its shape
    slope is restored and clipped to [0, 1].

    ``scale`` reproduces the source pyramid mapping (2 means half resolution,
    4 quarter resolution). Centers remain in full-resolution voxel units and
    default to half the declared image and template dimensions, as in the
    native-grid path.

    This ports the fixed-transform sampler, not atlas registration or
    registration-optimized shape bounds. A point outside the class-specific
    input bounds receives the source default (class 0 => 1, other classes =>
    0) and has no mapped template coordinate.
    """
    values = np.asarray(shape, dtype=np.float32)
    if values.ndim != 3 or any(size < 2 for size in values.shape):
        raise ValueError("TOADS shape prior must be 3D with dimensions >= 2")
    if not np.isfinite(values).all():
        raise ValueError("TOADS shape prior must contain only finite values")
    image_dims = _triple_int(image_dimensions_xyz, "image_dimensions_xyz", positive=True)
    point = _triple_int(image_point_xyz, "image_point_xyz")
    lower = _triple_int(bbox_low_xyz, "bbox_low_xyz")
    upper = _triple_int(bbox_high_xyz, "bbox_high_xyz")
    if any(point[a] < 0 or point[a] >= image_dims[a] for a in range(3)):
        raise ValueError("image_point_xyz must be inside the declared input grid")
    empty_support_bounds = lower == image_dims and upper == (0, 0, 0)
    if any(lower[a] < 0 or lower[a] > image_dims[a] for a in range(3)) or any(
        upper[a] < 0 or upper[a] > image_dims[a] for a in range(3)
    ) or (not empty_support_bounds and any(lower[a] > upper[a] for a in range(3))):
        raise ValueError("class bounds must be ordered and inside the declared input grid")
    if not isinstance(class_index, (int, np.integer)) or class_index < 0:
        raise ValueError("class_index must be a nonnegative integer")

    image_res = _triple_float(image_resolutions_mm, "image_resolutions_mm", positive=True)
    template_res = _triple_float(
        template_resolutions_mm, "template_resolutions_mm", positive=True
    )
    scale_value = np.float32(scale)
    if not np.isfinite(scale_value) or scale_value <= 0:
        raise ValueError("scale must be positive and finite")
    image_center = _triple_float(
        image_center_xyz
        if image_center_xyz is not None
        else tuple(np.float32(size / 2.0) for size in image_dims),
        "image_center_xyz",
    )
    template_center = _triple_float(
        template_center_xyz
        if template_center_xyz is not None
        else tuple(np.float32(size / 2.0) for size in values.shape),
        "template_center_xyz",
    )
    transform_values = np.asarray(transform, dtype=np.float32)
    if transform_values.shape != (6,) or not np.isfinite(transform_values).all():
        raise ValueError("transform must contain six finite rigid parameters")
    slope = np.float32(shape_slope)
    if not np.isfinite(slope):
        raise ValueError("shape_slope must be finite")

    if any(point[axis] < lower[axis] or point[axis] >= upper[axis] for axis in range(3)):
        default = np.float32(1.0 if class_index == 0 else 0.0)
        return ToadShapeSample(default, None)

    matrix = _rotation_matrix(transform_values[:3])
    physical = tuple(
        np.float32(
            np.float32(np.float32(point[axis]) * scale_value - image_center[axis])
            * image_res[axis]
        )
        for axis in range(3)
    )
    mapped: list[np.float32] = []
    for axis in range(3):
        rotated = np.float32(
            np.float32(
                np.float32(matrix[axis, 0] * physical[0])
                + np.float32(matrix[axis, 1] * physical[1])
            )
            + np.float32(matrix[axis, 2] * physical[2])
        )
        rotated = np.float32(rotated + transform_values[axis + 3])
        coordinate = np.float32(rotated / template_res[axis])
        coordinate = np.float32(coordinate + template_center[axis])
        coordinate = np.float32(coordinate / scale_value)
        mapped.append(coordinate)

    sampled = _linear_closest(values, mapped)
    recovered = np.float32(
        np.float32(slope * np.float32(sampled - np.float32(0.5))) + np.float32(0.5)
    )
    result = np.float32(np.clip(recovered, np.float32(0.0), np.float32(1.0)))
    return ToadShapeSample(result, tuple(mapped))


def _rotation_matrix(parameters: np.ndarray) -> np.ndarray:
    """Recreate ``RotationMatrix.setParameters`` from selected source."""
    a, b, c = (float(value) for value in parameters)
    norm = float(np.sqrt(a * a + b * b + c * c))
    if norm > 1.0:
        magnitude = norm % 2.0
        if magnitude > 1.0:
            a, b, c = (
                magnitude * value / norm - 2.0 * value / norm for value in (a, b, c)
            )
        else:
            a, b, c = (value * magnitude / norm for value in (a, b, c))
    a2, b2, c2 = a * a, b * b, c * c
    x = 1.0 - a2 - b2 - c2
    if x < 0.0 and x > -1e-10:
        x = 0.0
    scalar = float(np.sqrt(x))
    return np.asarray(
        [
            [1.0 - 2.0 * b2 - 2.0 * c2, 2.0 * a * b - 2.0 * scalar * c, 2.0 * a * c + 2.0 * scalar * b],
            [2.0 * a * b + 2.0 * scalar * c, 1.0 - 2.0 * a2 - 2.0 * c2, 2.0 * b * c - 2.0 * scalar * a],
            [2.0 * a * c - 2.0 * scalar * b, 2.0 * b * c + 2.0 * scalar * a, 1.0 - 2.0 * a2 - 2.0 * b2],
        ],
        dtype=np.float32,
    )


def _linear_closest(values: np.ndarray, point: list[np.float32]) -> np.float32:
    indices = [
        min(max(int(np.floor(point[axis])), 0), values.shape[axis] - 2)
        for axis in range(3)
    ]
    alpha = [
        np.float32(np.clip(np.float32(point[axis] - indices[axis]), 0.0, 1.0))
        for axis in range(3)
    ]
    complement = [np.float32(np.float32(1.0) - value) for value in alpha]
    x, y, z = indices
    ax, ay, az = alpha
    bx, by, bz = complement
    # Preserve the expression order in ImageFunctions.linearClosestInterpolation.
    terms = (
        np.float32(np.float32(np.float32(bx * by) * bz) * values[x, y, z]),
        np.float32(np.float32(np.float32(ax * by) * bz) * values[x + 1, y, z]),
        np.float32(np.float32(np.float32(bx * ay) * bz) * values[x, y + 1, z]),
        np.float32(np.float32(np.float32(bx * by) * az) * values[x, y, z + 1]),
        np.float32(np.float32(np.float32(ax * ay) * bz) * values[x + 1, y + 1, z]),
        np.float32(np.float32(np.float32(bx * ay) * az) * values[x, y + 1, z + 1]),
        np.float32(np.float32(np.float32(ax * by) * az) * values[x + 1, y, z + 1]),
        np.float32(np.float32(np.float32(ax * ay) * az) * values[x + 1, y + 1, z + 1]),
    )
    total = np.float32(terms[0] + terms[1])
    for term in terms[2:]:
        total = np.float32(total + term)
    return total


def _triple_int(
    values: tuple[int, int, int], name: str, *, positive: bool = False
) -> tuple[int, int, int]:
    array = np.asarray(values)
    if array.shape != (3,) or not np.issubdtype(array.dtype, np.integer):
        raise ValueError(f"{name} must contain three integers")
    result = tuple(int(value) for value in array)
    if positive and any(value <= 0 for value in result):
        raise ValueError(f"{name} values must be positive")
    return result


def _triple_float(
    values: tuple[float, float, float], name: str, *, positive: bool = False
) -> tuple[np.float32, np.float32, np.float32]:
    array = np.asarray(values, dtype=np.float32)
    if array.shape != (3,) or not np.isfinite(array).all():
        raise ValueError(f"{name} must contain three finite values")
    if positive and np.any(array <= 0):
        raise ValueError(f"{name} values must be positive")
    return tuple(np.float32(value) for value in array)
