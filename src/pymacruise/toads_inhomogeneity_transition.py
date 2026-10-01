"""Selected source-ordered, field-aware TOADS inhomogeneity update transition.

This isolates the previous-energy-ratio degree ramp and the per-channel
Chebyshev refits in ``AlgorithmMultipleToads``. It deliberately leaves
thinning, statistics, and membership evolution to the caller.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .toads_inhomogeneity import toads_chebyshev_basis
from .toads_inhomogeneity_fit import ToadsApproximateChebyshevFieldFit
from .toads_math import toads_powertable_lookup


@dataclass(frozen=True)
class ToadsInhomogeneityTransition:
    """Result of one degree-gated field update over all input channels."""

    degree: int
    fields: tuple[np.ndarray, ...]
    fits: tuple[ToadsApproximateChebyshevFieldFit, ...]


def toads_next_inhomogeneity_degree(
    degree: int,
    polynomial_degree: int,
    previous_energy_ratio: float,
    low_energy_ratio: float,
) -> int:
    """Apply the selected Java loop's one-step, strict-threshold degree ramp."""
    for name, value in (("degree", degree), ("polynomial_degree", polynomial_degree)):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise ValueError(f"{name} must be an integer")
    if degree < 0 or polynomial_degree < 0:
        raise ValueError("degrees must be nonnegative")
    ratio = np.float32(previous_energy_ratio)
    threshold = np.float32(low_energy_ratio)
    if not np.isfinite(ratio) or not np.isfinite(threshold):
        raise ValueError("energy ratios must be finite")
    if degree < polynomial_degree and ratio < threshold:
        return int(degree + 1)
    return int(degree)


def _fit_channel(
    image: np.ndarray,
    memberships: np.ndarray,
    best: np.ndarray,
    segmentation: np.ndarray,
    centroids: np.ndarray,
    field: np.ndarray,
    *,
    degree: int,
    background_class_index: int,
) -> ToadsApproximateChebyshevFieldFit:
    if degree not in (2, 3):
        raise ValueError("field-aware transition fits support degree 2 or 3")
    n_parameters = 10 if degree == 2 else 20
    mask = segmentation != background_class_index
    if not np.any(mask):
        raise ValueError("the selected source-mask fit requires nonempty correction support")

    nx, ny, nz = image.shape
    rhs = np.zeros(n_parameters, dtype=np.float64)
    normal = np.zeros((n_parameters, n_parameters), dtype=np.float64)
    for x in range(0, nx, 3):
        for y in range(0, ny, 3):
            for z in range(0, nz, 3):
                if not mask[x, y, z]:
                    continue
                basis = toads_chebyshev_basis((x, y, z), image.shape, degree=degree)
                voxel = image[x, y, z]
                corrected = np.float32(field[x, y, z] * voxel)
                for slot in range(memberships.shape[-1]):
                    klass = int(best[x, y, z, slot])
                    if klass == background_class_index:
                        continue
                    membership = memberships[x, y, z, slot]
                    membership_squared = np.float32(membership * membership)
                    residual = np.float32(corrected - centroids[klass])
                    distance = np.float32(residual * residual)
                    weight = np.float32(1.0)
                    if distance > np.float32(1e-30):
                        powered = toads_powertable_lookup(
                            np.asarray(distance, dtype=np.float32),
                            1.0,
                            cached_exponent=1.0,
                            minimum=0.0,
                            maximum=2.0,
                        )
                        weight = np.float32(float(powered) / float(distance))
                    vector_term = np.float32(
                        np.float32(np.float32(membership_squared * weight) * voxel)
                        * centroids[klass]
                    )
                    matrix_term = np.float32(
                        np.float32(np.float32(membership_squared * weight) * voxel)
                        * voxel
                    )
                    for row in range(n_parameters):
                        rhs[row] += np.float32(vector_term * basis[row])
                        row_term = np.float32(matrix_term * basis[row])
                        for column in range(n_parameters):
                            normal[row, column] += float(
                                np.float32(row_term * basis[column])
                            )

    volume = np.float32(nx * ny * nz)
    rhs /= float(volume)
    normal /= float(volume)
    if not np.isfinite(rhs).all() or not np.isfinite(normal).all():
        raise ValueError("source normal system contains non-finite values")
    try:
        coefficients = np.linalg.solve(normal, rhs)
    except np.linalg.LinAlgError as exc:
        raise ValueError("source normal system is singular for this fit") from exc
    if not np.isfinite(coefficients).all():
        raise ValueError("source normal system produced non-finite coefficients")

    mean = np.float32(0.0)
    count = np.float32(0.0)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                if not mask[x, y, z]:
                    continue
                basis = toads_chebyshev_basis((x, y, z), image.shape, degree=degree)
                value = np.float32(0.0)
                for index in range(n_parameters):
                    value = np.float32(
                        float(coefficients[index]) * float(basis[index]) + float(value)
                    )
                field[x, y, z] = value
                mean = np.float32(mean + value)
                count = np.float32(count + np.float32(1.0))
    mean = np.float32(mean / count)
    offset = np.float32(np.float32(1.0) - mean)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                if mask[x, y, z]:
                    field[x, y, z] = np.float32(field[x, y, z] + offset)

    return ToadsApproximateChebyshevFieldFit(
        field=field,
        coefficients=coefficients,
        normal_matrix=normal,
        right_hand_side=rhs,
        correction_mask=mask.copy(),
    )


def update_toads_approximate_field_transition(
    images: Any,
    memberships: Any,
    best_class_indices: Any,
    segmentation: Any,
    centroids_by_channel: Any,
    attached_fields: Any,
    *,
    degree: int,
    polynomial_degree: int,
    previous_energy_ratio: float,
    low_energy_ratio: float,
    background_class_index: int = 0,
) -> ToadsInhomogeneityTransition:
    """Gate the degree and update attached per-channel fields in place.

    ``images`` uses ``(x, y, z, channel)`` order and ``attached_fields`` is a
    sequence of per-channel ``(x, y, z)`` float32 arrays, matching the Java
    array references held by the segmentation. Centroids use
    ``(channel, class)``. Each supplied field object is retained and mutated.
    """
    image = np.asarray(images, dtype=np.float32)
    memberships_array = np.asarray(memberships, dtype=np.float32)
    best_raw = np.asarray(best_class_indices)
    seg_raw = np.asarray(segmentation)
    centers = np.asarray(centroids_by_channel, dtype=np.float32)
    try:
        fields = tuple(np.asarray(field) for field in attached_fields)
    except TypeError as exc:
        raise ValueError("attached_fields must be a sequence of channel fields") from exc
    if image.ndim != 4 or any(size < 2 for size in image.shape[:3]):
        raise ValueError("images must be an (x,y,z,channel) grid")
    if len(fields) != image.shape[3] or any(
        field.shape != image.shape[:3] or field.dtype != np.float32 for field in fields
    ):
        raise ValueError("attached_fields must contain one float32 grid per image channel")
    if memberships_array.ndim != 4 or memberships_array.shape[:3] != image.shape[:3]:
        raise ValueError("memberships must match the image spatial grid")
    if best_raw.shape != memberships_array.shape or not np.issubdtype(best_raw.dtype, np.integer):
        raise ValueError("best_class_indices must match memberships")
    if seg_raw.shape != image.shape[:3] or not np.issubdtype(seg_raw.dtype, np.integer):
        raise ValueError("segmentation must match the image spatial grid")
    if centers.ndim != 2 or centers.shape[0] != image.shape[3]:
        raise ValueError("centroids_by_channel must have one row per image channel")
    if not np.isfinite(image).all() or any(not np.isfinite(field).all() for field in fields):
        raise ValueError("images and attached_fields must be finite")
    if not np.isfinite(memberships_array).all() or not np.isfinite(centers).all():
        raise ValueError("memberships and centroids must be finite")
    if np.any(best_raw < 0) or np.any(best_raw >= centers.shape[1]):
        raise ValueError("best_class_indices contains a class outside the centroid vector")

    next_degree = toads_next_inhomogeneity_degree(
        degree, polynomial_degree, previous_energy_ratio, low_energy_ratio
    )
    fits: list[ToadsApproximateChebyshevFieldFit] = []
    if next_degree in (2, 3):
        best = best_raw.astype(np.intp, copy=False)
        for channel in range(image.shape[3]):
            channel_field = fields[channel]
            fits.append(
                _fit_channel(
                    image[..., channel],
                    memberships_array,
                    best,
                    seg_raw,
                    centers[channel],
                    channel_field,
                    degree=next_degree,
                    background_class_index=background_class_index,
                )
            )
    return ToadsInhomogeneityTransition(next_degree, fields, tuple(fits))


def toads_corrected_channel_samples(images: Any, fields: Any, coordinates: Any) -> np.ndarray:
    """Return source-consumer field×image samples without changing raw images."""
    image = np.asarray(images, dtype=np.float32)
    field = np.asarray(fields, dtype=np.float32)
    xyz = np.asarray(coordinates)
    if image.ndim != 4 or field.shape != image.shape:
        raise ValueError("images and fields must share (x,y,z,channel) shape")
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not np.issubdtype(xyz.dtype, np.integer):
        raise ValueError("coordinates must be an integer (n,3) array")
    if np.any(xyz < 0) or any(np.any(xyz[:, axis] >= image.shape[axis]) for axis in range(3)):
        raise ValueError("coordinates are outside the image grid")
    x, y, z = xyz.T
    return np.asarray(field[x, y, z, :] * image[x, y, z, :], dtype=np.float32)
