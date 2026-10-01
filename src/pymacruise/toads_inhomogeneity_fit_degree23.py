"""Selected fixed degree-two/three TOADS Chebyshev field fitting."""

from __future__ import annotations

from typing import Any

import numpy as np

from .toads_inhomogeneity import toads_chebyshev_basis
from .toads_inhomogeneity_fit import ToadsApproximateChebyshevFieldFit
from .toads_math import toads_powertable_lookup


def fit_toads_approximate_chebyshev_field_degree23(
    image_values: Any,
    memberships: Any,
    best_class_indices: Any,
    segmentation: Any,
    centroids: Any,
    *,
    degree: int,
    background_class_index: int = 0,
) -> ToadsApproximateChebyshevFieldFit:
    """Fit the selected fixed degree-2 or degree-3 global field.

    This follows `InhomogeneityCorrection2.computeApprox3DPolynomialField`
    for the wrapper's effective OPTIMIZED_IMAGE/GLOBAL/CHEBYSHEV/FIXED/
    EQUAL_VARIANCE/q=2 settings. It preserves source parameter ordering,
    stride-three sampling from the origin, p=1 PowerTable weights, background
    candidate exclusion, full-volume matrix/RHS scaling, all-grid evaluation,
    sequential float32 mask-mean accumulation, additive centering, and unit
    values outside the correction mask. As with the fixed degree-one port,
    singular systems raise ValueError instead of returning nonfinite fields.
    """
    if isinstance(degree, bool) or degree not in (2, 3):
        raise ValueError("degree must be 2 or 3")
    image = np.asarray(image_values, dtype=np.float32)
    mem = np.asarray(memberships, dtype=np.float32)
    best_raw = np.asarray(best_class_indices)
    seg_raw = np.asarray(segmentation)
    centers = np.asarray(centroids, dtype=np.float32)
    if image.ndim != 3 or any(size < 2 for size in image.shape):
        raise ValueError("image_values must be a 3D grid with dimensions >= 2")
    if not np.isfinite(image).all():
        raise ValueError("image_values must contain only finite values")
    if mem.ndim != 4 or mem.shape[:3] != image.shape or mem.shape[-1] < 2:
        raise ValueError("memberships must match the image grid and provide at least two slots")
    if not np.isfinite(mem).all() or np.any(mem < 0):
        raise ValueError("memberships must contain finite nonnegative values")
    if best_raw.shape != mem.shape or not np.issubdtype(best_raw.dtype, np.integer):
        raise ValueError("best_class_indices must be an integer grid matching memberships")
    if seg_raw.shape != image.shape or not np.issubdtype(seg_raw.dtype, np.integer):
        raise ValueError("segmentation must be an integer label grid matching the image")
    if centers.ndim != 1 or centers.size < 2 or not np.isfinite(centers).all():
        raise ValueError("centroids must contain one finite value for each class")
    if isinstance(background_class_index, bool) or not isinstance(
        background_class_index, (int, np.integer)
    ):
        raise ValueError("background_class_index must be an integer class id")
    background = int(background_class_index)
    if not 0 <= background < centers.size:
        raise ValueError("background_class_index is outside the centroid vector")
    if np.any(best_raw < 0) or np.any(best_raw >= centers.size):
        raise ValueError("best_class_indices contains a class outside the centroid vector")

    best = best_raw.astype(np.intp, copy=False)
    correction_mask = seg_raw != background
    if not np.any(correction_mask):
        raise ValueError("the selected source-mask fit requires nonempty correction support")

    n_parameters = 10 if degree == 2 else 20
    rhs = np.zeros(n_parameters, dtype=np.float64)
    normal = np.zeros((n_parameters, n_parameters), dtype=np.float64)
    nx, ny, nz = image.shape
    for x in range(0, nx, 3):
        for y in range(0, ny, 3):
            for z in range(0, nz, 3):
                if not correction_mask[x, y, z]:
                    continue
                basis = toads_chebyshev_basis((x, y, z), image.shape, degree=degree)
                for slot in range(mem.shape[-1]):
                    klass = int(best[x, y, z, slot])
                    if klass == background:
                        continue
                    membership = mem[x, y, z, slot]
                    membership_squared = np.float32(membership * membership)
                    voxel_value = image[x, y, z]
                    residual = np.float32(np.float32(voxel_value) - centers[klass])
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
                        np.float32(
                            np.float32(membership_squared * weight) * voxel_value
                        )
                        * centers[klass]
                    )
                    matrix_term = np.float32(
                        np.float32(
                            np.float32(membership_squared * weight) * voxel_value
                        )
                        * voxel_value
                    )
                    for row in range(n_parameters):
                        rhs[row] += np.float32(vector_term * basis[row])
                        row_term = np.float32(matrix_term * basis[row])
                        for column in range(n_parameters):
                            matrix_contribution = np.float32(row_term * basis[column])
                            normal[row, column] += float(matrix_contribution)

    full_volume = np.float32(nx * ny * nz)
    rhs /= float(full_volume)
    normal /= float(full_volume)
    if not np.isfinite(normal).all() or not np.isfinite(rhs).all():
        raise ValueError("source normal system contains non-finite values")
    try:
        coefficients = np.linalg.solve(normal, rhs)
    except np.linalg.LinAlgError as exc:
        raise ValueError("source normal system is singular for this fit") from exc
    if not np.isfinite(coefficients).all():
        raise ValueError("source normal system produced non-finite coefficients")

    field = np.ones(image.shape, dtype=np.float32)
    mean = np.float32(0.0)
    count = np.float32(0.0)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                if not correction_mask[x, y, z]:
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
                if correction_mask[x, y, z]:
                    field[x, y, z] = np.float32(field[x, y, z] + offset)

    return ToadsApproximateChebyshevFieldFit(
        field=field,
        coefficients=coefficients,
        normal_matrix=normal,
        right_hand_side=rhs,
        correction_mask=correction_mask.copy(),
    )
