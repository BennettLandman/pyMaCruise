"""Fixed-degree approximate Chebyshev TOADS field fitting.

This ports the selected `InhomogeneityCorrection2.computeApprox3DPolynomialField`
path at degree one and the wrapper's effective OPTIMIZED_IMAGE/GLOBAL/
CHEBYSHEV/FIXED/EQUAL_VARIANCE/q=2 settings. The next higher degree and the
TOADS iteration/segmentation integration remain separate work.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_inhomogeneity import toads_chebyshev_basis
from .toads_math import toads_powertable_lookup


@dataclass(frozen=True)
class ToadsApproximateChebyshevFieldFit:
    """Fit output and double-precision normal system for auditability."""

    field: NDArray[np.float32]
    coefficients: NDArray[np.float64]
    normal_matrix: NDArray[np.float64]
    right_hand_side: NDArray[np.float64]
    correction_mask: NDArray[np.bool_]


def fit_toads_approximate_chebyshev_field_degree1(
    image_values: Any,
    memberships: Any,
    best_class_indices: Any,
    segmentation: Any,
    centroids: Any,
    *,
    background_class_index: int = 0,
) -> ToadsApproximateChebyshevFieldFit:
    """Fit and mean-normalize one selected degree-one multiplicative field.

    Inputs use RAS-oriented `(x,y,z)` and `(x,y,z,slot)` order. This matches
    the wrapper's approximate correction constructor: membership exponent 2,
    effective optimized-image factor, global field, fixed coordinates, equal
    variance, and atlas class 0 as background. Since the default descriptor
    has `optimized_*` objects, the constructor changes IMAGE to
    OPTIMIZED_IMAGE. Source fitting samples at stride 3 from origin and
    divides the normal matrix and RHS by the full voxel count. The fitted
    field is evaluated over every correction-mask voxel; its mean over those
    voxels is shifted to one and values outside the mask remain one.

    A missing/singular fit is rejected instead of returning non-finite output.
    This source-safe finite subset does not define the historical solver's
    singular-matrix behavior.
    """
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
    if isinstance(background_class_index, bool) or not isinstance(background_class_index, (int, np.integer)):
        raise ValueError("background_class_index must be an integer class id")
    if not 0 <= int(background_class_index) < centers.size:
        raise ValueError("background_class_index is outside the centroid vector")
    if np.any(best_raw < 0) or np.any(best_raw >= centers.size):
        raise ValueError("best_class_indices contains a class outside the centroid vector")

    best = best_raw.astype(np.intp, copy=False)
    correction_mask = seg_raw != int(background_class_index)
    if not np.any(correction_mask):
        raise ValueError("the selected source-mask fit requires nonempty correction support")

    n_parameters = 4
    rhs = np.zeros(n_parameters, dtype=np.float64)
    normal = np.zeros((n_parameters, n_parameters), dtype=np.float64)
    nx, ny, nz = image.shape
    for x in range(0, nx, 3):
        for y in range(0, ny, 3):
            for z in range(0, nz, 3):
                if not correction_mask[x, y, z]:
                    continue
                basis = toads_chebyshev_basis((x, y, z), image.shape, degree=1)
                for slot in range(mem.shape[-1]):
                    klass = int(best[x, y, z, slot])
                    if klass == int(background_class_index):
                        continue
                    membership = mem[x, y, z, slot]
                    membership_squared = np.float32(membership * membership)
                    voxel_value = image[x, y, z]
                    residual = np.float32(np.float32(np.float32(1.0) * voxel_value) - centers[klass])
                    distance = np.float32(residual * residual)
                    weight = np.float32(1.0)
                    if distance > np.float32(1e-30):
                        powered = toads_powertable_lookup(
                            np.asarray(distance, dtype=np.float32), 1.0,
                            cached_exponent=1.0, minimum=0.0, maximum=2.0,
                        )
                        weight = np.float32(float(powered) / float(distance))
                    vector_term = np.float32(
                        np.float32(np.float32(np.float32(membership_squared * weight) * voxel_value) * centers[klass])
                    )
                    matrix_term = np.float32(
                        np.float32(np.float32(np.float32(membership_squared * weight) * voxel_value) * voxel_value)
                    )
                    for row in range(n_parameters):
                        rhs[row] += np.float32(vector_term * basis[row])
                        for column in range(n_parameters):
                            term = np.float32(np.float32(matrix_term * basis[row]) * basis[column])
                            normal[row, column] += term

    full_volume = np.float32(nx * ny * nz)
    rhs /= full_volume
    normal /= full_volume
    if not np.isfinite(normal).all() or not np.isfinite(rhs).all():
        raise ValueError("source normal system contains non-finite values")
    try:
        coefficients = np.linalg.solve(normal, rhs)
    except np.linalg.LinAlgError as exc:
        raise ValueError("source normal system is singular for this synthetic fit") from exc
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
                basis = toads_chebyshev_basis((x, y, z), image.shape, degree=1)
                value = np.float32(0.0)
                for index in range(n_parameters):
                    value = np.float32(float(coefficients[index]) * float(basis[index]) + float(value))
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
        field, coefficients, normal, rhs, correction_mask.copy(),
    )
