"""Selected exact-membership global TOADS polynomial field refit.

This implements the ``EXACT`` branch of VandyCRUISE
``InhomogeneityCorrection2.computeExact3DPolynomialField`` for the selected
TOADS defaults: stride three, q=2, optimized-image weights, a global
Chebyshev field, and the equal-variance mode. Degrees one through three use
the selected wrapper's coefficient counts and explicit basis order.
"""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Sequence
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_inhomogeneity import toads_chebyshev_basis
from .toads_math import toads_powertable_lookup


@dataclass(frozen=True)
class ToadsExactChebyshevFieldFit:
    """Diagnostics for one exact-membership global field fit."""

    field: NDArray[np.float32]
    coefficients: NDArray[np.float64]
    normal_matrix: NDArray[np.float64]
    right_hand_side: NDArray[np.float64]
    correction_mask: NDArray[np.bool_]


def _jama_lu_solve_square(
    matrix: NDArray[np.float64], right_hand_side: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Solve a square system in the selected JAMA loop and operation order.

    `Jama.Matrix.solve` dispatches square systems to `LUDecomposition.solve`.
    This ports that left-looking Crout/Doolittle factorization, strict-`>`
    partial-pivot selection, and forward/back substitutions. The order is
    material here: NumPy/LAPACK returns coefficients that differ by a few
    double-precision ULPs, which accumulate into later float32 field values.
    Source reference: `tests/java_reference/upstream/Jama/LUDecomposition.java`
    (SHA-256 pinned by the degree-1/2/3 Java fit harnesses).
    """
    a = np.asarray(matrix, dtype=np.float64)
    b = np.asarray(right_hand_side, dtype=np.float64)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("JAMA solve requires a square matrix")
    n = a.shape[0]
    if b.shape != (n,):
        raise ValueError("JAMA right-hand side must match the square matrix")
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("JAMA solve inputs must be finite")

    lu = a.copy()
    piv = np.arange(n, dtype=np.int64)
    column = np.empty(n, dtype=np.float64)
    for j in range(n):
        column[:] = lu[:, j]
        for i in range(n):
            total = np.float64(0.0)
            for k in range(min(i, j)):
                total = np.float64(total + np.float64(lu[i, k] * column[k]))
            lu[i, j] = np.float64(column[i] - total)
            column[i] = lu[i, j]

        pivot = j
        for i in range(j + 1, n):
            if abs(column[i]) > abs(column[pivot]):
                pivot = i
        if pivot != j:
            lu[[pivot, j], :] = lu[[j, pivot], :]
            piv[pivot], piv[j] = piv[j], piv[pivot]
        if lu[j, j] != 0.0:
            for i in range(j + 1, n):
                lu[i, j] = np.float64(lu[i, j] / lu[j, j])

    if np.any(np.diag(lu) == 0.0):
        raise ValueError("source normal system is singular for this fit")

    solution = b[piv].reshape(n, 1).copy()
    for k in range(n):
        for i in range(k + 1, n):
            solution[i, 0] = np.float64(
                solution[i, 0] - np.float64(solution[k, 0] * lu[i, k])
            )
    for k in range(n - 1, -1, -1):
        solution[k, 0] = np.float64(solution[k, 0] / lu[k, k])
        for i in range(k):
            solution[i, 0] = np.float64(
                solution[i, 0] - np.float64(solution[k, 0] * lu[i, k])
            )
    return solution[:, 0]


def fit_toads_exact_chebyshev_field(
    image_values: Any,
    memberships: Any,
    segmentation: Any,
    centroids: Any,
    field: Any,
    *,
    background_class_index: int = 0,
    degree: int = 1,
    class_types: Sequence[str] | None = None,
    backend: str = "numba",
) -> ToadsExactChebyshevFieldFit:
    """Refit an attached global field from the complete class-membership map.

    The selected short constructor sets ``useApprox=false``. Unlike the
    alternate approximate constructor, this path has no candidate-slot
    ``best`` array and iterates every non-background class. The default
    descriptor's optimized CSF/GM/WM power exponents are all one. The default
    class map is the selected 11-class MaCRUISE atlas. Other descriptor class
    layouts, including the 13-class LesionTOADS atlas, must pass ``class_types``
    in exact descriptor order. Source matrix/vector accumulation, the
    1e-30 residual cutoff, stride-three
    sampling, full-grid scaling, float32 field evaluation, and mask-mean
    centering are retained.

    ``degree`` selects the supported source basis size (4, 10, or 20 terms).
    ``field`` is mutated in place on the non-background mask, matching the
    source's attached-array update. Image, memberships, segmentation, and
    centroids are read-only inputs. A singular or nonfinite system raises
    ``ValueError`` instead of returning a nonfinite correction.
    """
    if backend not in ("numba", "python"):
        raise ValueError("backend must be 'numba' or 'python'")
    image = np.asarray(image_values, dtype=np.float32)
    mem = np.asarray(memberships, dtype=np.float32)
    seg = np.asarray(segmentation)
    centers = np.asarray(centroids, dtype=np.float32)
    field_values = np.asarray(field)
    if image.ndim != 3 or any(size < 2 for size in image.shape):
        raise ValueError("image_values must be a 3-D grid with dimensions >= 2")
    if mem.ndim != 4 or mem.shape[:3] != image.shape or mem.shape[-1] < 2:
        raise ValueError("memberships must match the image grid and provide class maps")
    if seg.shape != image.shape or not np.issubdtype(seg.dtype, np.integer):
        raise ValueError("segmentation must be an integer label grid matching image_values")
    if centers.shape != (mem.shape[-1],):
        raise ValueError("centroids must contain one value for each membership class")
    if field_values.shape != image.shape or field_values.dtype != np.float32:
        raise ValueError("field must be a float32 grid matching image_values")
    if not np.isfinite(image).all() or not np.isfinite(mem).all() or not np.isfinite(centers).all():
        raise ValueError("image_values, memberships, and centroids must be finite")
    if not np.isfinite(field_values).all() or np.any(mem < 0):
        raise ValueError("field must be finite and memberships must be nonnegative")
    if isinstance(background_class_index, bool) or not isinstance(background_class_index, (int, np.integer)):
        raise ValueError("background_class_index must be an integer")
    background = int(background_class_index)
    if not 0 <= background < mem.shape[-1]:
        raise ValueError("background_class_index is outside the class map")
    if np.any(seg < 0) or np.any(seg >= mem.shape[-1]):
        raise ValueError("segmentation contains a class outside the membership maps")
    if isinstance(degree, bool) or not isinstance(degree, (int, np.integer)) or int(degree) not in (1, 2, 3):
        raise ValueError("degree must be one of the selected supported values: 1, 2, or 3")
    degree = int(degree)

    mask = seg != background
    if not np.any(mask):
        raise ValueError("source correction mask is empty")

    dimensions = image.shape
    n_parameters = (4, 10, 20)[degree - 1]
    rhs = np.zeros(n_parameters, dtype=np.float64)
    normal = np.zeros((n_parameters, n_parameters), dtype=np.float64)
    power_table = {
        "optimized_csf": (0.0, 2.0, 1.0),
        "optimized_gm": (0.0, 2.0, 1.0),
        "optimized_wm": (0.0, 2.0, 1.0),
    }
    if class_types is None:
        if mem.shape[-1] != 11:
            raise ValueError("class_types is required unless using the selected 11-class atlas")
        # Base descriptor order: MASK, optimized CSF, optimized GM, optimized WM.
        groups = ("mask", "optimized_csf", "optimized_csf", *("optimized_gm",) * 5, *("optimized_wm",) * 3)
    else:
        if isinstance(class_types, (str, bytes)):
            raise ValueError("class_types must be a sequence of descriptor type names")
        groups = tuple(class_types)
        if len(groups) != mem.shape[-1] or any(not isinstance(group, str) or not group for group in groups):
            raise ValueError("class_types must contain one nonempty type name per membership class")
    if background != 0 and groups[background] != "mask":
        raise ValueError("background class type must be 'mask'")

    nx, ny, nz = dimensions
    if backend == "numba":
        from .toads_field_compiled import accumulate_field
        optimized_classes = np.asarray([group in power_table for group in groups], dtype=np.bool_)
        normal, rhs = accumulate_field(
            image, mem, mask, centers, field_values, degree, background, optimized_classes
        )
    else:
        nx, ny, nz = dimensions
        for x in range(0, nx, 3):
            for y in range(0, ny, 3):
                for z in range(0, nz, 3):
                    if not mask[x, y, z]:
                        continue
                    basis = toads_chebyshev_basis((x, y, z), dimensions, degree=degree)
                    for klass in range(mem.shape[-1]):
                        if klass == background:
                            continue
                        corrected = np.float32(field_values[x, y, z] * image[x, y, z])
                        residual = np.float32(corrected - centers[klass])
                        distance = np.float32(residual * residual)
                        optimized_weight = np.float32(1.0)
                        if distance > np.float32(1e-30) and groups[klass] in power_table:
                            low, high, exponent = power_table[groups[klass]]
                            powered = toads_powertable_lookup(
                                np.asarray(distance, dtype=np.float32),
                                exponent,
                                cached_exponent=exponent,
                                minimum=low,
                                maximum=high,
                            )
                            optimized_weight = np.float32(float(powered) / float(distance))

                        membership_squared = np.float32(mem[x, y, z, klass] * mem[x, y, z, klass])
                        vector_term = np.float32(membership_squared * optimized_weight)
                        vector_term = np.float32(vector_term * image[x, y, z])
                        vector_term = np.float32(vector_term * centers[klass])
                        matrix_term = np.float32(membership_squared * optimized_weight)
                        matrix_term = np.float32(matrix_term * image[x, y, z])
                        matrix_term = np.float32(matrix_term * image[x, y, z])
                        for row in range(n_parameters):
                            rhs[row] += float(np.float32(vector_term * basis[row]))
                            row_term = np.float32(matrix_term * basis[row])
                            for column in range(n_parameters):
                                normal[row, column] += float(np.float32(row_term * basis[column]))

    full_volume = np.float32(nx * ny * nz)
    rhs /= float(full_volume)
    normal /= float(full_volume)
    if not np.isfinite(rhs).all() or not np.isfinite(normal).all():
        raise ValueError("source normal system contains nonfinite values")
    try:
        coefficients = _jama_lu_solve_square(normal, rhs)
    except ValueError as exc:
        raise ValueError("source normal system is singular for this fit") from exc
    if not np.isfinite(coefficients).all():
        raise ValueError("source normal system produced nonfinite coefficients")

    if backend == "numba":
        from .toads_field_compiled import evaluate_field
        evaluate_field(field_values, mask, coefficients, degree)
    else:
        mean = np.float32(0.0)
        count = np.float32(0.0)
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    if not mask[x, y, z]:
                        continue
                    basis = toads_chebyshev_basis((x, y, z), dimensions, degree=degree)
                    value = np.float32(0.0)
                    for index in range(n_parameters):
                        value = np.float32(float(coefficients[index]) * float(basis[index]) + float(value))
                    field_values[x, y, z] = value
                    mean = np.float32(mean + value)
                    count = np.float32(count + np.float32(1.0))
        mean = np.float32(mean / count)
        offset = np.float32(np.float32(1.0) - mean)
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    if mask[x, y, z]:
                        field_values[x, y, z] = np.float32(field_values[x, y, z] + offset)

    return ToadsExactChebyshevFieldFit(field_values, coefficients, normal, rhs, mask.copy())


def fit_toads_exact_chebyshev_field_degree1(
    image_values: Any,
    memberships: Any,
    segmentation: Any,
    centroids: Any,
    field: Any,
    *,
    background_class_index: int = 0,
) -> ToadsExactChebyshevFieldFit:
    """Backward-compatible entry point for the exact degree-one field fit."""
    return fit_toads_exact_chebyshev_field(
        image_values,
        memberships,
        segmentation,
        centroids,
        field,
        background_class_index=background_class_index,
        degree=1,
    )
