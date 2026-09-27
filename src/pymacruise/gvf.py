"""Python port of the selected VandyCRUISE generalized GVF calculation."""

from __future__ import annotations

from dataclasses import dataclass
import warnings

import nibabel as nib
import numpy as np
from scipy.sparse.linalg import LinearOperator, cg

from pymacruise.image_io import normalize_image_to_ras


@dataclass(frozen=True)
class GVFImageResult:
    """RAS-oriented component-last NIfTI vector field."""

    vector_field: nib.Nifti1Image


def generalized_gvf(
    scalar_volume: np.ndarray,
    *,
    lambda_value: float = 1.0,
    max_iterations: int = 200,
    solver_tolerance: float = 1e-6,
) -> np.ndarray:
    """Compute generalized GVF using the source discretization and SciPy CG.

    The source solves each of three component equations with a multigrid
    V-cycle. This implementation solves the same matrix-free discrete system
    with conjugate gradients; returned components follow array axes 0, 1, 2.
    """
    values = np.asarray(scalar_volume, dtype=np.float32)
    if values.ndim != 3 or any(length < 2 for length in values.shape):
        raise ValueError("GVF input must be a 3D array with dimensions >= 2")
    if not np.isfinite(values).all():
        raise ValueError("GVF input must be finite")
    if not np.isfinite(lambda_value) or lambda_value <= 0:
        raise ValueError("lambda_value must be finite and positive")
    if max_iterations < 0:
        raise ValueError("max_iterations must be non-negative")
    if not np.isfinite(solver_tolerance) or solver_tolerance <= 0:
        raise ValueError("solver_tolerance must be finite and positive")

    minimum = float(values.min())
    maximum = float(values.max())
    if maximum == minimum:
        normalized = np.zeros(values.shape, dtype=np.float32)
    else:
        normalized = ((values - minimum) / (maximum - minimum)).astype(np.float32)
    gradient = _source_gradient(normalized)
    magnitude_squared = np.sum(gradient * gradient, axis=-1, dtype=np.float32)
    edge_weight = np.exp(-magnitude_squared / np.float32(lambda_value * lambda_value)).astype(np.float32)
    g = (edge_weight / 6.0).astype(np.float32)
    h = (1.0 - edge_weight).astype(np.float32)
    shape = values.shape
    size = values.size
    output = np.empty((*shape, 3), dtype=np.float32)

    for component in range(3):
        source = gradient[..., component].astype(np.float64)
        h64 = h.astype(np.float64)
        g64 = g.astype(np.float64)
        right_hand_side = h64 * source

        def apply(operator_values: np.ndarray) -> np.ndarray:
            field = np.asarray(operator_values, dtype=np.float64).reshape(shape)
            laplacian = _source_laplacian(field)
            return (h64 * field - g64 * laplacian).ravel()

        operator = LinearOperator((size, size), matvec=apply, dtype=np.float64)
        initial = source.ravel().copy()
        if max_iterations == 0:
            solution = initial
        else:
            try:
                solution, info = cg(
                    operator,
                    right_hand_side.ravel(),
                    x0=initial,
                    rtol=solver_tolerance,
                    atol=0.0,
                    maxiter=max_iterations,
                )
            except TypeError:  # SciPy 1.10/1.11 use `tol` instead of `rtol`.
                solution, info = cg(
                    operator,
                    right_hand_side.ravel(),
                    x0=initial,
                    tol=solver_tolerance,
                    maxiter=max_iterations,
                )
            if info < 0:
                raise RuntimeError(f"GVF conjugate-gradient solver failed with status {info}")
            if info > 0:
                warnings.warn(
                    f"GVF component {component} did not converge in {max_iterations} iterations",
                    RuntimeWarning,
                    stacklevel=2,
                )
        output[..., component] = solution.reshape(shape).astype(np.float32)
    return output


def generalized_gvf_image(
    scalar_volume: nib.spatialimages.SpatialImage,
    *,
    lambda_value: float = 1.0,
    max_iterations: int = 200,
    solver_tolerance: float = 1e-6,
) -> GVFImageResult:
    """Run generalized GVF on a scalar NIfTI image and return RAS vectors."""
    if not isinstance(scalar_volume, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError("GVF input must be a NIfTI image")
    image = normalize_image_to_ras(scalar_volume, source_name="GVF scalar input")
    if len(image.shape) != 3:
        raise ValueError("GVF scalar input must be 3D")
    vectors = generalized_gvf(
        np.asarray(image.dataobj),
        lambda_value=lambda_value,
        max_iterations=max_iterations,
        solver_tolerance=solver_tolerance,
    )
    header = image.header.copy()
    header.set_data_dtype(np.float32)
    header.set_slope_inter(None, None)
    header.set_intent("vector")
    return GVFImageResult(nib.Nifti1Image(vectors, image.affine, header=header))


def _source_gradient(values: np.ndarray) -> np.ndarray:
    gradient = []
    for axis, length in enumerate(values.shape):
        derivative = np.empty(values.shape, dtype=np.float32)
        low = [slice(None)] * 3
        high = [slice(None)] * 3
        if length == 2:
            low[axis], high[axis] = 0, 1
            difference = values[tuple(high)] - values[tuple(low)]
            derivative[tuple(low)] = difference
            derivative[tuple(high)] = difference
        else:
            low[axis], high[axis] = 0, 1
            derivative[tuple(low)] = values[tuple(high)] - values[tuple(low)]
            low[axis], high[axis] = -1, -2
            derivative[tuple(low)] = values[tuple(low)] - values[tuple(high)]
            center = [slice(None)] * 3
            before = [slice(None)] * 3
            after = [slice(None)] * 3
            center[axis] = slice(1, -1)
            before[axis] = slice(0, -2)
            after[axis] = slice(2, None)
            derivative[tuple(center)] = (values[tuple(after)] - values[tuple(before)]) * 0.5
        gradient.append(derivative)
    return np.stack(gradient, axis=-1)


def _source_laplacian(values: np.ndarray) -> np.ndarray:
    laplacian = np.zeros(values.shape, dtype=np.float64)
    for axis, length in enumerate(values.shape):
        indices_before = np.maximum(np.arange(length) - 1, 0)
        indices_after = np.minimum(np.arange(length) + 1, length - 1)
        laplacian += np.take(values, indices_before, axis=axis)
        laplacian += np.take(values, indices_after, axis=axis)
        laplacian -= 2.0 * values
    return laplacian
