"""Python implementation of the selected VandyCRUISE generalized GVF.

The model is the generalized gradient-vector-flow diffusion of Xu and Prince
(1998), DOI: 10.1016/S0165-1684(98)00140-6. The exact VandyCRUISE
parameterization and discrete operator are documented in ``docs/papers/gvf.md``;
they should not be confused with the paper's continuum equations.
"""

from __future__ import annotations

from dataclasses import dataclass

import nibabel as nib
import numpy as np
from numba import njit

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
    """Compute generalized GVF with the selected VandyCRUISE multigrid solver.

    At the model level this is the steady state of the GGVF diffusion equation
    in Xu and Prince (1998), applied independently to each vector component.
    The VandyCRUISE implementation supplies the actual weighting functions,
    gradient stencil, and boundary behavior used below (see
    ``docs/papers/gvf.md``). Its three component equations are solved with a
    multigrid V-cycle. This port preserves its averaged coefficient and
    residual downsampling, nearest-neighbor interpolation, and red/black
    Gauss-Seidel octant order. Numba compiles the voxel loops while retaining
    the source's single-precision operations. Returned components follow array
    axes 0, 1, 2.
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

    # Xu & Prince discuss [0, 1] edge-map normalization for the experiments
    # in §4; the selected VandyCRUISE implementation makes min/max
    # normalization an input step generally. Its separate `normalize`
    # option controls output vectors, not this step.
    minimum = float(values.min())
    maximum = float(values.max())
    if maximum == minimum:
        normalized = np.zeros(values.shape, dtype=np.float32)
    else:
        normalized = ((values - minimum) / (maximum - minimum)).astype(np.float32)
    # This central/one-sided voxel-grid gradient is the 3-D VandyCRUISE
    # Gradient3d stencil. No physical voxel-spacing factors are applied.
    gradient = _source_gradient(normalized)
    # ImageDataMath.vectorMagnitude4D returns sqrt(sum(component**2)); despite
    # the historical local name SqrMagf, the source uses that magnitude here.
    gradient_magnitude = np.sqrt(
        np.sum(gradient * gradient, axis=-1, dtype=np.float64)
    ).astype(np.float32)
    # Xu & Prince's GGVF equilibrium is §3, Eq. (5); their example spatial
    # weights are Eqs. (8)-(9). Preserve the selected Java implementation's exact
    # convention here: it computes |grad f| (despite its `SqrMagf` variable
    # name), places that magnitude rather than its square in the exponential,
    # and divides the smoothing weight by six for its six-neighbor stencil.
    # See docs/papers/gvf.md for this source-specific parameterization.
    edge_weight = np.exp(
        -gradient_magnitude / np.float32(lambda_value * lambda_value)
    ).astype(np.float32)
    g = (edge_weight / 6.0).astype(np.float32)
    h = (1.0 - edge_weight).astype(np.float32)
    shape = values.shape
    output = np.empty((*shape, 3), dtype=np.float32)

    # The Java wrapper supplies -h*gradient as Poisson3d's RHS; the source
    # smoother moves that term to the positive side of its update equation.
    for component in range(3):
        source = gradient[..., component]
        java_rhs = (-h * source).astype(np.float32)
        output[..., component] = _source_multigrid_solve(
            source, java_rhs, g, h, max_iterations, solver_tolerance
        )
    return output


@njit
def _mg_downsample(values: np.ndarray) -> np.ndarray:
    """Match DownSample3dByHalf's clamped 2x2x2 arithmetic mean."""
    x_size, y_size, z_size = values.shape
    result = np.empty(((x_size + 1) // 2, (y_size + 1) // 2, (z_size + 1) // 2), dtype=np.float32)
    for x in range(result.shape[0]):
        x0 = min(x * 2, x_size - 2)
        for y in range(result.shape[1]):
            y0 = min(y * 2, y_size - 2)
            for z in range(result.shape[2]):
                z0 = min(z * 2, z_size - 2)
                total = values[x0, y0, z0]
                total = np.float32(total + values[x0, y0, z0 + 1])
                total = np.float32(total + values[x0, y0 + 1, z0])
                total = np.float32(total + values[x0, y0 + 1, z0 + 1])
                total = np.float32(total + values[x0 + 1, y0, z0])
                total = np.float32(total + values[x0 + 1, y0, z0 + 1])
                total = np.float32(total + values[x0 + 1, y0 + 1, z0])
                total = np.float32(total + values[x0 + 1, y0 + 1, z0 + 1])
                result[x, y, z] = np.float32(0.125) * total
    return result


@njit
def _mg_upsample_nearest(values: np.ndarray, shape: tuple[int, int, int]) -> np.ndarray:
    """Match Interpolate3dNearest's floor(index/2) transfer."""
    result = np.empty(shape, dtype=np.float32)
    for x in range(shape[0]):
        for y in range(shape[1]):
            for z in range(shape[2]):
                result[x, y, z] = values[x // 2, y // 2, z // 2]
    return result


@njit
def _mg_residue(
    field: np.ndarray,
    right_hand_side: np.ndarray,
    weight: np.ndarray,
    data_weight: np.ndarray,
    spacing_squared: np.float32,
) -> np.ndarray:
    """Match Poisson3d.residue, including its mirrored edge indexing."""
    x_size, y_size, z_size = field.shape
    result = np.empty_like(field)
    scale = np.float32(1.0) / spacing_squared
    for x in range(x_size):
        xm = max(x - 1, 0)
        xp = min(x + 1, x_size - 1)
        for y in range(y_size):
            ym = max(y - 1, 0)
            yp = min(y + 1, y_size - 1)
            for z in range(z_size):
                zm = max(z - 1, 0)
                zp = min(z + 1, z_size - 1)
                neighbors = field[xm, y, z] + field[xp, y, z]
                neighbors = np.float32(neighbors + field[x, ym, z])
                neighbors = np.float32(neighbors + field[x, yp, z])
                neighbors = np.float32(neighbors + field[x, y, zm])
                neighbors = np.float32(neighbors + field[x, y, zp])
                laplacian = np.float32(neighbors - field[x, y, z] * np.float32(6.0))
                residual = data_weight[x, y, z] * field[x, y, z]
                residual = np.float32(residual + right_hand_side[x, y, z])
                residual = np.float32(residual - scale * weight[x, y, z] * laplacian)
                result[x, y, z] = residual
    return result


@njit
def _mg_octant(
    field: np.ndarray,
    right_hand_side: np.ndarray,
    weight: np.ndarray,
    data_weight: np.ndarray,
    spacing_squared: np.float32,
    x_parity: int,
    y_parity: int,
    z_parity: int,
) -> None:
    """In-place red/black sub-lattice update from GVFSolveGaussSeidel."""
    x_size, y_size, z_size = field.shape
    for x in range(x_parity, x_size, 2):
        xm = max(x - 1, 0)
        xp = min(x + 1, x_size - 1)
        for y in range(y_parity, y_size, 2):
            ym = max(y - 1, 0)
            yp = min(y + 1, y_size - 1)
            for z in range(z_parity, z_size, 2):
                zm = max(z - 1, 0)
                zp = min(z + 1, z_size - 1)
                neighbors = field[xm, y, z] + field[xp, y, z]
                neighbors = np.float32(neighbors + field[x, ym, z])
                neighbors = np.float32(neighbors + field[x, yp, z])
                neighbors = np.float32(neighbors + field[x, y, zm])
                neighbors = np.float32(neighbors + field[x, y, zp])
                weighted = np.float32(neighbors * weight[x, y, z])
                weighted_rhs = np.float32(right_hand_side[x, y, z] * spacing_squared)
                numerator = np.float32(weighted - weighted_rhs)
                denominator = np.float32(
                    np.float32(6.0) * weight[x, y, z]
                    + data_weight[x, y, z] * spacing_squared
                )
                field[x, y, z] = numerator / denominator


@njit
def _mg_smooth(
    field: np.ndarray,
    right_hand_side: np.ndarray,
    weight: np.ndarray,
    data_weight: np.ndarray,
    spacing: int,
    iterations: int,
) -> None:
    spacing_squared = np.float32(spacing * spacing)
    # Preserve the exact eight-octant update order selected by the Java code.
    parities = (
        (0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0),
        (1, 0, 0), (0, 0, 1), (0, 1, 0), (1, 1, 1),
    )
    for _ in range(iterations):
        for parity in parities:
            _mg_octant(
                field, right_hand_side, weight, data_weight, spacing_squared,
                parity[0], parity[1], parity[2],
            )


def _source_multigrid_solve(
    initial: np.ndarray,
    right_hand_side: np.ndarray,
    fine_weight: np.ndarray,
    fine_data_weight: np.ndarray,
    max_iterations: int,
    residual_threshold: float,
) -> np.ndarray:
    """Port Poisson3d.solve's fixed-order V-cycle for one GVF component."""
    fine_shape = initial.shape
    max_level = min((length.bit_length() - 1) for length in fine_shape)
    # Python levels are ordered coarse-to-fine as in Java's j=1..maxLevel.
    fine_to_coarse_shapes = [fine_shape]
    for _ in range(1, max_level):
        prior = fine_to_coarse_shapes[-1]
        fine_to_coarse_shapes.append(tuple((length + 1) // 2 for length in prior))
    shapes = list(reversed(fine_to_coarse_shapes))
    weights = [np.empty(shape, dtype=np.float32) for shape in shapes]
    data_weights = [np.empty(shape, dtype=np.float32) for shape in shapes]
    weights[-1] = np.ascontiguousarray(fine_weight)
    data_weights[-1] = np.ascontiguousarray(fine_data_weight)
    for level in range(max_level - 2, -1, -1):
        weights[level] = _mg_downsample(weights[level + 1])
        data_weights[level] = _mg_downsample(data_weights[level + 1])

    field = np.ascontiguousarray(initial).copy()
    rhs = np.ascontiguousarray(right_hand_side)
    threshold = np.float32(residual_threshold)
    for _ in range(max_iterations):
        residual_fine = _mg_residue(field, rhs, weights[-1], data_weights[-1], np.float32(1.0))
        maximum_error = np.max(np.abs(residual_fine))
        if maximum_error < threshold:
            break

        right_sides = [np.empty(shape, dtype=np.float32) for shape in shapes]
        right_sides[-1] = residual_fine
        for level in range(max_level - 2, -1, -1):
            right_sides[level] = _mg_downsample(right_sides[level + 1])

        corrections = [np.empty(shape, dtype=np.float32) for shape in shapes]
        corrections[0] = np.zeros(shapes[0], dtype=np.float32)
        coarse_spacing = 1 << (max_level - 1)
        _mg_smooth(
            corrections[0], right_sides[0], weights[0], data_weights[0],
            coarse_spacing, 2,
        )

        level_rhs = [np.empty(shape, dtype=np.float32) for shape in shapes]
        for level in range(1, max_level):
            corrections[level] = _mg_upsample_nearest(corrections[level - 1], shapes[level])
            level_rhs[level] = right_sides[level].copy()
            for current in range(level, 0, -1):
                java_level = current + 1
                spacing = 1 << (max_level - java_level)
                _mg_smooth(
                    corrections[current], level_rhs[current], weights[current],
                    data_weights[current], spacing, 2,
                )
                residual = _mg_residue(
                    corrections[current], level_rhs[current], weights[current],
                    data_weights[current], np.float32(spacing * spacing),
                )
                corrections[current - 1] = np.zeros(shapes[current - 1], dtype=np.float32)
                level_rhs[current - 1] = _mg_downsample(residual)

            _mg_smooth(
                corrections[0], level_rhs[0], weights[0], data_weights[0],
                coarse_spacing, 2,
            )
            for current in range(1, level + 1):
                interpolated = _mg_upsample_nearest(corrections[current - 1], shapes[current])
                np.add(corrections[current], interpolated, out=corrections[current])
                java_level = current + 1
                spacing = 1 << (max_level - java_level)
                _mg_smooth(
                    corrections[current], level_rhs[current], weights[current],
                    data_weights[current], spacing, 2,
                )

        np.add(field, corrections[-1], out=field)
    return field


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
    # The source returns vector components along its voxel axes. Reorienting
    # before calculating them makes those array axes RAS axes in the output.
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
    """Apply VandyCRUISE Gradient3d's central/one-sided finite differences."""
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
    """Apply its six-neighbor Laplacian with edge indices clamped in-bounds."""
    laplacian = np.zeros(values.shape, dtype=np.float64)
    for axis, length in enumerate(values.shape):
        indices_before = np.maximum(np.arange(length) - 1, 0)
        indices_after = np.minimum(np.arange(length) + 1, length - 1)
        laplacian += np.take(values, indices_before, axis=axis)
        laplacian += np.take(values, indices_after, axis=axis)
        laplacian -= 2.0 * values
    return laplacian
