"""Transform-dependent atlas fields for the TOADS rigid-registration kernels.

This module composes the already source-checked scalar transform sampler,
trilinear prior derivative, and rigid coordinate Jacobian. It deliberately
keeps those numerical primitives separate so their boundary behavior remains
visible and independently testable.

Paper context: the sampled atlas priors and rigid mapping provide the
registration component associated with the TOADS atlas framework (Bazin &
Pham, 2007, 2008). The selected VandyCRUISE implementation is the behavior
authority; these kernels' synthetic comparisons do not demonstrate full atlas
alignment or whole-method equivalence. See ``docs/papers/toads.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.toads_registration import (
    RegistrationCoefficientSums,
    RegistrationLineSearchResult,
    exact_registration_measure_coefficients,
    exact_registration_gradient_descent_step,
    rigid_image_to_template_jacobian,
    trilinear_prior_derivatives,
)
from pymacruise.toads_shapes import _rotation_matrix, sample_toad_shape_prior


@dataclass(frozen=True)
class ToadsRegistrationFields:
    """Sampled unsloped atlas values and derivatives on the input voxel grid.

    ``priors`` has shape ``(x, y, z, class)``; ``prior_gradients_xyz`` has
    shape ``(x, y, z, class, 3)`` and is derivative with respect to atlas
    template coordinates. ``template_coordinates_xyz`` is the rigidly mapped
    template coordinate for each input voxel, regardless of a class's support
    bounds. Arrays are float32 and use x/y/z index order.
    """

    priors: NDArray[np.float32]
    prior_gradients_xyz: NDArray[np.float32]
    template_coordinates_xyz: NDArray[np.float32]


@dataclass(frozen=True)
class ToadsRigidRegistrationIteration:
    """Summary of one source-shaped outer rigid-registration iteration."""

    transform: NDArray[np.float32]
    energy: np.float32
    lambda_value: np.float32
    relative_change: np.float32
    accepted: bool
    rejected_trials: int


@dataclass(frozen=True)
class ToadsRigidRegistrationResult:
    """Rigid registration result and auditable outer-loop history."""

    transform: NDArray[np.float32]
    energy: np.float32
    iterations: tuple[ToadsRigidRegistrationIteration, ...]
    stop_reason: str


@dataclass(frozen=True)
class ToadsRigidRegistrationPyramidLevelResult:
    """One source-scale registration result and its transform handoff."""

    scale: int
    initial_transform: NDArray[np.float32]
    result: ToadsRigidRegistrationResult


@dataclass(frozen=True)
class ToadsRigidRegistrationPyramidResult:
    """Coarse-to-fine registration result; source omits a scale-1 stage."""

    transform: NDArray[np.float32]
    levels: tuple[ToadsRigidRegistrationPyramidLevelResult, ...]


@dataclass(frozen=True)
class ToadsRegistrationPyramidLevel:
    """Source-compatible coarse TOADS memberships and atlas priors."""

    scale: int
    memberships: NDArray[np.float32]
    atlas_shapes: NDArray[np.float32]


def source_block_mean_subsample(image: ArrayLike, factor: int) -> NDArray[np.float32]:
    """Recreate VandyCRUISE `ImageFunctions.subsample` for 3D/4D arrays.

    Complete non-overlapping factor-cubed blocks are averaged in source
    `x/y/z` offset order using float32 additions. Trailing partial blocks are
    dropped; any trailing dimensions (normally classes) are kept unchanged.
    """
    values = np.asarray(image, dtype=np.float32)
    if values.ndim not in (3, 4) or any(size <= 0 for size in values.shape[:3]):
        raise ValueError("pyramid input must be a nonempty 3D or 4D array")
    if not np.isfinite(values).all():
        raise ValueError("pyramid input must be finite")
    if not isinstance(factor, (int, np.integer)) or factor < 1:
        raise ValueError("factor must be a positive integer")
    output_dims = tuple(size // factor for size in values.shape[:3])
    if any(size <= 0 for size in output_dims):
        raise ValueError("factor is too large for the spatial dimensions")
    output = np.zeros(output_dims + values.shape[3:], dtype=np.float32)
    for dx in range(factor):
        for dy in range(factor):
            for dz in range(factor):
                output += values[
                    dx : output_dims[0] * factor : factor,
                    dy : output_dims[1] * factor : factor,
                    dz : output_dims[2] * factor : factor,
                    ...,
                ]
    output /= np.float32(factor * factor * factor)
    return output


def build_toads_registration_pyramid(
    memberships: ArrayLike, atlas_shapes: ArrayLike, *, levels: int
) -> tuple[ToadsRegistrationPyramidLevel, ...]:
    """Create the selected exact-membership pyramid, coarsest level first.

    `registerShapesExactPyramid` starts at factor `2**(levels - 1)`, halves
    the factor each stage, and stops at factor 2; scale 1 is not optimized by
    that method. Each level is subsampled directly from the full-resolution
    source arrays.
    """
    mem = np.asarray(memberships, dtype=np.float32)
    shapes = np.asarray(atlas_shapes, dtype=np.float32)
    if mem.ndim != 4 or mem.shape[-1] < 2 or not np.isfinite(mem).all():
        raise ValueError("memberships must be finite 4D (x, y, z, class) values")
    if shapes.ndim != 4 or shapes.shape[-1] != mem.shape[-1] or not np.isfinite(shapes).all():
        raise ValueError("atlas_shapes must be finite 4D values with matching classes")
    if not isinstance(levels, (int, np.integer)) or levels < 1:
        raise ValueError("levels must be a positive integer")
    return tuple(
        ToadsRegistrationPyramidLevel(
            scale,
            source_block_mean_subsample(mem, scale),
            source_block_mean_subsample(shapes, scale),
        )
        for scale in (2**power for power in range(levels - 1, 0, -1))
    )


def sample_toads_registration_fields(
    atlas_shapes: ArrayLike,
    *,
    input_dimensions_xyz: tuple[int, int, int],
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    transform: tuple[float, float, float, float, float, float],
    class_bounds_low_xyz: ArrayLike,
    class_bounds_high_xyz: ArrayLike,
) -> ToadsRegistrationFields:
    """Sample all class priors and derivatives under one TOADS rigid map.

    Atlas arrays must be stored as ``(template_x, template_y, template_z,
    class)``. Bounds are the source's per-class input-grid lower and upper
    vectors; the upper vectors retain the selected source's historical
    exclusive-bound behavior. Outside a class bound, the source default is
    one for background class 0 and zero for every other class, with zero
    derivative. Registration priors are raw atlas values: unlike
    :func:`sample_toad_shape_prior`'s default segmentation path, this function
    does not apply shape-slope recovery or clipping.

    This composes validated kernels; it does not optimize the transform or
    claim registration/segmentation equivalence. The dense derivative output
    is intended for synthetic validation and small crops. A later production
    registration driver should stream voxel contributions to avoid retaining
    this five-dimensional array for full-resolution inputs.
    """
    shapes = np.asarray(atlas_shapes, dtype=np.float32)
    if shapes.ndim != 4 or any(size < 2 for size in shapes.shape[:3]):
        raise ValueError("atlas_shapes must be 4D (template x, y, z, class) with dimensions >= 2")
    if shapes.shape[3] < 2 or not np.isfinite(shapes).all():
        raise ValueError("atlas_shapes must contain at least two finite class priors")
    dimensions_array = np.asarray(input_dimensions_xyz)
    if (
        dimensions_array.shape != (3,)
        or not np.issubdtype(dimensions_array.dtype, np.integer)
        or np.any(dimensions_array <= 0)
    ):
        raise ValueError("input_dimensions_xyz must contain three positive integers")
    dimensions = tuple(int(value) for value in dimensions_array)
    low = np.asarray(class_bounds_low_xyz)
    high = np.asarray(class_bounds_high_xyz)
    expected_bounds_shape = (shapes.shape[3], 3)
    if low.shape != expected_bounds_shape or high.shape != expected_bounds_shape:
        raise ValueError("class bounds must each have shape (class, 3)")
    if not np.issubdtype(low.dtype, np.integer) or not np.issubdtype(high.dtype, np.integer):
        raise ValueError("class bounds must contain integers")

    priors = np.empty(dimensions + (shapes.shape[3],), dtype=np.float32)
    derivatives = np.zeros(dimensions + (shapes.shape[3], 3), dtype=np.float32)
    coordinates = np.empty(dimensions + (3,), dtype=np.float32)
    full_low = (0, 0, 0)
    full_high = dimensions

    for x in range(dimensions[0]):
        for y in range(dimensions[1]):
            for z in range(dimensions[2]):
                point = (x, y, z)
                # The full-grid bounds make the mapped coordinate available
                # even when every class-specific bound excludes this voxel.
                coordinate_sample = sample_toad_shape_prior(
                    shapes[..., 0],
                    point,
                    image_dimensions_xyz=dimensions,
                    image_resolutions_mm=image_resolutions_mm,
                    template_resolutions_mm=template_resolutions_mm,
                    transform=transform,
                    class_index=0,
                    bbox_low_xyz=full_low,
                    bbox_high_xyz=full_high,
                    shape_slope=1.0,
                )
                assert coordinate_sample.template_xyz is not None
                coordinates[x, y, z] = coordinate_sample.template_xyz
                for class_index in range(shapes.shape[3]):
                    sample = sample_toad_shape_prior(
                        shapes[..., class_index],
                        point,
                        image_dimensions_xyz=dimensions,
                        image_resolutions_mm=image_resolutions_mm,
                        template_resolutions_mm=template_resolutions_mm,
                        transform=transform,
                        class_index=class_index,
                        bbox_low_xyz=tuple(int(v) for v in low[class_index]),
                        bbox_high_xyz=tuple(int(v) for v in high[class_index]),
                        shape_slope=1.0,
                    )
                    priors[x, y, z, class_index] = sample.value
                    if sample.template_xyz is not None:
                        derivatives[x, y, z, class_index] = trilinear_prior_derivatives(
                            shapes[..., class_index], sample.template_xyz
                        )

    return ToadsRegistrationFields(priors, derivatives, coordinates)


def rigid_registration_coordinate_jacobians(
    input_dimensions_xyz: tuple[int, int, int],
    transform: tuple[float, float, float, float, float, float],
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
) -> NDArray[np.float32]:
    """Build the source rigid-map Jacobian field for a small input grid.

    The returned array is ``(x, y, z, template_coordinate, parameter)`` with
    final dimensions 3×6. Like the atlas derivative field, this dense helper
    supports validation and small test grids. Production integration should
    evaluate a single 3×6 Jacobian as each voxel is accumulated.
    """
    dimensions = np.asarray(input_dimensions_xyz)
    if (
        dimensions.shape != (3,)
        or not np.issubdtype(dimensions.dtype, np.integer)
        or np.any(dimensions <= 0)
    ):
        raise ValueError("input_dimensions_xyz must contain three positive integers")
    dims = tuple(int(value) for value in dimensions)
    center = tuple(np.float32(size / 2.0) for size in dims)
    result = np.empty(dims + (3, 6), dtype=np.float32)
    for x in range(dims[0]):
        for y in range(dims[1]):
            for z in range(dims[2]):
                result[x, y, z] = rigid_image_to_template_jacobian(
                    (x, y, z),
                    transform,
                    center,
                    image_resolutions_mm,
                    template_resolutions_mm,
                )
    return result


def stream_toads_registration_energy(
    memberships: ArrayLike,
    atlas_shapes: ArrayLike,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    transform: tuple[float, float, float, float, float, float],
    class_bounds_low_xyz: ArrayLike,
    class_bounds_high_xyz: ArrayLike,
    registered_classes: ArrayLike | None = None,
    offset: int = 0,
    subsample: int = 1,
) -> np.float32:
    """Evaluate the exact-membership objective without storing atlas fields."""
    mem, shapes, low, high, registered = _registration_inputs(
        memberships, atlas_shapes, class_bounds_low_xyz, class_bounds_high_xyz, registered_classes
    )
    _validate_scan(mem.shape[:3], offset, subsample)
    cost = np.float32(0.0)
    norm = np.float32(0.0)
    for x in range(offset, mem.shape[0], subsample):
        for y in range(offset, mem.shape[1], subsample):
            for z in range(offset, mem.shape[2], subsample):
                for class_index in range(1, mem.shape[-1]):
                    if not registered[class_index]:
                        continue
                    sample = sample_toad_shape_prior(
                        shapes[..., class_index],
                        (x, y, z),
                        image_dimensions_xyz=mem.shape[:3],
                        image_resolutions_mm=image_resolutions_mm,
                        template_resolutions_mm=template_resolutions_mm,
                        transform=transform,
                        class_index=class_index,
                        bbox_low_xyz=tuple(int(v) for v in low[class_index]),
                        bbox_high_xyz=tuple(int(v) for v in high[class_index]),
                        shape_slope=1.0,
                    )
                    p = sample.value
                    if p > np.float32(0.0):
                        m = mem[x, y, z, class_index]
                        m2 = np.float32(m * m)
                        p2 = np.float32(p * p)
                        cost = np.float32(cost + np.float32(m2 * p2))
                        norm = np.float32(norm + p2)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.float32(cost / norm)


def stream_toads_registration_coefficients(
    memberships: ArrayLike,
    atlas_shapes: ArrayLike,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    transform: tuple[float, float, float, float, float, float],
    class_bounds_low_xyz: ArrayLike,
    class_bounds_high_xyz: ArrayLike,
    registered_classes: ArrayLike | None = None,
    offset: int = 0,
    subsample: int = 1,
) -> RegistrationCoefficientSums:
    """Accumulate TOADS rigid coefficients voxel-at-a-time in source order.

    Only one class vector, its derivatives, and one 3×6 Jacobian are resident
    at a time. Float32 additions occur in the original x/y/z/class order, so
    the result can be compared exactly with both the Java-formula harness and
    the dense P18s/P18t path.
    """
    mem, shapes, low, high, registered = _registration_inputs(
        memberships, atlas_shapes, class_bounds_low_xyz, class_bounds_high_xyz, registered_classes
    )
    _validate_scan(mem.shape[:3], offset, subsample)
    cost = np.float32(0.0)
    norm = np.float32(0.0)
    gradient = np.zeros(6, dtype=np.float32)
    hessian = np.zeros(6, dtype=np.float32)
    norm_gradient = np.zeros(6, dtype=np.float32)
    norm_hessian = np.zeros(6, dtype=np.float32)
    center = tuple(np.float32(size / 2.0) for size in mem.shape[:3])
    parameters = tuple(float(v) for v in transform)

    for x in range(offset, mem.shape[0], subsample):
        for y in range(offset, mem.shape[1], subsample):
            for z in range(offset, mem.shape[2], subsample):
                point = (x, y, z)
                jacobian = rigid_image_to_template_jacobian(
                    point,
                    parameters,
                    center,
                    image_resolutions_mm,
                    template_resolutions_mm,
                )
                for class_index in range(1, mem.shape[-1]):
                    if not registered[class_index]:
                        continue
                    sample = sample_toad_shape_prior(
                        shapes[..., class_index],
                        point,
                        image_dimensions_xyz=mem.shape[:3],
                        image_resolutions_mm=image_resolutions_mm,
                        template_resolutions_mm=template_resolutions_mm,
                        transform=parameters,
                        class_index=class_index,
                        bbox_low_xyz=tuple(int(v) for v in low[class_index]),
                        bbox_high_xyz=tuple(int(v) for v in high[class_index]),
                        shape_slope=1.0,
                    )
                    p = sample.value
                    if p <= np.float32(0.0):
                        continue
                    assert sample.template_xyz is not None
                    prior_gradient = trilinear_prior_derivatives(
                        shapes[..., class_index], sample.template_xyz
                    )
                    m = mem[x, y, z, class_index]
                    m2 = np.float32(m * m)
                    p2 = np.float32(p * p)
                    cost = np.float32(cost + np.float32(np.float32(m2 * p) * p))
                    norm = np.float32(norm + p2)
                    dp = np.empty(6, dtype=np.float32)
                    for parameter in range(6):
                        value = np.float32(prior_gradient[0] * jacobian[0, parameter])
                        value = np.float32(value + np.float32(prior_gradient[1] * jacobian[1, parameter]))
                        dp[parameter] = np.float32(value + np.float32(prior_gradient[2] * jacobian[2, parameter]))
                    for parameter in range(6):
                        gradient_term = np.float32(np.float32(2.0) * m)
                        gradient_term = np.float32(gradient_term * m)
                        gradient_term = np.float32(gradient_term * p)
                        gradient_term = np.float32(gradient_term * dp[parameter])
                        gradient[parameter] = np.float32(gradient[parameter] + gradient_term)
                        hessian_term = np.float32(np.float32(2.0) * m)
                        hessian_term = np.float32(hessian_term * m)
                        hessian_term = np.float32(hessian_term * dp[parameter])
                        hessian_term = np.float32(hessian_term * dp[parameter])
                        hessian[parameter] = np.float32(hessian[parameter] + hessian_term)
                        norm_gradient_term = np.float32(np.float32(np.float32(2.0) * p) * dp[parameter])
                        norm_gradient[parameter] = np.float32(norm_gradient[parameter] + norm_gradient_term)
                        norm_hessian_term = np.float32(np.float32(np.float32(2.0) * dp[parameter]) * dp[parameter])
                        norm_hessian[parameter] = np.float32(norm_hessian[parameter] + norm_hessian_term)
    return RegistrationCoefficientSums(
        cost, norm, gradient, hessian, norm_gradient, norm_hessian
    )


def _registration_inputs(
    memberships: ArrayLike,
    atlas_shapes: ArrayLike,
    class_bounds_low_xyz: ArrayLike,
    class_bounds_high_xyz: ArrayLike,
    registered_classes: ArrayLike | None,
) -> tuple[NDArray[np.float32], NDArray[np.float32], NDArray[np.int64], NDArray[np.int64], NDArray[np.bool_]]:
    mem = np.asarray(memberships, dtype=np.float32)
    shapes = np.asarray(atlas_shapes, dtype=np.float32)
    if mem.ndim != 4 or any(size < 1 for size in mem.shape[:3]) or mem.shape[-1] < 2:
        raise ValueError("memberships must be a nonempty 4D (x, y, z, class) array")
    if shapes.ndim != 4 or shapes.shape[-1] != mem.shape[-1] or any(size < 2 for size in shapes.shape[:3]):
        raise ValueError("atlas_shapes must be 4D with one 3D prior per membership class")
    if not np.isfinite(mem).all() or not np.isfinite(shapes).all():
        raise ValueError("memberships and atlas shapes must be finite")
    low = np.asarray(class_bounds_low_xyz)
    high = np.asarray(class_bounds_high_xyz)
    expected_shape = (mem.shape[-1], 3)
    if low.shape != expected_shape or high.shape != expected_shape:
        raise ValueError("class bounds must each have shape (class, 3)")
    if not np.issubdtype(low.dtype, np.integer) or not np.issubdtype(high.dtype, np.integer):
        raise ValueError("class bounds must contain integers")
    registered = (
        np.ones(mem.shape[-1], dtype=bool)
        if registered_classes is None
        else np.asarray(registered_classes, dtype=bool)
    )
    if registered.shape != (mem.shape[-1],):
        raise ValueError("registered_classes must contain one flag per class")
    return mem, shapes, low.astype(np.int64, copy=False), high.astype(np.int64, copy=False), registered


def _validate_scan(
    dimensions: tuple[int, int, int], offset: int, subsample: int
) -> None:
    if not isinstance(offset, (int, np.integer)) or offset < 0:
        raise ValueError("offset must be a non-negative integer")
    if not isinstance(subsample, (int, np.integer)) or subsample < 1:
        raise ValueError("subsample must be a positive integer")
    if offset >= min(dimensions):
        raise ValueError("offset must be inside every spatial axis")


def vectorized_toads_registration_energy(
    memberships: ArrayLike,
    atlas_shapes: ArrayLike,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    transform: tuple[float, float, float, float, float, float],
    class_bounds_low_xyz: ArrayLike,
    class_bounds_high_xyz: ArrayLike,
    registered_classes: ArrayLike | None = None,
    offset: int = 0,
    subsample: int = 1,
    scale: float = 1.0,
    image_center_xyz: ArrayLike | None = None,
    template_center_xyz: ArrayLike | None = None,
) -> np.float32:
    """NumPy-slab exact-membership energy in source x/y/z/class order."""
    context = _vector_context(
        memberships, atlas_shapes, transform, image_resolutions_mm,
        template_resolutions_mm, class_bounds_low_xyz, class_bounds_high_xyz,
        registered_classes, offset, subsample, scale, image_center_xyz,
        template_center_xyz,
    )
    (mem, shapes, params, image_res, template_res, low, high, registered,
     rotation, ys, zs, scale_value, image_center, template_center) = context
    cost = np.float32(0.0)
    norm = np.float32(0.0)
    for x in range(offset, mem.shape[0], subsample):
        coordinates = _map_registration_slab(
            x, ys, zs, params, image_res, template_res, rotation,
            scale_value, image_center, template_center,
        )
        y_grid, z_grid = ys[:, None], zs[None, :]
        cost_terms = np.zeros((ys.size, zs.size, mem.shape[-1]), dtype=np.float32)
        norm_terms = np.zeros_like(cost_terms)
        for k in range(1, mem.shape[-1]):
            if not registered[k]:
                continue
            inside = (
                (x >= low[k, 0]) & (x < high[k, 0])
                & (y_grid >= low[k, 1]) & (y_grid < high[k, 1])
                & (z_grid >= low[k, 2]) & (z_grid < high[k, 2])
            )
            prior = _linear_closest_grid(shapes[..., k], coordinates)
            prior = np.where(inside, prior, np.float32(0.0)).astype(np.float32, copy=False)
            membership = mem[x, y_grid, z_grid, k]
            membership_square = np.float32(membership * membership)
            prior_square = np.float32(prior * prior)
            cost_terms[..., k] = np.float32(membership_square * prior_square)
            norm_terms[..., k] = prior_square
        cost = _accumulate_ordered(cost, cost_terms)
        norm = _accumulate_ordered(norm, norm_terms)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.float32(cost / norm)


def vectorized_toads_registration_coefficients(
    memberships: ArrayLike,
    atlas_shapes: ArrayLike,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    transform: tuple[float, float, float, float, float, float],
    class_bounds_low_xyz: ArrayLike,
    class_bounds_high_xyz: ArrayLike,
    registered_classes: ArrayLike | None = None,
    offset: int = 0,
    subsample: int = 1,
    scale: float = 1.0,
    image_center_xyz: ArrayLike | None = None,
    template_center_xyz: ArrayLike | None = None,
) -> RegistrationCoefficientSums:
    """Vectorize mapping/interpolation by x slab, preserving source sum order."""
    context = _vector_context(
        memberships, atlas_shapes, transform, image_resolutions_mm,
        template_resolutions_mm, class_bounds_low_xyz, class_bounds_high_xyz,
        registered_classes, offset, subsample, scale, image_center_xyz,
        template_center_xyz,
    )
    (mem, shapes, params, image_res, template_res, low, high, registered,
     rotation, ys, zs, scale_value, center, template_center) = context
    cost = np.float32(0.0)
    norm = np.float32(0.0)
    gradient = np.zeros(6, dtype=np.float32)
    hessian = np.zeros(6, dtype=np.float32)
    norm_gradient = np.zeros(6, dtype=np.float32)
    norm_hessian = np.zeros(6, dtype=np.float32)
    for x in range(offset, mem.shape[0], subsample):
        coordinates = _map_registration_slab(
            x, ys, zs, params, image_res, template_res, rotation,
            scale_value, center, template_center,
        )
        jacobian = _registration_jacobian_slab(
            x, ys, zs, center, image_res, template_res, params, scale_value
        )
        y_grid, z_grid = ys[:, None], zs[None, :]
        slab_shape = (ys.size, zs.size, mem.shape[-1])
        cost_terms = np.zeros(slab_shape, dtype=np.float32)
        norm_terms = np.zeros(slab_shape, dtype=np.float32)
        gradient_terms = np.zeros(slab_shape + (6,), dtype=np.float32)
        hessian_terms = np.zeros_like(gradient_terms)
        norm_gradient_terms = np.zeros_like(gradient_terms)
        norm_hessian_terms = np.zeros_like(gradient_terms)

        for k in range(1, mem.shape[-1]):
            if not registered[k]:
                continue
            inside = (
                (x >= low[k, 0]) & (x < high[k, 0])
                & (y_grid >= low[k, 1]) & (y_grid < high[k, 1])
                & (z_grid >= low[k, 2]) & (z_grid < high[k, 2])
            )
            prior = _linear_closest_grid(shapes[..., k], coordinates)
            prior = np.where(inside, prior, np.float32(0.0)).astype(np.float32, copy=False)
            contributing = prior > np.float32(0.0)
            if not np.any(contributing):
                continue
            prior_gradient = _trilinear_derivatives_grid(shapes[..., k], coordinates)
            membership = mem[x, y_grid, z_grid, k]
            membership_square = np.float32(membership * membership)
            prior_square = np.float32(prior * prior)
            cost_terms[..., k] = np.where(
                contributing,
                np.float32(np.float32(membership_square * prior) * prior),
                np.float32(0.0),
            )
            norm_terms[..., k] = np.where(contributing, prior_square, np.float32(0.0))
            for parameter in range(6):
                dp = np.float32(prior_gradient[..., 0] * jacobian[..., 0, parameter])
                dp = np.float32(dp + np.float32(prior_gradient[..., 1] * jacobian[..., 1, parameter]))
                dp = np.float32(dp + np.float32(prior_gradient[..., 2] * jacobian[..., 2, parameter]))

                gradient_term = np.float32(np.float32(np.float32(2.0) * membership) * membership)
                gradient_term = np.float32(gradient_term * prior)
                gradient_term = np.float32(gradient_term * dp)
                gradient_terms[..., k, parameter] = np.where(contributing, gradient_term, 0.0)

                hessian_term = np.float32(np.float32(np.float32(2.0) * membership) * membership)
                hessian_term = np.float32(hessian_term * dp)
                hessian_term = np.float32(hessian_term * dp)
                hessian_terms[..., k, parameter] = np.where(contributing, hessian_term, 0.0)

                norm_gradient_term = np.float32(np.float32(np.float32(2.0) * prior) * dp)
                norm_gradient_terms[..., k, parameter] = np.where(contributing, norm_gradient_term, 0.0)

                norm_hessian_term = np.float32(np.float32(np.float32(2.0) * dp) * dp)
                norm_hessian_terms[..., k, parameter] = np.where(contributing, norm_hessian_term, 0.0)

        cost = _accumulate_ordered(cost, cost_terms)
        norm = _accumulate_ordered(norm, norm_terms)
        for parameter in range(6):
            gradient[parameter] = _accumulate_ordered(gradient[parameter], gradient_terms[..., parameter])
            hessian[parameter] = _accumulate_ordered(hessian[parameter], hessian_terms[..., parameter])
            norm_gradient[parameter] = _accumulate_ordered(norm_gradient[parameter], norm_gradient_terms[..., parameter])
            norm_hessian[parameter] = _accumulate_ordered(norm_hessian[parameter], norm_hessian_terms[..., parameter])
    return RegistrationCoefficientSums(cost, norm, gradient, hessian, norm_gradient, norm_hessian)


def _vector_context(
    memberships: ArrayLike,
    atlas_shapes: ArrayLike,
    transform: ArrayLike,
    image_resolutions_mm: ArrayLike,
    template_resolutions_mm: ArrayLike,
    class_bounds_low_xyz: ArrayLike,
    class_bounds_high_xyz: ArrayLike,
    registered_classes: ArrayLike | None,
    offset: int,
    subsample: int,
    scale: float,
    image_center_xyz: ArrayLike | None,
    template_center_xyz: ArrayLike | None,
) -> tuple[
    NDArray[np.float32], NDArray[np.float32], NDArray[np.float32],
    NDArray[np.float32], NDArray[np.float32], NDArray[np.int64],
    NDArray[np.int64], NDArray[np.bool_], NDArray[np.float32], NDArray[np.int64], NDArray[np.int64],
    np.float32, tuple[np.float32, np.float32, np.float32], tuple[np.float32, np.float32, np.float32]
]:
    mem, shapes, low, high, registered = _registration_inputs(
        memberships, atlas_shapes, class_bounds_low_xyz, class_bounds_high_xyz, registered_classes
    )
    _validate_scan(mem.shape[:3], offset, subsample)
    dims = np.asarray(mem.shape[:3], dtype=np.int64)
    empty_bounds = np.all(low == dims, axis=1) & np.all(high == 0, axis=1)
    if np.any(low < 0) or np.any(high < 0) or np.any(low > dims) or np.any(high > dims):
        raise ValueError("class bounds must be inside the declared input grid")
    if np.any(np.any(low > high, axis=1) & ~empty_bounds):
        raise ValueError("class bounds must be ordered unless they encode empty support")
    params = np.asarray(transform, dtype=np.float32)
    image_res = np.asarray(image_resolutions_mm, dtype=np.float32)
    template_res = np.asarray(template_resolutions_mm, dtype=np.float32)
    if params.shape != (6,) or not np.isfinite(params).all():
        raise ValueError("transform must contain six finite values")
    if image_res.shape != (3,) or not np.isfinite(image_res).all() or np.any(image_res <= 0):
        raise ValueError("image_resolutions_mm must contain three positive finite values")
    if template_res.shape != (3,) or not np.isfinite(template_res).all() or np.any(template_res <= 0):
        raise ValueError("template_resolutions_mm must contain three positive finite values")
    scale_value = np.float32(scale)
    if not np.isfinite(scale_value) or scale_value <= 0:
        raise ValueError("scale must be positive and finite")
    if scale_value != np.float32(1.0) and (
        image_center_xyz is None or template_center_xyz is None
    ):
        raise ValueError(
            "coarse-scale registration requires full-resolution image and template centers"
        )
    center = np.asarray(
        image_center_xyz if image_center_xyz is not None else tuple(size / 2.0 for size in mem.shape[:3]),
        dtype=np.float32,
    )
    template_center = np.asarray(
        template_center_xyz if template_center_xyz is not None else tuple(size / 2.0 for size in shapes.shape[:3]),
        dtype=np.float32,
    )
    if center.shape != (3,) or not np.isfinite(center).all():
        raise ValueError("image_center_xyz must contain three finite values")
    if template_center.shape != (3,) or not np.isfinite(template_center).all():
        raise ValueError("template_center_xyz must contain three finite values")
    rotation = _rotation_matrix(params[:3])
    ys = np.arange(offset, mem.shape[1], subsample, dtype=np.int64)
    zs = np.arange(offset, mem.shape[2], subsample, dtype=np.int64)
    return (
        mem, shapes, params, image_res, template_res, low, high, registered,
        rotation, ys, zs, scale_value, tuple(center), tuple(template_center),
    )


def _map_registration_slab(
    x: int,
    y_indices: NDArray[np.int64],
    z_indices: NDArray[np.int64],
    transform: NDArray[np.float32],
    image_res: NDArray[np.float32],
    template_res: NDArray[np.float32],
    rotation: NDArray[np.float32],
    scale: np.float32,
    image_center: tuple[np.float32, np.float32, np.float32],
    template_center: tuple[np.float32, np.float32, np.float32],
) -> NDArray[np.float32]:
    y_grid = y_indices[:, None].astype(np.float32)
    z_grid = z_indices[None, :].astype(np.float32)
    physical = (
        np.float32(np.float32(np.float32(x) * scale - image_center[0]) * image_res[0]),
        np.float32(np.float32(y_grid * scale - image_center[1]) * image_res[1]),
        np.float32(np.float32(z_grid * scale - image_center[2]) * image_res[2]),
    )
    coordinates = np.empty((y_indices.size, z_indices.size, 3), dtype=np.float32)
    for axis in range(3):
        value = np.float32(rotation[axis, 0] * physical[0])
        value = np.float32(value + np.float32(rotation[axis, 1] * physical[1]))
        value = np.float32(value + np.float32(rotation[axis, 2] * physical[2]))
        value = np.float32(value + transform[axis + 3])
        value = np.float32(value / template_res[axis])
        value = np.float32(value + template_center[axis])
        coordinates[..., axis] = np.float32(value / scale)
    return coordinates


def _linear_closest_grid(values: NDArray[np.float32], coordinates: NDArray[np.float32]) -> NDArray[np.float32]:
    indices = [
        np.clip(np.floor(coordinates[..., axis]).astype(np.int64), 0, values.shape[axis] - 2)
        for axis in range(3)
    ]
    alpha = [
        np.clip(np.float32(coordinates[..., axis] - indices[axis].astype(np.float32)), 0.0, 1.0).astype(np.float32)
        for axis in range(3)
    ]
    beta = [np.float32(np.float32(1.0) - value) for value in alpha]
    x, y, z = indices
    ax, ay, az = alpha
    bx, by, bz = beta

    def weighted(weight_a: NDArray[np.float32], weight_b: NDArray[np.float32], weight_c: NDArray[np.float32], voxel: NDArray[np.float32]) -> NDArray[np.float32]:
        return np.float32(np.float32(np.float32(weight_a * weight_b) * weight_c) * voxel)

    terms = (
        weighted(bx, by, bz, values[x, y, z]),
        weighted(ax, by, bz, values[x + 1, y, z]),
        weighted(bx, ay, bz, values[x, y + 1, z]),
        weighted(bx, by, az, values[x, y, z + 1]),
        weighted(ax, ay, bz, values[x + 1, y + 1, z]),
        weighted(bx, ay, az, values[x, y + 1, z + 1]),
        weighted(ax, by, az, values[x + 1, y, z + 1]),
        weighted(ax, ay, az, values[x + 1, y + 1, z + 1]),
    )
    total = np.float32(terms[0] + terms[1])
    for term in terms[2:]:
        total = np.float32(total + term)
    return total


def _trilinear_derivatives_grid(
    values: NDArray[np.float32], coordinates: NDArray[np.float32]
) -> NDArray[np.float32]:
    valid = np.ones(coordinates.shape[:2], dtype=bool)
    for axis in range(3):
        valid &= (coordinates[..., axis] >= 0.0) & (coordinates[..., axis] <= values.shape[axis] - 2)
    result = np.zeros(coordinates.shape[:2] + (3,), dtype=np.float32)
    if not np.any(valid):
        return result
    indices = [
        np.clip(np.floor(coordinates[..., axis]).astype(np.int64), 0, values.shape[axis] - 2)
        for axis in range(3)
    ]
    alpha, beta, gamma = [
        np.float32(coordinates[..., axis] - indices[axis].astype(np.float32))
        for axis in range(3)
    ]
    nalpha, nbeta, ngamma = [np.float32(np.float32(1.0) - q) for q in (alpha, beta, gamma)]
    x, y, z = indices

    def product(*factors: NDArray[np.float32]) -> NDArray[np.float32]:
        value = np.asarray(factors[0], dtype=np.float32)
        for factor in factors[1:]:
            value = np.float32(value * factor)
        return value

    def difference(a: tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.int64]], b: tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.int64]]) -> NDArray[np.float32]:
        return np.float32(values[a] - values[b])

    dx_terms = (
        product(nbeta, ngamma, difference((x + 1, y, z), (x, y, z))),
        product(beta, ngamma, difference((x + 1, y + 1, z), (x, y + 1, z))),
        product(nbeta, gamma, difference((x + 1, y, z + 1), (x, y, z + 1))),
        product(beta, gamma, difference((x + 1, y + 1, z + 1), (x, y + 1, z + 1))),
    )
    dy_terms = (
        product(nalpha, ngamma, difference((x, y + 1, z), (x, y, z))),
        product(alpha, ngamma, difference((x + 1, y + 1, z), (x + 1, y, z))),
        product(nalpha, gamma, difference((x, y + 1, z + 1), (x, y, z + 1))),
        product(alpha, gamma, difference((x + 1, y + 1, z + 1), (x + 1, y, z + 1))),
    )
    dz_terms = (
        product(nalpha, nbeta, difference((x, y, z + 1), (x, y, z))),
        product(alpha, nbeta, difference((x + 1, y, z + 1), (x + 1, y, z))),
        product(nalpha, beta, difference((x, y + 1, z + 1), (x, y + 1, z))),
        product(alpha, beta, difference((x + 1, y + 1, z + 1), (x + 1, y + 1, z))),
    )
    for axis, terms in enumerate((dx_terms, dy_terms, dz_terms)):
        total = np.float32(terms[0] + terms[1])
        total = np.float32(total + terms[2])
        result[..., axis] = np.where(valid, np.float32(total + terms[3]), 0.0)
    return result


def _registration_jacobian_slab(
    x: int,
    y_indices: NDArray[np.int64],
    z_indices: NDArray[np.int64],
    center: tuple[np.float32, np.float32, np.float32],
    image_res: NDArray[np.float32],
    template_res: NDArray[np.float32],
    transform: NDArray[np.float32],
    scale: np.float32,
) -> NDArray[np.float32]:
    basis_jacobians = []
    for basis in ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)):
        basis_jacobians.append(
            rigid_image_to_template_jacobian(
                basis, transform, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0)
            )[:, :3]
        )
    y_grid = y_indices[:, None].astype(np.float32)
    z_grid = z_indices[None, :].astype(np.float32)
    physical = (
        np.float32(np.float32(np.float32(x) * scale - center[0]) * image_res[0]),
        np.float32(np.float32(y_grid * scale - center[1]) * image_res[1]),
        np.float32(np.float32(z_grid * scale - center[2]) * image_res[2]),
    )
    result = np.zeros((y_indices.size, z_indices.size, 3, 6), dtype=np.float32)
    for axis in range(3):
        for parameter in range(3):
            value = np.float32(basis_jacobians[0][axis, parameter] * physical[0])
            value = np.float32(value + np.float32(basis_jacobians[1][axis, parameter] * physical[1]))
            value = np.float32(value + np.float32(basis_jacobians[2][axis, parameter] * physical[2]))
            value = np.float32(value / template_res[axis])
            result[..., axis, parameter] = np.float32(value / scale)
        result[..., axis, axis + 3] = np.float32(
            np.float32(np.float32(1.0) / template_res[axis]) / scale
        )
    return result


def _accumulate_ordered(total: np.float32, terms: NDArray[np.float32]) -> np.float32:
    flat = terms.reshape(-1)
    if flat.size == 0:
        return total
    # Seed the first term with the prior-slab accumulator, preserving the
    # source's single sequential float32 addition stream across x slabs.
    flat[0] = np.float32(total + flat[0])
    return np.float32(np.cumsum(flat, dtype=np.float32)[-1])


def register_toads_atlas_rigid(
    memberships: ArrayLike,
    atlas_shapes: ArrayLike,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    initial_transform: tuple[float, float, float, float, float, float],
    class_bounds_low_xyz: ArrayLike,
    class_bounds_high_xyz: ArrayLike,
    registered_classes: ArrayLike | None = None,
    offset: int = 0,
    subsample: int = 1,
    max_iterations: int = 50,
    energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
    source_faithful: bool = True,
    scale: float = 1.0,
    image_center_xyz: ArrayLike | None = None,
    template_center_xyz: ArrayLike | None = None,
) -> ToadsRigidRegistrationResult:
    """Run the selected exact-membership rigid registration loop on a small grid.

    This connects source-matched coordinate mapping, prior interpolation and
    derivatives, rigid Jacobians, coefficient accumulation/normalization, and
    the existing line-search update. The source coefficient call-order
    behavior is preserved by default. The optimizer streams voxel samples and
    coefficients; the separate dense field-building utility remains useful
    only as a small-grid validation oracle.

    Atlas shapes and memberships use ``(x, y, z, class)`` array order. The
    class-bound upper values retain the source-exclusive interpretation.
    The function estimates only the six rigid parameters; deformable TOADS
    segmentation and alignment quality are outside its contract.
    """
    mem = np.asarray(memberships, dtype=np.float32)
    shapes = np.asarray(atlas_shapes, dtype=np.float32)
    if mem.ndim != 4 or any(size < 1 for size in mem.shape[:3]) or mem.shape[-1] < 2:
        raise ValueError("memberships must be a nonempty 4D (x, y, z, class) array")
    if shapes.ndim != 4 or shapes.shape[-1] != mem.shape[-1]:
        raise ValueError("atlas_shapes must be 4D with one template prior per membership class")
    if not np.isfinite(mem).all() or not np.isfinite(shapes).all():
        raise ValueError("memberships and atlas shapes must be finite")
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer")
    tolerance = np.float32(energy_tolerance)
    lambda_floor = np.float32(minimum_lambda)
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("energy_tolerance must be finite and non-negative")
    if not np.isfinite(lambda_floor) or lambda_floor <= 0:
        raise ValueError("minimum_lambda must be positive and finite")
    transform = np.asarray(initial_transform, dtype=np.float32)
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("initial_transform must contain six finite values")
    registered = (
        np.ones(mem.shape[-1], dtype=bool)
        if registered_classes is None
        else np.asarray(registered_classes, dtype=bool)
    )
    if registered.shape != (mem.shape[-1],):
        raise ValueError("registered_classes must contain one flag per class")

    def energy_at(candidate: NDArray[np.float32]) -> np.float32:
        return vectorized_toads_registration_energy(
            mem,
            shapes,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            transform=tuple(float(value) for value in candidate),
            class_bounds_low_xyz=class_bounds_low_xyz,
            class_bounds_high_xyz=class_bounds_high_xyz,
            registered_classes=registered,
            offset=offset,
            subsample=subsample,
            scale=scale,
            image_center_xyz=image_center_xyz,
            template_center_xyz=template_center_xyz,
        )

    history: list[ToadsRigidRegistrationIteration] = []
    current_lambda = np.float32(1.0)
    current_energy = energy_at(transform)
    if not np.isfinite(current_energy):
        raise ValueError("TOADS registration energy is not finite; no positive prior support was sampled")
    stop_reason = "iteration_limit"

    for _iteration in range(max_iterations):
        coefficients: RegistrationCoefficientSums = vectorized_toads_registration_coefficients(
            mem,
            shapes,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            transform=tuple(float(value) for value in transform),
            class_bounds_low_xyz=class_bounds_low_xyz,
            class_bounds_high_xyz=class_bounds_high_xyz,
            registered_classes=registered,
            offset=offset,
            subsample=subsample,
            scale=scale,
            image_center_xyz=image_center_xyz,
            template_center_xyz=template_center_xyz,
        )
        current_energy = energy_at(transform)
        if not np.isfinite(current_energy):
            raise ValueError("TOADS registration energy is not finite; no positive prior support was sampled")
        gradient, hessian = exact_registration_measure_coefficients(
            coefficients.gradient,
            coefficients.hessian,
            coefficients.norm_gradient,
            coefficients.norm_hessian,
            coefficients.norm,
            source_faithful=source_faithful,
        )
        update: RegistrationLineSearchResult = exact_registration_gradient_descent_step(
            transform,
            gradient,
            hessian,
            current_energy,
            energy_at,
            initial_lambda=current_lambda,
        )
        transform = update.transform.copy()
        current_energy = update.energy
        current_lambda = update.lambda_value
        history.append(
            ToadsRigidRegistrationIteration(
                transform.copy(),
                np.float32(current_energy),
                np.float32(current_lambda),
                np.float32(update.relative_change),
                update.accepted,
                update.rejected_trials,
            )
        )
        if update.relative_change <= tolerance:
            stop_reason = "energy_tolerance"
            break
        if current_lambda <= lambda_floor:
            stop_reason = "lambda_floor"
            break

    return ToadsRigidRegistrationResult(
        transform.copy(), np.float32(current_energy), tuple(history), stop_reason
    )


def register_toads_atlas_rigid_pyramid(
    memberships: ArrayLike,
    atlas_shapes: ArrayLike,
    *,
    levels: int,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    initial_transform: tuple[float, float, float, float, float, float],
    registered_classes: ArrayLike | None = None,
    max_iterations: int = 50,
    energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
    source_faithful: bool = True,
) -> ToadsRigidRegistrationPyramidResult:
    """Run exact-membership rigid registration on source pyramid levels.

    Membership and atlas arrays are independently downsampled from their
    native grids. Full-resolution image/template centers are retained while
    each coarse-grid level uses its source scale. The selected source method
    does not optimize a scale-1 level. It also evaluates every registered
    class across each coarse image grid, so this wrapper uses full-grid class
    bounds at each level rather than fixed-transform support boxes.
    """
    mem = np.asarray(memberships, dtype=np.float32)
    shapes = np.asarray(atlas_shapes, dtype=np.float32)
    if mem.ndim != 4 or shapes.ndim != 4:
        raise ValueError("pyramid inputs must be 4D (x, y, z, class) arrays")
    pyramid = build_toads_registration_pyramid(mem, shapes, levels=levels)
    transform = np.asarray(initial_transform, dtype=np.float32)
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("initial_transform must contain six finite values")
    image_center = tuple(np.float32(size / 2.0) for size in mem.shape[:3])
    template_center = tuple(np.float32(size / 2.0) for size in shapes.shape[:3])
    level_results: list[ToadsRigidRegistrationPyramidLevelResult] = []
    for level in pyramid:
        dimensions = level.memberships.shape[:3]
        low = np.zeros((mem.shape[-1], 3), dtype=np.int64)
        high = np.broadcast_to(np.asarray(dimensions, dtype=np.int64), low.shape).copy()
        level_initial = transform.copy()
        result = register_toads_atlas_rigid(
            level.memberships,
            level.atlas_shapes,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            initial_transform=tuple(float(value) for value in transform),
            class_bounds_low_xyz=low,
            class_bounds_high_xyz=high,
            registered_classes=registered_classes,
            max_iterations=max_iterations,
            energy_tolerance=energy_tolerance,
            minimum_lambda=minimum_lambda,
            source_faithful=source_faithful,
            scale=level.scale,
            image_center_xyz=image_center,
            template_center_xyz=template_center,
        )
        transform = result.transform.copy()
        level_results.append(
            ToadsRigidRegistrationPyramidLevelResult(
                level.scale, level_initial, result
            )
        )
    return ToadsRigidRegistrationPyramidResult(transform.copy(), tuple(level_results))
