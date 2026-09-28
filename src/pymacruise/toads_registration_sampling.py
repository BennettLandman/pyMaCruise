"""Transform-dependent atlas fields for the TOADS rigid-registration kernels.

This module composes the already source-checked scalar transform sampler,
trilinear prior derivative, and rigid coordinate Jacobian. It deliberately
keeps those numerical primitives separate so their boundary behavior remains
visible and independently testable.
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
from pymacruise.toads_shapes import sample_toad_shape_prior


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
        return stream_toads_registration_energy(
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
        )

    history: list[ToadsRigidRegistrationIteration] = []
    current_lambda = np.float32(1.0)
    current_energy = energy_at(transform)
    if not np.isfinite(current_energy):
        raise ValueError("TOADS registration energy is not finite; no positive prior support was sampled")
    stop_reason = "iteration_limit"

    for _iteration in range(max_iterations):
        coefficients: RegistrationCoefficientSums = stream_toads_registration_coefficients(
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
