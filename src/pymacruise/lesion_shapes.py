"""Selected LesionTOADS class-wise rigid shape-prior sampling."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.toads_registration import (
    RegistrationCoefficientSums,
    RegistrationLineSearchResult,
    exact_registration_gradient_descent_step,
    exact_registration_coefficient_contributions,
    rigid_image_to_template_jacobian,
    trilinear_prior_derivatives,
)
from pymacruise.toads_shapes import (
    _rotation_matrix as toads_rotation_matrix,
    sample_toad_shape_prior,
)
from pymacruise.toads_registration_sampling import (
    ToadsRegistrationPyramidLevel,
    _map_registration_slab,
    _registration_jacobian_slab,
    _trilinear_derivatives_grid,
    build_toads_registration_pyramid,
)


@dataclass(frozen=True)
class LesionShapeSample:
    """One sampled prior, mapped point, and optional template-space gradient."""

    value: np.float32
    template_xyz: tuple[np.float32, np.float32, np.float32] | None
    template_gradient_xyz: tuple[np.float32, np.float32, np.float32] | None = None


@dataclass(frozen=True)
class LesionRegistrationVoxelContributions:
    """One selected single-transform registration voxel before volume scaling."""

    shape_priors: NDArray[np.float32]
    prior_gradients_xyz: NDArray[np.float32]
    template_coordinates_xyz: tuple[np.float32, np.float32, np.float32] | None
    coordinate_jacobian: NDArray[np.float32]
    coefficients: RegistrationCoefficientSums


@dataclass(frozen=True)
class LesionRegistrationVolumeCoefficients:
    """Source-ordered selected LesionTOADS coefficient reduction."""

    cost: np.float32
    norm: np.float32
    gradient: NDArray[np.float32]
    hessian: NDArray[np.float32]
    energy: np.float32
    normalized_gradient: NDArray[np.float32]
    normalized_hessian: NDArray[np.float32]
    visited_voxels: int
    contributing_class_terms: int


@dataclass(frozen=True)
class LesionRegistrationOptimizationStep:
    """One fixed-scale LesionTOADS registration step and its trial trace."""

    result: RegistrationLineSearchResult
    trial_transforms: NDArray[np.float32]
    trial_energies: NDArray[np.float32]


@dataclass(frozen=True)
class LesionRegistrationOptimizationResult:
    """Fixed-scale registration iterations and the selected loop stop cause."""

    transform: NDArray[np.float32]
    energy: np.float32
    lambda_value: np.float32
    iterations: tuple[LesionRegistrationOptimizationStep, ...]
    stop_reason: str


@dataclass(frozen=True)
class LesionRegistrationPyramidLevelResult:
    """One source scale and its fixed-scale optimization result."""

    scale: int
    input_dimensions_xyz: tuple[int, int, int]
    template_dimensions_xyz: tuple[int, int, int]
    initial_transform: NDArray[np.float32]
    optimization: LesionRegistrationOptimizationResult


@dataclass(frozen=True)
class LesionRegistrationPyramidResult:
    """Coarse-to-fine LesionTOADS rigid registration result."""

    transform: NDArray[np.float32]
    levels: tuple[LesionRegistrationPyramidLevelResult, ...]


def sample_lesion_shape_prior(
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
    transform_mode: str = "rigid",
    scale: float = 1.0,
    image_center_xyz: tuple[float, float, float] | None = None,
    template_center_xyz: tuple[float, float, float] | None = None,
    _prevalidated_shape: bool = False,
    _precomputed_rotation_matrix: np.ndarray | None = None,
    _precomputed_template_xyz: tuple[float, float, float] | None = None,
) -> LesionShapeSample:
    """Sample one class as in `LesionToadAtlas.getRegularTransformedShape`.

    Rigid mapping reuses the source-compared `ParametricTransform` port in
    :mod:`pymacruise.toads_shapes`. The selected LesionTOADS wrapper differs
    from TOADS shape membership sampling: class 0 uses interpolation default
    1, other classes use 0, out-of-template coordinates return that complete
    default value, and no `shapeSlope` recovery or closest-edge clamping is
    applied. `transform_mode="none"` directly indexes the map as the Java
    method does. For registration pyramid scales greater than one, pass the
    coarse input/map arrays and the original full-resolution image/template
    centers; scaled coordinates follow `ParametricTransform.imageToTemplate`
    without class-bound optimization. Multiple/deformable transforms are
    outside this slice.
    """
    values = np.asarray(shape, dtype=np.float32)
    if values.ndim != 3 or any(size < 2 for size in values.shape):
        raise ValueError("LesionTOADS shape prior must be 3D with dimensions >= 2")
    if not isinstance(_prevalidated_shape, (bool, np.bool_)):
        raise TypeError("_prevalidated_shape must be boolean")
    if not _prevalidated_shape and not np.isfinite(values).all():
        raise ValueError("LesionTOADS shape prior must contain only finite values")
    image_dims = _triple_int(image_dimensions_xyz, "image_dimensions_xyz", positive=True)
    point = _triple_int(image_point_xyz, "image_point_xyz")
    if any(point[axis] < 0 or point[axis] >= image_dims[axis] for axis in range(3)):
        raise ValueError("image_point_xyz must be inside the declared input grid")
    if not isinstance(class_index, (int, np.integer)) or class_index < 0:
        raise ValueError("class_index must be a nonnegative integer")
    default = np.float32(1.0 if class_index == 0 else 0.0)

    mode = transform_mode.strip().lower()
    if mode == "none":
        if np.float32(scale) != np.float32(1.0):
            raise ValueError("untransformed shape sampling requires scale=1")
        if values.shape != image_dims:
            raise ValueError("untransformed shape dimensions must match the input grid")
        return LesionShapeSample(np.float32(values[point]), None)
    if mode != "rigid":
        raise ValueError("LesionTOADS shape sampling supports only 'none' or 'rigid' here")
    scale_value = np.float32(scale)
    if not np.isfinite(scale_value) or scale_value <= np.float32(0.0):
        raise ValueError("scale must be positive and finite")
    if scale_value != np.float32(1.0) and (
        image_center_xyz is None or template_center_xyz is None
    ):
        raise ValueError("pyramid sampling requires explicit full-resolution centers")

    if _precomputed_template_xyz is None:
        bounds = sample_toad_shape_prior(
            values,
            point,
            image_dimensions_xyz=image_dims,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            transform=transform,
            class_index=int(class_index),
            bbox_low_xyz=bbox_low_xyz,
            bbox_high_xyz=bbox_high_xyz,
            shape_slope=1.0,
            scale=float(scale_value),
            image_center_xyz=image_center_xyz,
            template_center_xyz=template_center_xyz,
            _java_float32_norm=True,
            _prevalidated_shape=_prevalidated_shape,
            _precomputed_rotation_matrix=_precomputed_rotation_matrix,
        )
        if bounds.template_xyz is None:
            return LesionShapeSample(default, None)
        mapped = bounds.template_xyz
    else:
        if len(_precomputed_template_xyz) != 3 or not np.isfinite(
            np.asarray(_precomputed_template_xyz, dtype=np.float32)
        ).all():
            raise ValueError("precomputed template coordinate must be a finite XYZ triple")
        lower = _triple_int(bbox_low_xyz, "bbox_low_xyz")
        upper = _triple_int(bbox_high_xyz, "bbox_high_xyz")
        if any(lower[axis] < 0 or lower[axis] > image_dims[axis] for axis in range(3)) or any(
            upper[axis] < 0 or upper[axis] > image_dims[axis] for axis in range(3)
        ) or any(lower[axis] > upper[axis] for axis in range(3)):
            raise ValueError("class bounds must be ordered and inside the declared input grid")
        if any(point[axis] < lower[axis] or point[axis] >= upper[axis] for axis in range(3)):
            return LesionShapeSample(default, None)
        mapped = tuple(np.float32(value) for value in _precomputed_template_xyz)
    sampled = _linear_interpolation_with_default(values, mapped, default)
    # LesionToadAtlas computes derivatives only after its priorT > 0 gate.
    # Reuse the derivative primitive already compared against the selected
    # ImageFunctions source, while retaining this class's constant-outside
    # value rule above (which differs from the TOADS shape sampler).
    gradient = None
    if sampled > np.float32(0.0):
        gradient = tuple(np.float32(v) for v in trilinear_prior_derivatives(values, mapped))
    return LesionShapeSample(sampled, mapped, gradient)


def lesion_registration_voxel_contributions(
    shape_priors: Sequence[ArrayLike],
    memberships: ArrayLike,
    image_point_xyz: tuple[int, int, int],
    *,
    image_dimensions_xyz: tuple[int, int, int],
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    transform: tuple[float, float, float, float, float, float],
    class_bounds_xyz: Sequence[
        tuple[tuple[int, int, int], tuple[int, int, int]]
    ],
    scale: float = 1.0,
    image_center_xyz: tuple[float, float, float] | None = None,
    template_center_xyz: tuple[float, float, float] | None = None,
) -> LesionRegistrationVoxelContributions:
    """Assemble one source-shaped LesionTOADS rigid registration voxel.

    The selected `LesionToadAtlas.computeRegistrationCoefficients` loop skips
    class 0, samples each remaining prior with default zero, and accumulates
    only where the interpolated prior is positive. This adapter composes that
    sampling/gating behavior with the already source-checked rigid Jacobian
    and common coefficient accumulator. Returned sums are unnormalized so
    callers can accumulate voxels in source scan order before the Java-style
    final division by total prior-square norm. For pyramid levels, ``scale``
    is passed to both source-matched samplers while centers stay in original
    resolution coordinates. Multiple-transform registration and volume
    iteration are outside this single-voxel slice.
    """
    membership_values = np.asarray(memberships, dtype=np.float32)
    class_count = len(shape_priors)
    if class_count < 2 or membership_values.shape != (class_count,):
        raise ValueError("one membership and shape prior are required for each class")
    if len(class_bounds_xyz) != class_count:
        raise ValueError("class_bounds_xyz must contain one bound pair per class")
    if not np.isfinite(membership_values).all():
        raise ValueError("memberships must contain only finite values")

    sampled_priors = np.zeros(class_count, dtype=np.float32)
    sampled_gradients = np.zeros((class_count, 3), dtype=np.float32)
    mapped_coordinates = None
    for class_index in range(1, class_count):
        bbox_low_xyz, bbox_high_xyz = class_bounds_xyz[class_index]
        sample = sample_lesion_shape_prior(
            np.asarray(shape_priors[class_index], dtype=np.float32),
            image_point_xyz,
            image_dimensions_xyz=image_dimensions_xyz,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            transform=transform,
            class_index=class_index,
            bbox_low_xyz=bbox_low_xyz,
            bbox_high_xyz=bbox_high_xyz,
            scale=scale,
            image_center_xyz=image_center_xyz,
            template_center_xyz=template_center_xyz,
        )
        sampled_priors[class_index] = sample.value
        if mapped_coordinates is None and sample.template_xyz is not None:
            mapped_coordinates = sample.template_xyz
        if sample.template_gradient_xyz is not None:
            sampled_gradients[class_index] = sample.template_gradient_xyz

    center = (
        image_center_xyz
        if image_center_xyz is not None
        else tuple(np.float32(size / 2.0) for size in image_dimensions_xyz)
    )
    jacobian = rigid_image_to_template_jacobian(
        image_point_xyz,
        transform,
        center,
        image_resolutions_mm,
        template_resolutions_mm,
        scale=scale,
    )
    coefficients = exact_registration_coefficient_contributions(
        membership_values,
        sampled_priors,
        sampled_gradients,
        jacobian,
    )
    return LesionRegistrationVoxelContributions(
        sampled_priors, sampled_gradients, mapped_coordinates, jacobian, coefficients
    )


def accumulate_lesion_registration_coefficients(
    memberships: ArrayLike,
    shape_priors: ArrayLike,
    prior_gradients_xyz: ArrayLike,
    coordinate_jacobians: ArrayLike,
    *,
    registered_classes: ArrayLike | None = None,
    offset: int = 0,
    subsample: int = 3,
) -> LesionRegistrationVolumeCoefficients:
    """Reduce pre-sampled LesionTOADS terms in selected Java loop order.

    Inputs are indexed ``[x, y, z, class]``; prior gradients add a trailing
    XYZ axis and rigid coordinate Jacobians add trailing ``(3, 6)`` axes.
    The selected ``computeScaledRegistrationCoefficients`` traversal is
    x-major, then y, z, and registered class, starting at ``offset`` and
    advancing by ``subsample``. Each class term is rounded into the running
    float32 volume sums before the next class/voxel, preserving source
    accumulation order rather than summing per-voxel aggregates.

    This reducer receives already sampled arrays so its acceptance tests can
    isolate scan/reduction semantics from atlas interpolation. Source
    normalization divides cost, gradient, and diagonal Hessian by the
    prior-square norm. As in Java, zero norm produces IEEE NaN/Inf outputs.
    """
    mem = np.asarray(memberships, dtype=np.float32)
    prior = np.asarray(shape_priors, dtype=np.float32)
    prior_grad = np.asarray(prior_gradients_xyz, dtype=np.float32)
    jacobian = np.asarray(coordinate_jacobians, dtype=np.float32)
    if mem.ndim != 4 or mem.shape[-1] < 2:
        raise ValueError("memberships must have shape (x, y, z, classes>=2)")
    if prior.shape != mem.shape:
        raise ValueError("shape_priors must match the membership grid and classes")
    if prior_grad.shape != mem.shape + (3,):
        raise ValueError("prior_gradients_xyz must have shape (x, y, z, classes, 3)")
    if jacobian.shape != mem.shape[:3] + (3, 6):
        raise ValueError("coordinate_jacobians must have shape (x, y, z, 3, 6)")
    if not all(np.isfinite(values).all() for values in (mem, prior, prior_grad, jacobian)):
        raise ValueError("registration volume terms must contain only finite values")
    if registered_classes is None:
        registered = np.ones(mem.shape[-1], dtype=bool)
    else:
        registered = np.asarray(registered_classes, dtype=bool)
        if registered.shape != (mem.shape[-1],):
            raise ValueError("registered_classes must have one value per class")
    if isinstance(offset, bool) or not isinstance(offset, (int, np.integer)) or offset < 0:
        raise ValueError("offset must be a nonnegative integer")
    if (
        isinstance(subsample, bool)
        or not isinstance(subsample, (int, np.integer))
        or subsample <= 0
    ):
        raise ValueError("subsample must be a positive integer")

    cost = np.float32(0.0)
    norm = np.float32(0.0)
    gradient = np.zeros(6, dtype=np.float32)
    hessian = np.zeros(6, dtype=np.float32)
    visited = 0
    terms = 0
    for x in range(int(offset), mem.shape[0], int(subsample)):
        for y in range(int(offset), mem.shape[1], int(subsample)):
            for z in range(int(offset), mem.shape[2], int(subsample)):
                visited += 1
                for class_index in range(1, mem.shape[-1]):
                    if not registered[class_index] or prior[x, y, z, class_index] <= np.float32(0.0):
                        continue
                    # The selected Java traversal adds one class contribution
                    # at a time. Keep its float32 rounding order directly here
                    # instead of allocating full class-sized temporary arrays
                    # and routing this single active class through the generic
                    # vector adapter on every contributing voxel.
                    p = prior[x, y, z, class_index]
                    m = mem[x, y, z, class_index]
                    m2 = np.float32(m * m)
                    p2 = np.float32(p * p)
                    cost = np.float32(
                        cost + np.float32(np.float32(m2 * p) * p)
                    )
                    norm = np.float32(norm + p2)
                    for parameter in range(6):
                        value = np.float32(
                            prior_grad[x, y, z, class_index, 0]
                            * jacobian[x, y, z, 0, parameter]
                        )
                        value = np.float32(value + np.float32(
                            prior_grad[x, y, z, class_index, 1]
                            * jacobian[x, y, z, 1, parameter]
                        ))
                        dp = np.float32(value + np.float32(
                            prior_grad[x, y, z, class_index, 2]
                            * jacobian[x, y, z, 2, parameter]
                        ))
                        gradient_term = np.float32(np.float32(2.0) * m)
                        gradient_term = np.float32(gradient_term * m)
                        gradient_term = np.float32(gradient_term * p)
                        gradient_term = np.float32(gradient_term * dp)
                        gradient[parameter] = np.float32(
                            gradient[parameter] + gradient_term
                        )
                        hessian_term = np.float32(np.float32(2.0) * m)
                        hessian_term = np.float32(hessian_term * m)
                        hessian_term = np.float32(hessian_term * dp)
                        hessian_term = np.float32(hessian_term * dp)
                        hessian[parameter] = np.float32(
                            hessian[parameter] + hessian_term
                        )
                    terms += 1
    with np.errstate(divide="ignore", invalid="ignore"):
        energy = np.divide(cost, norm, dtype=np.float32)
        normalized_gradient = np.divide(gradient, norm, dtype=np.float32)
        normalized_hessian = np.divide(hessian, norm, dtype=np.float32)
    return LesionRegistrationVolumeCoefficients(
        cost,
        norm,
        gradient,
        hessian,
        energy,
        normalized_gradient,
        normalized_hessian,
        visited,
        terms,
    )


def lesion_registration_volume_coefficients(
    shape_priors: Sequence[ArrayLike],
    memberships: ArrayLike,
    *,
    image_dimensions_xyz: tuple[int, int, int],
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    transform: tuple[float, float, float, float, float, float],
    scale: float = 1.0,
    image_center_xyz: tuple[float, float, float] | None = None,
    template_center_xyz: tuple[float, float, float] | None = None,
    offset: int = 0,
    subsample: int = 3,
) -> LesionRegistrationVolumeCoefficients:
    """Sample coarse LesionTOADS maps and reduce fixed-transform coefficients.

    ``memberships`` is indexed ``[x, y, z, class]`` on the coarse input grid;
    each shape map is already loaded as an XYZ array on the corresponding
    registration-pyramid atlas grid. Pyramid centers, when provided, remain
    in full-resolution coordinates. The selected scaled method does not use
    optimized per-class input bounds, so each registered class uses the full
    coarse input extent. Class 0 is not sampled or registered.

    This composes fixed-transform sampling and coefficient reduction only; it
    does not perform a line search or update the transform.
    """
    mem = np.asarray(memberships, dtype=np.float32)
    if mem.ndim != 4 or mem.shape[-1] < 2:
        raise ValueError("memberships must have shape (x, y, z, classes>=2)")
    image_dims = _triple_int(image_dimensions_xyz, "image_dimensions_xyz", positive=True)
    if mem.shape[:3] != image_dims:
        raise ValueError("membership grid must match image_dimensions_xyz")
    if len(shape_priors) != mem.shape[-1]:
        raise ValueError("one shape prior is required for each membership class")
    if not np.isfinite(mem).all():
        raise ValueError("memberships must contain only finite values")
    maps = [np.asarray(shape, dtype=np.float32) for shape in shape_priors]
    for values in maps:
        if values.ndim != 3 or any(size < 2 for size in values.shape):
            raise ValueError("each LesionTOADS shape prior must be 3D with dimensions >= 2")
        if not np.isfinite(values).all():
            raise ValueError("LesionTOADS shape priors must contain only finite values")
    if (
        isinstance(offset, bool)
        or not isinstance(offset, (int, np.integer))
        or offset < 0
    ):
        raise ValueError("offset must be a nonnegative integer")
    if (
        isinstance(subsample, bool)
        or not isinstance(subsample, (int, np.integer))
        or subsample <= 0
    ):
        raise ValueError("subsample must be a positive integer")

    priors = np.zeros(mem.shape, dtype=np.float32)
    gradients = np.zeros(mem.shape + (3,), dtype=np.float32)
    jacobians = np.zeros(mem.shape[:3] + (3, 6), dtype=np.float32)
    image_res = np.asarray(image_resolutions_mm, dtype=np.float32)
    template_res = np.asarray(template_resolutions_mm, dtype=np.float32)
    if image_res.shape != (3,) or not np.isfinite(image_res).all() or np.any(image_res == 0):
        raise ValueError("image_resolutions_mm must contain three finite nonzero signed values")
    if template_res.shape != (3,) or not np.isfinite(template_res).all() or np.any(template_res <= 0):
        raise ValueError("template_resolutions_mm must contain three positive finite values")
    transform_values = np.asarray(transform, dtype=np.float32)
    if transform_values.shape != (6,) or not np.isfinite(transform_values).all():
        raise ValueError("transform must contain six finite rigid parameters")
    rotation_matrix = toads_rotation_matrix(
        transform_values[:3], java_float32_norm=True
    )
    scale_value = np.float32(scale)
    if not np.isfinite(scale_value) or scale_value <= 0:
        raise ValueError("scale must be positive and finite")
    if scale_value != np.float32(1.0) and (
        image_center_xyz is None or template_center_xyz is None
    ):
        raise ValueError("pyramid sampling requires explicit full-resolution centers")
    center_values = np.asarray(
        image_center_xyz
        if image_center_xyz is not None
        else tuple(np.float32(size / 2.0) for size in image_dims),
        dtype=np.float32,
    )
    if center_values.shape != (3,) or not np.isfinite(center_values).all():
        raise ValueError("image_center_xyz must contain three finite values")
    center = tuple(np.float32(value) for value in center_values)
    template_center_values = np.asarray(
        template_center_xyz
        if template_center_xyz is not None
        else tuple(np.float32(size / 2.0) for size in maps[0].shape),
        dtype=np.float32,
    )
    if template_center_values.shape != (3,) or not np.isfinite(template_center_values).all():
        raise ValueError("template_center_xyz must contain three finite values")
    template_center = tuple(np.float32(value) for value in template_center_values)
    image_indices_y = np.arange(int(offset), image_dims[1], int(subsample), dtype=np.int64)
    image_indices_z = np.arange(int(offset), image_dims[2], int(subsample), dtype=np.int64)
    shared_template_grid = len({values.shape for values in maps}) == 1
    full_input_bounds = (tuple(0 for _ in range(3)), image_dims)
    y_grid, z_grid = np.meshgrid(image_indices_y, image_indices_z, indexing="ij")
    for x in range(int(offset), image_dims[0], int(subsample)):
        jacobian_slab = _registration_jacobian_slab(
            x,
            image_indices_y,
            image_indices_z,
            center,
            image_res,
            template_res,
            transform_values,
            scale_value,
        )
        if shared_template_grid:
            mapped_grid = _map_registration_slab(
                x,
                image_indices_y,
                image_indices_z,
                transform_values,
                image_res,
                template_res,
                rotation_matrix,
                scale_value,
                center,
                template_center,
            )
            for class_index in range(1, mem.shape[-1]):
                shape = maps[class_index]
                priors[x, y_grid, z_grid, class_index] = _linear_interpolation_default_grid(
                    shape, mapped_grid, np.float32(0.0)
                )
                gradients[x, y_grid, z_grid, class_index] = _trilinear_derivatives_grid(
                    shape, mapped_grid
                )
            jacobians[x, y_grid, z_grid] = jacobian_slab
            continue
        for y_index, y in enumerate(image_indices_y):
            for z_index, z in enumerate(image_indices_z):
                point = (x, y, z)
                for class_index in range(1, mem.shape[-1]):
                    sample = sample_lesion_shape_prior(
                        maps[class_index],
                        point,
                        image_dimensions_xyz=image_dims,
                        image_resolutions_mm=image_resolutions_mm,
                        template_resolutions_mm=template_resolutions_mm,
                        transform=transform,
                        class_index=class_index,
                        bbox_low_xyz=full_input_bounds[0],
                        bbox_high_xyz=full_input_bounds[1],
                        scale=scale,
                        image_center_xyz=image_center_xyz,
                        template_center_xyz=template_center_xyz,
                        _prevalidated_shape=True,
                        _precomputed_rotation_matrix=rotation_matrix,
                    )
                    priors[x, y, z, class_index] = sample.value
                    if sample.template_gradient_xyz is not None:
                        gradients[x, y, z, class_index] = sample.template_gradient_xyz
                jacobians[x, y, z] = jacobian_slab[y_index, z_index]
    return accumulate_lesion_registration_coefficients(
        mem,
        priors,
        gradients,
        jacobians,
        offset=offset,
        subsample=subsample,
    )


def register_lesion_shapes_scaled_step(
    shape_priors: Sequence[ArrayLike],
    memberships: ArrayLike,
    transform: ArrayLike,
    *,
    image_dimensions_xyz: tuple[int, int, int],
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    scale: float = 1.0,
    image_center_xyz: tuple[float, float, float] | None = None,
    template_center_xyz: tuple[float, float, float] | None = None,
    offset: int = 0,
    subsample: int = 3,
    initial_lambda: float = 1.0,
    lambda_factor: float = 1.5,
    max_rejections: int = 10,
) -> LesionRegistrationOptimizationStep:
    """Take one scaled rigid-registration step for the fixed atlas fields.

    Coefficients come from :func:`lesion_registration_volume_coefficients`;
    each line-search trial recomputes the sampled volume energy at that
    candidate transform. The source-compatible diagonal update, strict
    maximization rule, lambda adaptation, and rejection limit are delegated to
    the existing checked registration helper. No outer iteration is run.
    """
    current = np.asarray(transform, dtype=np.float32)
    if current.shape != (6,) or not np.isfinite(current).all():
        raise ValueError("transform must be a finite six-value rigid parameter vector")
    coefficients = lesion_registration_volume_coefficients(
        shape_priors,
        memberships,
        image_dimensions_xyz=image_dimensions_xyz,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        transform=tuple(current),
        scale=scale,
        image_center_xyz=image_center_xyz,
        template_center_xyz=template_center_xyz,
        offset=offset,
        subsample=subsample,
    )
    trials: list[NDArray[np.float32]] = []
    energies: list[np.float32] = []

    def evaluate(candidate: NDArray[np.float32]) -> np.float32:
        trial = candidate.copy()
        energy = lesion_registration_volume_coefficients(
            shape_priors,
            memberships,
            image_dimensions_xyz=image_dimensions_xyz,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            transform=tuple(trial),
            scale=scale,
            image_center_xyz=image_center_xyz,
            template_center_xyz=template_center_xyz,
            offset=offset,
            subsample=subsample,
        ).energy
        trials.append(trial)
        energies.append(energy)
        return energy

    result = exact_registration_gradient_descent_step(
        current,
        coefficients.normalized_gradient,
        coefficients.normalized_hessian,
        coefficients.energy,
        evaluate,
        initial_lambda=initial_lambda,
        lambda_factor=lambda_factor,
        max_rejections=max_rejections,
    )
    trial_array = np.asarray(trials, dtype=np.float32).reshape((-1, 6))
    energy_array = np.asarray(energies, dtype=np.float32)
    return LesionRegistrationOptimizationStep(result, trial_array, energy_array)


def optimize_lesion_shapes_scaled(
    shape_priors: Sequence[ArrayLike],
    memberships: ArrayLike,
    initial_transform: ArrayLike,
    *,
    image_dimensions_xyz: tuple[int, int, int],
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    scale: float = 1.0,
    image_center_xyz: tuple[float, float, float] | None = None,
    template_center_xyz: tuple[float, float, float] | None = None,
    offset: int = 0,
    subsample: int = 3,
    max_iterations: int = 50,
    min_energy_change: float = 1e-6,
    min_lambda: float = 0.001,
    initial_lambda: float = 1.0,
    lambda_factor: float = 1.5,
    max_rejections: int = 10,
) -> LesionRegistrationOptimizationResult:
    """Repeat fixed-scale LesionTOADS updates until a source loop stop.

    This follows the outer `registerShapesPyramid` scale loop body only. The
    caller supplies already selected coarse maps/memberships and this
    function does not construct a pyramid or advance to another scale.
    """
    if (
        isinstance(max_iterations, bool)
        or not isinstance(max_iterations, (int, np.integer))
        or max_iterations < 0
    ):
        raise ValueError("max_iterations must be a nonnegative integer")
    min_change = np.float32(min_energy_change)
    lambda_floor = np.float32(min_lambda)
    if not np.isfinite(min_change) or min_change < 0:
        raise ValueError("min_energy_change must be finite and nonnegative")
    if not np.isfinite(lambda_floor) or lambda_floor < 0:
        raise ValueError("min_lambda must be finite and nonnegative")
    transform = np.asarray(initial_transform, dtype=np.float32)
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("initial_transform must be a finite six-value rigid vector")
    current = lesion_registration_volume_coefficients(
        shape_priors,
        memberships,
        image_dimensions_xyz=image_dimensions_xyz,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
        transform=tuple(transform),
        scale=scale,
        image_center_xyz=image_center_xyz,
        template_center_xyz=template_center_xyz,
        offset=offset,
        subsample=subsample,
    )
    energy = current.energy
    lam = np.float32(initial_lambda)
    difference = np.float32(1.0)
    history: list[LesionRegistrationOptimizationStep] = []
    while (
        len(history) < int(max_iterations)
        and difference > min_change
        and lam > lambda_floor
    ):
        step = register_lesion_shapes_scaled_step(
            shape_priors,
            memberships,
            transform,
            image_dimensions_xyz=image_dimensions_xyz,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            scale=scale,
            image_center_xyz=image_center_xyz,
            template_center_xyz=template_center_xyz,
            offset=offset,
            subsample=subsample,
            initial_lambda=lam,
            lambda_factor=lambda_factor,
            max_rejections=max_rejections,
        )
        history.append(step)
        transform = step.result.transform.copy()
        energy = step.result.energy
        lam = step.result.lambda_value
        difference = step.result.relative_change
    if len(history) >= int(max_iterations):
        reason = "iteration_limit"
    elif difference <= min_change:
        reason = "min_energy_change"
    else:
        reason = "lambda_floor"
    return LesionRegistrationOptimizationResult(
        transform.copy(), energy, lam, tuple(history), reason
    )


def register_lesion_shapes_pyramid(
    shape_priors: Sequence[ArrayLike],
    memberships: ArrayLike,
    initial_transform: ArrayLike,
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    levels: int,
    image_center_xyz: tuple[float, float, float] | None = None,
    template_center_xyz: tuple[float, float, float] | None = None,
    offset: int = 0,
    subsample: int = 2,
    max_iterations: int = 50,
    min_energy_change: float = 1e-6,
    min_lambda: float = 0.001,
    lambda_factor: float = 1.5,
    max_rejections: int = 10,
) -> LesionRegistrationPyramidResult:
    """Run selected coarse-to-fine fixed-transform LesionTOADS registration.

    Memberships and atlas priors are downsampled directly from the supplied
    full-resolution arrays with the source-checked block-mean helper. The
    selected `registerShapesPyramid` visits scales ``2**(levels-1)`` down to
    2 and has no scale-1 stage. Each scale resets lambda to its initial value,
    while the transform from the coarser level is carried to the next finer
    level. Atlas interpolation and per-level optimization are delegated to
    the previously validated helpers.
    """
    if isinstance(levels, bool) or not isinstance(levels, (int, np.integer)) or levels < 1:
        raise ValueError("levels must be a positive integer")
    mem = np.asarray(memberships, dtype=np.float32)
    if mem.ndim != 4 or mem.shape[-1] < 2:
        raise ValueError("memberships must have shape (x, y, z, classes>=2)")
    if len(shape_priors) != mem.shape[-1]:
        raise ValueError("one atlas shape prior is required for each membership class")
    shape_arrays = [np.asarray(shape, dtype=np.float32) for shape in shape_priors]
    if any(shape.ndim != 3 for shape in shape_arrays):
        raise ValueError("each atlas shape prior must be a 3D array")
    if len({shape.shape for shape in shape_arrays}) != 1:
        raise ValueError("all atlas shape priors must share one spatial grid")
    atlas_stack = np.stack(shape_arrays, axis=-1)
    pyramid: tuple[ToadsRegistrationPyramidLevel, ...] = build_toads_registration_pyramid(
        mem, atlas_stack, levels=int(levels)
    )
    transform = np.asarray(initial_transform, dtype=np.float32)
    if transform.shape != (6,) or not np.isfinite(transform).all():
        raise ValueError("initial_transform must be a finite six-value rigid vector")
    full_image_center = (
        image_center_xyz
        if image_center_xyz is not None
        else tuple(np.float32(size / 2.0) for size in mem.shape[:3])
    )
    full_template_center = (
        template_center_xyz
        if template_center_xyz is not None
        else tuple(np.float32(size / 2.0) for size in atlas_stack.shape[:3])
    )
    history: list[LesionRegistrationPyramidLevelResult] = []
    for level in pyramid:
        input_dims = tuple(int(size) for size in level.memberships.shape[:3])
        template_dims = tuple(int(size) for size in level.atlas_shapes.shape[:3])
        if any(size < 2 for size in template_dims):
            raise ValueError("pyramid level leaves atlas map dimensions below two voxels")
        coarse_shapes = tuple(
            level.atlas_shapes[..., class_index]
            for class_index in range(level.atlas_shapes.shape[-1])
        )
        stage_initial_transform = transform.copy()
        result = optimize_lesion_shapes_scaled(
            coarse_shapes,
            level.memberships,
            transform,
            image_dimensions_xyz=input_dims,
            image_resolutions_mm=image_resolutions_mm,
            template_resolutions_mm=template_resolutions_mm,
            scale=float(level.scale),
            image_center_xyz=full_image_center,
            template_center_xyz=full_template_center,
            offset=offset,
            subsample=subsample,
            max_iterations=max_iterations,
            min_energy_change=min_energy_change,
            min_lambda=min_lambda,
            initial_lambda=1.0,
            lambda_factor=lambda_factor,
            max_rejections=max_rejections,
        )
        transform = result.transform
        history.append(
            LesionRegistrationPyramidLevelResult(
                level.scale,
                input_dims,
                template_dims,
                stage_initial_transform,
                result,
            )
        )
    return LesionRegistrationPyramidResult(transform.copy(), tuple(history))


def _linear_interpolation_with_default(
    image: np.ndarray,
    point_xyz: tuple[np.float32, np.float32, np.float32],
    default: np.float32,
) -> np.float32:
    """Match ImageFunctions.linearInterpolation(value,...), including bounds."""
    if any(
        point_xyz[axis] < np.float32(0.0)
        or point_xyz[axis] > np.float32(image.shape[axis] - 2)
        for axis in range(3)
    ):
        return np.float32(default)

    lower = tuple(int(np.floor(point_xyz[axis])) for axis in range(3))
    fraction = tuple(
        np.float32(point_xyz[axis] - np.float32(lower[axis])) for axis in range(3)
    )
    complement = tuple(np.float32(np.float32(1.0) - fraction[axis]) for axis in range(3))
    x, y, z = lower
    a, b, c = fraction
    na, nb, nc = complement
    terms = (
        np.float32(np.float32(np.float32(na * nb) * nc) * image[x, y, z]),
        np.float32(np.float32(np.float32(a * nb) * nc) * image[x + 1, y, z]),
        np.float32(np.float32(np.float32(na * b) * nc) * image[x, y + 1, z]),
        np.float32(np.float32(np.float32(na * nb) * c) * image[x, y, z + 1]),
        np.float32(np.float32(np.float32(a * b) * nc) * image[x + 1, y + 1, z]),
        np.float32(np.float32(np.float32(na * b) * c) * image[x, y + 1, z + 1]),
        np.float32(np.float32(np.float32(a * nb) * c) * image[x + 1, y, z + 1]),
        np.float32(np.float32(np.float32(a * b) * c) * image[x + 1, y + 1, z + 1]),
    )
    result = np.float32(terms[0] + terms[1])
    for term in terms[2:]:
        result = np.float32(result + term)
    return result


def _linear_interpolation_default_grid(
    image: np.ndarray,
    points_xyz: np.ndarray,
    default: np.float32,
) -> NDArray[np.float32]:
    """Apply LesionTOADS' zero-outside interpolation to a coordinate slab.

    This retains ``ImageFunctions.linearInterpolation`` bounds and term order
    while evaluating independent points with NumPy float32 operations. It is
    used only by the shared-grid registration path; exact selected Java volume
    fixtures guard its values and accumulation.
    """
    values = np.asarray(image, dtype=np.float32)
    points = np.asarray(points_xyz, dtype=np.float32)
    if values.ndim != 3 or points.ndim != 3 or points.shape[-1] != 3:
        raise ValueError("image and coordinate slab must be 3D with an XYZ trailing axis")
    valid = np.ones(points.shape[:2], dtype=bool)
    for axis in range(3):
        valid &= (points[..., axis] >= np.float32(0.0)) & (
            points[..., axis] <= np.float32(values.shape[axis] - 2)
        )
    indices = []
    alpha = []
    for axis in range(3):
        safe_coordinate = np.where(valid, points[..., axis], np.float32(0.0))
        index = np.floor(safe_coordinate).astype(np.int64)
        fraction = np.float32(safe_coordinate - index.astype(np.float32))
        indices.append(index)
        alpha.append(fraction)
    complement = [np.float32(np.float32(1.0) - fraction) for fraction in alpha]
    x, y, z = indices
    a, b, c = alpha
    na, nb, nc = complement
    terms = (
        np.float32(np.float32(np.float32(na * nb) * nc) * values[x, y, z]),
        np.float32(np.float32(np.float32(a * nb) * nc) * values[x + 1, y, z]),
        np.float32(np.float32(np.float32(na * b) * nc) * values[x, y + 1, z]),
        np.float32(np.float32(np.float32(na * nb) * c) * values[x, y, z + 1]),
        np.float32(np.float32(np.float32(a * b) * nc) * values[x + 1, y + 1, z]),
        np.float32(np.float32(np.float32(na * b) * c) * values[x, y + 1, z + 1]),
        np.float32(np.float32(np.float32(a * nb) * c) * values[x + 1, y, z + 1]),
        np.float32(np.float32(np.float32(a * b) * c) * values[x + 1, y + 1, z + 1]),
    )
    result = np.float32(terms[0] + terms[1])
    for term in terms[2:]:
        result = np.float32(result + term)
    return np.where(valid, result, np.float32(default)).astype(np.float32, copy=False)


def _triple_int(
    values: tuple[int, ...], name: str, *, positive: bool = False
) -> tuple[int, int, int]:
    if len(values) != 3 or any(not isinstance(item, (int, np.integer)) for item in values):
        raise ValueError(f"{name} must contain three integers")
    result = tuple(int(item) for item in values)
    if positive and any(item <= 0 for item in result):
        raise ValueError(f"{name} values must be positive")
    return result
