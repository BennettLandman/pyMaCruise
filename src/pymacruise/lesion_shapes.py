"""Selected LesionTOADS class-wise rigid shape-prior sampling."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.toads_registration import (
    RegistrationCoefficientSums,
    exact_registration_coefficient_contributions,
    rigid_image_to_template_jacobian,
    trilinear_prior_derivatives,
)
from pymacruise.toads_shapes import sample_toad_shape_prior


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
    if not np.isfinite(values).all():
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
    )
    if bounds.template_xyz is None:
        return LesionShapeSample(default, None)
    mapped = bounds.template_xyz
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
                    selected = np.zeros(mem.shape[-1], dtype=bool)
                    selected[class_index] = True
                    class_priors = np.zeros(mem.shape[-1], dtype=np.float32)
                    class_priors[class_index] = prior[x, y, z, class_index]
                    class_gradients = np.zeros((mem.shape[-1], 3), dtype=np.float32)
                    class_gradients[class_index] = prior_grad[x, y, z, class_index]
                    one_term = exact_registration_coefficient_contributions(
                        mem[x, y, z],
                        class_priors,
                        class_gradients,
                        jacobian[x, y, z],
                        registered_classes=selected,
                    )
                    cost = np.float32(cost + one_term.cost)
                    norm = np.float32(norm + one_term.norm)
                    gradient = np.asarray(gradient + one_term.gradient, dtype=np.float32)
                    hessian = np.asarray(hessian + one_term.hessian, dtype=np.float32)
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


def _triple_int(
    values: tuple[int, ...], name: str, *, positive: bool = False
) -> tuple[int, int, int]:
    if len(values) != 3 or any(not isinstance(item, (int, np.integer)) for item in values):
        raise ValueError(f"{name} must contain three integers")
    result = tuple(int(item) for item in values)
    if positive and any(item <= 0 for item in result):
        raise ValueError(f"{name} values must be positive")
    return result
