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
    coordinate_jacobian: NDArray[np.float32]
    coefficients: RegistrationCoefficientSums


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
    method does. Multiple/deformable transforms are outside this slice.
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
        if values.shape != image_dims:
            raise ValueError("untransformed shape dimensions must match the input grid")
        return LesionShapeSample(np.float32(values[point]), None)
    if mode != "rigid":
        raise ValueError("LesionTOADS shape sampling supports only 'none' or 'rigid' here")

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
        scale=1.0,
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
    final division by total prior-square norm. Multiple-transform registration
    and volume iteration are outside this single-voxel slice.
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
            image_center_xyz=image_center_xyz,
            template_center_xyz=template_center_xyz,
        )
        sampled_priors[class_index] = sample.value
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
    )
    coefficients = exact_registration_coefficient_contributions(
        membership_values,
        sampled_priors,
        sampled_gradients,
        jacobian,
    )
    return LesionRegistrationVoxelContributions(
        sampled_priors, sampled_gradients, jacobian, coefficients
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
