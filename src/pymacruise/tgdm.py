"""Source-derived local geometry and force terms for VandyCRUISE TGDM."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

TGDMProfile = Literal["inner", "central", "outer"]


@dataclass(frozen=True)
class TGDMGeometry:
    """Finite-difference quantities computed by `GenericTGDM.updatePoint`."""

    gradient: np.ndarray
    mean_curvature_numerator: np.ndarray
    gaussian_curvature_numerator: np.ndarray
    backward: np.ndarray
    forward: np.ndarray


def tgdm_geometry(level_set: np.ndarray) -> TGDMGeometry:
    """Compute source mirror-boundary first/second derivatives and curvature."""
    phi = np.asarray(level_set, dtype=np.float32)
    if phi.ndim != 3 or any(length < 2 for length in phi.shape):
        raise ValueError("TGDM geometry requires a 3D field with dimensions >= 2")
    if not np.isfinite(phi).all():
        raise ValueError("TGDM level set must be finite")

    before = tuple(_shift(phi, axis, -1) for axis in range(3))
    after = tuple(_shift(phi, axis, 1) for axis in range(3))
    backward = np.stack(tuple(phi - n for n in before), axis=-1)
    forward = np.stack(tuple(n - phi for n in after), axis=-1)
    gradient = (np.stack(after, -1) - np.stack(before, -1)) * np.float32(0.5)
    dxx, dyy, dzz = tuple(before[a] + after[a] - 2.0 * phi for a in range(3))
    dxy = _cross_difference(phi, 0, 1)
    dxz = _cross_difference(phi, 0, 2)
    dyz = _cross_difference(phi, 1, 2)
    gx, gy, gz = (gradient[..., axis] for axis in range(3))
    sx, sy, sz = gx * gx, gy * gy, gz * gz

    mean = (
        (dyy + dzz) * sx
        + (dxx + dzz) * sy
        + (dxx + dyy) * sz
        - 2.0 * (gx * gy * dxy + gx * gz * dxz + gy * gz * dyz)
    )
    gaussian = (
        (dyy * dzz - dyz * dyz) * sx
        + (dxx * dzz - dxz * dxz) * sy
        + (dxx * dyy - dxy * dxy) * sz
        + 2.0
        * (
            gx * gy * (dxz * dyz - dxy * dzz)
            + gx * gz * (dxy * dyz - dxz * dyy)
            + gy * gz * (dxy * dxz - dyz * dxx)
        )
    )
    return TGDMGeometry(
        np.asarray(gradient, dtype=np.float32),
        np.asarray(mean, dtype=np.float32),
        np.asarray(gaussian, dtype=np.float32),
        np.asarray(backward, dtype=np.float32),
        np.asarray(forward, dtype=np.float32),
    )


def tgdm_force(
    geometry: TGDMGeometry,
    *,
    profile: TGDMProfile,
    curvature_force: float,
    pressure_force: float,
    external_force: float = 0.0,
    region_force: np.ndarray | None = None,
    advection: np.ndarray | None = None,
) -> np.ndarray:
    """Evaluate one source profile's `getEffectiveForce` over a whole volume."""
    gradient = np.asarray(geometry.gradient, dtype=np.float64)
    mean = np.asarray(geometry.mean_curvature_numerator, dtype=np.float64)
    gaussian = np.asarray(geometry.gaussian_curvature_numerator, dtype=np.float64)
    backward = np.asarray(geometry.backward, dtype=np.float64)
    forward = np.asarray(geometry.forward, dtype=np.float64)
    if gradient.ndim != 4 or gradient.shape[-1] != 3:
        raise ValueError("TGDM geometry must contain three derivative components")
    if profile not in ("inner", "central", "outer"):
        raise ValueError("profile must be 'inner', 'central', or 'outer'")
    if not np.isfinite([curvature_force, pressure_force, external_force]).all():
        raise ValueError("TGDM force weights must be finite")

    region_coefficient = (
        np.full(gradient.shape[:-1], pressure_force, dtype=np.float64)
        if region_force is None
        else pressure_force * np.asarray(region_force, dtype=np.float64)
    )
    if region_coefficient.shape != gradient.shape[:-1] or not np.isfinite(region_coefficient).all():
        raise ValueError("region_force must be a finite scalar field on the TGDM grid")
    vector = (
        np.zeros_like(gradient)
        if advection is None
        else np.asarray(advection, dtype=np.float64)
    )
    if vector.shape != gradient.shape or not np.isfinite(vector).all():
        raise ValueError("advection must be a finite three-component TGDM field")

    grad_magnitude = np.linalg.norm(gradient, axis=-1)
    if profile == "inner":
        denominator = grad_magnitude + 0.001
        normalized_mean = mean / denominator**3
        normalized_gaussian = gaussian / denominator**4
        transformed_mean = normalized_mean.copy()
        discriminant = normalized_mean**2 - 2.0 * normalized_gaussian
        positive = discriminant > 0
        transformed_mean[positive] = (
            normalized_mean[positive] * normalized_gaussian[positive] / discriminant[positive]
        )
        curvature_term = curvature_force * transformed_mean * denominator
    elif profile == "outer":
        denominator = grad_magnitude + 0.001
        special_case = (mean < -1.5) & (gaussian < 0)
        region_coefficient = np.where(special_case, 0.0, region_coefficient)
        curvature_multiplier = np.where(special_case, 3.0 * curvature_force, curvature_force)
        transformed_mean = mean.copy()
        discriminant = mean**2 - 2.0 * gaussian
        positive = discriminant > 0
        transformed_mean[positive] = mean[positive] * gaussian[positive] / discriminant[positive]
        transformed_mean = np.clip(transformed_mean, -1.5, 1.5)
        curvature_term = curvature_multiplier * transformed_mean * denominator
    else:
        transformed_mean = np.zeros_like(mean)
        nonzero_gradient = grad_magnitude > 1e-5
        normalized_mean = np.zeros_like(mean)
        normalized_gaussian = np.zeros_like(gaussian)
        normalized_mean[nonzero_gradient] = (
            mean[nonzero_gradient] / grad_magnitude[nonzero_gradient] ** 3
        )
        normalized_gaussian[nonzero_gradient] = (
            gaussian[nonzero_gradient] / grad_magnitude[nonzero_gradient] ** 4
        )
        discriminant = normalized_mean**2 - 2.0 * normalized_gaussian
        positive = nonzero_gradient & (discriminant > 0)
        transformed_mean[positive] = (
            normalized_mean[positive]
            * normalized_gaussian[positive]
            / discriminant[positive]
        )
        curvature_term = curvature_force * transformed_mean * grad_magnitude

    upwind = np.where(
        region_coefficient[..., None] > 0,
        np.maximum(backward, 0.0) + np.minimum(forward, 0.0),
        np.minimum(backward, 0.0) + np.maximum(forward, 0.0),
    )
    region_term = region_coefficient * np.linalg.norm(upwind, axis=-1)
    scaled_advection = vector * external_force
    advection_gradient = np.where(scaled_advection > 0, backward, forward)
    advection_term = np.sum(scaled_advection * advection_gradient, axis=-1)
    return np.asarray(curvature_term + region_term + advection_term, dtype=np.float32)


def _shift(values: np.ndarray, axis: int, offset: int) -> np.ndarray:
    indices = np.clip(np.arange(values.shape[axis]) + offset, 0, values.shape[axis] - 1)
    return np.take(values, indices, axis=axis)


def _cross_difference(values: np.ndarray, axis_a: int, axis_b: int) -> np.ndarray:
    plus_plus = _shift(_shift(values, axis_a, 1), axis_b, 1)
    minus_minus = _shift(_shift(values, axis_a, -1), axis_b, -1)
    minus_plus = _shift(_shift(values, axis_a, -1), axis_b, 1)
    plus_minus = _shift(_shift(values, axis_a, 1), axis_b, -1)
    return (plus_plus + minus_minus - minus_plus - plus_minus) * 0.25
