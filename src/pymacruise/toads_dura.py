"""Source-faithful distance fields used by the MaCRUISE dura stage.

The selected `MedicAlgorithmRemoveDura09` implementation delegates signed
distance construction to JIST `DistanceField`. The implementation is shared
with the verified six-neighbor fast-marching port used elsewhere in the
package; Java stores every distance update as float32, so this wrapper selects
that arithmetic explicitly. Distances and caps are in voxel units, matching
the source class (which does not use image voxel spacing).
"""

from __future__ import annotations

import numpy as np

from pymacruise.ace import _signed_distance_fmm
from pymacruise.tgdm import tgdm_geometry


def dura_signed_distance_field(
    level_set: np.ndarray,
    max_distance_voxels: float,
) -> np.ndarray:
    """Return the selected JIST signed-distance field for a 3D level set.

    Inputs are expected in the already RAS-normalized array order used by the
    workflow. The algorithm operates on voxel indices and deliberately does
    not apply physical image spacing. Values are cast to float32, and the
    result is a float32 array with the source class's finite-distance cap.
    """
    return _signed_distance_fmm(
        level_set,
        max_distance=max_distance_voxels,
        round_updates_to_float32=True,
    )


def dynamic_pressure_tgdm_point_update(
    level_set: np.ndarray,
    region_forces: np.ndarray,
    reference_level_set: np.ndarray,
    point: tuple[int, int, int],
    *,
    time_step: float,
    alpha: float,
    alignment_force: float,
    curvature_force: float,
    pressure_profile_force: float,
    curvature_field: np.ndarray | None = None,
    curvature_bias_force: float = 0.0,
) -> float:
    """Update one voxel with VandyCRUISE ``DynamicPressureTGDM`` forces.

    This source-level primitive reproduces the selected subclass' point update
    and ``DynamicPressureCortexProfile`` effective force. It does not run the
    GenericTGDM narrow-band scheduler, topology veto, or reinitialization.
    Level-set and pressure arrays are interpreted in source x/y/z order and
    must already share a grid.
    """
    phi = np.asarray(level_set, dtype=np.float32)
    pressure = np.asarray(region_forces, dtype=np.float32)
    reference = np.asarray(reference_level_set, dtype=np.float32)
    if phi.ndim != 3 or any(size < 2 for size in phi.shape):
        raise ValueError("dura TGDM fields must be 3D with dimensions >= 2")
    if pressure.shape != phi.shape or reference.shape != phi.shape:
        raise ValueError("dura TGDM fields must have identical shapes")
    if not np.isfinite(phi).all() or not np.isfinite(pressure).all() or not np.isfinite(reference).all():
        raise ValueError("dura TGDM fields must contain only finite values")
    if curvature_field is not None:
        curvature = np.asarray(curvature_field, dtype=np.float32)
        if curvature.shape != phi.shape or not np.isfinite(curvature).all():
            raise ValueError("curvature_field must be finite and match the dura TGDM grid")
    else:
        curvature = None

    x, y, z = point
    if len(point) != 3 or any(
        not isinstance(index, (int, np.integer)) or index < 0 or index >= phi.shape[axis]
        for axis, index in enumerate(point)
    ):
        raise ValueError("point must be an in-bounds integer x/y/z coordinate")
    weights = np.asarray(
        [time_step, alpha, alignment_force, curvature_force, pressure_profile_force, curvature_bias_force],
        dtype=np.float64,
    )
    if not np.isfinite(weights).all():
        raise ValueError("dura TGDM force and time-step parameters must be finite")

    # GenericTGDM.updatePoint uses mirror-boundary neighbor indices. Its
    # second-order and centered first-order terms are shared with the ported
    # source geometry calculation.
    geometry = tgdm_geometry(phi)
    grad = geometry.gradient[x, y, z]
    mean = float(geometry.mean_curvature_numerator[x, y, z])
    gaussian = float(geometry.gaussian_curvature_numerator[x, y, z])
    backward = geometry.backward[x, y, z]
    forward = geometry.forward[x, y, z]

    # DynamicPressureTGDM's region gradient uses Godunov upwinding selected by
    # the signed, profile-weighted regional force.
    region_force = float(np.float64(pressure_profile_force) * np.float64(pressure[x, y, z]))
    if region_force > 0:
        region_grad = np.maximum(backward, 0.0) + np.minimum(forward, 0.0)
    else:
        region_grad = np.minimum(backward, 0.0) + np.maximum(forward, 0.0)

    reference_grad = _fast_marching_unit_gradient_at(reference, (x, y, z))
    region_length = float(np.linalg.norm(region_grad))
    angular_force = 0.5 - 0.5 * float(np.dot(reference_grad, region_grad)) / (
        region_length if region_length > 0 else 1.0
    )
    alpha32 = np.float32(alpha)
    if alpha32 != np.float32(1.0):
        angular_force = abs(angular_force) ** float(alpha32)
    region_force -= float(np.float32(alignment_force)) * angular_force
    if curvature is not None:
        region_force += float(curvature[x, y, z]) * float(np.float32(curvature_bias_force))

    grad_length = float(np.linalg.norm(grad))
    # DynamicPressureCortexProfile.getEffectiveForce's float literal is
    # promoted to double, as are the curvature numerators and profile force.
    profile_gradient_length = grad_length + float(np.float32(0.001))
    normalized_mean = mean / (profile_gradient_length * profile_gradient_length * profile_gradient_length)
    normalized_gaussian = gaussian / (profile_gradient_length**4)
    discriminant = normalized_mean * normalized_mean - 2.0 * normalized_gaussian
    transformed_mean = normalized_mean
    if discriminant > 0:
        transformed_mean = normalized_mean * normalized_gaussian / discriminant
    effective_force = (
        float(curvature_force) * transformed_mean * profile_gradient_length
        + region_force * region_length
    )
    return float(phi[x, y, z]) - float(time_step) * effective_force


def _fast_marching_unit_gradient_at(
    level_set: np.ndarray, point: tuple[int, int, int]
) -> np.ndarray:
    """Match ``FastMarchingGradient.gradient(..., MAGNITUDE)`` at one voxel."""
    x, y, z = point
    shape = level_set.shape
    low = (max(x - 1, 0), max(y - 1, 0), max(z - 1, 0))
    high = (min(x + 1, shape[0] - 1), min(y + 1, shape[1] - 1), min(z + 1, shape[2] - 1))
    center = level_set[x, y, z]
    neighbors = (
        (level_set[low[0], y, z], level_set[high[0], y, z]),
        (level_set[x, low[1], z], level_set[x, high[1], z]),
        (level_set[x, y, low[2]], level_set[x, y, high[2]]),
    )
    components = []
    for before, after in neighbors:
        if center > 0:
            component = center - before if before < center and before < after else (
                after - center if after < center else 0.0
            )
        else:
            component = center - before if before > center and before > after else (
                after - center if after > center else 0.0
            )
        components.append(np.float32(component))
    vector = np.asarray(components, dtype=np.float32)
    length = np.float32(np.sqrt(np.sum(vector * vector, dtype=np.float32)))
    if length > 0:
        vector = np.asarray(vector / length, dtype=np.float32)
    return vector.astype(np.float64)
