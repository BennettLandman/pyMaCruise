"""Source-faithful distance fields used by the MaCRUISE dura stage.

The selected `MedicAlgorithmRemoveDura09` implementation delegates signed
distance construction to JIST `DistanceField`. The implementation is shared
with the verified six-neighbor fast-marching port used elsewhere in the
package; Java stores every distance update as float32, so this wrapper selects
that arithmetic explicitly. Distances and caps are in voxel units, matching
the source class (which does not use image voxel spacing).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pymacruise.ace import _signed_distance_fmm
from pymacruise.tgdm import (
    TGDMGeometry,
    TGDMProfileSettings,
    TGDMEvolutionResult,
    tgdm_evolve,
    tgdm_geometry,
)


@dataclass(frozen=True)
class DuraOutputMaskResult:
    """Corrected volumes and mask products from the Dura wrapper's final loop."""

    stripped_volume: np.ndarray
    corrected_white_matter: np.ndarray
    corrected_gray_matter: np.ndarray
    brain_mask_level_set: np.ndarray
    corrected_brain_mask: np.ndarray


def apply_dura_output_mask(
    original_volume: np.ndarray,
    white_matter: np.ndarray,
    gray_matter: np.ndarray,
    evolved_level_set: np.ndarray,
    *,
    mask_translation: float = 2.0,
    level_set_translation: float = 0.0,
) -> DuraOutputMaskResult:
    """Apply the final `MedicAlgorithmRemoveDura09` voxelwise output loop.

    Voxels whose evolved level set is greater than or equal to the mask
    threshold are zeroed from all three corrected volumes. Remaining voxels
    are marked 1 in the binary mask. The level-set output is translated after
    this decision, matching the order in the source wrapper.
    """
    source_arrays = tuple(
        np.asarray(volume)
        for volume in (original_volume, white_matter, gray_matter, evolved_level_set)
    )
    original, wm, gm, level_set = source_arrays
    if level_set.ndim != 3 or any(size < 1 for size in level_set.shape):
        raise ValueError("Dura output fields must be non-empty 3D arrays")
    if any(volume.shape != level_set.shape for volume in source_arrays):
        raise ValueError("Dura output fields must have identical shapes")
    if any(not np.issubdtype(volume.dtype, np.number) for volume in source_arrays):
        raise ValueError("Dura output fields must be numeric")
    if any(not np.isfinite(volume).all() for volume in source_arrays):
        raise ValueError("Dura output fields must be finite")
    controls = np.asarray([mask_translation, level_set_translation], dtype=np.float32)
    if not np.isfinite(controls).all():
        raise ValueError("Dura mask translations must be finite")

    inside_mask = np.asarray(level_set < controls[0], dtype=np.uint8)
    outside = inside_mask == 0
    stripped = np.array(original, copy=True)
    corrected_wm = np.array(wm, copy=True)
    corrected_gm = np.array(gm, copy=True)
    stripped[outside] = 0
    corrected_wm[outside] = 0
    corrected_gm[outside] = 0
    translated_level_set = np.asarray(level_set, dtype=np.float32) - controls[1]
    return DuraOutputMaskResult(
        stripped,
        corrected_wm,
        corrected_gm,
        translated_level_set,
        inside_mask,
    )


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


def dura_tgdm_evolve(
    level_set: np.ndarray,
    gray_matter_volume: np.ndarray,
    reference_level_set: np.ndarray,
    *,
    membership_level: float = 0.4,
    alpha: float = 0.5,
    alignment_force: float = 1.0,
    curvature_force: float = 0.5,
    pressure_force: float = 0.5,
    outer_iterations: int = 5,
    record_passes: bool = False,
) -> TGDMEvolutionResult:
    """Run Dura's selected DynamicPressureTGDM narrow-band evolution.

    This applies the `MedicAlgorithmRemoveDura09` membership normalization
    and sign convention, the VandyCRUISE dynamic-pressure point force,
    GenericTGDM's selected 18/6 scheduler/topology gate, and original-JIST
    float32-rounded distance reinitialization after each outer pass. The rest
    of the removal wrapper's mask cleanup and membership edits is separate.
    """
    phi = np.asarray(level_set, dtype=np.float32)
    gm = np.asarray(gray_matter_volume, dtype=np.float32)
    reference = np.asarray(reference_level_set, dtype=np.float32)
    if phi.ndim != 3 or any(size < 2 for size in phi.shape):
        raise ValueError("dura TGDM fields must be 3D with dimensions >= 2")
    if gm.shape != phi.shape or reference.shape != phi.shape:
        raise ValueError("dura TGDM fields must share identical shapes")
    if not np.isfinite(phi).all() or not np.isfinite(gm).all() or not np.isfinite(reference).all():
        raise ValueError("dura TGDM fields must contain only finite values")
    controls = np.asarray(
        [membership_level, alpha, alignment_force, curvature_force, pressure_force],
        dtype=np.float32,
    )
    if not np.isfinite(controls).all():
        raise ValueError("dura TGDM controls must be finite")
    if membership_level <= 0:
        raise ValueError("membership_level must be positive")
    if outer_iterations < 0:
        raise ValueError("outer_iterations must be nonnegative")

    # MedicAlgorithmRemoveDura09 -> DynamicPressureTGDM wrapper arithmetic.
    membership = np.float32(membership_level)
    denominator = np.float32(membership * np.float32(255.0))
    normalized_region = np.asarray((gm - denominator) / denominator, dtype=np.float32)
    reference_gradient = _fast_marching_unit_gradient_field(reference)
    settings = TGDMProfileSettings(
        max_distance=4.8,
        how_close=3.6,
        max_signed_distance=6.0,
        final_max_signed_distance=5.0,
        internal_iterations=30,
        final_internal_iterations=30,
        delta_t=0.1,
        final_delta_t=0.045,
    )

    def evaluate(
        geometry: TGDMGeometry, roi: tuple[slice, slice, slice]
    ) -> np.ndarray:
        return _dynamic_pressure_force_field(
            geometry,
            normalized_region[roi],
            reference_gradient[roi],
            alpha=float(np.float32(alpha)),
            alignment_force=float(np.float32(alignment_force)),
            curvature_force=-float(np.float32(curvature_force)),
            pressure_force=float(np.float32(pressure_force)),
        )

    return tgdm_evolve(
        phi,
        profile="inner",
        curvature_force=0.0,
        pressure_force=0.0,
        outer_iterations=outer_iterations,
        settings=settings,
        force_evaluator=evaluate,
        round_reinitialization_updates_to_float32=True,
        record_passes=record_passes,
    )


def _fast_marching_unit_gradient_field(level_set: np.ndarray) -> np.ndarray:
    result = np.empty((*level_set.shape, 3), dtype=np.float32)
    for point in np.ndindex(level_set.shape):
        result[point] = _fast_marching_unit_gradient_at(level_set, point).astype(np.float32)
    return result


def _dynamic_pressure_force_field(
    geometry: TGDMGeometry,
    normalized_region: np.ndarray,
    reference_gradient: np.ndarray,
    *,
    alpha: float,
    alignment_force: float,
    curvature_force: float,
    pressure_force: float,
) -> np.ndarray:
    """Vectorized point force, preserving DynamicPressureCortexProfile math."""
    region_force = pressure_force * normalized_region.astype(np.float64)
    backward = geometry.backward
    forward = geometry.forward
    positive_pressure = region_force > 0.0
    region_gradient = np.where(
        positive_pressure[..., None],
        np.maximum(backward, 0.0) + np.minimum(forward, 0.0),
        np.minimum(backward, 0.0) + np.maximum(forward, 0.0),
    )
    region_length = np.linalg.norm(region_gradient, axis=-1)
    dot = np.sum(reference_gradient.astype(np.float64) * region_gradient, axis=-1)
    angular_force = 0.5 - 0.5 * dot / np.where(region_length > 0.0, region_length, 1.0)
    if np.float32(alpha) != np.float32(1.0):
        angular_force = np.abs(angular_force) ** float(np.float32(alpha))
    region_force = region_force - float(np.float32(alignment_force)) * angular_force

    gradient_length = np.linalg.norm(geometry.gradient, axis=-1)
    profile_gradient_length = gradient_length + float(np.float32(0.001))
    normalized_mean = geometry.mean_curvature_numerator / profile_gradient_length**3
    normalized_gaussian = geometry.gaussian_curvature_numerator / profile_gradient_length**4
    discriminant = normalized_mean**2 - 2.0 * normalized_gaussian
    transformed_mean = normalized_mean.copy()
    positive = discriminant > 0.0
    transformed_mean[positive] = (
        normalized_mean[positive] * normalized_gaussian[positive] / discriminant[positive]
    )
    return (
        curvature_force * transformed_mean * profile_gradient_length
        + region_force * region_length
    )
