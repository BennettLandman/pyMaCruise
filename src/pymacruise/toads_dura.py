"""Source-faithful operations for the MaCRUISE second-pass dura-removal stage.

Landman et al. (2013), §“Dura Removal” (printed p. 95; DOI
10.1007/s12021-012-9159-9), is the primary published description identified
for the specific “Remove Dura
2009” operation. It calls the stage a previously unpublished contribution and
describes a second pass after tissue segmentation: initialize from the
topologically corrected WM membership, evolve a topology-preserving geometric
deformable model (TGDM) between WM and WM+GM memberships, add a cortical-
thickness regularization force to prune dura attachments using fast-marching
topology correction, then derive a brain mask and apply it to the original
image. The paper provides the method-level account, not the selected JIST
revision’s complete parameter and arithmetic contract.

The implementation authority for source-specific details is the selected
VandyCRUISE ``MedicAlgorithmRemoveDura09`` and its called JIST classes. The
article's software listing names “Remove Dura 2009 (1.6 RC)”; this is a
procedure-level mapping and does not establish that the selected source
revision reproduces every detail of that release. Han, Xu, and Prince
(2003; DOI 10.1109/TPAMI.2003.1201824) describe the underlying topology-
preserving TGDM method; Bazin and Pham (2007; DOI
10.1016/j.cmpb.2007.08.006) describe the fast-marching topology-correction
method. These component papers do not, by themselves, specify the complete
dura-removal procedure. The Python port composes selected-source topology,
distance-field, and TGDM operations and does not claim independent clinical
validation or end-to-end equivalence. See ``docs/papers/dura-removal.md`` for
the operation-to-code crosswalk and source/validation limits.

The signed-distance adapter delegates the algorithm to JIST ``DistanceField``.
Its six-neighbor fast-marching implementation is shared with the separately
validated package primitive; Java stores each distance update as float32, so
this wrapper selects that arithmetic explicitly. Distances and caps are in
voxel units, matching the source class (which does not use image voxel
spacing).
"""

from __future__ import annotations

from dataclasses import dataclass

import nibabel as nib
import numpy as np

from pymacruise.ace import _signed_distance_fmm
from pymacruise.image_io import normalize_image_to_ras
from pymacruise.tgdm import (
    TGDMGeometry,
    TGDMProfileSettings,
    TGDMEvolutionResult,
    tgdm_evolve,
    tgdm_geometry,
)
from pymacruise.topology import topology_correct
from pymacruise.volume_ops import calculate_volume


@dataclass(frozen=True)
class DuraOutputMaskResult:
    """Corrected volumes and mask products from the Dura wrapper's final loop."""

    stripped_volume: np.ndarray
    corrected_white_matter: np.ndarray
    corrected_gray_matter: np.ndarray
    brain_mask_level_set: np.ndarray
    corrected_brain_mask: np.ndarray


@dataclass(frozen=True)
class DuraInitialLevelSets:
    """Membership level sets and capped distance fields used before Dura TGDM."""

    inner_level_set: np.ndarray
    outer_level_set: np.ndarray
    inner_signed_distance: np.ndarray
    outer_signed_distance: np.ndarray


@dataclass(frozen=True)
class DuraRemovalResult:
    """RAS NIfTI products and stage state from a composed Dura removal run."""

    stripped_t1: nib.Nifti1Image
    corrected_white_matter: nib.Nifti1Image
    corrected_gray_matter: nib.Nifti1Image
    brain_mask_level_set: nib.Nifti1Image
    brain_mask: nib.Nifti1Image
    wm_topology: nib.Nifti1Image
    gm_topology: nib.Nifti1Image
    evolution: TGDMEvolutionResult


def run_dura_removal(
    original_t1: nib.spatialimages.SpatialImage,
    white_matter_membership: nib.spatialimages.SpatialImage,
    gray_matter_membership: nib.spatialimages.SpatialImage,
    original_white_matter_mask: nib.spatialimages.SpatialImage,
    *,
    membership_level: float = 0.4,
    curvature_force: float = 0.5,
    alignment_force: float = 1.0,
    pressure_force: float = 0.5,
    alignment_alpha: float = 0.5,
    max_iterations: int = 5,
    mask_translation: float = 2.0,
    level_set_translation: float = 0.0,
    record_passes: bool = False,
) -> DuraRemovalResult:
    """Compose the selected Dura wrapper on a shared native T1 RAS grid.

    This wires `MedicAlgorithmRemoveDura09.execute()` in source order:
    correct WM topology, form WM-corrected + raw-GM and correct GM topology,
    build capped inner/outer signed-distance fields, evolve the pial mask
    level set, then threshold it before zeroing the three returned images.
    Inputs are axis-reoriented to RAS with warnings and must share one physical
    grid; outputs retain that grid. Topology, distance, and TGDM numerics are
    supplied by the separately validated stage implementations.
    """
    named_inputs = (
        ("original T1", original_t1),
        ("white-matter membership", white_matter_membership),
        ("gray-matter membership", gray_matter_membership),
        ("original white-matter mask", original_white_matter_mask),
    )
    normalized: list[nib.spatialimages.SpatialImage] = []
    for name, image in named_inputs:
        if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
            raise ValueError(f"{name} must be a NIfTI image")
        ras = normalize_image_to_ras(image, source_name=name)
        data = np.asanyarray(ras.dataobj)
        if data.ndim != 3 or any(size < 3 for size in data.shape):
            raise ValueError(f"{name} must be a 3D image with dimensions >= 3")
        if not np.issubdtype(data.dtype, np.number) or not np.isfinite(data).all():
            raise ValueError(f"{name} must contain finite numeric voxels")
        normalized.append(ras)
    t1, raw_wm, raw_gm, wm_mask = normalized
    reference_shape = t1.shape
    reference_affine = t1.affine
    for name, image in zip(
        ("white-matter membership", "gray-matter membership", "original white-matter mask"),
        normalized[1:],
    ):
        if image.shape != reference_shape or not np.allclose(
            image.affine, reference_affine, rtol=0, atol=1e-4
        ):
            raise ValueError(f"{name} must match the original T1 RAS grid")

    wm_tc = dura_correct_white_matter_topology(raw_wm, wm_mask)
    gm_tc = dura_correct_gray_matter_topology(wm_tc, raw_gm, wm_mask)
    wm_values = np.asarray(wm_tc.dataobj, dtype=np.float32)
    gm_values = np.asarray(gm_tc.dataobj, dtype=np.float32)
    levels = dura_initial_level_sets(
        wm_values,
        gm_values,
        membership_level=membership_level,
    )
    evolution = dura_tgdm_evolve(
        levels.outer_signed_distance,
        gm_values,
        levels.inner_signed_distance,
        membership_level=membership_level,
        alpha=alignment_alpha,
        alignment_force=alignment_force,
        curvature_force=curvature_force,
        pressure_force=pressure_force,
        outer_iterations=max_iterations,
        record_passes=record_passes,
    )
    final = apply_dura_output_mask(
        np.asanyarray(t1.dataobj),
        wm_values,
        gm_values,
        evolution.level_set,
        mask_translation=mask_translation,
        level_set_translation=level_set_translation,
    )
    def output_image(values: np.ndarray, name: str) -> nib.Nifti1Image:
        header = t1.header.copy()
        header.set_data_dtype(values.dtype)
        image = nib.Nifti1Image(values, reference_affine, header=header)
        _, qcode = t1.get_qform(coded=True)
        _, scode = t1.get_sform(coded=True)
        if qcode:
            image.set_qform(reference_affine, code=int(qcode), strip_shears=True)
        else:
            image.set_qform(None, code=0)
        image.set_sform(reference_affine, code=int(scode) if scode else 1)
        image.set_filename(f"{name}.nii.gz")
        return image

    return DuraRemovalResult(
        stripped_t1=output_image(final.stripped_volume, "t1_durastripped"),
        corrected_white_matter=output_image(final.corrected_white_matter, "wm_fixed"),
        corrected_gray_matter=output_image(final.corrected_gray_matter, "gm_fixed"),
        brain_mask_level_set=output_image(final.brain_mask_level_set, "brain_mask_levelset"),
        brain_mask=output_image(final.corrected_brain_mask, "brain_mask"),
        wm_topology=wm_tc,
        gm_topology=gm_tc,
        evolution=evolution,
    )


def dura_correct_white_matter_topology(
    white_matter: nib.spatialimages.SpatialImage,
    original_white_matter_mask: nib.spatialimages.SpatialImage,
) -> nib.Nifti1Image:
    """Apply the first Dura topology adapter call in the current Python port.

    The selected Java wrapper passes its source WM membership field and the
    original WM mask to ``MedicAlgorithmTopologyCorrection`` with foreground/
    background connectivities 6/18, object-to-background propagation,
    normalization enabled over levels 0–1, paint threshold 0.5, and regularize
    amount 0.05. The shared Python topology implementation handles RAS
    normalization and grid checks; this adapter supplies the Dura-specific
    parameter that differs from its general default.

    The returned image is the first corrected field (``wmTC``). This function
    deliberately does not construct ``wmTC + raw GM`` or perform the separate
    second ``tcGM`` call. Five selected-Java synthetic fixtures match
    voxelwise exactly; see
    ``docs/reports/topology-correction-java-comparison-2026-09-29.md``.
    """
    return topology_correct(
        white_matter,
        original_white_matter_mask,
        paint_threshold=0.5,
        lowest_level=0.0,
        highest_level=1.0,
        regularize=0.05,
    ).corrected


def dura_correct_gray_matter_topology(
    topologically_corrected_white_matter: nib.spatialimages.SpatialImage,
    raw_gray_matter: nib.spatialimages.SpatialImage,
    original_white_matter_mask: nib.spatialimages.SpatialImage,
) -> nib.Nifti1Image:
    """Build the Dura WM+raw-GM union and run the second topology adapter.

    ``MedicAlgorithmRemoveDura09.execute()`` adds the already corrected WM
    field to the original, uncorrected GM field in float32, then passes that
    union and the same original WM mask to a second topology adapter. Generic
    NIfTI addition and physical-grid checks reuse ``calculate_volume``; the
    Dura-specific 6/18, object-to-background settings and regularization are
    applied by ``topology_correct``. The returned image is ``gmTC``. The
    selected algorithm-pair call sequence and float32 union now match one
    deterministic selected-Java harness fixture exactly; the harness follows
    the wrapper sequence but does not instantiate the full JIST wrapper. This
    slice does not evolve TGDM or assemble the complete dura-removal workflow.
    """
    union = calculate_volume(
        topologically_corrected_white_matter,
        raw_gray_matter,
        "ADD",
    )
    return topology_correct(
        union,
        original_white_matter_mask,
        paint_threshold=0.5,
        lowest_level=0.0,
        highest_level=1.0,
        regularize=0.05,
    ).corrected


def dura_initial_level_sets(
    topologically_corrected_white_matter: np.ndarray,
    topologically_corrected_gray_matter: np.ndarray,
    *,
    membership_level: float = 0.4,
) -> DuraInitialLevelSets:
    """Build the source wrapper's inner and outer Dura TGDM initial fields."""
    wm = np.asarray(topologically_corrected_white_matter, dtype=np.float32)
    gm = np.asarray(topologically_corrected_gray_matter, dtype=np.float32)
    if wm.ndim != 3 or any(size < 1 for size in wm.shape):
        raise ValueError("Dura membership fields must be non-empty 3D arrays")
    if gm.shape != wm.shape:
        raise ValueError("Dura membership fields must have identical shapes")
    if not np.isfinite(wm).all() or not np.isfinite(gm).all():
        raise ValueError("Dura membership fields must be finite")
    membership = np.float32(membership_level)
    if not np.isfinite(membership) or not 0.0 <= membership <= 1.0:
        raise ValueError("membership_level must be finite and in [0, 1]")

    inner = np.float32(0.5) - wm / np.float32(255.0)
    outer = membership - gm / np.float32(255.0)
    return DuraInitialLevelSets(
        inner,
        outer,
        dura_signed_distance_field(inner, max_distance_voxels=15.0),
        dura_signed_distance_field(outer, max_distance_voxels=10.0),
    )


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
