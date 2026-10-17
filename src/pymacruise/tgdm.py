"""Source-derived local geometry and force terms for VandyCRUISE TGDM.

Published-method context: the nested cortical-surface workflow follows the
CRUISE method of Han et al. (NeuroImage, 2004; DOI
10.1016/j.neuroimage.2004.06.043), whose geometric-deformable-model evolution
is grounded in the topology-preserving level-set method of Han, Xu, and Prince
(IEEE TPAMI, 2003; DOI 10.1109/TPAMI.2003.1201824). The selected
MaCRUISE/VandyCRUISE source and wrapper remain the authority for exact force
profiles, iteration settings, connectivity, visitation order, and adapter
quirks. These citations provide method context; they do not imply that every
step or setting in either paper is reproduced by this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Literal

import numpy as np
from scipy import ndimage
import nibabel as nib

from pymacruise.image_io import normalize_image_to_ras

TGDMProfile = Literal["inner", "central", "outer"]
_FACE = ndimage.generate_binary_structure(3, 1)
_EDGE = ndimage.generate_binary_structure(3, 2)


@dataclass(frozen=True)
class TGDMGeometry:
    """Finite-difference quantities computed by `GenericTGDM.updatePoint`."""

    gradient: np.ndarray
    mean_curvature_numerator: np.ndarray
    gaussian_curvature_numerator: np.ndarray
    backward: np.ndarray
    forward: np.ndarray


@dataclass(frozen=True)
class TGDMEvolutionResult:
    """Level-set output and execution accounting for one TGDM profile."""

    level_set: np.ndarray
    completed_steps: int
    completed_outer_passes: int
    landmine_restarts: int
    pass_level_sets: tuple[np.ndarray, ...] = ()
    pre_reinitialization_level_sets: tuple[np.ndarray, ...] = ()
    point_trace: tuple[tuple[float, ...], ...] = ()
    point_event_trace: tuple["TGDMPointUpdateTrace", ...] = ()


@dataclass(frozen=True)
class TGDMPointUpdateTrace:
    """Opt-in decision details for a selected narrow-band point update."""

    iteration: int
    point: tuple[int, int, int]
    old_value: float
    proposed_value: float
    accepted_value: float
    reason: str
    simple_point: bool | None
    status_neighborhood: tuple[int, ...]
    topology_cube: tuple[int, ...]


@dataclass(frozen=True)
class TGDMProfileSettings:
    """Source profile time and narrow-band constants for a surface stage."""

    max_distance: float
    how_close: float
    max_signed_distance: float
    final_max_signed_distance: float
    internal_iterations: int
    final_internal_iterations: int
    delta_t: float
    final_delta_t: float


@dataclass(frozen=True)
class TGDMBaselineParameters:
    """Effective values passed by the selected MaCRUISE JIST wrapper.

    The MATLAB command wrapper sends the central iteration option to Java's
    outer iteration parameter. Both defaults are seven; retain this mapping
    when callers customize the source wrapper values.
    """

    initial_iso_level: float = 126.9
    white_matter_iso: float = 0.5
    gray_matter_iso: float = 0.4
    initial_curvature: float = 1.0
    initial_iterations: int = 6
    inner_curvature: float = 0.2
    inner_pressure: float = 1.0
    inner_external: float = 0.0
    inner_iterations: int = 4
    central_shrink_iterations: int = 1
    central_curvature: float = 0.15
    central_pressure: float = 1.5
    central_external: float = 1.0
    central_iterations: int = 7
    outer_curvature: float = 1.0
    outer_pressure: float = 1.0
    outer_external: float = 0.0
    outer_iterations: int = 7


@dataclass(frozen=True)
class NestedTGDMResult:
    """Array outputs from the nested CRUISE level-set stages."""

    initial_level_set: np.ndarray
    smoothed_initial: TGDMEvolutionResult
    tissue_sign: np.ndarray
    initial_inner: TGDMEvolutionResult
    central: TGDMEvolutionResult
    central_shrink: TGDMEvolutionResult
    inner: TGDMEvolutionResult
    outer: TGDMEvolutionResult


@dataclass(frozen=True)
class TGDMMaCRUISEImageResult:
    """RAS-oriented NIfTI images from the MaCRUISE TGDM array workflow."""

    initial_level_set: nib.Nifti1Image
    tissue_sign: nib.Nifti1Image
    initial_inner: nib.Nifti1Image
    inner: nib.Nifti1Image
    central: nib.Nifti1Image
    outer: nib.Nifti1Image
    execution: NestedTGDMResult


_PROFILE_SETTINGS = {
    # The Java subclasses pass different final caps, but the selected
    # ProfileTGDM constructor never assigns its finalMaxSignDist argument.
    # Its field initializer therefore remains 20 for all these source calls.
    "inner": TGDMProfileSettings(4.8, 3.6, 6.0, 20.0, 30, 10, 0.1, 0.045),
    "central": TGDMProfileSettings(4.7, 3.6, 6.2, 20.0, 30, 10, 0.1, 0.02),
    "outer": TGDMProfileSettings(4.7, 3.6, 6.2, 20.0, 30, 30, 0.1, 0.045),
}
_SMOOTHING_PROFILE_SETTINGS = TGDMProfileSettings(4.7, 3.6, 6.2, 20.0, 30, 30, 0.1, 0.01)
_INNER_FROM_CENTRAL_PROFILE_SETTINGS = TGDMProfileSettings(
    4.8, 3.6, 6.0, 20.0, 30, 10, 0.1, 0.045
)


def tgdm_evolve(
    level_set: np.ndarray,
    *,
    profile: TGDMProfile,
    curvature_force: float,
    pressure_force: float,
    external_force: float = 0.0,
    region_force: np.ndarray | None = None,
    advection: np.ndarray | None = None,
    inner_bound: np.ndarray | None = None,
    outer_bound: np.ndarray | None = None,
    skeletal_bound: np.ndarray | None = None,
    outer_iterations: int = 7,
    settings: TGDMProfileSettings | None = None,
    force_evaluator: Callable[[TGDMGeometry, tuple[slice, slice, slice]], np.ndarray]
    | None = None,
    round_reinitialization_updates_to_float32: bool = False,
    record_passes: bool = False,
    trace_points: tuple[tuple[int, int, int], ...] = (),
) -> TGDMEvolutionResult:
    """Evolve one source-style TGDM surface from an already signed level set.

    The loop follows GenericTGDM's synchronous narrow-band updates and
    landmine-triggered outer-pass reinitialization. Topology flips use the
    selected `(18,6)` wrapper's observed Java behavior. This low-level API does
    not create the source tissue-force fields or derive initial membership
    distances; those belong to the nested-surface adapters.

    Paper connection: the narrow-band level-set evolution and local topology
    admissibility checks are the TGDM concepts described by Han, Xu, and
    Prince (2003, §3). The exact `CONNECT_18_6` simple-point test and
    reinitialization schedule here follow selected VandyCRUISE `GenericTGDM`
    source behavior, not a parameter prescription from the paper.
    """
    from pymacruise.ace import _tgdm_signed_distance_fmm

    phi = np.asarray(level_set, dtype=np.float32).copy()
    if phi.ndim != 3 or any(size < 2 for size in phi.shape):
        raise ValueError("TGDM evolution requires a 3D level set with dimensions >= 2")
    if not np.isfinite(phi).all():
        raise ValueError("TGDM level set must be finite")
    for point in trace_points:
        if len(point) != 3 or any(index < 0 or index >= phi.shape[axis]
                                   for axis, index in enumerate(point)):
            raise ValueError("trace_points must contain in-bounds (x, y, z) indices")
    if profile not in _PROFILE_SETTINGS:
        raise ValueError("profile must be 'inner', 'central', or 'outer'")
    if outer_iterations < 0:
        raise ValueError("outer_iterations must be non-negative")
    cfg = _PROFILE_SETTINGS[profile] if settings is None else settings
    if (
        cfg.max_distance <= 0
        or cfg.how_close < 0
        or cfg.max_signed_distance <= 0
        or cfg.final_max_signed_distance <= 0
        or cfg.internal_iterations < 0
        or cfg.final_internal_iterations < 0
        or cfg.delta_t < 0
        or cfg.final_delta_t < 0
        or not np.isfinite(
            [
                cfg.max_distance,
                cfg.how_close,
                cfg.max_signed_distance,
                cfg.final_max_signed_distance,
                cfg.delta_t,
                cfg.final_delta_t,
            ]
        ).all()
    ):
        raise ValueError("TGDM profile settings must be finite and non-negative")
    if not np.isfinite([curvature_force, pressure_force, external_force]).all():
        raise ValueError("TGDM force weights must be finite")

    region = _optional_field(region_force, phi.shape, "region_force")
    vector = _optional_vector(advection, phi.shape)
    inner = _optional_field(inner_bound, phi.shape, "inner_bound")
    outer = _optional_field(outer_bound, phi.shape, "outer_bound")
    skeleton = _optional_field(skeletal_bound, phi.shape, "skeletal_bound")
    completed_steps = completed_passes = restarts = 0
    pass_level_sets: list[np.ndarray] = []
    pre_reinitialization_level_sets: list[np.ndarray] = []
    point_trace: list[tuple[float, ...]] = []
    point_event_trace: list[TGDMPointUpdateTrace] = []
    trace_point_set = set(trace_points)

    for pass_index in range(outer_iterations):
        final_pass = pass_index == outer_iterations - 1
        steps = cfg.final_internal_iterations if final_pass else cfg.internal_iterations
        time_step = cfg.final_delta_t if final_pass else cfg.delta_t
        max_signed_distance = (
            cfg.final_max_signed_distance if final_pass else cfg.max_signed_distance
        )
        signs = np.where(phi > 0, -1, 1).astype(np.int8)
        narrow = np.abs(phi) <= cfg.max_distance
        points = np.argwhere(narrow)
        landmine_points = np.argwhere(narrow & (np.abs(phi) > cfg.how_close))
        if landmine_points.size:
            landmine_values = phi[tuple(landmine_points.T)]
            landmine_signs = np.where(landmine_values > 0, -1, 1).astype(np.int8)
        else:
            landmine_signs = np.empty(0, dtype=np.int8)
        mine_hit = False

        for _ in range(steps):
            if mine_hit:
                restarts += 1
                break
            if points.size:
                lower = np.maximum(points.min(axis=0) - 1, 0)
                upper = np.minimum(points.max(axis=0) + 2, phi.shape)
                roi = tuple(slice(int(a), int(b)) for a, b in zip(lower, upper))
                geometry = tgdm_geometry(phi[roi])
                if force_evaluator is None:
                    force = tgdm_force(
                        geometry,
                        profile=profile,
                        curvature_force=curvature_force,
                        pressure_force=pressure_force,
                        external_force=external_force,
                        region_force=None if region is None else region[roi],
                        advection=None if vector is None else vector[roi],
                    )
                else:
                    force = np.asarray(force_evaluator(geometry, roi), dtype=np.float64)
                    expected_shape = tuple(int(part.stop) - int(part.start) for part in roi)
                    if force.shape != expected_shape or not np.isfinite(force).all():
                        raise ValueError("custom TGDM force must be finite and match the update ROI")
            else:
                lower = np.zeros(3, dtype=np.int64)
                force = None
            next_phi = phi.copy()
            # GenericTGDM computes the candidate for every narrow-band voxel,
            # but topology decisions are needed only when its sign changes.
            # Compute the unchanged-sign majority in NumPy, then process sign
            # flips in the same x/y/z order as the Java loop so accepted flips
            # update the topology state seen by later candidates.
            flip_positions = np.empty(0, dtype=np.int64)
            flip_flags = np.empty(0, dtype=np.bool_)
            proposed_values = np.empty(0, dtype=np.float64)
            old_values = np.empty(0, dtype=np.float64)
            trace_positions = np.empty(0, dtype=np.int64)
            if points.size:
                point_coordinates = tuple(points.T)
                old_values = np.asarray(phi[point_coordinates], dtype=np.float64)
                local_points = points - lower
                local_coordinates = tuple(local_points.T)
                proposed_values = old_values - time_step * np.asarray(
                    force[local_coordinates], dtype=np.float64
                )
                flips = ((old_values > 0) & (proposed_values <= 0)) | (
                    (old_values <= 0) & (proposed_values > 0)
                )
                flip_flags = flips
                flip_positions = np.flatnonzero(flips)
                next_phi[point_coordinates] = proposed_values.astype(np.float32)

                if trace_point_set:
                    point_flat = np.ravel_multi_index(points.T, phi.shape)
                    trace_flat = np.asarray(
                        [np.ravel_multi_index(point, phi.shape) for point in trace_point_set],
                        dtype=np.int64,
                    )
                    candidates = np.searchsorted(point_flat, trace_flat)
                    in_bounds = candidates < point_flat.size
                    candidates = candidates[in_bounds]
                    trace_positions = candidates[
                        point_flat[candidates] == trace_flat[in_bounds]
                    ]
            selected_positions = np.union1d(flip_positions, trace_positions)
            # Java visits points in x/y/z array order and stores candidates
            # until the entire narrow-band update has been validated.
            for position in selected_positions:
                point_array = points[position]
                point = tuple(int(v) for v in point_array)
                old_value = float(old_values[position])
                proposed = float(proposed_values[position])
                flips = bool(flip_flags[position])
                event_reason = "no_sign_change"
                event_simple: bool | None = None
                event_neighborhood: tuple[int, ...] = ()
                event_topology_cube: tuple[int, ...] = ()
                if flips:
                    if skeleton is not None and old_value > 0 and skeleton[point] > 0:
                        proposed = 0.1
                        event_reason = "skeleton_bound"
                    elif outer is not None and old_value > 0 and outer[point] >= 0:
                        proposed = 0.1
                        event_reason = "outer_bound"
                    elif inner is not None and old_value <= 0 and inner[point] <= 0:
                        proposed = -0.1
                        event_reason = "inner_bound"
                    else:
                        neighborhood = _status_neighborhood(signs, point)
                        topology_cube = _tgdm_simple_point_cube_618(neighborhood)
                        simple_point = (
                            _source_component_count(topology_cube == 1, _FACE) == 1
                            and _source_component_count(topology_cube == -1, _EDGE) == 1
                        )
                        event_simple = simple_point
                        if point in trace_point_set:
                            event_neighborhood = tuple(int(value) for value in signs[
                                point[0] - 1 : point[0] + 2,
                                point[1] - 1 : point[1] + 2,
                                point[2] - 1 : point[2] + 2,
                            ].ravel())
                            event_topology_cube = tuple(int(value) for value in topology_cube.ravel())
                        if not simple_point:
                            proposed = 0.1 if old_value > 0 else -0.1
                            event_reason = "topology_rejected"
                        else:
                            signs[point] = -signs[point]
                            event_reason = "topology_accepted"
                    next_phi[point] = np.float32(proposed)
                if point in trace_point_set:
                    point_event_trace.append(
                        TGDMPointUpdateTrace(
                            iteration=completed_steps + 1,
                            point=point,
                            old_value=old_value,
                            proposed_value=proposed,
                            accepted_value=float(next_phi[point]),
                            reason=event_reason,
                            simple_point=event_simple,
                            status_neighborhood=event_neighborhood,
                            topology_cube=event_topology_cube,
                        )
                    )

            # GenericTGDM applies bound clamps after all candidate updates.
            next_phi = _apply_level_set_bounds(next_phi, inner, outer)
            phi = next_phi
            completed_steps += 1
            if trace_points:
                point_trace.append(tuple(float(phi[point]) for point in trace_points))
            if landmine_points.size:
                current_values = phi[tuple(landmine_points.T)]
                current_signs = np.where(current_values > 0, -1, 1).astype(np.int8)
                mine_hit = bool(np.any(current_signs != landmine_signs))

        if record_passes:
            pre_reinitialization_level_sets.append(phi.copy())
        phi = _tgdm_signed_distance_fmm(
            phi,
            max_distance=max_signed_distance,
        )
        phi = _apply_level_set_bounds(
            phi, inner, outer, narrow_only=True, band=max_signed_distance
        )
        completed_passes += 1
        if record_passes:
            pass_level_sets.append(phi.copy())

    return TGDMEvolutionResult(
        phi,
        completed_steps,
        completed_passes,
        restarts,
        tuple(pass_level_sets),
        tuple(pre_reinitialization_level_sets),
        tuple(point_trace),
        tuple(point_event_trace),
    )


def tgdm_initial_membership_level_set(
    membership: np.ndarray,
    *,
    iso_level: float = 126.9,
    max_distance: float = 8.0,
) -> np.ndarray:
    """Create the source signed-distance seed from a membership image."""
    from pymacruise.ace import _tgdm_signed_distance_fmm

    values = np.asarray(membership, dtype=np.float32)
    if values.ndim != 3 or any(size < 2 for size in values.shape):
        raise ValueError("membership must be a 3D volume with dimensions >= 2")
    if not np.isfinite(values).all() or not np.isfinite([iso_level, max_distance]).all():
        raise ValueError("membership and distance parameters must be finite")
    if max_distance <= 0:
        raise ValueError("max_distance must be positive")
    return _tgdm_signed_distance_fmm(np.float32(iso_level) - values, max_distance)


def tgdm_tissue_sign_field(
    white_matter: np.ndarray,
    ace_gray_matter: np.ndarray,
    *,
    white_matter_iso: float = 0.1,
    gray_matter_iso: float = 0.4,
) -> np.ndarray:
    """Recreate MedicAlgorithmNestedTGDM's quantized 0..254 tissue sign map."""
    wm = np.asarray(white_matter, dtype=np.float32)
    gm = np.asarray(ace_gray_matter, dtype=np.float32)
    if wm.ndim != 3 or wm.shape != gm.shape:
        raise ValueError("WM and ACE-GM must share one 3D grid")
    if not np.isfinite(wm).all() or not np.isfinite(gm).all():
        raise ValueError("WM and ACE-GM fields must be finite")
    if not 0 < white_matter_iso < 1 or not 0 < gray_matter_iso <= 1:
        raise ValueError("tissue iso levels must satisfy 0 < WM < 1 and 0 < GM <= 1")

    wm_probability = wm / np.float32(255.0)
    gm_probability = gm / np.float32(255.0)
    csf_probability = np.float32(1.0) - wm_probability - gm_probability
    wm_dominant = wm_probability > csf_probability
    sign = np.empty(wm.shape, dtype=np.float32)
    wm_denominator = np.float32(np.float32(1.0) - np.float32(white_matter_iso))
    sign[wm_dominant] = (
        wm_probability[wm_dominant] - np.float32(white_matter_iso)
    ) / wm_denominator
    sign[~wm_dominant] = (
        np.float32(1.0)
        - csf_probability[~wm_dominant]
        - np.float32(gray_matter_iso)
    ) / np.float32(gray_matter_iso)
    sign[wm_dominant] = np.maximum(sign[wm_dominant], 0.0)
    sign[~wm_dominant] = np.minimum(sign[~wm_dominant], 0.0)
    quantized = np.trunc(sign * np.float32(127.0) + np.float32(127.0))
    return np.minimum(quantized, 254).astype(np.uint8)


def tgdm_inner_pressure_field(
    white_matter: np.ndarray, *, white_matter_iso: float = 0.1
) -> np.ndarray:
    """Build the selected inner-surface white-matter pressure map.

    Han et al. (2004, Eq. 9) describe inner-surface pressure around the 0.5
    WM-membership interface. This adapter scales source 0..255 WM values
    against the selected ``white_matter_iso`` (default 0.1); that threshold and
    normalization are wrapper behavior, not the paper's exact formula.
    """
    wm = np.asarray(white_matter, dtype=np.float32)
    if wm.ndim != 3 or not np.isfinite(wm).all():
        raise ValueError("white_matter must be a finite 3D array")
    if not 0 < white_matter_iso <= 1:
        raise ValueError("white_matter_iso must be in (0, 1]")
    denominator = np.float32(white_matter_iso * 255.0)
    return (wm - denominator) / denominator


def tgdm_central_pressure_field(sign_field: np.ndarray) -> np.ndarray:
    """Build the central surface's signed tissue pressure from sign volume.

    Han et al. (2004, Eq. 11) describe a central-surface regional pressure
    that pushes outward from WM, inward from CSF, and is inactive within GM,
    alongside the GM-derived GVF force. This function decodes the selected
    wrapper's quantized 0..254 sign image; its conversion and scaling are
    source-specific and are not a direct implementation of Eq. 11.
    """
    values = np.asarray(sign_field, dtype=np.float32)
    if values.ndim != 3 or not np.isfinite(values).all():
        raise ValueError("sign_field must be a finite 3D array")
    return (values - np.float32(127.0)) / np.float32(127.0)


def tgdm_outer_pressure_field(
    central_level_set: np.ndarray,
    gray_matter: np.ndarray,
    ace_skeleton: np.ndarray,
    *,
    gray_matter_iso: float = 0.4,
) -> np.ndarray:
    """Recreate the selected three-branch GM pressure map and ACE bias.

    Han et al. (2004, Eq. 12) describe outer-surface pressure as
    ``2 * (GM + WM) - 1``. The selected wrapper operation here is different:
    it builds three GM-intensity branches from central-level-set regions and
    assigns a strong inward force on ACE skeleton voxels. Those are
    VandyCRUISE/MaCRUISE source choices, not the paper's Eq. 12 formula.
    """
    central = np.asarray(central_level_set, dtype=np.float32)
    gray = np.asarray(gray_matter, dtype=np.float32)
    skeleton = np.asarray(ace_skeleton)
    if central.ndim != 3 or central.shape != gray.shape or central.shape != skeleton.shape:
        raise ValueError("central level set, GM, and ACE skeleton must share a 3D grid")
    if not np.isfinite(central).all() or not np.isfinite(gray).all():
        raise ValueError("central level set and GM must be finite")
    if not 0 < gray_matter_iso <= 1:
        raise ValueError("gray_matter_iso must be in (0, 1]")
    gray_probability = np.float32(gray_matter_iso * 254.0)
    result = np.empty(central.shape, dtype=np.float32)
    inner_side = central < 0.45
    outer_side = central > 4.5
    middle = ~(inner_side | outer_side)
    result[inner_side] = (gray[inner_side] - 32.0) / 32.0
    result[outer_side] = (gray[outer_side] - 220.0) / 34.0
    result[middle] = (gray[middle] - gray_probability) / gray_probability
    result[middle] = np.where(result[middle] < 0, result[middle] * 1.5, result[middle])
    result = np.clip(result, -1.5, 1.5)
    result[skeleton > 0] = -1.5
    return result


def tgdm_nested_arrays(
    initial_membership: np.ndarray,
    white_matter: np.ndarray,
    gray_matter: np.ndarray,
    ace_gray_matter: np.ndarray,
    *,
    brain_mask: np.ndarray | None = None,
    brain_mask_level_set: np.ndarray | None = None,
    ace_skeleton: np.ndarray | None = None,
    gvf: np.ndarray | None = None,
    parameters: TGDMBaselineParameters = TGDMBaselineParameters(),
) -> NestedTGDMResult:
    """Run the selected nested inner/central/outer CRUISE level-set stages.

    Inputs are arrays already aligned on a RAS-ordered voxel grid. NIfTI
    wrappers that enforce/reorient physical image geometry are a separate
    adapter task; this function never changes grids or resamples inputs.

    Paper connection: CRUISE (Han et al., 2004, “Nested cortical surface
    reconstruction,” Eqs. 9–12) describes a WM-constrained inner surface,
    central evolution from the inner surface using ACE-GM GVF and tissue
    pressure with an inner nesting bound, then an outer surface initialized
    from the central surface with a combined GM/WM pressure and nesting bound.
    This adapter uses selected MaCRUISE/VandyCRUISE tissue fields and force
    profiles: `gvf` is caller-supplied rather than computed here, and the
    wrapper's three-branch outer map is not Eq. 12. Its extra central shrink
    and inner-surface recomputation are source-wrapper behavior, not steps
    attributed to the 2004 paper.
    """
    membership = np.asarray(initial_membership, dtype=np.float32)
    wm = np.asarray(white_matter, dtype=np.float32)
    gm = np.asarray(gray_matter, dtype=np.float32)
    ace_gm = np.asarray(ace_gray_matter, dtype=np.float32)
    if membership.ndim != 3:
        raise ValueError("nested TGDM inputs must be 3D arrays")
    shape = membership.shape
    for name, values in (("white_matter", wm), ("gray_matter", gm), ("ace_gray_matter", ace_gm)):
        if values.shape != shape or not np.isfinite(values).all():
            raise ValueError(f"{name} must be finite and match the initial grid")
    if not np.isfinite(membership).all():
        raise ValueError("initial membership must be finite")
    if brain_mask is not None and brain_mask_level_set is not None:
        raise ValueError("provide brain_mask or brain_mask_level_set, not both")
    mask_values = _optional_field(brain_mask, shape, "brain_mask")
    mask_level_values = _optional_field(
        brain_mask_level_set, shape, "brain_mask_level_set"
    )
    # The public workflow contract supplies a binary brain mask, while
    # GenericTGDM's outer-bound argument is a signed field (inside < 0,
    # outside >= 0). Passing a 0/255 mask directly made every voxel satisfy
    # the positive-bound test and blocked legitimate surface evolution.
    # Only the sign is needed for this containment gate; use unit magnitude
    # to avoid inventing physical-distance units in the array-only adapter.
    # Public workflows normally supply a binary mask; in that case make the
    # inside-negative unit bound used by the adapter. Source wrappers can
    # instead supply their actual signed Brain Mask Level Set. Its magnitude
    # matters because GenericTGDM clamps to the bound value, not just its sign.
    mask = (
        mask_level_values
        if mask_level_values is not None
        else None if mask_values is None else _binary_brain_mask_bound(mask_values)
    )
    skeleton = _optional_field(ace_skeleton, shape, "ace_skeleton")
    velocity = _optional_vector(gvf, shape)

    initial = tgdm_initial_membership_level_set(
        membership, iso_level=parameters.initial_iso_level, max_distance=8.0
    )
    # MaGenericTGDM.solveSmoothSurface runs DistanceField once more before
    # entering its profile loop, even though this wrapper seed is already a
    # signed-distance image.
    from pymacruise.ace import _tgdm_signed_distance_fmm

    smoothing_input = _tgdm_signed_distance_fmm(initial, max_distance=8.0)
    # SmoothCortexProfile passes 5f, but ProfileTGDM never stores that argument.
    smoothed = tgdm_evolve(
        smoothing_input,
        profile="central",
        curvature_force=-parameters.initial_curvature,
        pressure_force=0.0,
        outer_iterations=parameters.initial_iterations,
        settings=_SMOOTHING_PROFILE_SETTINGS,
    )
    sign = tgdm_tissue_sign_field(
        wm,
        ace_gm,
        white_matter_iso=parameters.white_matter_iso,
        gray_matter_iso=parameters.gray_matter_iso,
    )

    inner_force = tgdm_inner_pressure_field(wm, white_matter_iso=parameters.white_matter_iso)
    initial_inner = tgdm_evolve(
        smoothed.level_set,
        profile="inner",
        curvature_force=-parameters.inner_curvature,
        pressure_force=parameters.inner_pressure,
        external_force=parameters.inner_external,
        region_force=inner_force,
        outer_bound=mask,
        outer_iterations=parameters.inner_iterations,
    )

    central_force = tgdm_central_pressure_field(sign)
    central = tgdm_evolve(
        initial_inner.level_set,
        profile="central",
        curvature_force=-parameters.central_curvature,
        pressure_force=parameters.central_pressure,
        external_force=parameters.central_external,
        region_force=central_force,
        advection=velocity,
        inner_bound=initial_inner.level_set,
        outer_bound=mask,
        outer_iterations=parameters.central_iterations,
    )

    # The selected MaTGDMcrop adapter first pushes the central surface inward
    # by one iteration, then recomputes the delivered WM-constrained inner
    # surface from that shrunken central level set.
    # solveInnerSurfaceFromCentral uses InnerCortexProfile; its 1000f argument
    # is likewise ignored by ProfileTGDM, leaving the active value at 20.
    central_shrink = tgdm_evolve(
        central.level_set,
        profile="inner",
        curvature_force=-parameters.central_curvature,
        pressure_force=-1.0,
        outer_iterations=parameters.central_shrink_iterations,
        settings=_INNER_FROM_CENTRAL_PROFILE_SETTINGS,
    )
    inner = tgdm_evolve(
        central_shrink.level_set,
        profile="inner",
        curvature_force=-parameters.inner_curvature,
        pressure_force=-parameters.inner_pressure,
        external_force=parameters.inner_external,
        region_force=-inner_force,
        outer_iterations=parameters.inner_iterations,
    )

    outer_force = tgdm_outer_pressure_field(
        central.level_set,
        gm,
        np.zeros(shape, dtype=np.uint8) if skeleton is None else skeleton,
        gray_matter_iso=parameters.gray_matter_iso,
    )
    outer = tgdm_evolve(
        central.level_set,
        profile="outer",
        curvature_force=-parameters.outer_curvature,
        pressure_force=parameters.outer_pressure,
        external_force=parameters.outer_external,
        region_force=outer_force,
        inner_bound=central.level_set,
        outer_bound=mask,
        skeletal_bound=skeleton,
        # The selected MedicAlgorithmNestedTGDM path passes null GVF here.
        advection=None,
        outer_iterations=parameters.central_iterations,
    )
    return NestedTGDMResult(
        initial,
        smoothed,
        sign,
        initial_inner,
        central,
        central_shrink,
        inner,
        outer,
    )


def tgdm_nested_images(
    initial_membership: nib.spatialimages.SpatialImage,
    white_matter: nib.spatialimages.SpatialImage,
    gray_matter: nib.spatialimages.SpatialImage,
    ace_gray_matter: nib.spatialimages.SpatialImage,
    *,
    brain_mask: nib.spatialimages.SpatialImage | None = None,
    brain_mask_level_set: nib.spatialimages.SpatialImage | None = None,
    ace_skeleton: nib.spatialimages.SpatialImage | None = None,
    gvf: nib.spatialimages.SpatialImage | None = None,
    parameters: TGDMBaselineParameters = TGDMBaselineParameters(),
) -> TGDMMaCRUISEImageResult:
    """Run nested TGDM after RAS reorientation and strict grid validation.

    Paper context: Han et al. (2004), ``docs/papers/cruise.md``, describes
    coupled nested cortical surfaces; the level-set/TGDM method context is
    discussed by Han, Xu, and Prince (2003), ``docs/papers/tgdm.md``. This
    image wrapper normalizes inputs to RAS and delegates evolution to the
    selected source-mapped nested-array adapter. Its effective fields and
    profile settings are MaCRUISE/VandyCRUISE wrapper behavior, not a claim to
    reproduce every stage or parameter in either paper.
    """
    reference, ref_data = _ras_scalar_data(initial_membership, "initial membership")
    aligned: dict[str, np.ndarray] = {}
    for name, image in (
        ("white matter", white_matter),
        ("gray matter", gray_matter),
        ("ACE gray matter", ace_gray_matter),
        ("brain mask", brain_mask),
        ("brain mask level set", brain_mask_level_set),
        ("ACE skeleton", ace_skeleton),
    ):
        if image is None:
            continue
        candidate, data = _ras_scalar_data(image, name)
        _require_ras_same_grid(reference, candidate, name)
        aligned[name] = data

    vector = None
    if gvf is not None:
        vector_image = _require_nifti(gvf, "GVF")
        source_affine = vector_image.affine.copy()
        vector_image = normalize_image_to_ras(vector_image, source_name="GVF")
        _require_ras_same_grid(reference, vector_image, "GVF")
        if len(vector_image.shape) != 4 or vector_image.shape[-1] != 3:
            raise ValueError("GVF must be a component-last 3D vector NIfTI with 3 components")
        vector = _ras_vector_components(
            source_affine, np.asarray(vector_image.dataobj, dtype=np.float32)
        )
        if not np.isfinite(vector).all():
            raise ValueError("GVF vector field must be finite")

    execution = tgdm_nested_arrays(
        ref_data,
        aligned["white matter"],
        aligned["gray matter"],
        aligned["ACE gray matter"],
        brain_mask=aligned.get("brain mask"),
        brain_mask_level_set=aligned.get("brain mask level set"),
        ace_skeleton=aligned.get("ACE skeleton"),
        gvf=vector,
        parameters=parameters,
    )
    return TGDMMaCRUISEImageResult(
        initial_level_set=_nifti_like(execution.initial_level_set, reference),
        tissue_sign=_nifti_like(execution.tissue_sign, reference, np.uint8),
        initial_inner=_nifti_like(execution.initial_inner.level_set, reference),
        inner=_nifti_like(execution.inner.level_set, reference),
        central=_nifti_like(execution.central.level_set, reference),
        outer=_nifti_like(execution.outer.level_set, reference),
        execution=execution,
    )


def _optional_field(field: np.ndarray | None, shape: tuple[int, ...], name: str) -> np.ndarray | None:
    if field is None:
        return None
    values = np.asarray(field, dtype=np.float32)
    if values.shape != shape or not np.isfinite(values).all():
        raise ValueError(f"{name} must be finite and match the TGDM grid")
    return values


def _binary_brain_mask_bound(mask: np.ndarray) -> np.ndarray:
    """Convert a binary mask to TGDM's inside-negative containment field."""
    values = np.asarray(mask, dtype=np.float32)
    if values.ndim != 3 or not np.isfinite(values).all():
        raise ValueError("brain mask must be a finite 3D field")
    # GenericTGDM tests the sign of its bound and clamps same-sign values.
    # Keep only that sign here: array adapters do not carry voxel spacing with
    # which to define a physical signed-distance magnitude.
    return np.where(values > 0, -1.0, 1.0).astype(np.float32)


def _optional_vector(field: np.ndarray | None, shape: tuple[int, ...]) -> np.ndarray | None:
    if field is None:
        return None
    values = np.asarray(field, dtype=np.float32)
    if values.shape != (*shape, 3) or not np.isfinite(values).all():
        raise ValueError("advection must be a finite three-component TGDM field")
    return values


def _level_sign(value: float) -> int:
    return 1 if value > 0 else -1


def _status_neighborhood(signs: np.ndarray, point: tuple[int, int, int]) -> np.ndarray:
    x, y, z = point
    block = np.zeros((3, 3, 3), dtype=np.int8)
    bounds = tuple(
        slice(max(0, p - 1), min(length, p + 2))
        for p, length in zip(point, signs.shape)
    )
    destination = tuple(
        slice(max(0, 1 - p), min(3, length + 1 - p))
        for p, length in zip(point, signs.shape)
    )
    block[destination] = signs[bounds]
    block[1, 1, 1] = 0
    return block


def _apply_level_set_bounds(
    phi: np.ndarray,
    inner: np.ndarray | None,
    outer: np.ndarray | None,
    *,
    narrow_only: bool = False,
    band: float = 0,
) -> np.ndarray:
    result = phi.copy()
    selected = np.ones(phi.shape, dtype=bool) if not narrow_only else np.abs(phi) < band
    if inner is not None:
        same_sign = result * inner >= 0
        clamp = selected & same_sign
        result[clamp] = np.minimum(result[clamp], inner[clamp])
        selected &= same_sign
    if outer is not None:
        same_sign = result * outer >= 0
        clamp = selected & same_sign
        result[clamp] = np.maximum(result[clamp], outer[clamp])
    return result


def tgdm_geometry(level_set: np.ndarray) -> TGDMGeometry:
    """Compute source mirror-boundary first/second derivatives and curvature."""
    # Java reads float32 level-set samples but performs finite differences and
    # curvature algebra in double precision.
    phi = np.asarray(level_set, dtype=np.float64)
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
        np.asarray(gradient, dtype=np.float64),
        np.asarray(mean, dtype=np.float64),
        np.asarray(gaussian, dtype=np.float64),
        np.asarray(backward, dtype=np.float64),
        np.asarray(forward, dtype=np.float64),
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

    # GenericTGDM stores three scalar gradient fields in a Vecmath Vector3d and
    # calls Vector3d.length(), whose source expression is sqrt(x*x + y*y + z*z).
    # Keep that operation order explicit; np.linalg.norm can use a different
    # reduction path, and tiny changes can flip a topology-protected update at
    # a near-zero level-set value.
    grad_magnitude = np.sqrt(
        gradient[..., 0] * gradient[..., 0]
        + gradient[..., 1] * gradient[..., 1]
        + gradient[..., 2] * gradient[..., 2]
    )
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
        grad_squared = grad_magnitude * grad_magnitude
        normalized_mean[nonzero_gradient] = (
            mean[nonzero_gradient]
            / (grad_squared[nonzero_gradient] * grad_magnitude[nonzero_gradient])
        )
        normalized_gaussian[nonzero_gradient] = (
            gaussian[nonzero_gradient] / (grad_squared[nonzero_gradient] * grad_squared[nonzero_gradient])
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
    # Vecmath Vector3d.length() evaluates x*x + y*y + z*z in component order.
    # Preserve that order for the upwind region gradient as well as the level
    # set gradient above; np.linalg.norm may use a different reduction path.
    upwind_magnitude = np.sqrt(
        upwind[..., 0] * upwind[..., 0]
        + upwind[..., 1] * upwind[..., 1]
        + upwind[..., 2] * upwind[..., 2]
    )
    region_term = region_coefficient * upwind_magnitude
    scaled_advection = vector * external_force
    advection_gradient = np.where(scaled_advection > 0, backward, forward)
    advection_term = (
        scaled_advection[..., 0] * advection_gradient[..., 0]
        + scaled_advection[..., 1] * advection_gradient[..., 1]
        + scaled_advection[..., 2] * advection_gradient[..., 2]
    )
    return np.asarray(curvature_term + region_term + advection_term, dtype=np.float64)


def _tgdm_simple_point_cube_618(neighborhood: np.ndarray) -> np.ndarray:
    """Build MaGenericTGDM.SimplePointTest18's geodesic 3x3x3 coding."""
    raw_status = np.asarray(neighborhood)
    if raw_status.dtype == np.bool_:
        status = np.where(raw_status, 1, -1).astype(np.int8)
    else:
        status = np.asarray(raw_status, dtype=np.int8)
    if status.shape != (3, 3, 3):
        raise ValueError("TGDM simple-point neighborhood must be 3x3x3")
    if not np.isin(status, (-1, 0, 1)).all():
        raise ValueError("TGDM status values must be -1, 0, or 1")
    cube = np.zeros((3, 3, 3), dtype=np.int8)
    cube[1, 1, 1] = 0

    # Face neighbors (the six axial positions) are copied directly.
    for point in _offsets(1):
        cube[tuple(v + 1 for v in point)] = status[tuple(v + 1 for v in point)]

    # Edge neighbors are admitted to the geodesic neighborhood only if they
    # are interior or touch an already admitted interior face neighbor.
    for axis in range(3):
        other = [value for value in range(3) if value != axis]
        for first in (-1, 1):
            for second in (-1, 1):
                point = [0, 0, 0]
                point[other[0]], point[other[1]] = first, second
                index = tuple(v + 1 for v in point)
                if status[index] == -1:
                    cube[index] = -1
                # A zero in this padded block denotes an out-of-volume
                # neighbor. MaGenericTGDM.SimplePointTest18 checks bounds
                # before its interior-fill branch, so only explicit +1
                # statuses may enter the geodesic neighborhood.
                elif status[index] == 1 and any(
                    cube[tuple(v + 1 for v in face)] == 1
                    for face in _adjacent_faces(point)
                ):
                    cube[index] = 1

    # Corner voxels are admitted only when interior and touching one of their
    # three edge neighbors, matching the source's ordered geodesic fill.
    for x in (-1, 1):
        for y in (-1, 1):
            for z in (-1, 1):
                point = (x, y, z)
                index = tuple(v + 1 for v in point)
                if status[index] == -1:
                    cube[index] = -1
                # Match the source's explicit bounds check for the corner
                # neighbors too; padded zero cells remain excluded.
                elif status[index] == 1 and any(
                    cube[tuple(v + 1 for v in edge)] == 1
                    for edge in ((x, y, 0), (x, 0, z), (0, y, z))
                ):
                    cube[index] = 1

    return cube


def tgdm_simple_point_618(neighborhood: np.ndarray) -> bool:
    """Apply the selected `(18,6)` source simple-point check to a 3x3x3 field.

    The geodesic cube follows `MaGenericTGDM.SimplePointTest18`; source
    `SixConnected` and `EighteenConnected` component rules then test interior
    and exterior connectivity, respectively.
    """
    cube = _tgdm_simple_point_cube_618(neighborhood)
    interior = cube == 1
    exterior = cube == -1
    interior[1, 1, 1] = False
    exterior[1, 1, 1] = False
    return (
        _source_component_count(interior, _FACE) == 1
        and _source_component_count(exterior, _EDGE) == 1
    )


def _offsets(manhattan_distance: int) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        (x, y, z)
        for x in (-1, 0, 1)
        for y in (-1, 0, 1)
        for z in (-1, 0, 1)
        if abs(x) + abs(y) + abs(z) == manhattan_distance
    )


def _adjacent_faces(point: tuple[int, int, int] | list[int]) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        tuple(0 if axis == changed else point[axis] for axis in range(3))
        for changed in range(3)
        if point[changed] != 0
        and sum(value != 0 for value in point) == 2
    )


def _source_component_count(mask: np.ndarray, structure: np.ndarray) -> int:
    """Count components reached from the center's neighbors as Java does.

    The Java connectivity classes seed their flood fill from 6- or 18-adjacent
    voxels around the center. Components consisting only of unseeded corners
    are therefore ignored, while SciPy's full-volume label count includes them.
    """
    labels, _ = ndimage.label(mask, structure=structure)
    seeds = np.asarray(structure, dtype=bool).copy()
    seeds[tuple(size // 2 for size in seeds.shape)] = False
    neighbor_labels = labels[seeds]
    return int(np.unique(neighbor_labels[neighbor_labels > 0]).size)


def _require_nifti(
    image: nib.spatialimages.SpatialImage, name: str
) -> nib.Nifti1Image | nib.Nifti2Image:
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError(f"{name} must be a NIfTI image")
    return image


def _ras_scalar_data(
    image: nib.spatialimages.SpatialImage, name: str
) -> tuple[nib.Nifti1Image | nib.Nifti2Image, np.ndarray]:
    ras = normalize_image_to_ras(_require_nifti(image, name), source_name=name)
    if len(ras.shape) != 3:
        raise ValueError(f"{name} must be a 3D scalar image")
    data = np.asarray(ras.dataobj, dtype=np.float32)
    if not np.isfinite(data).all():
        raise ValueError(f"{name} must be finite")
    return ras, data


def _require_ras_same_grid(
    reference: nib.spatialimages.SpatialImage,
    candidate: nib.spatialimages.SpatialImage,
    name: str,
) -> None:
    if candidate.shape[:3] != reference.shape[:3] or not np.allclose(
        candidate.affine, reference.affine, rtol=0.0, atol=1e-4
    ):
        raise ValueError(f"{name} must match the initial-membership RAS grid")


def _nifti_like(
    data: np.ndarray,
    reference: nib.Nifti1Image | nib.Nifti2Image,
    dtype: np.dtype | type | None = np.float32,
) -> nib.Nifti1Image:
    values = np.asarray(data, dtype=dtype)
    header = reference.header.copy()
    header.set_data_dtype(values.dtype)
    header.set_slope_inter(None, None)
    output = nib.Nifti1Image(values, reference.affine, header=header)
    qform, qcode = reference.get_qform(coded=True)
    _, scode = reference.get_sform(coded=True)
    if qcode:
        output.set_qform(qform, code=int(qcode), strip_shears=True)
    else:
        output.set_qform(None, code=0)
    output.set_sform(reference.affine, code=int(scode) if scode else 1)
    return output


def _ras_vector_components(source_affine: np.ndarray, vectors: np.ndarray) -> np.ndarray:
    """Convert voxel-axis vector components along with a spatial RAS reorientation."""
    orientation = nib.orientations.ornt_transform(
        nib.orientations.io_orientation(source_affine),
        nib.orientations.axcodes2ornt(("R", "A", "S")),
    )
    converted = np.empty_like(vectors)
    for source_axis, (ras_axis, direction) in enumerate(orientation):
        converted[..., int(ras_axis)] = vectors[..., source_axis] * direction
    return converted


def _shift(values: np.ndarray, axis: int, offset: int) -> np.ndarray:
    indices = np.clip(np.arange(values.shape[axis]) + offset, 0, values.shape[axis] - 1)
    return np.take(values, indices, axis=axis)


def _cross_difference(values: np.ndarray, axis_a: int, axis_b: int) -> np.ndarray:
    plus_plus = _shift(_shift(values, axis_a, 1), axis_b, 1)
    minus_minus = _shift(_shift(values, axis_a, -1), axis_b, -1)
    minus_plus = _shift(_shift(values, axis_a, -1), axis_b, 1)
    plus_minus = _shift(_shift(values, axis_a, 1), axis_b, -1)
    return (plus_plus + minus_minus - minus_plus - plus_minus) * 0.25
