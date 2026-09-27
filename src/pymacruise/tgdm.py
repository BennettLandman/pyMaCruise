"""Source-derived local geometry and force terms for VandyCRUISE TGDM."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

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
    "inner": TGDMProfileSettings(4.8, 3.6, 6.0, 20.0, 30, 10, 0.1, 0.045),
    "central": TGDMProfileSettings(4.7, 3.6, 6.2, 20.0, 30, 10, 0.1, 0.02),
    "outer": TGDMProfileSettings(4.7, 3.6, 6.2, 20.0, 30, 30, 0.1, 0.045),
}


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
) -> TGDMEvolutionResult:
    """Evolve one source-style TGDM surface from an already signed level set.

    The loop follows GenericTGDM's synchronous narrow-band updates and
    landmine-triggered outer-pass reinitialization. Topology flips use the
    selected `(18,6)` wrapper's observed Java behavior. This low-level API does
    not create the source tissue-force fields or derive initial membership
    distances; those belong to the nested-surface adapters.
    """
    from pymacruise.ace import _signed_distance_fmm

    phi = np.asarray(level_set, dtype=np.float32).copy()
    if phi.ndim != 3 or any(size < 2 for size in phi.shape):
        raise ValueError("TGDM evolution requires a 3D level set with dimensions >= 2")
    if not np.isfinite(phi).all():
        raise ValueError("TGDM level set must be finite")
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
        landmine_signs = np.asarray([_level_sign(phi[tuple(p)]) for p in landmine_points])
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
                lower = np.zeros(3, dtype=np.int64)
                force = None
            next_phi = phi.copy()
            # Java visits points in x/y/z array order and stores candidates
            # until the entire narrow-band update has been validated.
            for point_array in points:
                point = tuple(int(v) for v in point_array)
                old_value = float(phi[point])
                local_point = tuple(int(v - start) for v, start in zip(point, lower))
                proposed = float(old_value - time_step * float(force[local_point]))
                flips = (old_value > 0 and proposed <= 0) or (old_value <= 0 and proposed > 0)
                if flips:
                    if skeleton is not None and old_value > 0 and skeleton[point] > 0:
                        proposed = 0.1
                    elif outer is not None and old_value > 0 and outer[point] >= 0:
                        proposed = 0.1
                    elif inner is not None and old_value <= 0 and inner[point] <= 0:
                        proposed = -0.1
                    elif not tgdm_simple_point_618(_status_neighborhood(signs, point)):
                        proposed = 0.1 if old_value > 0 else -0.1
                    else:
                        signs[point] = -signs[point]
                next_phi[point] = np.float32(proposed)

            # GenericTGDM applies bound clamps after all candidate updates.
            next_phi = _apply_level_set_bounds(next_phi, inner, outer)
            phi = next_phi
            completed_steps += 1
            if landmine_points.size:
                current_signs = np.asarray([_level_sign(phi[tuple(p)]) for p in landmine_points])
                mine_hit = bool(np.any(current_signs != landmine_signs))

        phi = _signed_distance_fmm(phi, max_distance=max_signed_distance)
        phi = _apply_level_set_bounds(
            phi, inner, outer, narrow_only=True, band=max_signed_distance
        )
        completed_passes += 1

    return TGDMEvolutionResult(phi, completed_steps, completed_passes, restarts)


def tgdm_initial_membership_level_set(
    membership: np.ndarray,
    *,
    iso_level: float = 126.9,
    max_distance: float = 8.0,
) -> np.ndarray:
    """Create the source signed-distance seed from a membership image."""
    from pymacruise.ace import _signed_distance_fmm

    values = np.asarray(membership, dtype=np.float32)
    if values.ndim != 3 or any(size < 2 for size in values.shape):
        raise ValueError("membership must be a 3D volume with dimensions >= 2")
    if not np.isfinite(values).all() or not np.isfinite([iso_level, max_distance]).all():
        raise ValueError("membership and distance parameters must be finite")
    if max_distance <= 0:
        raise ValueError("max_distance must be positive")
    return _signed_distance_fmm(np.float32(iso_level) - values, max_distance)


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
    sign[wm_dominant] = (
        wm_probability[wm_dominant] - np.float32(white_matter_iso)
    ) / np.float32(1.0 - white_matter_iso)
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
    """Build `solveInnerSurfaceFromOuter`'s normalized white-matter force."""
    wm = np.asarray(white_matter, dtype=np.float32)
    if wm.ndim != 3 or not np.isfinite(wm).all():
        raise ValueError("white_matter must be a finite 3D array")
    if not 0 < white_matter_iso <= 1:
        raise ValueError("white_matter_iso must be in (0, 1]")
    denominator = np.float32(white_matter_iso * 255.0)
    return (wm - denominator) / denominator


def tgdm_central_pressure_field(sign_field: np.ndarray) -> np.ndarray:
    """Build the central surface's signed tissue pressure from sign volume."""
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
    """Recreate the three-branch GM pressure map and ACE skeleton bias."""
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
    ace_skeleton: np.ndarray | None = None,
    gvf: np.ndarray | None = None,
    parameters: TGDMBaselineParameters = TGDMBaselineParameters(),
) -> NestedTGDMResult:
    """Run the selected nested inner/central/outer CRUISE level-set stages.

    Inputs are arrays already aligned on a RAS-ordered voxel grid. NIfTI
    wrappers that enforce/reorient physical image geometry are a separate
    adapter task; this function never changes grids or resamples inputs.
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
    mask = _optional_field(brain_mask, shape, "brain_mask")
    skeleton = _optional_field(ace_skeleton, shape, "ace_skeleton")
    velocity = _optional_vector(gvf, shape)

    initial = tgdm_initial_membership_level_set(
        membership, iso_level=parameters.initial_iso_level, max_distance=8.0
    )
    # MaGenericTGDM.solveSmoothSurface runs DistanceField once more before
    # entering its profile loop, even though this wrapper seed is already a
    # signed-distance image.
    from pymacruise.ace import _signed_distance_fmm

    smoothing_input = _signed_distance_fmm(initial, max_distance=8.0)
    smoothing_settings = TGDMProfileSettings(4.7, 3.6, 6.2, 20.0, 30, 30, 0.1, 0.01)
    smoothed = tgdm_evolve(
        smoothing_input,
        profile="central",
        curvature_force=-parameters.initial_curvature,
        pressure_force=0.0,
        outer_iterations=parameters.initial_iterations,
        settings=smoothing_settings,
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
    central_shrink_settings = TGDMProfileSettings(4.8, 3.6, 6.0, 20.0, 30, 10, 0.1, 0.045)
    central_shrink = tgdm_evolve(
        central.level_set,
        profile="inner",
        curvature_force=-parameters.central_curvature,
        pressure_force=-1.0,
        outer_iterations=parameters.central_shrink_iterations,
        settings=central_shrink_settings,
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
    ace_skeleton: nib.spatialimages.SpatialImage | None = None,
    gvf: nib.spatialimages.SpatialImage | None = None,
    parameters: TGDMBaselineParameters = TGDMBaselineParameters(),
) -> TGDMMaCRUISEImageResult:
    """Run nested TGDM after RAS reorientation and strict grid validation."""
    reference, ref_data = _ras_scalar_data(initial_membership, "initial membership")
    aligned: dict[str, np.ndarray] = {}
    for name, image in (
        ("white matter", white_matter),
        ("gray matter", gray_matter),
        ("ACE gray matter", ace_gray_matter),
        ("brain mask", brain_mask),
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
    return np.asarray(curvature_term + region_term + advection_term, dtype=np.float64)


def tgdm_simple_point_618(neighborhood: np.ndarray) -> bool:
    """Apply GenericTGDM's geodesic 18-neighborhood simple-point check.

    ``neighborhood`` is a 3x3x3 boolean array: True is the source status +1
    (the interior for the selected `(18,6)` wrapper setting). This reproduces
    the source `SimplePointTest18` geodesic coding followed by its
    `SixConnected` + `EighteenConnected` component checks. The apparent
    constructor/connectivity naming mismatch is documented in the TGDM source
    map and intentionally represented here as source behavior.
    """
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
                elif any(
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
                elif any(
                    cube[tuple(v + 1 for v in edge)] == 1
                    for edge in ((x, y, 0), (x, 0, z), (0, y, z))
                ):
                    cube[index] = 1

    interior = cube == 1
    exterior = cube == -1
    interior[1, 1, 1] = False
    exterior[1, 1, 1] = False
    return (
        _component_count(interior, _FACE) == 1
        and _component_count(exterior, _EDGE) == 1
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


def _component_count(mask: np.ndarray, structure: np.ndarray) -> int:
    return int(ndimage.label(mask, structure=structure)[1])


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
