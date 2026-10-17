"""Selected MaCRUISEplus entry-point options, separate from base defaults.

This module records the option values passed by the MATLAB plus entry point.
It does not implement ACE, topology correction, CRUISE, or LesionTOADS stages.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.ace import ACEImageResult, ace_enhance
from pymacruise.boundary import (
    clean_csf_boundary_image,
    combine_multi_atlas_toads_images,
    multi_atlas_csf_boundary_image,
)
from pymacruise.cropping import CropParameters, uncrop_image
from pymacruise.gvf import GVFImageResult, generalized_gvf_image
from pymacruise.tgdm import (
    TGDMBaselineParameters,
    TGDMMaCRUISEImageResult,
    tgdm_nested_images,
)
from pymacruise.surfaces import SurfaceMesh, extract_levelset_surface
from pymacruise.thickness import (
    SurfaceThicknessResult,
    ThicknessVolumeResult,
    closest_point_thickness,
    levelset_subtraction_thickness,
)
from pymacruise.topology import TopologyCorrectionResult, topology_correct


@dataclass(frozen=True)
class PlusPipelineOptions:
    """Configuration overrides selected by ``maCRUISEpluspipeline.m``."""

    lesion_toads_output: str = "dura removal inputs"
    lesion_toads_correct: bool = False
    wm_iso_level: float = 0.5
    gvf_lambda: float = 0.2
    update_enhance: int = 0
    generate_intermediate: bool = True
    use_updated_csf_in_ace: int = 3
    boundary_erosion_method: str = "3D"
    boundary_overlook_labels: tuple[int, ...] = (120, 121, 202, 203, 156, 157)
    refine_topology_by_maseg: bool = True

    def __post_init__(self) -> None:
        if self.lesion_toads_output != "dura removal inputs":
            raise ValueError("selected plus pipeline requires LesionTOADS dura removal inputs")
        if not isinstance(self.lesion_toads_correct, bool):
            raise ValueError("lesion_toads_correct must be boolean")
        if not 0.0 < self.wm_iso_level < 1.0:
            raise ValueError("wm_iso_level must be strictly between zero and one")
        if self.gvf_lambda < 0.0:
            raise ValueError("gvf_lambda must be nonnegative")
        if self.update_enhance not in (0, 1):
            raise ValueError("update_enhance must be 0 or 1")
        if not isinstance(self.generate_intermediate, bool):
            raise ValueError("generate_intermediate must be boolean")
        if self.use_updated_csf_in_ace not in (0, 1, 2, 3):
            raise ValueError("use_updated_csf_in_ace must be one of 0, 1, 2, or 3")
        if self.boundary_erosion_method not in ("3D", "3Dmatlab"):
            raise ValueError("unsupported boundary erosion method")
        if any(not isinstance(label, int) for label in self.boundary_overlook_labels):
            raise ValueError("boundary_overlook_labels must contain integers")
        if not isinstance(self.refine_topology_by_maseg, bool):
            raise ValueError("refine_topology_by_maseg must be boolean")

    def as_dict(self) -> dict[str, Any]:
        """Return JSON-ready Python names/values for run provenance."""
        result = asdict(self)
        result["boundary_overlook_labels"] = list(self.boundary_overlook_labels)
        return result


@dataclass(frozen=True)
class PlusBoundaryResult:
    """Raw and TOADS-WM-cleaned RAS boundary products for MaCRUISEplus."""

    raw_boundary: nib.spatialimages.SpatialImage
    cleaned_boundary: nib.spatialimages.SpatialImage


@dataclass(frozen=True)
class PlusTissueUpdateResult:
    """MaCRUISEplus atlas/TOADS-combined tissue products on the RAS grid."""

    gm: nib.spatialimages.SpatialImage
    wm: nib.spatialimages.SpatialImage
    mask: nib.spatialimages.SpatialImage


@dataclass(frozen=True)
class PlusPreTGDMResult:
    """Selected boundary, tissue, ACE, and GVF products before TGDM."""

    boundary: PlusBoundaryResult
    tissue: PlusTissueUpdateResult
    ace: ACEImageResult
    gvf: GVFImageResult


@dataclass(frozen=True)
class PlusUncroppedTGDMResult:
    """Nested TGDM scalar products restored to the full RAS reference grid."""

    initial_level_set: nib.spatialimages.SpatialImage
    tissue_sign: nib.spatialimages.SpatialImage
    initial_inner: nib.spatialimages.SpatialImage
    inner: nib.spatialimages.SpatialImage
    central: nib.spatialimages.SpatialImage
    outer: nib.spatialimages.SpatialImage


@dataclass(frozen=True)
class PlusSurfaceResult:
    """RAS+ millimetre meshes for the three nested cortical surfaces."""

    inner: SurfaceMesh
    central: SurfaceMesh
    outer: SurfaceMesh


@dataclass(frozen=True)
class PlusThicknessResult:
    """Level-set and closest-point thickness products for plus outputs."""

    volume: ThicknessVolumeResult
    surface: SurfaceThicknessResult


@dataclass(frozen=True)
class PlusPostTGDMResult:
    """Selected plus outputs after cropped TGDM, uncropping, and measurements."""

    level_sets: PlusUncroppedTGDMResult
    surfaces: PlusSurfaceResult
    thickness: PlusThicknessResult


def run_plus_gvf_stage(
    scalar_volume: nib.spatialimages.SpatialImage,
    options: PlusPipelineOptions = PlusPipelineOptions(),
) -> GVFImageResult:
    """Run the existing RAS GVF adapter with selected MaCRUISEplus settings.

    The source wrapper passes ``GVF_Lambda`` and fixes maximum iterations at
    200 with output normalization disabled. Solver implementation remains in
    :mod:`pymacruise.gvf`; this function only maps plus configuration to it.
    """
    if not isinstance(options, PlusPipelineOptions):
        raise ValueError("options must be PlusPipelineOptions")
    return generalized_gvf_image(
        scalar_volume,
        lambda_value=options.gvf_lambda,
        max_iterations=200,
    )


def plus_tgdm_parameters(
    options: PlusPipelineOptions = PlusPipelineOptions(),
) -> TGDMBaselineParameters:
    """Map the plus white-matter isolevel into the existing TGDM settings.

    The MATLAB plus runner overrides `WMIsoLevel` to 0.5 before calling its
    surface pipeline. All other selected TGDM values retain the baseline
    Python defaults; this adapter does not implement surface evolution.
    """
    if not isinstance(options, PlusPipelineOptions):
        raise ValueError("options must be PlusPipelineOptions")
    return TGDMBaselineParameters(
        white_matter_iso=options.wm_iso_level,
    )


def run_plus_ace_stage(
    gray_matter: nib.spatialimages.SpatialImage,
    white_matter: nib.spatialimages.SpatialImage,
    multi_atlas_csf: nib.spatialimages.SpatialImage,
    options: PlusPipelineOptions = PlusPipelineOptions(),
) -> ACEImageResult:
    """Run the selected MaCRUISEplus ACE mode from typed plus settings.

    The selected entry point uses a positive-CSF merge after the speed
    Gaussian (`UseUpdatedCSFinACE=3`) and disables the later enhanced-GM /
    skeleton update branch. The boundary input is supplied on the shared grid.

    Huo et al. (2016), ``docs/papers/ace-maace.md``, describe a distinct
    multi-atlas ``mu_boundary`` speed term. This adapter accepts a precomputed
    atlas-CSF/boundary image and applies the selected Java mode-3 merge; it
    does not compute or validate the paper's ``mu_boundary`` term or its
    max-based speed equation. See ``run_plus_boundary_stage`` for the separate
    selected ``MinFilter3D`` label-boundary operation.
    """
    if not isinstance(options, PlusPipelineOptions):
        raise ValueError("options must be PlusPipelineOptions")
    if options.update_enhance != 0:
        raise ValueError("plus ACE currently supports only update_enhance=0")
    if options.use_updated_csf_in_ace != 3:
        raise ValueError("selected plus ACE requires use_updated_csf_in_ace=3")
    return ace_enhance(
        gray_matter,
        white_matter,
        multi_atlas_csf=multi_atlas_csf,
        use_updated_csf_in_ace=options.use_updated_csf_in_ace,
        merge_method="Merge",
    )


def run_plus_topology_stage(
    volume: nib.spatialimages.SpatialImage,
    paint_mask: nib.spatialimages.SpatialImage,
    ma_wm: nib.spatialimages.SpatialImage | None,
    options: PlusPipelineOptions = PlusPipelineOptions(),
) -> TopologyCorrectionResult:
    """Run selected plus topology correction and optional MaWM refinement."""
    if not isinstance(options, PlusPipelineOptions):
        raise ValueError("options must be PlusPipelineOptions")
    if options.refine_topology_by_maseg and ma_wm is None:
        raise ValueError("plus topology refinement requires the selected MaWM input")
    return topology_correct(
        volume,
        paint_mask,
        multi_atlas_wm=ma_wm if options.refine_topology_by_maseg else None,
    )


def run_plus_boundary_stage(
    segmentation: nib.spatialimages.SpatialImage,
    toads_wm_mask: nib.spatialimages.SpatialImage,
    options: PlusPipelineOptions = PlusPipelineOptions(),
) -> PlusBoundaryResult:
    """Create and clean the selected plus atlas-label boundary.

    Huo et al. (2016), ``docs/papers/ace-maace.md``, describe a multi-atlas
    cortical-boundary term related to MaACE. This function creates the
    selected V3 ``MinFilter3D`` boundary and applies its TOADS-WM cleanup. Its
    relation to the paper's exact ``mu_boundary`` definition is not asserted;
    the neighborhood and cleanup semantics here are source-derived.
    """
    if not isinstance(options, PlusPipelineOptions):
        raise ValueError("options must be PlusPipelineOptions")
    if options.boundary_erosion_method != "3D":
        raise ValueError("selected plus boundary stage supports only ErodeMethod='3D'")
    raw = multi_atlas_csf_boundary_image(
        segmentation,
        overlooked_labels=options.boundary_overlook_labels,
    )
    cleaned = clean_csf_boundary_image(raw, toads_wm_mask, mode="toads_wm")
    return PlusBoundaryResult(raw_boundary=raw, cleaned_boundary=cleaned)


def run_plus_tissue_update_stage(
    multi_atlas_gm: nib.spatialimages.SpatialImage,
    multi_atlas_wm: nib.spatialimages.SpatialImage,
    multi_atlas_brain_mask: nib.spatialimages.SpatialImage,
    toads_gm: nib.spatialimages.SpatialImage,
    toads_wm: nib.spatialimages.SpatialImage,
    toads_mask: nib.spatialimages.SpatialImage,
) -> PlusTissueUpdateResult:
    """Combine cropped atlas and TOADS tissue maps as selected by MaCRUISEplus.

    The selected wrapper calls ``CombineMAandTOADS`` with 3D erosion and
    ``GMconfident=0.001``. The source-checked image helper owns the six-face
    erosion, confident-WM priority, strict GM threshold, fallback, and
    outside-brain behavior; this adapter only routes the six stage products.
    """
    gm, wm, mask = combine_multi_atlas_toads_images(
        multi_atlas_gm,
        multi_atlas_wm,
        multi_atlas_brain_mask,
        toads_gm,
        toads_wm,
        toads_mask,
        gm_confident_level=0.001,
    )
    return PlusTissueUpdateResult(gm=gm, wm=wm, mask=mask)


def run_plus_pre_tgdm_stages(
    segmentation: nib.spatialimages.SpatialImage,
    toads_wm_for_boundary: nib.spatialimages.SpatialImage,
    multi_atlas_gm: nib.spatialimages.SpatialImage,
    multi_atlas_wm: nib.spatialimages.SpatialImage,
    multi_atlas_brain_mask: nib.spatialimages.SpatialImage,
    toads_gm: nib.spatialimages.SpatialImage,
    toads_wm: nib.spatialimages.SpatialImage,
    toads_mask: nib.spatialimages.SpatialImage,
    options: PlusPipelineOptions = PlusPipelineOptions(),
) -> PlusPreTGDMResult:
    """Compose selected MaCRUISEplus boundary -> tissue -> ACE -> GVF stages.

    The selected source passes the fused TOADS/atlas GM and WM into ACE, then
    sends ACE's enhanced gray-matter output to GVF. Atlas maps and TOADS stage
    products are supplied on their selected RAS grid, so this adapter begins
    after TOADS/atlas registration and does not claim those upstream handoffs.
    The topology-corrected initial level set is a separate TGDM input; TGDM
    remains a separate call because its crop provenance is not wired here.
    """
    boundary = run_plus_boundary_stage(
        segmentation, toads_wm_for_boundary, options
    )
    tissue = run_plus_tissue_update_stage(
        multi_atlas_gm,
        multi_atlas_wm,
        multi_atlas_brain_mask,
        toads_gm,
        toads_wm,
        toads_mask,
    )
    ace = run_plus_ace_stage(
        tissue.gm, tissue.wm, boundary.cleaned_boundary, options
    )
    gvf = run_plus_gvf_stage(ace.enhanced_gray_matter, options)
    return PlusPreTGDMResult(boundary=boundary, tissue=tissue, ace=ace, gvf=gvf)


def run_plus_tgdm_stage(
    corrected_level_set: nib.spatialimages.SpatialImage,
    white_matter: nib.spatialimages.SpatialImage,
    gray_matter: nib.spatialimages.SpatialImage,
    ace_gray_matter: nib.spatialimages.SpatialImage,
    ace_thin_skeleton: nib.spatialimages.SpatialImage,
    gvf: nib.spatialimages.SpatialImage,
    options: PlusPipelineOptions = PlusPipelineOptions(),
    *,
    brain_mask: nib.spatialimages.SpatialImage | None = None,
) -> TGDMMaCRUISEImageResult:
    """Run the plus nested-surface stage on its selected cropped RAS products.

    The MaCRUISEplus wrapper supplies corrected level-set, updated tissue,
    MaACE-enhanced GM and thin skeleton, the plus GVF field, and an optional
    brain-mask level set to its TGDM wrapper. This adapter routes those inputs and selected
    `WMIsoLevel` settings; output inverse-cropping belongs to the runner.

    Paper context: Han et al. (2004), ``docs/papers/cruise.md``, describes the
    coupled nested-surface role. MaCRUISEplus supplies its own MaSeg/TOADS
    fields and wrapper settings; this adapter does not reproduce the paper's
    full upstream segmentation pipeline.
    """
    if not isinstance(options, PlusPipelineOptions):
        raise ValueError("options must be PlusPipelineOptions")
    return tgdm_nested_images(
        corrected_level_set,
        white_matter,
        gray_matter,
        ace_gray_matter,
        ace_skeleton=ace_thin_skeleton,
        gvf=gvf,
        brain_mask=brain_mask,
        parameters=plus_tgdm_parameters(options),
    )


def uncrop_plus_tgdm_outputs(
    result: TGDMMaCRUISEImageResult,
    crop_parameters: CropParameters,
) -> PlusUncroppedTGDMResult:
    """Restore nested TGDM scalar volumes using the selected JIST crop rule.

    ``MaTGDMcrop`` uncrops inner, central, outer, and tissue-sign volumes to
    the original GM reference grid. The Python adapter also restores the
    computed initial level set for consistent intermediate provenance. The
    source MaTGDMcrop adapter constructs its inverse-crop parameters with a
    20-valued exterior. Keep that level-set background separate from the
    zero-valued crop background used for tissue and label inputs: surface
    extraction intentionally maps exact zero to -0.1, so uncrop zeros would
    otherwise create false surfaces outside a partial crop. The crop helper
    preserves the selected compiled-JIST exclusive-maximum behavior.
    """
    if not isinstance(result, TGDMMaCRUISEImageResult):
        raise ValueError("result must be TGDMMaCRUISEImageResult")
    if not isinstance(crop_parameters, CropParameters):
        raise ValueError("crop_parameters must be CropParameters")
    if nib.aff2axcodes(crop_parameters.original_affine) != ("R", "A", "S"):
        raise ValueError("plus TGDM crop parameters must describe an RAS reference grid")

    expected_crop_affine = np.asarray(crop_parameters.original_affine, dtype=np.float64).copy()
    translation = np.eye(4)
    translation[:3, 3] = crop_parameters.offset
    expected_crop_affine = expected_crop_affine @ translation
    # MaTGDMcrop.java constructs CropParameters(..., 20) for all four
    # uncropped outputs. Cropping inputs still uses its independent background
    # (normally zero); only restored level-set products need this exterior.
    levelset_uncrop_parameters = replace(crop_parameters, background=20.0)
    restored = {}
    for name in (
        "initial_level_set",
        "tissue_sign",
        "initial_inner",
        "inner",
        "central",
        "outer",
    ):
        image = getattr(result, name)
        if image.shape != crop_parameters.cropped_shape or not np.allclose(
            image.affine, expected_crop_affine, rtol=0.0, atol=1e-4
        ):
            raise ValueError(
                f"{name} must match the cropped RAS grid described by crop_parameters"
            )
        restored[name] = uncrop_image(
            image,
            levelset_uncrop_parameters,
            legacy_exclusive_max=True,
        )
    return PlusUncroppedTGDMResult(**restored)


def extract_plus_surfaces(
    result: PlusUncroppedTGDMResult,
) -> PlusSurfaceResult:
    """Extract nested meshes from full-grid TGDM level sets.

    This uses the Python-native case table generated from selected JIST
    ``IsoSurfaceOnGrid`` source and the complete-affine RAS conversion. The
    case table and direct selected-source comparisons are documented in
    ``docs/reports/jist-connectivity18-case-template-port-2026-10-08.md``.
    """
    if not isinstance(result, PlusUncroppedTGDMResult):
        raise ValueError("result must be PlusUncroppedTGDMResult")
    return PlusSurfaceResult(
        inner=extract_levelset_surface(result.inner, level=0.0),
        central=extract_levelset_surface(result.central, level=0.0),
        outer=extract_levelset_surface(result.outer, level=0.0),
    )


def run_plus_thickness_stage(
    level_sets: PlusUncroppedTGDMResult,
    surfaces: PlusSurfaceResult,
) -> PlusThicknessResult:
    """Run selected MaCRUISEplus volume and closest-point thickness methods."""
    if not isinstance(level_sets, PlusUncroppedTGDMResult):
        raise ValueError("level_sets must be PlusUncroppedTGDMResult")
    if not isinstance(surfaces, PlusSurfaceResult):
        raise ValueError("surfaces must be PlusSurfaceResult")
    volume = levelset_subtraction_thickness(
        level_sets.inner,
        level_sets.outer,
        iso_level=0.0,
    )
    surface = closest_point_thickness(
        surfaces.inner,
        surfaces.outer,
    )
    return PlusThicknessResult(volume=volume, surface=surface)


def run_plus_post_tgdm_stages(
    cropped_tgdm: TGDMMaCRUISEImageResult,
    crop_parameters: CropParameters,
) -> PlusPostTGDMResult:
    """Compose the selected post-TGDM MaCRUISEplus stages in source order.

    This bounded continuation starts with already-produced cropped TGDM
    images. It restores them with the selected JIST-compatible crop rule,
    extracts inner/central/outer meshes from the full-grid RAS level sets, and
    computes the selected level-set and closest-point thickness products.
    It does not run TOADS/LesionTOADS, MaACE/GVF, or full MaCRUISEplus.

    The mesh stage uses the source-derived JIST connectivity-18 case table;
    see :func:`extract_plus_surfaces` and the linked source comparison report.
    """
    level_sets = uncrop_plus_tgdm_outputs(cropped_tgdm, crop_parameters)
    surfaces = extract_plus_surfaces(level_sets)
    thickness = run_plus_thickness_stage(level_sets, surfaces)
    return PlusPostTGDMResult(
        level_sets=level_sets,
        surfaces=surfaces,
        thickness=thickness,
    )
