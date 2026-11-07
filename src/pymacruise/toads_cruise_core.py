"""Compose the selected MaCRUISE base tissue-to-TGDM core.

This follows the base wrapper after second-pass TOADS: prepare/crop/toplevel
WM and its cleaned MaSeg CSF boundary; run base ACE; compute generalized GVF
from enhanced GM; run nested TGDM on the crop; and restore its three level sets
to the full GM reference grid with the source cropper's background value.
Visualization surfaces are intentionally omitted.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import nibabel as nib
import numpy as np

from .ace import ACEImageResult, ace_enhance
from .cropping import CropParameters, crop_image, uncrop_image
from .image_io import normalize_image_to_ras
from .gvf import GVFImageResult, generalized_gvf_image
from .tgdm import (
    TGDMBaselineParameters,
    TGDMMaCRUISEImageResult,
    tgdm_nested_images,
)
from .toads_cruise_preparation import (
    ToadsCruisePreparationResult,
    prepare_toads_cruise_inputs,
)
from .topology import topology_correct


@dataclass(frozen=True)
class ToadsCruiseCoreResult:
    """Intermediate selected ACE/GVF/TGDM products and native-grid level sets."""

    preparation: ToadsCruisePreparationResult
    ace: ACEImageResult
    gvf: GVFImageResult
    tgdm_crop: TGDMMaCRUISEImageResult
    inner_native: nib.spatialimages.SpatialImage
    central_native: nib.spatialimages.SpatialImage
    outer_native: nib.spatialimages.SpatialImage


def run_toads_cruise_core(
    gray_matter: nib.spatialimages.SpatialImage,
    white_matter_fill: nib.spatialimages.SpatialImage,
    white_matter_mask: nib.spatialimages.SpatialImage,
    *,
    segmentation: nib.spatialimages.SpatialImage | None = None,
    boundary_csf: nib.spatialimages.SpatialImage | None = None,
    brain_mask: nib.spatialimages.SpatialImage | None = None,
    brain_mask_level_set: nib.spatialimages.SpatialImage | None = None,
    remove_cerebellum: bool = False,
    crop_edge_width: int = 5,
    ace_diffusion_iterations: int = 25,
    ace_use_updated_csf_in_ace: int = 0,
    wm_iso_level: float = 0.6,
    gvf_lambda: float = 1.0,
    gvf_iterations: int = 200,
    gvf_solver_tolerance: float = 1e-6,
    tgdm_parameters: TGDMBaselineParameters | None = None,
    prepared_crop_parameters: CropParameters | None = None,
    topology_propagation: str = "object->background",
    topology_source_axis_order: str = "RAS",
) -> ToadsCruiseCoreResult:
    """Run preparation, base ACE/GVF, nested TGDM and selected inverse crop.

    ACE updated-CSF placement and TGDM WM iso are explicit parameters because
    the selected V3 MATLAB pipeline overrides the older wrapper defaults:
    it uses ACE placement mode 2 and WM iso 0.5. The generic function defaults
    remain available for compatibility callers; the selected base workflow
    passes the V3 values explicitly. Native-grid TGDM outputs use the JIST
    uncrop's 20-valued exterior and legacy exclusive high-bound behavior.
    """
    if ace_use_updated_csf_in_ace == 2 and boundary_csf is None:
        raise ValueError("ACE mode 2 requires the cleaned MaSeg boundary_csf image")
    if brain_mask is not None and brain_mask_level_set is not None:
        raise ValueError("provide brain_mask or brain_mask_level_set, not both")
    if prepared_crop_parameters is None:
        prepared = prepare_toads_cruise_inputs(
            gray_matter,
            white_matter_fill,
            white_matter_mask,
            segmentation=segmentation,
            remove_cerebellum=remove_cerebellum,
            crop_edge_width=crop_edge_width,
            topology_propagation=topology_propagation,
            topology_source_axis_order=topology_source_axis_order,
        )
    else:
        if remove_cerebellum or segmentation is not None:
            raise ValueError(
                "pre-cropped fused maps cannot repeat cerebellum removal; use the base preparation stage"
            )
        gm_ras, wm_ras, mask_ras = (
            normalize_image_to_ras(image, source_name=name)
            for image, name in (
                (gray_matter, "fused base GM"),
                (white_matter_fill, "fused base WM-fill"),
                (white_matter_mask, "fused base WM-mask"),
            )
        )
        if any(
            image.shape != gm_ras.shape
            or not np.allclose(image.affine, gm_ras.affine, rtol=0, atol=1e-5)
            for image in (wm_ras, mask_ras)
        ):
            raise ValueError("pre-cropped fused base maps must share one RAS grid")
        if gm_ras.shape != prepared_crop_parameters.cropped_shape:
            raise ValueError("fused base map shape differs from the recorded MaSeg crop")
        translation = np.eye(4)
        translation[:3, 3] = prepared_crop_parameters.offset
        cropped_affine = prepared_crop_parameters.original_affine @ translation
        if not np.allclose(gm_ras.affine, cropped_affine, rtol=0, atol=1e-5):
            raise ValueError("fused base map affine differs from the recorded MaSeg crop")
        corrected = topology_correct(
            wm_ras,
            mask_ras,
            paint_threshold=0.5,
            lowest_level=0.0,
            highest_level=1.0,
            regularize=1e-4,
            propagation_direction=topology_propagation,
            source_axis_order=topology_source_axis_order,
        )
        prepared = ToadsCruisePreparationResult(
            crop_reference_gray_matter=gm_ras,
            gray_matter=gm_ras,
            white_matter_fill=wm_ras,
            white_matter_mask=mask_ras,
            topology_corrected_white_matter_fill=corrected.corrected,
            crop_parameters=prepared_crop_parameters,
            removed_voxel_counts=None,
            stage_order=("useFusedMaSegToadsMaps", "topologyCorrectWhiteMatterFill"),
        )
    ace_boundary = None
    if boundary_csf is not None:
        boundary_ras = normalize_image_to_ras(
            boundary_csf, source_name="base ACE multi-atlas CSF boundary"
        )
        same_cropped_grid = (
            boundary_ras.shape == prepared.gray_matter.shape
            and np.allclose(boundary_ras.affine, prepared.gray_matter.affine, rtol=0, atol=1e-5)
        )
        if same_cropped_grid:
            ace_boundary = boundary_ras
        elif (
            boundary_ras.shape == prepared.crop_parameters.original_shape
            and np.allclose(
                boundary_ras.affine,
                prepared.crop_parameters.original_affine,
                rtol=0,
                atol=1e-5,
            )
        ):
            ace_boundary = crop_image(boundary_ras, prepared.crop_parameters)
        else:
            raise ValueError(
                "ACE multi-atlas CSF boundary must match the prepared crop or its full reference grid"
            )
    ace = ace_enhance(
        prepared.gray_matter,
        prepared.white_matter_fill,
        multi_atlas_csf=ace_boundary,
        intensity_threshold=126.9,
        tensor_alpha=1.6,
        gaussian_alpha=1.8,
        diffusion_iterations=ace_diffusion_iterations,
        use_updated_csf_in_ace=ace_use_updated_csf_in_ace,
        merge_method="Merge",
    )
    gvf = generalized_gvf_image(
        ace.enhanced_gray_matter,
        lambda_value=gvf_lambda,
        max_iterations=gvf_iterations,
        solver_tolerance=gvf_solver_tolerance,
    )
    baseline = (
        replace(TGDMBaselineParameters(), white_matter_iso=wm_iso_level)
        if tgdm_parameters is None
        else tgdm_parameters
    )
    tgdm_brain_mask = None
    tgdm_brain_mask_level_set = None
    mask_input = brain_mask_level_set if brain_mask_level_set is not None else brain_mask
    if mask_input is not None:
        mask_name = (
            "base brain mask level set"
            if brain_mask_level_set is not None
            else "base brain mask"
        )
        mask_ras = normalize_image_to_ras(mask_input, source_name=mask_name)
        crop = prepared.crop_parameters
        crop_affine = crop.original_affine.copy()
        crop_affine[:3, 3] += crop.original_affine[:3, :3] @ np.asarray(crop.offset)
        same_full_grid = (
            mask_ras.shape == crop.original_shape
            and np.allclose(mask_ras.affine, crop.original_affine, rtol=0, atol=1e-5)
        )
        same_cropped_grid = (
            mask_ras.shape == crop.cropped_shape
            and np.allclose(mask_ras.affine, crop_affine, rtol=0, atol=1e-5)
        )
        if same_full_grid:
            mask_ras = crop_image(mask_ras, crop)
        elif not same_cropped_grid:
            raise ValueError(
                f"{mask_name} must match the source T1 grid or the prepared TGDM crop"
            )
        if brain_mask_level_set is not None:
            tgdm_brain_mask_level_set = mask_ras
        else:
            tgdm_brain_mask = mask_ras
    tgdm_crop = tgdm_nested_images(
        prepared.topology_corrected_white_matter_fill,
        prepared.white_matter_fill,
        prepared.gray_matter,
        ace.enhanced_gray_matter,
        brain_mask=tgdm_brain_mask,
        brain_mask_level_set=tgdm_brain_mask_level_set,
        ace_skeleton=ace.thinned_skeleton,
        gvf=gvf.vector_field,
        parameters=baseline,
    )
    native_crop = replace(prepared.crop_parameters, background=20.0)
    inner_native = uncrop_image(tgdm_crop.inner, native_crop)
    central_native = uncrop_image(tgdm_crop.central, native_crop)
    outer_native = uncrop_image(tgdm_crop.outer, native_crop)
    return ToadsCruiseCoreResult(
        preparation=prepared,
        ace=ace,
        gvf=gvf,
        tgdm_crop=tgdm_crop,
        inner_native=inner_native,
        central_native=central_native,
        outer_native=outer_native,
    )
