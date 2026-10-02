"""Compose the selected MaCRUISE base tissue-to-TGDM core.

This follows the base wrapper after second-pass TOADS: prepare/crop/toplevel
WM; run base ACE (without updated multi-atlas CSF); compute generalized GVF
from enhanced GM; run nested TGDM on the crop; and restore its three level sets
to the full GM reference grid with the source cropper's background value.
Visualization surfaces are intentionally omitted.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import nibabel as nib

from .ace import ACEImageResult, ace_enhance
from .cropping import CropParameters, uncrop_image
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
    remove_cerebellum: bool = False,
    crop_edge_width: int = 5,
    ace_diffusion_iterations: int = 25,
    gvf_lambda: float = 1.0,
    gvf_iterations: int = 200,
    gvf_solver_tolerance: float = 1e-6,
    tgdm_parameters: TGDMBaselineParameters | None = None,
) -> ToadsCruiseCoreResult:
    """Run preparation, base ACE/GVF, nested TGDM and selected inverse crop.

    The MaCRUISE base MATLAB route leaves updated multi-atlas CSF disabled in
    `run_Ace_jist`; therefore this composition explicitly selects ACE merge
    placement 0. The TGDM wrapper selects WM iso 0.6 after `run_ToadsCruise`
    assigns `opts.WMIsoLevel=0.6`; other baseline values follow the selected
    wrapper defaults. Native-grid TGDM outputs use the JIST uncrop's 20-valued
    exterior and legacy exclusive high-bound behavior.
    """
    prepared = prepare_toads_cruise_inputs(
        gray_matter,
        white_matter_fill,
        white_matter_mask,
        segmentation=segmentation,
        remove_cerebellum=remove_cerebellum,
        crop_edge_width=crop_edge_width,
    )
    ace = ace_enhance(
        prepared.gray_matter,
        prepared.white_matter_fill,
        intensity_threshold=126.9,
        tensor_alpha=1.6,
        gaussian_alpha=1.8,
        diffusion_iterations=ace_diffusion_iterations,
        use_updated_csf_in_ace=0,
        merge_method="Merge",
    )
    gvf = generalized_gvf_image(
        ace.enhanced_gray_matter,
        lambda_value=gvf_lambda,
        max_iterations=gvf_iterations,
        solver_tolerance=gvf_solver_tolerance,
    )
    baseline = (
        replace(TGDMBaselineParameters(), white_matter_iso=0.6)
        if tgdm_parameters is None
        else tgdm_parameters
    )
    tgdm_crop = tgdm_nested_images(
        prepared.topology_corrected_white_matter_fill,
        prepared.white_matter_fill,
        prepared.gray_matter,
        ace.enhanced_gray_matter,
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
