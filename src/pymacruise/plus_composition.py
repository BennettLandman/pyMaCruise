"""Narrow, explicit composition of the selected MaCRUISEplus stages.

This module joins existing boundary/tissue/ACE/GVF and TGDM adapters when a
caller supplies all crop-ready stage inputs on one compatible RAS grid. It
does not infer or implement upstream atlas registration, TOADS, skull
stripping, topology correction, or crop construction.
"""

from __future__ import annotations

from dataclasses import dataclass

import nibabel as nib
import numpy as np

from pymacruise.image_io import normalize_image_to_ras
from pymacruise.plus_pipeline import (
    PlusPipelineOptions,
    PlusPreTGDMResult,
    run_plus_pre_tgdm_stages,
    run_plus_tgdm_stage,
)
from pymacruise.tgdm import TGDMMaCRUISEImageResult


@dataclass(frozen=True)
class PlusPreTGDMAndTGDMResult:
    """Outputs from the selected pre-TGDM slice followed by nested TGDM."""

    pre_tgdm: PlusPreTGDMResult
    tgdm: TGDMMaCRUISEImageResult


def run_plus_pre_tgdm_and_tgdm(
    segmentation: nib.spatialimages.SpatialImage,
    toads_wm_for_boundary: nib.spatialimages.SpatialImage,
    multi_atlas_gm: nib.spatialimages.SpatialImage,
    multi_atlas_wm: nib.spatialimages.SpatialImage,
    multi_atlas_brain_mask: nib.spatialimages.SpatialImage,
    toads_gm: nib.spatialimages.SpatialImage,
    toads_wm: nib.spatialimages.SpatialImage,
    toads_mask: nib.spatialimages.SpatialImage,
    topology_corrected_initial_membership: nib.spatialimages.SpatialImage,
    options: PlusPipelineOptions = PlusPipelineOptions(),
) -> PlusPreTGDMAndTGDMResult:
    """Compose supplied crop-ready inputs through selected plus TGDM.

    ``topology_corrected_initial_membership`` is an explicit upstream product:
    this function neither computes nor validates its topology. The nested
    TGDM adapter thresholds that membership to construct its initial signed
    distance field. All nine supplied images must already occupy the same
    RAS-ordered cropped grid; crop construction and provenance remain the
    caller's responsibility.

    The order follows the selected plus wrapper: boundary and tissue update,
    ACE, GVF from ACE enhanced GM, then nested TGDM. This is a composition of
    selected Python adapters on supplied stage products, not a claim of
    whole MaCRUISEplus equivalence. In particular, TOADS, multi-atlas
    registration, topology correction, and crop setup are outside this call.
    """
    if not isinstance(options, PlusPipelineOptions):
        raise ValueError("options must be PlusPipelineOptions")
    inputs = (
        segmentation,
        toads_wm_for_boundary,
        multi_atlas_gm,
        multi_atlas_wm,
        multi_atlas_brain_mask,
        toads_gm,
        toads_wm,
        toads_mask,
        topology_corrected_initial_membership,
    )
    names = (
        "segmentation",
        "TOADS WM for boundary",
        "multi-atlas GM",
        "multi-atlas WM",
        "multi-atlas brain mask",
        "TOADS GM",
        "TOADS WM",
        "TOADS mask",
        "topology-corrected initial membership",
    )
    ras_inputs: list[nib.spatialimages.SpatialImage] = []
    for name, image in zip(names, inputs, strict=True):
        if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
            raise ValueError(f"{name} must be a NIfTI image")
        if len(image.shape) != 3:
            raise ValueError(f"{name} must be a 3D scalar image")
        ras_inputs.append(normalize_image_to_ras(image, source_name=name))

    reference: nib.spatialimages.SpatialImage | None = None
    for name, image in zip(names, ras_inputs, strict=True):
        if reference is None:
            reference = image
        elif image.shape != reference.shape or not np.allclose(
            image.affine, reference.affine, rtol=0.0, atol=1e-4
        ):
            raise ValueError(f"{name} must match the supplied crop-ready RAS grid after reorientation")

    pre = run_plus_pre_tgdm_stages(*ras_inputs[:8], options=options)
    nested = run_plus_tgdm_stage(
        ras_inputs[8],
        pre.tissue.wm,
        pre.tissue.gm,
        pre.ace.enhanced_gray_matter,
        pre.ace.thinned_skeleton,
        pre.gvf.vector_field,
        options=options,
    )
    return PlusPreTGDMAndTGDMResult(pre_tgdm=pre, tgdm=nested)
