"""Source-ordered MaCRUISEplus composition with topology correction.

This bounded composition starts from products already prepared on a common
cropped RAS grid. It joins the selected boundary and tissue update, plus
topology correction/refinement, ACE, GVF, and TGDM adapters. It does not build
the crop or run upstream registration, TOADS, or LesionTOADS.
"""

from __future__ import annotations

from dataclasses import dataclass

import nibabel as nib
import numpy as np

from pymacruise.ace import ACEImageResult
from pymacruise.image_io import normalize_image_to_ras
from pymacruise.gvf import GVFImageResult
from pymacruise.plus_pipeline import (
    PlusBoundaryResult,
    PlusPipelineOptions,
    PlusTissueUpdateResult,
    run_plus_ace_stage,
    run_plus_boundary_stage,
    run_plus_gvf_stage,
    run_plus_tgdm_stage,
    run_plus_tissue_update_stage,
    run_plus_topology_stage,
)
from pymacruise.tgdm import TGDMMaCRUISEImageResult
from pymacruise.topology import TopologyCorrectionResult


@dataclass(frozen=True)
class PlusTopologyTGDMResult:
    """Intermediate products from boundary through the selected TGDM stage."""

    boundary: PlusBoundaryResult
    tissue: PlusTissueUpdateResult
    topology: TopologyCorrectionResult
    tgdm_initial_membership: nib.spatialimages.SpatialImage
    ace: ACEImageResult
    gvf: GVFImageResult
    tgdm: TGDMMaCRUISEImageResult


def run_plus_topology_and_tgdm(
    segmentation: nib.spatialimages.SpatialImage,
    toads_wm_for_boundary: nib.spatialimages.SpatialImage,
    multi_atlas_gm: nib.spatialimages.SpatialImage,
    multi_atlas_wm: nib.spatialimages.SpatialImage,
    multi_atlas_brain_mask: nib.spatialimages.SpatialImage,
    toads_gm: nib.spatialimages.SpatialImage,
    toads_wm: nib.spatialimages.SpatialImage,
    toads_mask: nib.spatialimages.SpatialImage,
    options: PlusPipelineOptions = PlusPipelineOptions(),
) -> PlusTopologyTGDMResult:
    """Run the selected crop-ready MaCRUISEplus stages through nested TGDM.

    The eight inputs are already cropped and registered onto one grid.
    ``multi_atlas_wm`` is the raw cropped MaSeg-derived WM map: it supplies
    both atlas tissue fusion and the separate topology-refinement argument.
    Stage order follows ``maCRUISEpluspipeline.m`` and
    ``run_CruisePlusSurf_pipeline.m``: boundary/tissue preparation, topology
    correction of fused WM using the fused mask and raw atlas WM, ACE, GVF,
    then TGDM. When atlas refinement is enabled, the source wrapper captures
    the first ``run_TopoCorrect`` return value, which is the MaWM-refined
    product; that is the initial membership passed here to TGDM.

    This is not full-workflow equivalence. It does not perform cropping,
    atlas/TOADS registration, TOADS/LesionTOADS, or upstream image preparation.
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
    )
    names = (
        "segmentation",
        "TOADS WM for boundary cleanup",
        "multi-atlas GM",
        "raw multi-atlas WM",
        "multi-atlas brain mask",
        "TOADS GM",
        "TOADS WM",
        "TOADS mask",
    )
    ras_inputs: list[nib.spatialimages.SpatialImage] = []
    for name, image in zip(names, inputs, strict=True):
        if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
            raise ValueError(f"{name} must be a NIfTI image")
        if len(image.shape) != 3:
            raise ValueError(f"{name} must be a 3D scalar image")
        ras_inputs.append(normalize_image_to_ras(image, source_name=name))

    reference = ras_inputs[0]
    for name, image in zip(names[1:], ras_inputs[1:], strict=True):
        if image.shape != reference.shape or not np.allclose(
            image.affine, reference.affine, rtol=0.0, atol=1e-4
        ):
            raise ValueError(f"{name} must match the crop-ready RAS grid")

    boundary = run_plus_boundary_stage(ras_inputs[0], ras_inputs[1], options)
    tissue = run_plus_tissue_update_stage(*ras_inputs[2:])
    topology = run_plus_topology_stage(
        tissue.wm,
        tissue.mask,
        ras_inputs[3] if options.refine_topology_by_maseg else None,
        options,
    )
    initial_membership = (
        topology.atlas_refined
        if options.refine_topology_by_maseg
        else topology.corrected
    )
    if initial_membership is None:
        raise RuntimeError("selected topology stage did not produce its TGDM seed")

    ace = run_plus_ace_stage(tissue.gm, tissue.wm, boundary.cleaned_boundary, options)
    gvf = run_plus_gvf_stage(ace.enhanced_gray_matter, options)
    tgdm = run_plus_tgdm_stage(
        initial_membership,
        tissue.wm,
        tissue.gm,
        ace.enhanced_gray_matter,
        ace.thinned_skeleton,
        gvf.vector_field,
        options,
    )
    return PlusTopologyTGDMResult(
        boundary=boundary,
        tissue=tissue,
        topology=topology,
        tgdm_initial_membership=initial_membership,
        ace=ace,
        gvf=gvf,
        tgdm=tgdm,
    )
