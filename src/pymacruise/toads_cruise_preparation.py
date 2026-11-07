"""Prepare the selected second-pass TOADS maps for ACE/CRUISE stages.

The MaCRUISE base route optionally removes brainstem/cerebellum labels,
derives one crop from second-pass GM and applies it to GM/WM-fill/WM-mask,
then topology-corrects cropped WM-fill using cropped WM-mask as the paint.
Each operation delegates to its independently tested selected-stage adapter.
"""

from __future__ import annotations

from dataclasses import dataclass

import nibabel as nib

from .cropping import CropParameters
from .image_io import normalize_image_to_ras
from .toads_cerebellum import remove_brainstem_and_cerebellum
from .toads_tissue_crop import crop_toads_tissue_maps
from .topology import topology_correct


@dataclass(frozen=True)
class ToadsCruisePreparationResult:
    """Crop-ready tissue maps plus the selected topology-corrected WM field."""

    crop_reference_gray_matter: nib.spatialimages.SpatialImage
    gray_matter: nib.spatialimages.SpatialImage
    white_matter_fill: nib.spatialimages.SpatialImage
    white_matter_mask: nib.spatialimages.SpatialImage
    topology_corrected_white_matter_fill: nib.spatialimages.SpatialImage
    crop_parameters: CropParameters
    removed_voxel_counts: tuple[int, int, int] | None
    stage_order: tuple[str, ...]


def prepare_toads_cruise_inputs(
    gray_matter: nib.spatialimages.SpatialImage,
    white_matter_fill: nib.spatialimages.SpatialImage,
    white_matter_mask: nib.spatialimages.SpatialImage,
    *,
    segmentation: nib.spatialimages.SpatialImage | None = None,
    remove_cerebellum: bool = False,
    crop_edge_width: int = 5,
    topology_propagation: str = "object->background",
    topology_source_axis_order: str = "RAS",
) -> ToadsCruisePreparationResult:
    """Compose optional label cleanup, GM crop, and selected topology stage.

    The `MaCRUISEplus`-specific multi-atlas topology refinement is deliberately
    not applied; the selected base wrapper defaults `refineTopoCorrectByMa` to
    false. This helper ends with the WM level set that feeds subsequent ACE.
    """
    if not isinstance(remove_cerebellum, bool):
        raise TypeError("remove_cerebellum must be boolean")
    if remove_cerebellum and segmentation is None:
        raise ValueError("segmentation is required when remove_cerebellum is enabled")

    gm = gray_matter
    wm_fill = white_matter_fill
    wm_mask = white_matter_mask
    stages: list[str] = []
    removed_counts: tuple[int, int, int] | None = None
    if remove_cerebellum:
        assert segmentation is not None
        cleaned = tuple(
            remove_brainstem_and_cerebellum(image, segmentation)
            for image in (gm, wm_fill, wm_mask)
        )
        gm, wm_fill, wm_mask = (result.image for result in cleaned)
        removed_counts = tuple(result.removed_voxels for result in cleaned)
        stages.append("removeBrainstemAndCerebellum")

    cropped = crop_toads_tissue_maps(
        gm,
        wm_fill,
        wm_mask,
        edge_width=crop_edge_width,
    )
    stages.append("cropByPositiveGrayMatter")
    crop_reference = normalize_image_to_ras(gm, source_name="TOADS crop-reference GM")
    topology = topology_correct(
        cropped.white_matter_fill,
        cropped.white_matter_mask,
        paint_threshold=0.5,
        lowest_level=0.0,
        highest_level=1.0,
        regularize=1e-4,
        propagation_direction=topology_propagation,
        source_axis_order=topology_source_axis_order,
    )
    stages.append("topologyCorrectWhiteMatterFill")
    return ToadsCruisePreparationResult(
        crop_reference_gray_matter=crop_reference,
        gray_matter=cropped.gray_matter,
        white_matter_fill=cropped.white_matter_fill,
        white_matter_mask=cropped.white_matter_mask,
        topology_corrected_white_matter_fill=topology.corrected,
        crop_parameters=cropped.crop_parameters,
        removed_voxel_counts=removed_counts,
        stage_order=tuple(stages),
    )
