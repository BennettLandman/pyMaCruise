"""Source-selected MaSeg and TOADS preparation for the base MaCRUISE route.

The historical V3 base wrapper derives cortical masks from MaSeg, transfers
the second-pass TOADS maps to that reference grid, crops every map from the
MaSeg GM support, and fuses TOADS memberships with confident atlas tissue.
This module exposes that block independently and wires it into the restartable
base workflow so the composition can be validated on synthetic and human inputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import nibabel as nib
import numpy as np

from .boundary import (
    clean_csf_boundary_image,
    combine_multi_atlas_toads_images,
    multi_atlas_csf_boundary_image,
)
from .cortical_separation import cortical_separation
from .cropping import CropParameters, crop_image, find_crop_parameters
from .image_io import normalize_image_to_ras
from .spatial_transfer import transfer_to_reference_space
from .toads_cerebellum import remove_brainstem_and_cerebellum


@dataclass(frozen=True)
class BaseTissuePreparationResult:
    """MaSeg-derived crop, aligned source maps, and fused base tissue maps."""

    atlas_gm: nib.spatialimages.SpatialImage
    atlas_wm: nib.spatialimages.SpatialImage
    atlas_mask: nib.spatialimages.SpatialImage
    atlas_brain_tissue_mask: nib.spatialimages.SpatialImage
    transferred_toads_gm: nib.spatialimages.SpatialImage
    transferred_toads_wm_fill: nib.spatialimages.SpatialImage
    transferred_toads_mask: nib.spatialimages.SpatialImage
    cropped_atlas_gm: nib.spatialimages.SpatialImage
    cropped_atlas_wm: nib.spatialimages.SpatialImage
    cropped_atlas_mask: nib.spatialimages.SpatialImage
    cropped_atlas_brain_tissue_mask: nib.spatialimages.SpatialImage
    cropped_toads_gm: nib.spatialimages.SpatialImage
    cropped_toads_wm_fill: nib.spatialimages.SpatialImage
    cropped_toads_mask: nib.spatialimages.SpatialImage
    raw_boundary_csf: nib.spatialimages.SpatialImage
    toads_cleaned_boundary_csf: nib.spatialimages.SpatialImage
    boundary_csf: nib.spatialimages.SpatialImage
    gray_matter: nib.spatialimages.SpatialImage
    white_matter: nib.spatialimages.SpatialImage
    mask: nib.spatialimages.SpatialImage
    crop_parameters: CropParameters
    stage_order: tuple[str, ...]
    provenance: dict[str, Any]
    removed_voxel_counts: tuple[int, int, int] | None


def prepare_base_tissue_maps(
    maseg: nib.spatialimages.SpatialImage,
    toads_gray_matter: nib.spatialimages.SpatialImage,
    toads_white_matter_fill: nib.spatialimages.SpatialImage,
    toads_mask: nib.spatialimages.SpatialImage,
    *,
    crop_edge_width: int = 5,
    gm_confident_level: float = 0.001,
    remove_cerebellum: bool = False,
    separate_hemispheres: bool = False,
) -> BaseTissuePreparationResult:
    """Prepare source-selected MaSeg/TOADS maps for base CRUISE.

    MaSeg is the reference grid because the original ``maCRUISEpipeline``
    passes its prepared ``Seg.nii.gz`` to ``ToadsCruiseSpaceToOriginalSpace``.
    Transfer is nearest-neighbor with unchanged physical resolution. A shared
    crop is determined by positive MaSeg-derived GM and applied to the three
    cortical masks plus the three transferred TOADS maps. The source's
    ``CombineMAandTOADS`` receives atlas GM, atlas WM, atlas BrainTissueMask,
    and TOADS GM/WM-fill/Mask; its default selected options are 3D erosion and
    GM confidence level 0.001.

    The user-provided skull-stripping mask is intentionally not accepted here:
    it controls TOADS skull stripping, while ``Atlas_BrainTissueMask`` is
    independently derived from MaSeg labels and clips this fusion stage.
    """
    if crop_edge_width < 0:
        raise ValueError("crop_edge_width must be non-negative")
    if not np.isfinite(gm_confident_level) or gm_confident_level < 0:
        raise ValueError("gm_confident_level must be finite and non-negative")
    if not isinstance(remove_cerebellum, bool):
        raise TypeError("remove_cerebellum must be boolean")

    maseg_ras = normalize_image_to_ras(maseg, source_name="base preparation MaSeg")
    atlas = cortical_separation(maseg_ras, separate_hemispheres=separate_hemispheres)
    toads_inputs = (toads_gray_matter, toads_white_matter_fill, toads_mask)
    removed_voxel_counts: tuple[int, int, int] | None = None
    stages = ["corticalSeparationFromMaSeg"]
    if separate_hemispheres:
        stages.append("v3.2LeftRightCorticalBoundaryExclusion")
    if remove_cerebellum:
        label_values = np.asarray(maseg_ras.dataobj)
        if not np.isfinite(label_values).all() or not np.equal(
            label_values, np.floor(label_values)
        ).all():
            raise ValueError("MaSeg cerebellum-removal labels must be finite integers")
        if not np.issubdtype(label_values.dtype, np.integer):
            label_header = maseg_ras.header.copy()
            label_header.set_data_dtype(np.int32)
            maseg_labels = nib.Nifti1Image(
                label_values.astype(np.int32), maseg_ras.affine, header=label_header
            )
        else:
            maseg_labels = maseg_ras
        cleaned = tuple(
            remove_brainstem_and_cerebellum(image, maseg_labels)
            for image in toads_inputs
        )
        toads_inputs = tuple(item.image for item in cleaned)
        removed_voxel_counts = tuple(item.removed_voxels for item in cleaned)
        stages.append("removeBrainstemAndCerebellumBeforeTransfer")
    transferred = tuple(
        transfer_to_reference_space(image, maseg_ras, cval=0.0)
        for image in toads_inputs
    )
    stages.extend((
        "nearestNeighborTransferToMaSegReference",
        "sharedCropFromPositiveMaSegGM",
        "deriveMaSegCSFBoundary",
        "cleanBoundaryBySecondPassToadsWM",
        "combineMaSegAtlasAndTOADS",
        "cleanBoundaryByUpdatedWM",
    ))

    crop_parameters = find_crop_parameters(
        atlas["Atlas_GM"], threshold=0.0, edge_width=crop_edge_width, background=0.0
    )
    atlas_gm = crop_image(atlas["Atlas_GM"], crop_parameters)
    atlas_wm = crop_image(atlas["Atlas_WM"], crop_parameters)
    atlas_mask = crop_image(atlas["Atlas_Mask"], crop_parameters)
    atlas_brain = crop_image(atlas["Atlas_BrainTissueMask"], crop_parameters)
    toads_gm, toads_wm, toads_mask_crop = (
        crop_image(image, crop_parameters) for image in transferred
    )
    cropped_maseg = crop_image(maseg_ras, crop_parameters)

    # The V3 wrapper extracts and TOADS-cleans this boundary before tissue-map
    # update; only the final updated-WM cleanup follows the fusion stage.
    raw_boundary = multi_atlas_csf_boundary_image(cropped_maseg)
    toads_cleaned_boundary = clean_csf_boundary_image(
        raw_boundary, toads_mask_crop, mode="toads_wm"
    )
    fused = combine_multi_atlas_toads_images(
        atlas_gm,
        atlas_wm,
        atlas_brain,
        toads_gm,
        toads_wm,
        toads_mask_crop,
        gm_confident_level=gm_confident_level,
    )
    # V3's ACE mode 2 uses the twice-cleaned MaSeg boundary after diffusion.
    cleaned_boundary = clean_csf_boundary_image(
        toads_cleaned_boundary, fused[1], mode="updated_wm"
    )
    order = tuple(stages)
    provenance: dict[str, Any] = {
        "source_profile": "MaCRUISE V3 base selected wrapper",
        "stage_order": list(order),
        "reference_grid": "MaSeg (prepared Seg.nii.gz)",
        "remove_cerebellum": bool(remove_cerebellum),
        "separate_hemispheres": bool(separate_hemispheres),
        "removed_voxel_counts": (
            None if removed_voxel_counts is None else list(removed_voxel_counts)
        ),
        "transfer": {"interpolation": "nearest_neighbor", "outside_value": 0.0},
        "crop": {
            "support": "Atlas_GM > 0",
            "edge_width_voxels": int(crop_edge_width),
            "bounds_inclusive": [list(bound) for bound in crop_parameters.bounds],
            "original_shape": list(crop_parameters.original_shape),
            "original_affine": crop_parameters.original_affine.tolist(),
        },
        "fusion": {
            "operator": "CombineMAandTOADS",
            "erosion": "3D center plus six neighbors",
            "gm_confident_level": float(gm_confident_level),
        },
        "boundary_csf": {
            "source": "cropped MaSeg labels",
            "overlooked_labels": [120, 121, 202, 203, 156, 157],
            "cleanup_order": ["toads_wm_gt_0.01", "updated_wm_gt_25.5"],
            "ace_updated_csf_mode": 2,
            "merge_location": "after diffusion, before Gaussian smoothing",
        },
        "mask_distinction": {
            "fusion_mask": "MaSeg-derived Atlas_BrainTissueMask",
            "skull_stripping_mask": "external matched skull-stripping mask; not consumed by this stage",
        },
    }
    return BaseTissuePreparationResult(
        atlas_gm=atlas["Atlas_GM"],
        atlas_wm=atlas["Atlas_WM"],
        atlas_mask=atlas["Atlas_Mask"],
        atlas_brain_tissue_mask=atlas["Atlas_BrainTissueMask"],
        transferred_toads_gm=transferred[0],
        transferred_toads_wm_fill=transferred[1],
        transferred_toads_mask=transferred[2],
        cropped_atlas_gm=atlas_gm,
        cropped_atlas_wm=atlas_wm,
        cropped_atlas_mask=atlas_mask,
        cropped_atlas_brain_tissue_mask=atlas_brain,
        cropped_toads_gm=toads_gm,
        cropped_toads_wm_fill=toads_wm,
        cropped_toads_mask=toads_mask_crop,
        raw_boundary_csf=raw_boundary,
        toads_cleaned_boundary_csf=toads_cleaned_boundary,
        boundary_csf=cleaned_boundary,
        gray_matter=fused[0],
        white_matter=fused[1],
        mask=fused[2],
        crop_parameters=crop_parameters,
        stage_order=order,
        provenance=provenance,
        removed_voxel_counts=removed_voxel_counts,
    )
