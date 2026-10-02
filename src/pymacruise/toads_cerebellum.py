"""Selected MaCRUISE brainstem/cerebellum exclusion for TOADS tissue maps.

The MATLAB wrapper transfers the supplied segmentation onto each second-pass
TOADS map with ``ToadsCruiseSpaceToOriginalSpace`` (unchanged resolution,
nearest-neighbor), then zeros labels 35, 38–41, and 71–73. The transfer helper
uses physical-affine nearest-neighbor sampling; exact agreement with JIST's
orientation-based operation is established only for shared physical grids.
"""

from __future__ import annotations

from dataclasses import dataclass

import nibabel as nib
import numpy as np

from .image_io import nifti_ras_on_reference, normalize_image_to_ras
from .spatial_transfer import transfer_to_reference_space


CEREBELLUM_AND_BRAINSTEM_LABELS = frozenset({35, 38, 39, 40, 41, 71, 72, 73})


@dataclass(frozen=True)
class CerebellumRemovalResult:
    """One cleaned map plus the transferred segmentation and changed count."""

    image: nib.spatialimages.SpatialImage
    transferred_segmentation: nib.spatialimages.SpatialImage
    removed_voxels: int


def remove_brainstem_and_cerebellum(
    tissue_map: nib.spatialimages.SpatialImage,
    segmentation: nib.spatialimages.SpatialImage,
) -> CerebellumRemovalResult:
    """Zero selected brainstem/cerebellum voxels in a scalar TOADS map.

    The source applies this operation independently to second-pass GM,
    WM-fill, and WM-mask. Segmentation is nearest-neighbor transferred onto
    the map grid first. The tissue-map dtype and RAS reference geometry are
    retained; `removed_voxels` counts selected segmentation voxels on that
    target grid, regardless of whether their tissue value was already zero.
    """
    for name, image in (("tissue_map", tissue_map), ("segmentation", segmentation)):
        if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
            raise TypeError(f"{name} must be a NIfTI image")
        if len(image.shape) != 3:
            raise ValueError(f"{name} must be a 3D scalar image")

    target = normalize_image_to_ras(tissue_map, source_name="TOADS tissue map")
    values = np.asarray(target.dataobj)
    labels = np.asarray(segmentation.dataobj)
    if not np.issubdtype(values.dtype, np.number) or not np.isfinite(values).all():
        raise ValueError("TOADS tissue map must contain finite numeric voxels")
    if not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("MaSeg cerebellum labels must use an integer data type")

    transferred = transfer_to_reference_space(segmentation, target, cval=0.0)
    transferred_labels = np.asarray(transferred.dataobj)
    remove = np.isin(transferred_labels, tuple(CEREBELLUM_AND_BRAINSTEM_LABELS))
    cleaned = np.array(values, copy=True)
    cleaned[remove] = 0
    result = nifti_ras_on_reference(target, cleaned, dtype=values.dtype)
    return CerebellumRemovalResult(
        image=result,
        transferred_segmentation=transferred,
        removed_voxels=int(np.count_nonzero(remove)),
    )
