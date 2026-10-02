"""Shared GM-derived crop for selected second-pass TOADS tissue maps.

This composes the source `CropBasedOnGM` wrapper with the repository's
JIST-checked `CubicVolumeCropper` adapter. GM support determines one crop box;
the same inclusive bounds and edge padding are applied to WM-fill and WM-mask.
"""

from __future__ import annotations

from dataclasses import dataclass

import nibabel as nib

from .cropping import CropParameters, crop_image, find_crop_parameters


@dataclass(frozen=True)
class ToadsCroppedTissueMaps:
    """Second-pass GM/WM maps cropped on a shared GM-derived RAS grid."""

    gray_matter: nib.spatialimages.SpatialImage
    white_matter_fill: nib.spatialimages.SpatialImage
    white_matter_mask: nib.spatialimages.SpatialImage
    crop_parameters: CropParameters


def crop_toads_tissue_maps(
    gray_matter: nib.spatialimages.SpatialImage,
    white_matter_fill: nib.spatialimages.SpatialImage,
    white_matter_mask: nib.spatialimages.SpatialImage,
    *,
    edge_width: int = 5,
) -> ToadsCroppedTissueMaps:
    """Crop second-pass GM, WM-fill, and WM-mask with one GM support box."""
    parameters = find_crop_parameters(
        gray_matter,
        threshold=0.0,
        edge_width=edge_width,
        background=0.0,
    )
    return ToadsCroppedTissueMaps(
        gray_matter=crop_image(gray_matter, parameters),
        white_matter_fill=crop_image(white_matter_fill, parameters),
        white_matter_mask=crop_image(white_matter_mask, parameters),
        crop_parameters=parameters,
    )
