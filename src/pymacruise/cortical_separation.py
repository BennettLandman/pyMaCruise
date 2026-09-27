"""MaCRUISE MaSeg label-to-mask preparation."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import normalize_image_to_ras, save_nifti_ras

# Exact label lists from MaCRUISE `pipeline/cortical_seperation.m` (SVN 9360).
GM_LABELS = frozenset(
    [
        *range(100, 110),
        *range(112, 126),
        128,
        129,
        *range(132, 158),
        *range(160, 188),
        *range(190, 208),
    ]
)
WM_LABELS = frozenset([4, 23, 30, 36, 37, 44, 45, 51, 52, *range(55, 63)])
ADDITIONAL_BRAIN_TISSUE_LABELS = frozenset([75, 76, 47, 48, 31, 32, 49, 50])

MASK_LABELS: dict[str, frozenset[int]] = {
    "Atlas_GM": GM_LABELS,
    "Atlas_WM": WM_LABELS,
    # The selected MATLAB source intentionally uses the WM list again here.
    "Atlas_Mask": WM_LABELS,
    "Atlas_BrainTissueMask": GM_LABELS | WM_LABELS | ADDITIONAL_BRAIN_TISSUE_LABELS,
}


def cortical_separation(
    maseg: nib.spatialimages.SpatialImage,
    output_dir: str | Path | None = None,
) -> dict[str, nib.spatialimages.SpatialImage]:
    """Derive exact MaCRUISE cortical masks from an integer-valued MaSeg map.

    Selected label voxels are encoded as uint8 value 255; every other voxel is
    zero. The map is axis-reoriented to RAS with a warning and no interpolation.
    """
    if not isinstance(maseg, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError("MaSeg input must be a NIfTI image")
    if len(maseg.shape) != 3:
        raise ValueError(f"MaSeg must be a 3D single-volume image; got {maseg.shape}")
    normalized = normalize_image_to_ras(maseg, source_name="MaSeg")
    data = np.asanyarray(normalized.dataobj)
    if not np.isfinite(data).all():
        raise ValueError("MaSeg contains NaN or infinite label values")
    if not np.equal(data, np.floor(data)).all():
        raise ValueError("MaSeg labels must have integer-valued voxels")

    header = normalized.header.copy()
    header.set_data_dtype(np.uint8)
    header.set_slope_inter(None, None)
    outputs: dict[str, nib.spatialimages.SpatialImage] = {}
    destination = Path(output_dir) if output_dir is not None else None
    for name, selected_labels in MASK_LABELS.items():
        mask = np.isin(data, tuple(selected_labels)).astype(np.uint8) * np.uint8(255)
        image = nib.Nifti1Image(mask, normalized.affine, header=header.copy())
        qform, qcode = normalized.get_qform(coded=True)
        _, scode = normalized.get_sform(coded=True)
        if qcode:
            image.set_qform(qform, code=int(qcode), strip_shears=True)
        else:
            image.set_qform(None, code=0)
        image.set_sform(normalized.affine, code=int(scode) if scode else 1)
        outputs[name] = image
        if destination is not None:
            save_nifti_ras(image, destination / f"{name}.nii.gz")
    return outputs
