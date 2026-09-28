"""Shared RAS NIfTI to ANTsPy image conversion."""

from __future__ import annotations

import nibabel as nib
import numpy as np


class ANTsDependencyError(RuntimeError):
    """Raised when the optional ANTsPy dependency cannot be loaded."""


def nifti_to_ants(
    image: nib.spatialimages.SpatialImage,
    *,
    dtype: np.dtype | str = np.float32,
):
    """Convert a RAS-oriented NIfTI image to ANTsPy's physical LPS convention."""
    try:
        import ants
    except (ImportError, OSError) as exc:
        raise ANTsDependencyError(
            "ANTsPy is required for registration; install PyMaCRUISE with its 'ants' extra"
        ) from exc

    affine = np.asarray(image.affine, dtype=np.float64)
    spacing = nib.affines.voxel_sizes(affine).astype(np.float64)
    ras_to_lps = np.diag([-1.0, -1.0, 1.0])
    direction_lps = ras_to_lps @ (affine[:3, :3] / spacing[np.newaxis, :])
    origin_lps = ras_to_lps @ affine[:3, 3]
    values = np.asarray(image.dataobj, dtype=np.dtype(dtype))
    ants_image = ants.from_numpy(
        values,
        origin=tuple(float(value) for value in origin_lps),
        spacing=tuple(float(value) for value in spacing),
        direction=direction_lps,
    )
    return ants, ants_image
