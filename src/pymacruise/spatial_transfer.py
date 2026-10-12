"""Geometry-aware TOADS-to-reference volume transfer.

Implements the selected adapter behavior of
``ToadsCruiseSpaceToOriginalSpace`` (Yuankai Huo, CVS revision 1.3): nearest
neighbor interpolation, with the supplied reference defining output geometry.
The implementation uses NiBabel's physical affine mapping and outputs RAS.
"""

from __future__ import annotations

import nibabel as nib
import numpy as np
from nibabel.processing import resample_from_to

from pymacruise.image_io import normalize_image_to_ras


def transfer_to_reference_space(
    source: nib.spatialimages.SpatialImage,
    reference: nib.spatialimages.SpatialImage,
    *,
    cval: float = 0.0,
) -> nib.spatialimages.SpatialImage:
    """Resample a source volume to a reference grid using nearest neighbor.

    The legacy adapter uses ``ReorientVolume.solve`` with unchanged source
    resolution and nearest-neighbor interpolation, then assigns the reference
    header. The Python implementation resolves that contract in physical
    coordinates: output shape/affine come from the reference, and source data
    are sampled through both RAS affines without applying an estimated
    registration. Use a separate ANTs registration stage when anatomy is not
    already aligned.
    """
    if not isinstance(source, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError("TOADS source must be a NIfTI image")
    if not isinstance(reference, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError("Reference must be a NIfTI image")
    if len(source.shape) != 3 or len(reference.shape) != 3:
        raise ValueError("TOADS transfer accepts only 3D scalar NIfTI images")
    if not np.isfinite(cval):
        raise ValueError("Outside-grid fill value must be finite")

    source_ras = normalize_image_to_ras(source, source_name="TOADS source")
    reference_ras = normalize_image_to_ras(reference, source_name="TOADS reference")
    source_data = np.asanyarray(source_ras.dataobj)
    if not np.isfinite(source_data).all():
        raise ValueError("TOADS source contains NaN or infinite voxel values")

    sampled = resample_from_to(
        source_ras,
        (reference_ras.shape, reference_ras.affine),
        order=0,
        mode="constant",
        cval=float(cval),
    )
    data = np.asarray(sampled.dataobj).astype(source_ras.get_data_dtype(), copy=False)
    header = reference_ras.header.copy()
    header.set_data_dtype(data.dtype)
    header.set_slope_inter(None, None)
    output = nib.Nifti1Image(data, reference_ras.affine, header=header)
    qform, qcode = reference_ras.get_qform(coded=True)
    _, scode = reference_ras.get_sform(coded=True)
    output.set_qform(qform if qcode else None, code=int(qcode), strip_shears=True)
    output.set_sform(reference_ras.affine, code=int(scode) if scode else 1)
    return output


def assign_reference_geometry_without_resampling(
    source: nib.spatialimages.SpatialImage,
    reference: nib.spatialimages.SpatialImage,
) -> nib.Nifti1Image:
    """Reproduce selected JIST ``UNCHANGED`` extent plus reference-header behavior.

    Source values and extents are retained after RAS reorientation; only the
    reference geometry is assigned. This intentionally does not map samples
    through physical coordinates and is for legacy compatibility comparisons.
    """
    if not isinstance(source, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError("TOADS source must be a NIfTI image")
    if not isinstance(reference, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError("Reference must be a NIfTI image")
    if len(source.shape) != 3 or len(reference.shape) != 3:
        raise ValueError("TOADS transfer accepts only 3D scalar NIfTI images")
    source_ras = normalize_image_to_ras(source, source_name="TOADS source")
    reference_ras = normalize_image_to_ras(reference, source_name="TOADS reference")
    data = np.asarray(source_ras.dataobj)
    if not np.issubdtype(data.dtype, np.number) or not np.isfinite(data).all():
        raise ValueError("TOADS source must contain finite numeric values")
    header = reference_ras.header.copy()
    header.set_data_dtype(data.dtype)
    header.set_slope_inter(None, None)
    output = nib.Nifti1Image(data, reference_ras.affine, header=header)
    qform, qcode = reference_ras.get_qform(coded=True)
    _, scode = reference_ras.get_sform(coded=True)
    output.set_qform(qform if qcode else None, code=int(qcode), strip_shears=True)
    output.set_sform(reference_ras.affine, code=int(scode) if scode else 1)
    return output
