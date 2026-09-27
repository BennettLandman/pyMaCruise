"""MaCRUISE multi-atlas CSF boundary and boundary-merge operations.

Ports the selected V3 workflow's default ``MinFilter3D`` path and the used
``merge_csf_method`` values (1 and 3). Thin-skeleton generation is a separate
CRUISE skeleton operation and is intentionally handled by P27.
"""

from __future__ import annotations

from typing import Literal, Sequence

import nibabel as nib
import numpy as np
from scipy import ndimage

from pymacruise.image_io import normalize_image_to_ras

DEFAULT_OVERLOOKED_LABELS = (120, 121, 202, 203, 156, 157)
MergeMode = Literal["boundary_as_one", "speed_function"]
CleanupMode = Literal["toads_wm", "updated_wm"]


def multi_atlas_csf_boundary(
    segmentation: np.ndarray,
    *,
    overlooked_labels: Sequence[int] = DEFAULT_OVERLOOKED_LABELS,
) -> np.ndarray:
    """Return the exact selected ``MinFilter3D`` boundary as uint8 0/1.

    The source removes overlooked labels, computes a 3-D erosion using the
    center plus six face-neighbors, then accepts only candidates that touch a
    different nonzero label in the 26-neighborhood. Any neighboring label in
    ``1..99`` suppresses that boundary voxel. Volume-edge voxels are omitted,
    as in the Java implementation's interior-only loops.
    """
    values = np.asarray(segmentation)
    if values.ndim != 3:
        raise ValueError(f"Segmentation must be 3D; got {values.shape}")
    if not np.issubdtype(values.dtype, np.number) or not np.isfinite(values).all():
        raise ValueError("Segmentation must contain finite numeric labels")
    if not np.equal(values, np.floor(values)).all():
        raise ValueError("Segmentation labels must be integer-valued")

    source = values.astype(np.int64, copy=True)
    if overlooked_labels:
        source[np.isin(source, np.asarray(overlooked_labels, dtype=np.int64))] = 0

    result = np.zeros(source.shape, dtype=np.uint8)
    if any(length < 3 for length in source.shape):
        return result

    six_neighbor = ndimage.generate_binary_structure(3, 1)
    interior = tuple(slice(1, -1) for _ in range(3))
    candidate = np.zeros(source.shape, dtype=bool)
    for label in np.unique(source):
        if label < 100:
            continue
        binary = source == label
        eroded = ndimage.binary_erosion(binary, structure=six_neighbor, border_value=0)
        candidate[interior] |= (binary[interior] != eroded[interior])

    low_contact_any = np.zeros(source.shape, dtype=bool)
    distinct_contact_any = np.zeros(source.shape, dtype=bool)
    center = source[interior]
    for offset in np.ndindex((3, 3, 3)):
        delta = tuple(component - 1 for component in offset)
        if delta == (0, 0, 0):
            continue
        neighbor = source[
            1 + delta[0] : source.shape[0] - 1 + delta[0],
            1 + delta[1] : source.shape[1] - 1 + delta[1],
            1 + delta[2] : source.shape[2] - 1 + delta[2],
        ]
        low_contact_any[interior] |= (neighbor != 0) & (neighbor < 100)
        distinct_contact_any[interior] |= (neighbor != 0) & (neighbor != center)
    result[interior] = (
        candidate[interior]
        & distinct_contact_any[interior]
        & ~low_contact_any[interior]
    ).astype(np.uint8)
    return result


def merge_csf_boundary(
    boundary: np.ndarray,
    csf: np.ndarray,
    *,
    mode: MergeMode,
) -> np.ndarray:
    """Merge a boundary with a CSF speed volume using selected V3 rules.

    ``boundary_as_one`` implements source mode 1: ``min(boundary + csf, 1)``.
    ``speed_function`` implements source mode 3: set every positive-boundary
    voxel to 0.1 in the CSF volume. Mode 2 is unused by the selected runner.
    """
    edge = np.asarray(boundary)
    values = np.asarray(csf)
    if edge.shape != values.shape:
        raise ValueError(f"Boundary and CSF shapes differ: {edge.shape} vs {values.shape}")
    if edge.ndim != 3:
        raise ValueError(f"Boundary and CSF volumes must be 3D; got {edge.shape}")
    if not np.isfinite(edge).all() or not np.isfinite(values).all():
        raise ValueError("Boundary and CSF volumes must contain finite values")
    if mode == "boundary_as_one":
        merged = edge.astype(np.float64) + values.astype(np.float64)
        merged[merged > 1] = 1
        return merged
    if mode == "speed_function":
        merged = values.astype(np.float64, copy=True)
        merged[edge > 0] = 0.1
        return merged
    raise ValueError(f"Unsupported selected CSF merge mode: {mode!r}")


def multi_atlas_csf_boundary_image(
    segmentation: nib.spatialimages.SpatialImage,
    *,
    overlooked_labels: Sequence[int] = DEFAULT_OVERLOOKED_LABELS,
) -> nib.spatialimages.SpatialImage:
    """Create the RAS-oriented NIfTI boundary image from a MaSeg image."""
    normalized = normalize_image_to_ras(segmentation, source_name="MaSeg")
    boundary = multi_atlas_csf_boundary(
        np.asanyarray(normalized.dataobj), overlooked_labels=overlooked_labels
    )
    header = normalized.header.copy()
    header.set_data_dtype(np.uint8)
    header.set_slope_inter(None, None)
    output = nib.Nifti1Image(boundary, normalized.affine, header=header)
    _, qcode = normalized.get_qform(coded=True)
    _, scode = normalized.get_sform(coded=True)
    output.set_qform(normalized.affine if qcode else None, code=int(qcode))
    output.set_sform(normalized.affine, code=int(scode) if scode else 1)
    return output


def merge_csf_boundary_images(
    boundary: nib.spatialimages.SpatialImage,
    csf: nib.spatialimages.SpatialImage,
    *,
    mode: MergeMode,
) -> nib.spatialimages.SpatialImage:
    """Merge RAS-oriented boundary and CSF images on exactly the same grid."""
    boundary_ras = normalize_image_to_ras(boundary, source_name="CSF boundary")
    csf_ras = normalize_image_to_ras(csf, source_name="CSF")
    if boundary_ras.shape != csf_ras.shape or not np.allclose(
        boundary_ras.affine, csf_ras.affine, rtol=0, atol=1e-4
    ):
        raise ValueError("CSF boundary and CSF images must share the same physical grid")
    result = merge_csf_boundary(
        np.asanyarray(boundary_ras.dataobj),
        np.asanyarray(csf_ras.dataobj),
        mode=mode,
    ).astype(np.float32)
    header = csf_ras.header.copy()
    header.set_data_dtype(np.float32)
    header.set_slope_inter(None, None)
    output = nib.Nifti1Image(result, csf_ras.affine, header=header)
    _, qcode = csf_ras.get_qform(coded=True)
    _, scode = csf_ras.get_sform(coded=True)
    output.set_qform(csf_ras.affine if qcode else None, code=int(qcode))
    output.set_sform(csf_ras.affine, code=int(scode) if scode else 1)
    return output


def clean_csf_boundary(
    boundary: np.ndarray,
    wm: np.ndarray,
    *,
    mode: CleanupMode,
) -> np.ndarray:
    """Apply one of the source-compatible V3 WM boundary cleanups.

    ``toads_wm`` removes boundary wherever TOADS WM is strictly greater than
    0.01. ``updated_wm`` uses the source's byte-scaled threshold of 25.5.
    """
    edge = np.asarray(boundary)
    white_matter = np.asarray(wm)
    if edge.shape != white_matter.shape:
        raise ValueError(f"Boundary and WM shapes differ: {edge.shape} vs {white_matter.shape}")
    if edge.ndim != 3:
        raise ValueError(f"Boundary and WM volumes must be 3D; got {edge.shape}")
    if not np.isfinite(edge).all() or not np.isfinite(white_matter).all():
        raise ValueError("Boundary and WM volumes must contain finite values")
    if mode == "toads_wm":
        threshold = 0.01
    elif mode == "updated_wm":
        threshold = 255 * 0.1
    else:
        raise ValueError(f"Unsupported boundary cleanup mode: {mode!r}")
    cleaned = edge.copy()
    cleaned[white_matter > threshold] = 0
    return cleaned


def clean_csf_boundary_image(
    boundary: nib.spatialimages.SpatialImage,
    wm: nib.spatialimages.SpatialImage,
    *,
    mode: CleanupMode,
) -> nib.spatialimages.SpatialImage:
    """Apply boundary cleanup to aligned NIfTI images in RAS orientation."""
    boundary_ras = normalize_image_to_ras(boundary, source_name="CSF boundary")
    wm_ras = normalize_image_to_ras(wm, source_name="white matter mask")
    if boundary_ras.shape != wm_ras.shape or not np.allclose(
        boundary_ras.affine, wm_ras.affine, rtol=0, atol=1e-4
    ):
        raise ValueError("CSF boundary and WM images must share the same physical grid")
    values = clean_csf_boundary(
        np.asanyarray(boundary_ras.dataobj), np.asanyarray(wm_ras.dataobj), mode=mode
    )
    header = boundary_ras.header.copy()
    header.set_data_dtype(values.dtype)
    header.set_slope_inter(None, None)
    output = nib.Nifti1Image(values, boundary_ras.affine, header=header)
    _, qcode = boundary_ras.get_qform(coded=True)
    _, scode = boundary_ras.get_sform(coded=True)
    output.set_qform(boundary_ras.affine if qcode else None, code=int(qcode))
    output.set_sform(boundary_ras.affine, code=int(scode) if scode else 1)
    return output
