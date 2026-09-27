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


def combine_multi_atlas_toads(
    gm_ma: np.ndarray,
    wm_ma: np.ndarray,
    brain_mask_ma: np.ndarray,
    gm_toads: np.ndarray,
    wm_toads: np.ndarray,
    mask_toads: np.ndarray,
    *,
    gm_confident_level: float = 0.001,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Combine atlas and TOADS tissue maps using the selected V3 Java rules.

    Preserves the source priority: confident eroded atlas WM first, confident
    TOADS/atlas GM second, then TOADS values; voxels outside the atlas brain
    mask are zeroed. ``MASImorph.Minimize3DFilter`` uses integer inputs, a
    center-plus-six-neighbor minimum, and leaves the output volume border at
    zero. The MATLAB-only alternative is not the selected V3 default.
    """
    arrays = [
        np.asarray(value)
        for value in (gm_ma, wm_ma, brain_mask_ma, gm_toads, wm_toads, mask_toads)
    ]
    shapes = {value.shape for value in arrays}
    if len(shapes) != 1 or arrays[0].ndim != 3:
        raise ValueError("All atlas and TOADS inputs must be aligned 3D arrays of the same shape")
    if any(not np.issubdtype(value.dtype, np.number) or not np.isfinite(value).all() for value in arrays):
        raise ValueError("All atlas and TOADS inputs must contain finite numeric values")
    if not np.isfinite(gm_confident_level) or gm_confident_level < 0:
        raise ValueError("GM confidence level must be finite and non-negative")

    gm_atlas, wm_atlas, brain_mask, gm_toads_data, wm_toads_data, mask_toads_data = arrays
    # ImageDataInt conversion in the source precedes integer erosion and mask
    # reads. NumPy's truncation toward zero matches Java's numeric cast here.
    wm_int = np.trunc(wm_atlas).astype(np.int32)
    brain_int = np.trunc(brain_mask).astype(np.int32)
    gm_atlas_int = np.trunc(gm_atlas).astype(np.int32)
    wm_toads_int = np.trunc(wm_toads_data).astype(np.int32)
    mask_toads_int = np.trunc(mask_toads_data).astype(np.int32)

    eroded_wm = np.zeros(wm_int.shape, dtype=np.int32)
    if all(length >= 3 for length in wm_int.shape):
        six_neighbor = ndimage.generate_binary_structure(3, 1)
        min_neighbors = ndimage.minimum_filter(
            wm_int, footprint=six_neighbor, mode="constant", cval=0
        )
        interior = tuple(slice(1, -1) for _ in range(3))
        eroded_wm[interior] = min_neighbors[interior]

    inside = brain_int > 0
    confident_wm = eroded_wm > 0
    confident_gm = (
        (gm_toads_data > 0)
        & (gm_toads_data < (255 * gm_confident_level))
        & (gm_atlas_int > 0)
    )

    gm_out = np.where(inside, gm_toads_data, 0).astype(np.float32)
    wm_out = np.where(inside, wm_toads_data, 0).astype(np.float32)
    mask_out = np.where(inside, mask_toads_int, 0).astype(np.int32)

    gm_out[inside & confident_wm] = 0
    wm_out[inside & confident_wm] = 255
    mask_out[inside & confident_wm] = 255

    gm_only = inside & ~confident_wm & confident_gm
    gm_out[gm_only] = 255
    wm_out[gm_only] = 0
    mask_out[gm_only] = 0
    return gm_out, wm_out, mask_out


def combine_multi_atlas_toads_images(
    gm_ma: nib.spatialimages.SpatialImage,
    wm_ma: nib.spatialimages.SpatialImage,
    brain_mask_ma: nib.spatialimages.SpatialImage,
    gm_toads: nib.spatialimages.SpatialImage,
    wm_toads: nib.spatialimages.SpatialImage,
    mask_toads: nib.spatialimages.SpatialImage,
    *,
    gm_confident_level: float = 0.001,
) -> tuple[nib.spatialimages.SpatialImage, nib.spatialimages.SpatialImage, nib.spatialimages.SpatialImage]:
    """Run V3 atlas/TOADS combination with aligned RAS NIfTI inputs/outputs."""
    inputs = [
        normalize_image_to_ras(image, source_name=name)
        for image, name in zip(
            (gm_ma, wm_ma, brain_mask_ma, gm_toads, wm_toads, mask_toads),
            ("atlas GM", "atlas WM", "atlas brain mask", "TOADS GM", "TOADS WM", "TOADS mask"),
            strict=True,
        )
    ]
    reference = inputs[0]
    if any(
        image.shape != reference.shape
        or not np.allclose(image.affine, reference.affine, rtol=0, atol=1e-4)
        for image in inputs[1:]
    ):
        raise ValueError("All atlas and TOADS images must share the same physical grid")
    arrays = [np.asanyarray(image.dataobj) for image in inputs]
    outputs = combine_multi_atlas_toads(*arrays, gm_confident_level=gm_confident_level)
    result_images = []
    for output_data, template, dtype in zip(
        outputs, (inputs[3], inputs[4], inputs[5]), (np.float32, np.float32, np.int32), strict=True
    ):
        output_data = output_data.astype(dtype, copy=False)
        header = template.header.copy()
        header.set_data_dtype(dtype)
        header.set_slope_inter(None, None)
        image = nib.Nifti1Image(output_data, template.affine, header=header)
        _, qcode = template.get_qform(coded=True)
        _, scode = template.get_sform(coded=True)
        image.set_qform(template.affine if qcode else None, code=int(qcode))
        image.set_sform(template.affine, code=int(scode) if scode else 1)
        result_images.append(image)
    return tuple(result_images)  # type: ignore[return-value]
