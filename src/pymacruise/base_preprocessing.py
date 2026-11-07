"""Selected BASE mask-to-axial-image handoff before the first TOADS pass.

Behavior references: MaCRUISE SVN r9360 ``run_MultiAtlasStripping.m``;
``MedicAlgorithmSPECTRE``; Blake Lucas' ``ReorientVolume`` revision 1.6;
MIPAV ``AlgorithmTransform`` and ``AlgorithmWSinc``. The source retains
the axial WINDOWED_SINC transform even when skull stripping is already done.
Its unnormalized lookup table and float coordinate arithmetic change values
even on an unchanged lattice. This module ports that selected axial operation,
not general oblique resampling or the SPECTRE skull-stripping algorithm.
"""

from __future__ import annotations

import math

import nibabel as nib
import numpy as np
from numba import njit
from scipy.ndimage import maximum_filter

from pymacruise.image_io import normalize_image_to_ras, nifti_ras_on_reference


def _axis_weights(size: int, spacing: float, flip: bool, table: np.ndarray):
    resolution = np.float32(spacing)
    translation = float(np.float32(resolution * np.float32(size - 1))) if flip else 0.0
    indices = np.empty((size, 8), dtype=np.int64)
    weights = np.empty((size, 8), dtype=np.float64)
    valid = np.empty(size, dtype=np.bool_)
    for index in range(size):
        physical = float(np.float32(np.float32(index) * resolution))
        position = (translation - physical if flip else physical) / float(resolution)
        valid[index] = 0 <= position < size
        base = int(position)
        for slot, sample in enumerate(range(base - 3, base + 5)):
            indices[index, slot] = max(0, min(size - 1, sample))
            weights[index, slot] = table[int(abs(249.75 * (position - sample)))]
    return indices, weights, valid


@njit(cache=True, fastmath=False)
def _evaluate_axial_sinc(values, active, ix, wx, vx, iy, wy, vy, iz, wz, vz, fill):
    nx, ny, nz = values.shape
    result = np.zeros(values.shape, dtype=np.float32)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                if not (vx[x] and vy[y] and vz[z]):
                    result[x, y, z] = fill
                    continue
                # Skip only neighborhoods whose complete source support is zero.
                if not active[nx - 1 - x, y, z]:
                    continue
                total = 0.0
                for qx in range(8):
                    sx = ix[x, qx]
                    for qy in range(8):
                        sy = iy[y, qy]
                        wxy = wx[x, qx] * wy[y, qy]
                        for qz in range(8):
                            total += wxy * wz[z, qz] * float(values[sx, sy, iz[z, qz]])
                # The selected transform sets clip=false: retain ringing.
                result[x, y, z] = total
    return result


def axial_windowed_sinc_ras(values: np.ndarray, spacing) -> np.ndarray:
    """Apply the selected RAS-to-axial source operation and return RAS values.

    This is an axis transform with unchanged dimensions/resolutions, not an
    affine registration. FileIO reverses NIfTI y in its buffer; ReorientVolume
    then reverses x. SPECTRE marks that buffer LPS. Undoing its LPS storage
    order returns values on the original RAS lattice.
    """
    image = np.asarray(values, dtype=np.float32)
    resolution = np.asarray(spacing, dtype=np.float32)
    if image.ndim != 3 or not image.size or not np.isfinite(image).all():
        raise ValueError("Axial sinc input must be a nonempty finite 3D array")
    if resolution.shape != (3,) or not np.isfinite(resolution).all() or np.any(resolution <= 0):
        raise ValueError("Axial sinc spacing must have three finite positive values")
    source = np.ascontiguousarray(image[:, ::-1, :])
    table = np.empty(1000, dtype=np.float64)
    table[0] = 1.0
    inverse_scale = 1.0 / 249.75
    for index in range(1, 1000):
        argument = math.pi * (index * inverse_scale)
        table[index] = (math.sin(argument) / argument) * 0.5 * (1.0 + math.cos(math.pi * index / 999.0))
    x = _axis_weights(source.shape[0], resolution[0], True, table)
    y = _axis_weights(source.shape[1], resolution[1], False, table)
    z = _axis_weights(source.shape[2], resolution[2], False, table)
    active = maximum_filter(source != 0, size=9, mode="nearest")
    lps = _evaluate_axial_sinc(source, active, *x, *y, *z, float(source.min()))
    return np.ascontiguousarray(lps[::-1, ::-1, :])


def prepare_base_stripped_t1(t1_n4, brain_mask) -> tuple[nib.Nifti1Image, nib.Nifti1Image]:
    """Mask before source interpolation; supply an all-voxel TOADS mask.

    Remasking after interpolation would erase source ringing outside the
    initial mask and change TOADS' histogram and crop. Both returned images
    retain the corrected T1's RAS lattice. The second output deliberately
    enables all voxels for the existing TOADS normalization adapter.
    """
    t1 = normalize_image_to_ras(t1_n4, source_name="BASE corrected T1")
    mask = normalize_image_to_ras(brain_mask, source_name="BASE brain mask")
    if len(t1.shape) != 3 or mask.shape != t1.shape or not np.allclose(mask.affine, t1.affine, rtol=0, atol=1e-4):
        raise ValueError("BASE corrected T1 and brain mask must share one 3D grid")
    mask_values = np.asanyarray(mask.dataobj)
    if not np.isin(mask_values, (0, 1, 255)).all() or not np.any(mask_values):
        raise ValueError("BASE brain mask must be nonempty and binary")
    values = np.asarray(t1.dataobj, dtype=np.float32)
    if not np.isfinite(values).all():
        raise ValueError("BASE corrected T1 must contain finite intensities")
    masked = np.where(mask_values > 0, values, np.float32(0))
    prepared = axial_windowed_sinc_ras(masked, t1.header.get_zooms()[:3])
    return (nifti_ras_on_reference(t1, prepared, dtype=np.float32),
            nifti_ras_on_reference(t1, np.ones(t1.shape, dtype=np.uint8), dtype=np.uint8))
