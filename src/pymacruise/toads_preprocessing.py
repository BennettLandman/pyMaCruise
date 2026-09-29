"""Source-derived intensity preparation for the selected TOADS workflow.

Paper context: robust intensity normalization supports the tissue-classification
stage described by Bazin & Pham (2007, 2008); the concrete histogram bins and
iteration order here follow the selected Java ``ImageFunctions`` source, not a
formula inferred from those papers. See ``docs/papers/toads.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import nibabel as nib
import numpy as np

from pymacruise.image_io import normalize_image_to_ras


@dataclass(frozen=True)
class ToadsIntensityResult:
    image: nib.Nifti1Image
    robust_minimum: float
    robust_maximum: float
    masked_voxels: int


def robust_intensity_bounds(
    values: np.ndarray,
    *,
    fraction: float = 0.01,
    scales: int = 2,
) -> tuple[np.float32, np.float32]:
    """Return the source `ImageFunctions` robust min/max histogram estimates."""
    image = np.asarray(values, dtype=np.float32)
    if image.ndim != 3 or any(size == 0 for size in image.shape):
        raise ValueError(f"TOADS intensity input must be a nonempty 3D array; got {image.shape}")
    if not np.isfinite(image).all():
        raise ValueError("TOADS intensity input contains NaN or infinite values")
    if not np.isfinite(fraction) or fraction <= 0 or fraction > 1:
        raise ValueError("Robust intensity fraction must be finite and within (0, 1]")
    if not isinstance(scales, int) or scales < 1:
        raise ValueError("Robust intensity scales must be a positive integer")
    flat = image.ravel(order="F")
    return (
        _robust_minimum(flat, fraction, scales),
        _robust_maximum(flat, fraction, scales),
    )


def prepare_toads_intensity(
    t1_n4: nib.spatialimages.SpatialImage,
    brain_mask: nib.spatialimages.SpatialImage,
    *,
    robust_fraction: float = 0.01,
    robust_scales: int = 2,
) -> ToadsIntensityResult:
    """Mask corrected T1 and apply TOADS' robust-span normalization.

    Both images are axis-reoriented to RAS without interpolation. The mask
    must be binary and occupy the exact T1 physical grid. Matching the source,
    intensities are divided by ``robust_maximum - robust_minimum`` without
    subtracting the robust minimum; zeros outside the mask remain zero.
    """
    source_orientation = nib.aff2axcodes(t1_n4.affine)
    normalized_t1 = normalize_image_to_ras(t1_n4, source_name="TOADS N4 T1")
    normalized_mask = normalize_image_to_ras(brain_mask, source_name="TOADS brain mask")
    if len(normalized_t1.shape) != 3 or len(normalized_mask.shape) != 3:
        raise ValueError("TOADS T1 and brain mask must be 3D single-volume images")
    if normalized_t1.shape != normalized_mask.shape:
        raise ValueError("TOADS T1 and brain mask shapes differ after RAS normalization")
    if not np.allclose(normalized_t1.affine, normalized_mask.affine, rtol=0, atol=1e-4):
        raise ValueError("TOADS brain mask and T1 are not aligned on the same physical grid")

    t1_values = np.asarray(normalized_t1.dataobj, dtype=np.float32)
    mask_values = np.asanyarray(normalized_mask.dataobj)
    if not np.isfinite(t1_values).all():
        raise ValueError("TOADS N4 T1 contains NaN or infinite values")
    if not np.isin(mask_values, (0, 1, 255)).all() or not np.any(mask_values):
        raise ValueError("TOADS brain mask must be nonempty and binary (0/1 or 0/255)")
    masked = np.where(mask_values > 0, t1_values, np.float32(0.0)).astype(np.float32)
    low, high = robust_intensity_bounds(
        masked,
        fraction=robust_fraction,
        scales=robust_scales,
    )
    span = np.float32(high - low)
    if span <= 0:
        raise ValueError(
            "TOADS robust intensity range is zero; source normalization would divide by zero"
        )
    normalized_values = np.asarray(masked / span, dtype=np.float32)
    result = nib.Nifti1Image(normalized_values, normalized_t1.affine)
    return ToadsIntensityResult(
        result,
        float(low),
        float(high),
        int(mask_values.size - np.count_nonzero(mask_values)),
    )


def _scaled_boundary(low: np.float32, high: np.float32, numerator: int) -> np.float32:
    fraction = np.float32(np.float32(numerator) / np.float32(10.0))
    return np.float32(low + np.float32(fraction * np.float32(high - low)))


def _histogram(values: np.ndarray, low: np.float32, high: np.float32) -> np.ndarray:
    counts = np.zeros(10, dtype=np.float32)
    first_upper = np.float32(
        low + np.float32(np.float32(1.0) / np.float32(10.0) * np.float32(high - low))
    )
    counts[0] = np.float32(np.count_nonzero((values >= low) & (values <= first_upper)))
    for index in range(1, 10):
        lower = _scaled_boundary(low, high, index)
        upper = _scaled_boundary(low, high, index + 1)
        counts[index] = np.float32(np.count_nonzero((values > lower) & (values <= upper)))
    return counts


def _robust_minimum(values: np.ndarray, fraction: float, scales: int) -> np.float32:
    low = np.float32(values.min())
    high = np.float32(values.max())
    ratio = np.float32(np.float32(fraction) * np.float32(values.size))
    result = low
    for _ in range(scales):
        counts = _histogram(values, low, high)
        count = np.float32(0.0)
        index = 0
        while count < ratio and index < 10:
            count = np.float32(count + counts[index])
            index += 1
        position = np.float32(np.float32(index) - np.float32(0.5))
        fraction_of_span = np.float32(position / np.float32(10.0))
        result = np.float32(
            low + np.float32(fraction_of_span * np.float32(high - low))
        )
        old_low = low
        low = _scaled_boundary(old_low, high, index - 1)
        high = _scaled_boundary(old_low, high, index)
        ratio = np.float32(ratio - np.float32(count - counts[index - 1]))
    return result


def _robust_maximum(values: np.ndarray, fraction: float, scales: int) -> np.float32:
    low = np.float32(values.min())
    high = np.float32(values.max())
    ratio = np.float32(np.float32(fraction) * np.float32(values.size))
    result = high
    for _ in range(scales):
        counts = _histogram(values, low, high)
        count = np.float32(0.0)
        index = 10
        while count < ratio and index > 0:
            index -= 1
            count = np.float32(count + counts[index])
        position = np.float32(np.float32(index) + np.float32(0.5))
        fraction_of_span = np.float32(position / np.float32(10.0))
        result = np.float32(
            low + np.float32(fraction_of_span * np.float32(high - low))
        )
        old_low = low
        low = _scaled_boundary(old_low, high, index)
        high = _scaled_boundary(old_low, high, index + 1)
        ratio = np.float32(ratio - np.float32(count - counts[index]))
    return result
