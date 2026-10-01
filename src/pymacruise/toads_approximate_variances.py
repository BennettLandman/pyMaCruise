"""Selected default-atlas TOADS approximate class-variance update."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np


@dataclass(frozen=True)
class ToadsApproximateVarianceUpdate:
    """Raw selected source accumulators and per-channel class deviations."""

    numerators: np.ndarray
    denominators: np.ndarray
    standard_deviations: np.ndarray


def compute_default_toads_approximate_variances(
    images: Any,
    memberships: Any,
    best_class_indices: Any,
    segmentation: Any,
    centroids: Any,
    *,
    intensity_field: Any | None = None,
) -> ToadsApproximateVarianceUpdate:
    """Compute selected 11-class approximate `computeVariances` output.

    This implements the non-LESION path selected by the default TOADS atlas.
    For each channel it scans x/y/z, then the four candidate slots, and counts
    a slot only when its best class equals the current segmentation label.
    Arrays use ``(x,y,z,channel)`` image/field, ``(x,y,z,slot)`` memberships
    and best labels, and ``(channel,class)`` centroids. A missing field models
    the source's optional uncorrected branch; the selected P18 path supplies
    the field from the inhomogeneity fit.
    """
    image = np.asarray(images, dtype=np.float32)
    mems = np.asarray(memberships, dtype=np.float32)
    best_raw = np.asarray(best_class_indices)
    labels = np.asarray(segmentation)
    centers = np.asarray(centroids, dtype=np.float32)
    field = None if intensity_field is None else np.asarray(intensity_field, dtype=np.float32)
    if image.ndim != 4 or image.shape[3] < 1 or any(n < 1 for n in image.shape[:3]):
        raise ValueError("images must be a nonempty (x,y,z,channel) volume")
    spatial = image.shape[:3]
    if not np.isfinite(image).all():
        raise ValueError("images must be finite")
    if mems.shape != spatial + (4,) or not np.isfinite(mems).all() or np.any(mems < 0):
        raise ValueError("memberships must be finite, nonnegative, and have four candidate slots")
    if best_raw.shape != mems.shape or not np.issubdtype(best_raw.dtype, np.integer):
        raise ValueError("best_class_indices must be an integer grid matching memberships")
    if labels.shape != spatial or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be an integer label grid matching images")
    if centers.shape != (image.shape[3], 11) or not np.isfinite(centers).all():
        raise ValueError("centroids must be finite with shape (channel,11)")
    if field is not None and (field.shape != image.shape or not np.isfinite(field).all()):
        raise ValueError("intensity_field must be finite and match images")
    if np.any(best_raw < 0) or np.any(best_raw >= 11):
        raise ValueError("best_class_indices must refer to the selected 11 classes")
    if np.any(labels < 0) or np.any(labels >= 11):
        raise ValueError("segmentation labels must refer to the selected 11 classes")

    best = best_raw.astype(np.intp, copy=False)
    numerators = np.zeros((image.shape[3], 11), dtype=np.float32)
    denominators = np.zeros_like(numerators)
    standard_deviations = np.zeros_like(numerators)
    epsilon = np.float32(1e-30)
    nx, ny, nz, channels = image.shape
    for channel in range(channels):
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    label = int(labels[x, y, z])
                    for slot in range(4):
                        klass = int(best[x, y, z, slot])
                        if label != klass:
                            continue
                        membership = mems[x, y, z, slot]
                        weight = np.float32(np.float32(membership * membership) + epsilon)
                        sample = (
                            np.float32(field[x, y, z, channel] * image[x, y, z, channel])
                            if field is not None
                            else image[x, y, z, channel]
                        )
                        residual = np.float32(sample - centers[channel, klass])
                        weighted_residual = np.float32(weight * residual)
                        contribution = np.float32(weighted_residual * residual)
                        numerators[channel, klass] = np.float32(
                            numerators[channel, klass] + contribution
                        )
                        denominators[channel, klass] = np.float32(
                            denominators[channel, klass] + weight
                        )
        for klass in range(11):
            denominator = denominators[channel, klass]
            if denominator > np.float32(0.0):
                quotient = np.float32(numerators[channel, klass] / denominator)
                standard_deviations[channel, klass] = np.float32(math.sqrt(float(quotient)))
            else:
                standard_deviations[channel, klass] = np.float32(0.0)

    return ToadsApproximateVarianceUpdate(numerators, denominators, standard_deviations)
