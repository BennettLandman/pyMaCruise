"""Selected default-atlas TOADS approximate centroid update.

This ports only the 11-class/four-candidate `computeApproxCentroids` branch.
The fixed class order and source traversal are kept explicit; it is not a
general replacement for all TOADS centroid modes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .toads_math import toads_powertable_lookup


DEFAULT_TOADS_11_CLASS_TYPES = (
    "MASK",
    "OPTIMIZED_CSF",
    "OPTIMIZED_CSF",
    "OPTIMIZED_GM",
    "OPTIMIZED_GM",
    "OPTIMIZED_GM",
    "OPTIMIZED_GM",
    "OPTIMIZED_GM",
    "OPTIMIZED_WM",
    "OPTIMIZED_WM",
    "OPTIMIZED_WM",
)
_OPTIMIZED_TYPES = ("OPTIMIZED_CSF", "OPTIMIZED_GM", "OPTIMIZED_WM")
_ZERO = np.float32(1e-10)
_POWER_STEP = float(np.float32(1e-6))


@dataclass(frozen=True)
class ToadsApproximateCentroidUpdate:
    """Accumulators and updated channel-by-class centroids."""

    numerators: np.ndarray
    denominators: np.ndarray
    centroids: np.ndarray


def update_default_toads_approximate_centroids(
    images: Any,
    memberships: Any,
    best_class_indices: Any,
    previous_centroids: Any,
    intensity_field: Any,
    optimized_exponents: Any,
    *,
    centroid_smoothness: float = 0.1,
) -> ToadsApproximateCentroidUpdate:
    """Update selected default-atlas centroids from four approximate slots.

    The Java traversal is channel, full-grid x/y/z, then approximate slot.
    Inputs are never modified. Arrays use ``(x,y,z,channel)`` for images and
    field, ``(x,y,z,slot)`` for memberships and best-class indices, and
    ``(channel,class)`` for centroids and exponents.
    """
    image = np.asarray(images, dtype=np.float32)
    mems = np.asarray(memberships, dtype=np.float32)
    best_raw = np.asarray(best_class_indices)
    old = np.asarray(previous_centroids, dtype=np.float32)
    field = np.asarray(intensity_field, dtype=np.float32)
    exponents = np.asarray(optimized_exponents, dtype=np.float32)
    if image.ndim != 4 or image.shape[3] < 1 or any(n < 1 for n in image.shape[:3]):
        raise ValueError("images must be a nonempty (x,y,z,channel) volume")
    if field.shape != image.shape or not np.isfinite(field).all():
        raise ValueError("intensity_field must be finite and match images")
    if mems.shape != image.shape[:3] + (4,) or not np.isfinite(mems).all():
        raise ValueError("memberships must be finite and match the image with four candidate slots")
    if np.any(mems < np.float32(0.0)):
        raise ValueError("memberships must be nonnegative")
    if best_raw.shape != mems.shape or not np.issubdtype(best_raw.dtype, np.integer):
        raise ValueError("best_class_indices must be an integer grid matching memberships")
    if old.shape != (image.shape[3], 11) or not np.isfinite(old).all():
        raise ValueError("previous_centroids must be finite with shape (channel,11)")
    if exponents.shape != (image.shape[3], 3) or not np.isfinite(exponents).all():
        raise ValueError("optimized_exponents must be finite with shape (channel,3)")
    if not np.isfinite(image).all():
        raise ValueError("images must be finite")
    if np.any(best_raw < 0) or np.any(best_raw >= len(DEFAULT_TOADS_11_CLASS_TYPES)):
        raise ValueError("best_class_indices contains a class outside the selected 11-class atlas")
    smoothness = np.float32(centroid_smoothness)
    if not np.isfinite(smoothness):
        raise ValueError("centroid_smoothness must be finite")

    best = best_raw.astype(np.intp, copy=False)
    numerators = np.zeros(old.shape, dtype=np.float32)
    denominators = np.zeros(old.shape, dtype=np.float32)
    nx, ny, nz, channels = image.shape
    for channel in range(channels):
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    for slot in range(4):
                        klass = int(best[x, y, z, slot])
                        kind = DEFAULT_TOADS_11_CLASS_TYPES[klass]
                        membership = mems[x, y, z, slot]
                        corrected = np.float32(field[x, y, z, channel] * image[x, y, z, channel])
                        factor = np.float32(1.0)
                        if kind in _OPTIMIZED_TYPES:
                            residual = np.float32(corrected - old[channel, klass])
                            distance = np.float32(residual * residual)
                            if distance > _ZERO:
                                exponent_index = _OPTIMIZED_TYPES.index(kind)
                                exponent = exponents[channel, exponent_index]
                                powered = toads_powertable_lookup(
                                    np.asarray(distance, dtype=np.float32),
                                    float(exponent),
                                    cached_exponent=float(exponent),
                                    minimum=0.0,
                                    maximum=4.0,
                                    step=_POWER_STEP,
                                )
                                factor = np.float32(float(powered) / float(distance))
                        weighted = np.float32(np.float32(factor * membership) * membership)
                        numerator_term = np.float32(weighted * corrected)
                        numerators[channel, klass] = np.float32(
                            numerators[channel, klass] + numerator_term
                        )
                        denominators[channel, klass] = np.float32(
                            denominators[channel, klass] + weighted
                        )

    updated = old.copy()
    one_minus_smoothness = np.float32(np.float32(1.0) - smoothness)
    for channel in range(channels):
        for klass in range(11):
            denominator = denominators[channel, klass]
            if denominator > np.float32(0.0):
                numerator_scaled = np.float32(one_minus_smoothness * numerators[channel, klass])
                weighted_mean = np.float32(numerator_scaled / denominator)
                retained = np.float32(smoothness * old[channel, klass])
                updated[channel, klass] = np.float32(weighted_mean + retained)
    return ToadsApproximateCentroidUpdate(numerators, denominators, updated)
