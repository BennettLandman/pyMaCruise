"""Source-derived intensity initialization for the LesionTOADS algorithm.

Source anchor: ``toads-cruise/src/edu/jhmi/rad/medic/methods/LesionToadPreprocess.java``
(March 2007; Pierre-Louis Bazin). The adapter and algorithm revision details
are recorded in ``docs/toads-port-plan.md`` and ``docs/assets-and-dependencies.md``.

The source ``LesionToadPreprocess`` constructor clamps negative channel values
to zero, estimates robust bounds independently per channel, and
``normalizeImages`` divides values by the robust span without subtracting the
robust minimum. Lesion and black-hole initial centroids use the same
``low / span + atlas_prior`` expression. The calculations here retain Java
single-precision operation order.

Inputs are channel-first voxel arrays with RAS voxel-axis ordering established
by the caller. NIfTI orientation and shared-grid validation belong to the I/O
workflow, before arrays enter this numerical helper.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from pymacruise.toads_preprocessing import robust_intensity_bounds


@dataclass(frozen=True)
class LesionToadsInitialization:
    """Normalized channels and source-initialized lesion model values."""

    normalized_images: np.ndarray
    robust_minima: np.ndarray
    robust_maxima: np.ndarray
    lesion_centroid: np.ndarray
    black_hole_centroid: np.ndarray


def initialize_lesion_toads_intensities(
    images: Sequence[np.ndarray] | np.ndarray,
    lesion_priors: Sequence[float] | np.ndarray,
    black_hole_priors: Sequence[float] | np.ndarray,
    *,
    robust_fraction: float = 0.01,
    robust_scales: int = 2,
) -> LesionToadsInitialization:
    """Apply LesionTOADS robust-span normalization and initialize centroids.

    ``images`` has shape ``(channels, x, y, z)``. Negative voxel values are
    replaced with zero before bounds are estimated, matching
    ``LesionToadPreprocess``. Each channel is independently normalized by its
    robust maximum minus robust minimum, retaining zero as zero. Initial
    centroids use the source expression ``low / span + prior`` for each
    channel.

    A zero robust span is rejected with ``ValueError``. The Java source does
    not guard this case and consequently produces NaN for zero-valued input or
    infinity for a nonzero constant image. This is a deliberate, documented
    invalid-input deviation; normal finite-span results are source-faithful.
    """
    channel_images = np.asarray(images, dtype=np.float32)
    if channel_images.ndim != 4 or channel_images.shape[0] == 0:
        raise ValueError(
            "LesionTOADS intensities must have shape (channels, x, y, z) with at least one channel"
        )
    if any(size == 0 for size in channel_images.shape[1:]):
        raise ValueError("LesionTOADS intensity channels must be nonempty")
    if not np.isfinite(channel_images).all():
        raise ValueError("LesionTOADS intensities contain NaN or infinite values")

    channels = channel_images.shape[0]
    lesion_prior_values = _validated_priors(lesion_priors, channels, "lesion")
    black_hole_prior_values = _validated_priors(
        black_hole_priors, channels, "black-hole"
    )

    prepared = np.maximum(channel_images, np.float32(0.0)).astype(np.float32)
    robust_minima = np.empty(channels, dtype=np.float32)
    robust_maxima = np.empty(channels, dtype=np.float32)
    normalized = np.empty_like(prepared, dtype=np.float32)
    lesion_centroid = np.empty(channels, dtype=np.float32)
    black_hole_centroid = np.empty(channels, dtype=np.float32)

    for channel in range(channels):
        low, high = robust_intensity_bounds(
            prepared[channel],
            fraction=robust_fraction,
            scales=robust_scales,
        )
        span = np.float32(high - low)
        if not np.isfinite(span) or span <= np.float32(0.0):
            raise ValueError(
                f"LesionTOADS channel {channel} has zero or invalid robust intensity span; "
                "the Java source would divide by zero"
            )

        robust_minima[channel] = low
        robust_maxima[channel] = high
        normalized[channel] = np.asarray(prepared[channel] / span, dtype=np.float32)

        low_over_span = np.float32(low / span)
        lesion_centroid[channel] = np.float32(
            low_over_span + lesion_prior_values[channel]
        )
        black_hole_centroid[channel] = np.float32(
            low_over_span + black_hole_prior_values[channel]
        )

    return LesionToadsInitialization(
        normalized_images=normalized,
        robust_minima=robust_minima,
        robust_maxima=robust_maxima,
        lesion_centroid=lesion_centroid,
        black_hole_centroid=black_hole_centroid,
    )


def _validated_priors(
    values: Sequence[float] | np.ndarray, channels: int, name: str
) -> np.ndarray:
    priors = np.asarray(values, dtype=np.float32)
    if priors.shape != (channels,):
        raise ValueError(f"LesionTOADS {name} priors must have shape ({channels},)")
    if not np.isfinite(priors).all():
        raise ValueError(f"LesionTOADS {name} priors contain NaN or infinite values")
    return priors
