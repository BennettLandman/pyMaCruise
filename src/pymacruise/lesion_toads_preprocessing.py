"""Source-derived intensity normalization and centroid initialization.

Source anchor: ``toads-cruise/src/edu/jhmi/rad/medic/methods/LesionToadPreprocess.java``
(March 2007; Pierre-Louis Bazin). The adapter and algorithm revision details
are recorded in ``docs/toads-port-plan.md`` and ``docs/assets-and-dependencies.md``.

Paper context: Shiee et al. (2010) describe Lesion-TOADS as multichannel,
statistical-atlas segmentation with a lesion class. This helper covers only
the source-derived intensity normalization and initial lesion/black-hole
centroid values used to prepare that model; it does not implement the paper's
iterative segmentation. See ``docs/papers/lesion-toads.md`` for the paper-to-
code boundary and the selected-source-version caveat.

The selected ``AlgorithmLesionToads`` pins ``robustMaxMin=false``. The source
constructor clamps negative channel values and records ordinary minima and
maxima. The wrapper then calls ``normalizeImages`` before cropping, which
applies ``(image - Imin) / (Imax - Imin)``. Later, after crop, separate
``initialLesionCentroid`` and ``initialBlackHoleCentroid`` calls estimate
robust bounds from those normalized images and use
``low + prior * (high - low)``. Robust-span normalization and
``low / span + prior`` are only the inactive ``robustMaxMin=true`` branch.
The calculations here retain Java single-precision operation order.

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
class LesionToadsImageNormalization:
    """Active wrapper normalization outputs before crop."""

    normalized_images: np.ndarray
    intensity_minima: np.ndarray
    intensity_maxima: np.ndarray


@dataclass(frozen=True)
class LesionToadsCentroids:
    """Active wrapper's post-normalization, post-crop centroid outputs."""

    robust_minima: np.ndarray
    robust_maxima: np.ndarray
    lesion_centroid: np.ndarray
    black_hole_centroid: np.ndarray


@dataclass(frozen=True)
class LesionToadsInitialization:
    """Convenience no-crop composition of the selected wrapper substeps."""

    normalized_images: np.ndarray
    intensity_minima: np.ndarray
    intensity_maxima: np.ndarray
    centroid_robust_minima: np.ndarray
    centroid_robust_maxima: np.ndarray
    lesion_centroid: np.ndarray
    black_hole_centroid: np.ndarray


def normalize_lesion_toads_images(
    images: Sequence[np.ndarray] | np.ndarray,
) -> LesionToadsImageNormalization:
    """Apply the selected wrapper's active min/max normalization before crop.

    This ports ``LesionToadPreprocess`` constructor clamp/min/max followed by
    ``normalizeImages`` with ``robustMaxMin=false``. Every channel is clamped
    at zero, then normalized by ``(image - Imin) / (Imax - Imin)``. It does not
    calculate robust intensity bounds or centroids; the wrapper computes those
    later from the cropped normalized images.

    Zero or nonfinite ordinary spans raise ``ValueError``. Java divides by
    zero in this case and returns NaN/infinity; rejecting it is a deliberate
    invalid-input deviation.
    """
    channel_images = _validated_channel_images(images)
    channels = channel_images.shape[0]
    prepared = channel_images.copy()
    minima = np.empty(channels, dtype=np.float32)
    maxima = np.empty(channels, dtype=np.float32)
    normalized = np.empty_like(prepared, dtype=np.float32)

    for channel in range(channels):
        source = prepared[channel]
        minimum = np.float32(source[0, 0, 0])
        maximum = np.float32(source[0, 0, 0])
        for x in range(source.shape[0]):
            for y in range(source.shape[1]):
                for z in range(source.shape[2]):
                    value = np.float32(source[x, y, z])
                    if value < np.float32(0.0):
                        source[x, y, z] = np.float32(0.0)
                        minimum = np.float32(0.0)
                    else:
                        if value > maximum:
                            maximum = value
                        if value < minimum:
                            minimum = value

        span = np.float32(maximum - minimum)
        if not np.isfinite(span) or span <= np.float32(0.0):
            raise ValueError(
                f"LesionTOADS channel {channel} has zero or invalid ordinary intensity span; "
                "the Java source would divide by zero"
            )
        minima[channel] = minimum
        maxima[channel] = maximum
        normalized[channel] = np.asarray(
            (source - minimum) / span, dtype=np.float32
        )

    return LesionToadsImageNormalization(normalized, minima, maxima)


def initialize_lesion_toads_centroids(
    normalized_cropped_images: Sequence[np.ndarray] | np.ndarray,
    lesion_priors: Sequence[float] | np.ndarray,
    black_hole_priors: Sequence[float] | np.ndarray,
    *,
    robust_fraction: float = 0.01,
    robust_scales: int = 2,
) -> LesionToadsCentroids:
    """Initialize lesion/black-hole centroids from the wrapper's later inputs.

    The selected wrapper calls this source operation after image normalization
    and ``cropImages``. ``robustMinimum``/``robustMaximum`` estimate bounds on
    each supplied normalized cropped channel. Centroids are then
    ``low + prior * (high - low)`` in Java float32 operation order. This does
    not perform normalization or cropping; callers must provide the current
    cropped normalized arrays to match the wrapper exactly.
    """
    channel_images = _validated_channel_images(normalized_cropped_images)
    channels = channel_images.shape[0]
    lesion_prior_values = _validated_priors(lesion_priors, channels, "lesion")
    black_hole_prior_values = _validated_priors(
        black_hole_priors, channels, "black-hole"
    )
    robust_minima = np.empty(channels, dtype=np.float32)
    robust_maxima = np.empty(channels, dtype=np.float32)
    lesion_centroid = np.empty(channels, dtype=np.float32)
    black_hole_centroid = np.empty(channels, dtype=np.float32)

    for channel in range(channels):
        low, high = robust_intensity_bounds(
            channel_images[channel],
            fraction=robust_fraction,
            scales=robust_scales,
        )
        span = np.float32(high - low)
        if not np.isfinite(span) or span <= np.float32(0.0):
            raise ValueError(
                f"LesionTOADS channel {channel} has zero or invalid robust centroid span; "
                "the Java source would divide by zero"
            )
        robust_minima[channel] = low
        robust_maxima[channel] = high
        low_span = np.float32(low + np.float32(lesion_prior_values[channel] * span))
        bh_span = np.float32(
            low + np.float32(black_hole_prior_values[channel] * span)
        )
        lesion_centroid[channel] = low_span
        black_hole_centroid[channel] = bh_span

    return LesionToadsCentroids(
        robust_minima,
        robust_maxima,
        lesion_centroid,
        black_hole_centroid,
    )


def initialize_lesion_toads_intensities(
    images: Sequence[np.ndarray] | np.ndarray,
    lesion_priors: Sequence[float] | np.ndarray,
    black_hole_priors: Sequence[float] | np.ndarray,
    *,
    robust_fraction: float = 0.01,
    robust_scales: int = 2,
) -> LesionToadsInitialization:
    """Compose active normalization and centroid initialization without crop.

    This convenience function runs the source normalization stage and then
    computes the later source centroids on the same whole-volume arrays. The
    selected wrapper crops between those stages, so use
    :func:`normalize_lesion_toads_images` and
    :func:`initialize_lesion_toads_centroids` separately when a cropped image
    is available. ``robust_fraction`` and ``robust_scales`` apply only to the
    centroid bounds, never to active image normalization.
    """
    normalized = normalize_lesion_toads_images(images)
    centroids = initialize_lesion_toads_centroids(
        normalized.normalized_images,
        lesion_priors,
        black_hole_priors,
        robust_fraction=robust_fraction,
        robust_scales=robust_scales,
    )
    return LesionToadsInitialization(
        normalized_images=normalized.normalized_images,
        intensity_minima=normalized.intensity_minima,
        intensity_maxima=normalized.intensity_maxima,
        centroid_robust_minima=centroids.robust_minima,
        centroid_robust_maxima=centroids.robust_maxima,
        lesion_centroid=centroids.lesion_centroid,
        black_hole_centroid=centroids.black_hole_centroid,
    )


def _validated_channel_images(images: Sequence[np.ndarray] | np.ndarray) -> np.ndarray:
    channel_images = np.asarray(images, dtype=np.float32)
    if channel_images.ndim != 4 or channel_images.shape[0] == 0:
        raise ValueError(
            "LesionTOADS intensities must have shape (channels, x, y, z) with at least one channel"
        )
    if any(size == 0 for size in channel_images.shape[1:]):
        raise ValueError("LesionTOADS intensity channels must be nonempty")
    if not np.isfinite(channel_images).all():
        raise ValueError("LesionTOADS intensities contain NaN or infinite values")
    return channel_images


def _validated_priors(
    values: Sequence[float] | np.ndarray, channels: int, name: str
) -> np.ndarray:
    priors = np.asarray(values, dtype=np.float32)
    if priors.shape != (channels,):
        raise ValueError(f"LesionTOADS {name} priors must have shape ({channels},)")
    if not np.isfinite(priors).all():
        raise ValueError(f"LesionTOADS {name} priors contain NaN or infinite values")
    return priors
