"""Compose the selected TOADS segmentation initialization handoff.

This adapter follows the selected ``AlgorithmMultipleToads`` call sequence
after atlas registration: initialize centroid/prior-centroid arrays, retain
robust intensity max/scale, set the atlas-derived lesion centroid, and compute
the fresh initial membership volume. Generic image handling remains in the
existing NumPy/NIfTI adapters; these fields and call order follow VandyCRUISE.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .toads_initialization import initial_exact_memberships
from .toads_lesion_centroid_initialization import initial_toads_lesion_centroid


@dataclass(frozen=True)
class ToadsSegmentationInitialization:
    """Source-shaped state handed to the next TOADS initialization stages."""

    centroids: NDArray[np.float32]
    prior_centroids: NDArray[np.float32]
    intensity_max: NDArray[np.float32]
    intensity_scale: NDArray[np.float32]
    lesion_centroid: NDArray[np.float32]
    memberships: NDArray[np.float32]
    approximate_class_count: int
    use_approximate_memberships: bool
    call_order: tuple[str, ...]


TOADS_SEGMENTATION_INITIALIZATION_CALL_ORDER = (
    "segmentation.initCentroids",
    "segmentation.setIntensityMax",
    "segmentation.setIntensityScale",
    "segmentation.setLesionCentroid",
    "preprocess.initialMemberships",
)


def initialize_toads_segmentation_state(
    *,
    images: Any,
    centroids: Any,
    robust_minimum: Any,
    robust_maximum: Any,
    lesion_prior_by_channel: Any,
    crop_border: int = 10,
) -> ToadsSegmentationInitialization:
    """Build the wrapper's post-registration initial segmentation state.

    The selected base wrapper sets ``useApprox=false`` and ``approx=classes``;
    its atlas-aware membership overload therefore calls exact memberships.
    Arrays are RAS voxel ordered, with channel-first images and channel/class
    centroid matrices. Source array setters retain references; this value
    object instead snapshots arrays so later pipeline stages cannot mutate the
    recorded evidence accidentally.
    """
    image_values = np.asarray(images, dtype=np.float32)
    center_values = np.asarray(centroids, dtype=np.float32)
    low = np.asarray(robust_minimum, dtype=np.float32)
    high = np.asarray(robust_maximum, dtype=np.float32)
    lesion_prior = np.asarray(lesion_prior_by_channel, dtype=np.float32)
    if image_values.ndim != 4 or any(size <= 0 for size in image_values.shape):
        raise ValueError("images must have shape (channels,x,y,z) with nonempty axes")
    channels = image_values.shape[0]
    if center_values.ndim != 2 or center_values.shape[0] != channels or center_values.shape[1] < 1:
        raise ValueError("centroids must have shape (channels,classes)")
    if low.shape != (channels,) or high.shape != (channels,) or lesion_prior.shape != (channels,):
        raise ValueError("intensity bounds and lesion priors must have one value per channel")
    if not all(np.isfinite(value).all() for value in (image_values, center_values, low, high, lesion_prior)):
        raise ValueError("segmentation initialization arrays must be finite")
    span = np.subtract(high, low, dtype=np.float32)
    if np.any(span <= 0):
        raise ValueError("robust intensity spans must be positive")

    memberships = initial_exact_memberships(
        image_values, center_values, crop_border=crop_border
    )
    lesion_centroid = initial_toads_lesion_centroid(low, high, lesion_prior)
    return ToadsSegmentationInitialization(
        center_values.copy(),
        center_values.copy(),
        high.copy(),
        span.copy(),
        lesion_centroid,
        memberships,
        int(center_values.shape[1]),
        False,
        TOADS_SEGMENTATION_INITIALIZATION_CALL_ORDER,
    )
