"""Compose exact TOADS membership initialization with the first thinning pass.

This narrow adapter follows the selected ``MedicAlgorithmToads08`` default
constructor path: full 11-class exact memberships are generated first, then
passed unchanged into the existing source-matched ``initThinningLabels`` and
``propagateThinning`` ports. Earlier registration/statistics stages are
represented by the supplied centroids and segmentation state; this is not a
complete TOADS startup or iteration.

Paper context: Bazin and Pham (2008), §2.3, describe a membership-guided,
topology-preserving fast-marching stage that thins structures before a later
growing pass. This composition corresponds only to one selected-source first
thinning pass; its heap/queue order, topology gates, and update details come
from VandyCRUISE, not from the paper. See ``docs/papers/toads.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_initialization import initial_exact_memberships
from .toads_thinning import ToadsThinningPass, toads_propagate_thinning


@dataclass(frozen=True)
class ToadsFirstThinningPass:
    """Exact initialized memberships and the selected first thinning result."""

    memberships: NDArray[np.float32]
    thinning: ToadsThinningPass


def toads_first_exact_thinning_pass(
    images: ArrayLike,
    centroids: ArrayLike,
    segmentation: ArrayLike,
    relation_matrix: ArrayLike,
    transformed_shape_values: ArrayLike,
    *,
    crop_border: int = 10,
    volume_ratio: float = 0.5,
    limit: float = 3.0,
    initial_skeleton: ArrayLike | None = None,
) -> ToadsFirstThinningPass:
    """Run the selected exact initializer into the first TOADS thinning pass.

    Inputs are already in true-RAS ``(channel,x,y,z)`` / ``(x,y,z)`` order.
    The selected VandyCRUISE wrapper uses the 11-class default atlas and its
    constructor disables approximate memberships. The wrapper defaults used
    here are a ten-voxel preprocessing border, ``volumeRatio=0.5`` and
    ``limit=3``. Callers may explicitly set those scalar parameters for
    controlled source comparisons.
    """
    image_values = np.asarray(images, dtype=np.float32)
    center_values = np.asarray(centroids, dtype=np.float32)
    labels = np.asarray(segmentation)
    if center_values.ndim != 2 or center_values.shape[1] != 11:
        raise ValueError("selected default TOADS centroids must have shape (channels, 11)")
    if labels.ndim != 3 or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be a 3-D integer label map")
    if labels.shape != image_values.shape[1:]:
        raise ValueError("segmentation must match the image spatial grid")

    memberships = initial_exact_memberships(
        image_values,
        center_values,
        crop_border=crop_border,
    )
    thinning = toads_propagate_thinning(
        labels,
        memberships,
        relation_matrix,
        transformed_shape_values,
        volume_ratio=volume_ratio,
        limit=limit,
        # The production constructor starts with a zero-filled Java byte
        # array. Accept an explicit seed only for callers replaying a harness
        # or other already-initialized Java state.
        initial_skeleton=initial_skeleton,
    )
    return ToadsFirstThinningPass(memberships=memberships, thinning=thinning)
