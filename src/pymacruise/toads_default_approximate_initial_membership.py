"""Intensity-only top-k memberships used for initial default TOADS alignment.

This implements the no-atlas overload of
``MultipleToadPreprocess.initialApproxMemberships`` called by the selected
wrapper before ``registerShapesApproxPyramid``. The separate atlas-aware
overload is exposed by ``toads_default_approximate_membership`` for the later
reinitialization path.

Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``, describe
statistical tissue classification as part of atlas-driven TOADS. This
intensity-only seed supplies class candidates for the selected rigid
registration prepass; its duplicate-centroid scaling, strict top-k ties, and
call order are VandyCRUISE details, not equations asserted by those papers.
It is an initialization slice, not the paper's iterative segmentation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES


@dataclass(frozen=True)
class ToadsDefaultApproximateInitialMemberships:
    """Source-ordered intensity-only top-k memberships and class slots."""

    memberships: NDArray[np.float32]
    best_class_indices: NDArray[np.uint8]
    class_types: tuple[str, ...]


def initialize_toads_default_approximate_initial_memberships(
    image_values: ArrayLike,
    centroids: ArrayLike,
    *,
    approx: int,
) -> ToadsDefaultApproximateInitialMemberships:
    """Compute the selected no-atlas initialization used before the pyramid.

    Each class cost is the float32 sum of squared channel residuals. The
    source inverts each cost with its ``1e-20`` cutoff, sums inverse costs,
    multiplies by the count of classes sharing the exact first-channel
    centroid, and divides by the common denominator. Java's strict top-k scan
    then returns candidates in descending membership order, retaining lower
    class indices on ties. There is no shape penalty and no forced slot.
    """

    image = np.asarray(image_values, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    classes = len(DEFAULT_ATLAS_CLASS_TYPES)
    if image.ndim != 4 or any(size <= 0 for size in image.shape[:3]) or image.shape[-1] <= 0:
        raise ValueError("image_values must have shape (x,y,z,channel) with nonempty axes")
    if centers.shape != (image.shape[-1], classes):
        raise ValueError("centroids must have shape (channel,11)")
    if not np.isfinite(image).all() or not np.isfinite(centers).all():
        raise ValueError("image and centroid values must be finite")
    if isinstance(approx, bool) or not isinstance(approx, (int, np.integer)) or not (1 < approx < classes):
        raise ValueError("approx must select between two and class_count-1 candidates")

    memberships = np.zeros(image.shape[:3] + (int(approx),), dtype=np.float32)
    best = np.zeros(image.shape[:3] + (int(approx),), dtype=np.uint8)
    grouping = np.zeros(classes, dtype=np.float32)
    for klass in range(classes):
        count = np.float32(0.0)
        for other in range(classes):
            if centers[0, other] == centers[0, klass]:
                count = np.float32(count + np.float32(1.0))
        grouping[klass] = count

    for x in range(image.shape[0]):
        for y in range(image.shape[1]):
            for z in range(image.shape[2]):
                reciprocal = np.empty(classes, dtype=np.float32)
                denominator = np.float32(0.0)
                for klass in range(classes):
                    distance = np.float32(0.0)
                    for channel in range(image.shape[-1]):
                        residual = np.float32(image[x, y, z, channel] - centers[channel, klass])
                        distance = np.float32(distance + np.float32(residual * residual))
                    reciprocal[klass] = (
                        np.float32(1.0 / distance)
                        if distance > np.float32(1e-20)
                        else np.float32(1e20)
                    )
                    denominator = np.float32(denominator + reciprocal[klass])
                values = np.empty(classes, dtype=np.float32)
                if denominator > np.float32(0.0):
                    for klass in range(classes):
                        scaled = np.float32(grouping[klass] * reciprocal[klass])
                        values[klass] = np.float32(scaled / denominator)
                else:
                    values[:] = reciprocal

                scratch = values.copy()
                for slot in range(int(approx)):
                    winner = 0
                    for candidate in range(1, classes):
                        if scratch[candidate] > scratch[winner]:
                            winner = candidate
                    best[x, y, z, slot] = winner
                    memberships[x, y, z, slot] = scratch[winner]
                    scratch[winner] = np.float32(-1e30)
    return ToadsDefaultApproximateInitialMemberships(
        memberships, best, DEFAULT_ATLAS_CLASS_TYPES,
    )
