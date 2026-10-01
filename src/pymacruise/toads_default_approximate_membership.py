"""Optional top-k initial memberships for the selected 11-class TOADS atlas.

This ports only the atlas-aware overload of
``MultipleToadPreprocess.initialApproxMemberships``. The selected wrapper's
first rigid atlas pyramid instead calls the no-atlas ``initialMemberships``
overload; this atlas-aware state belongs to a later reinitialization after
segmentation state exists. Neither overload implies that the subsequent
exact-membership segmentation sweep is sparse.

Paper context: Bazin and Pham (2007, 2008), ``docs/papers/toads.md``, describe
statistical class memberships and atlas priors in TOADS. This atlas-aware
top-k overload is a selected registration-initialization stage. Its class
grouping, candidate ordering, and forced-slot details come from VandyCRUISE
source; it is not the papers' full membership update or iterative
segmentation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES


@dataclass(frozen=True)
class ToadsDefaultApproximateMembershipInitialization:
    """RAS-grid top-k values and source class indices for alignment."""

    memberships: NDArray[np.float32]
    best_class_indices: NDArray[np.uint8]
    class_types: tuple[str, ...]


def _best_shape_index(values: NDArray[np.float32]) -> int:
    winner = 0
    for index in range(1, values.size):
        if values[index] > values[winner]:
            winner = index
    return winner


def _top_k(
    values: NDArray[np.float32], count: int, *, first: int,
) -> tuple[NDArray[np.uint8], NDArray[np.float32]]:
    """Match selected ``Numerics.bestIndex(..., first)`` scan/tie order."""
    scratch = values.copy()
    indices = np.empty(count, dtype=np.uint8)
    selected = np.empty(count, dtype=np.float32)
    for slot in range(count):
        winner = int(first) if slot == 0 else 0
        for candidate in range(scratch.size):
            if candidate != winner and scratch[candidate] > scratch[winner]:
                winner = candidate
        indices[slot] = winner
        selected[slot] = scratch[winner]
        scratch[winner] = np.float32(-1e30)
    return indices, selected


def initialize_toads_default_approximate_memberships(
    image_values: ArrayLike,
    centroids: ArrayLike,
    transformed_shapes: ArrayLike,
    *,
    approx: int,
    process_mask: ArrayLike | None = None,
) -> ToadsDefaultApproximateMembershipInitialization:
    """Compute selected 11-class atlas-seeded initial top-k memberships.

    For each processed voxel and class, the selected source adds squared
    channel residuals and every other transformed class shape squared, then
    normalizes reciprocal costs and retains ``approx`` candidates. Candidate
    slot zero is seeded from the maximum transformed shape before the strict
    membership scan, exactly as the Java ``bestIndex(..., first)`` overload
    does. Values are float32 and the input arrays are treated as RAS ordered.
    """

    image = np.asarray(image_values, dtype=np.float32)
    centers = np.asarray(centroids, dtype=np.float32)
    shapes = np.asarray(transformed_shapes, dtype=np.float32)
    classes = len(DEFAULT_ATLAS_CLASS_TYPES)
    if image.ndim != 4 or any(size <= 0 for size in image.shape[:3]) or image.shape[-1] <= 0:
        raise ValueError("image_values must have shape (x,y,z,channel) with nonempty axes")
    if shapes.shape != image.shape[:3] + (classes,):
        raise ValueError("transformed_shapes must match the subject grid and default class order")
    if centers.shape != (image.shape[-1], classes):
        raise ValueError("centroids must have shape (channel,11)")
    if not np.isfinite(image).all() or not np.isfinite(shapes).all() or not np.isfinite(centers).all():
        raise ValueError("image, atlas, and centroid values must be finite")
    if isinstance(approx, bool) or not isinstance(approx, (int, np.integer)) or not (1 < approx < classes):
        raise ValueError("approx must select between two and class_count-1 candidates")
    if process_mask is None:
        process = np.ones(image.shape[:3], dtype=bool)
    else:
        process = np.asarray(process_mask, dtype=bool)
        if process.shape != image.shape[:3]:
            raise ValueError("process_mask must match the subject grid")

    memberships = np.zeros(image.shape[:3] + (int(approx),), dtype=np.float32)
    best = np.zeros(image.shape[:3] + (int(approx),), dtype=np.uint8)
    for x in range(image.shape[0]):
        for y in range(image.shape[1]):
            for z in range(image.shape[2]):
                if not process[x, y, z]:
                    continue
                costs = np.zeros(classes, dtype=np.float32)
                for klass in range(classes):
                    total = np.float32(0.0)
                    for channel in range(image.shape[-1]):
                        residual = np.float32(image[x, y, z, channel] - centers[channel, klass])
                        total = np.float32(total + np.float32(residual * residual))
                    for other in range(classes):
                        if other != klass:
                            prior = shapes[x, y, z, other]
                            total = np.float32(total + np.float32(prior * prior))
                    costs[klass] = total
                reciprocal = np.empty(classes, dtype=np.float32)
                denominator = np.float32(0.0)
                for klass, cost in enumerate(costs):
                    reciprocal[klass] = (
                        np.float32(1.0 / cost) if cost > np.float32(1e-20) else np.float32(1e20)
                    )
                    denominator = np.float32(denominator + reciprocal[klass])
                if denominator > np.float32(0.0):
                    probabilities = np.divide(reciprocal, denominator, dtype=np.float32)
                else:
                    probabilities = reciprocal
                first = _best_shape_index(shapes[x, y, z])
                ids, values = _top_k(probabilities, int(approx), first=first)
                best[x, y, z] = ids
                memberships[x, y, z] = values
    return ToadsDefaultApproximateMembershipInitialization(
        memberships, best, DEFAULT_ATLAS_CLASS_TYPES,
    )
