"""Compose the selected LesionTOADS post-thinning statistics update."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_statistics import (
    LesionToadsAuxiliaryCentroids,
    LesionToadsVariances,
    compute_lesion_toads_centroids,
    compute_lesion_toads_channel_weights,
    compute_lesion_toads_variances,
)
from .toads_math import toads_powertable_lookup


@dataclass(frozen=True)
class LesionToadsStatisticsRefresh:
    """Results of one source-ordered statistics refresh."""

    centroids: LesionToadsAuxiliaryCentroids
    variances: LesionToadsVariances
    channel_weights: NDArray[np.float32]


def refresh_lesion_toads_statistics(
    images: ArrayLike,
    memberships: ArrayLike,
    segmentation: ArrayLike,
    class_centroids: ArrayLike,
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    outlier_membership: ArrayLike,
    lesion_centroids: ArrayLike,
    black_hole_centroids: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    *,
    fields: ArrayLike | None = None,
    centroid_options: Mapping[str, Any] | None = None,
    channel_weight_options: Mapping[str, Any] | None = None,
    use_lesion_weight: bool = False,
    flair_channel: int | None = None,
) -> LesionToadsStatisticsRefresh:
    """Apply ``computeCentroids``, ``computeVariances``, and
    ``computeChannelWeights`` in the selected wrapper order.

    Inputs are fixed maps/state after the wrapper's thinning and optional
    correction-field update. No correction field is estimated and no input
    arrays are mutated. Centroid options are forwarded to the already
    validated mixed-category centroid port; channel-weight options must
    include the source ``outlier_distance`` setting. Optional ``fields`` are
    supplied to all three stages consistently.
    """
    centroid_kwargs = dict(centroid_options or {})
    weight_kwargs = dict(channel_weight_options or {})
    if "fields" in centroid_kwargs or "fields" in weight_kwargs:
        raise ValueError("pass fields only through the fields argument")

    updated_centroids = compute_lesion_toads_centroids(
        images,
        memberships,
        segmentation,
        class_centroids,
        object_types,
        fields=fields,
        lesion_membership=lesion_membership,
        black_hole_membership=black_hole_membership,
        lesion_centroids=lesion_centroids,
        black_hole_centroids=black_hole_centroids,
        **centroid_kwargs,
    )
    if updated_centroids.lesion_centroids is None:
        raise ValueError("statistics refresh requires a LESION class and lesion centroid")
    if updated_centroids.black_hole_centroids is None:
        raise ValueError("statistics refresh requires a BLACKHOLE class and black-hole centroid")
    if flair_channel is None:
        modality_names = centroid_kwargs.get("modalities") or ()
        flair_channel = max(
            (index for index, name in enumerate(modality_names) if name == "FLAIR"),
            default=-1,
        )

    # The selected wrapper's active computeCentroids() branch has
    # linkedBlackHoleCentroid hard-coded false, so it writes the direct
    # lesion weighted mean. The standalone mixed-centroid API intentionally
    # preserves the separately validated 20% damping behavior; replace only
    # this auxiliary result locally to match the active wrapper branch.
    direct_lesion_centroid = _source_direct_lesion_centroid(
        images,
        memberships,
        segmentation,
        lesion_membership,
        black_hole_membership,
        lesion_centroids,
        object_types,
        fields=fields,
        lesion_power_exponents=centroid_kwargs.get("lesion_power_exponents"),
        optimized_auxiliary_power_exponents=centroid_kwargs.get(
            "optimized_auxiliary_power_exponents"
        ),
        use_lesion_weight=use_lesion_weight,
        flair_channel=flair_channel,
    )
    updated_centroids = LesionToadsAuxiliaryCentroids(
        updated_centroids.class_centroids,
        direct_lesion_centroid,
        updated_centroids.black_hole_centroids,
    )

    updated_variances = compute_lesion_toads_variances(
        images,
        memberships,
        segmentation,
        updated_centroids.class_centroids,
        lesion_membership,
        black_hole_membership,
        updated_centroids.lesion_centroids,
        updated_centroids.black_hole_centroids,
        object_types,
        fields=fields,
    )
    updated_weights = compute_lesion_toads_channel_weights(
        images,
        memberships,
        updated_centroids.class_centroids,
        lesion_membership,
        black_hole_membership,
        outlier_membership,
        updated_centroids.lesion_centroids,
        updated_centroids.black_hole_centroids,
        object_types,
        fields=fields,
        **weight_kwargs,
    )
    return LesionToadsStatisticsRefresh(
        updated_centroids, updated_variances, updated_weights
    )


def _source_direct_lesion_centroid(
    images: ArrayLike,
    memberships: ArrayLike,
    segmentation: ArrayLike,
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    previous_centroid: ArrayLike,
    object_types: tuple[str, ...] | list[str],
    *,
    fields: ArrayLike | None,
    lesion_power_exponents: ArrayLike | None,
    optimized_auxiliary_power_exponents: ArrayLike | None,
    use_lesion_weight: bool,
    flair_channel: int,
) -> NDArray[np.float32]:
    """Reproduce the wrapper-active direct LESION centroid branch only."""
    values = np.asarray(images, dtype=np.float32)
    weights = np.asarray(memberships, dtype=np.float32)
    labels = np.asarray(segmentation)
    lesion_map = np.asarray(lesion_membership, dtype=np.float32)
    black_map = np.asarray(black_hole_membership, dtype=np.float32)
    previous = np.asarray(previous_centroid, dtype=np.float32)
    field_values = None if fields is None else np.asarray(fields, dtype=np.float32)
    types = tuple(item.lower() for item in object_types)
    if not isinstance(use_lesion_weight, bool):
        raise ValueError("use_lesion_weight must be a boolean")
    channels = values.shape[3]
    if not isinstance(flair_channel, int) or flair_channel < -1 or flair_channel >= channels:
        raise ValueError("flair_channel must be -1 or index an image channel")
    exponents = (
        None if lesion_power_exponents is None
        else np.asarray(lesion_power_exponents, dtype=np.float32)
    )
    optimized_exponents = (
        None if optimized_auxiliary_power_exponents is None
        else np.asarray(optimized_auxiliary_power_exponents, dtype=np.float32)
    )
    if exponents is not None and exponents.shape != (channels,):
        raise ValueError("lesion_power_exponents must have one value per channel")
    optimized_auxiliary = any(
        kind in {"optimized_lesion", "optimized_blackhole"} for kind in types
    )
    if optimized_auxiliary and (
        optimized_exponents is None or optimized_exponents.shape != (channels,)
    ):
        raise ValueError("optimized_auxiliary_power_exponents must match image channels")
    if exponents is None:
        # Selected unfielded FLAIR always uses the lesion table; fielded uses
        # it only when the wrapper's explicit lesion-weight switch is enabled.
        if flair_channel >= 0 and (fields is None or use_lesion_weight):
            raise ValueError("lesion_power_exponents are required for the active FLAIR lesion term")

    result = np.zeros(channels, dtype=np.float32)
    shape = values.shape[:3]
    for channel in range(channels):
        numerator = np.float32(0.0)
        denominator = np.float32(0.0)
        for class_index, kind in enumerate(types):
            if kind not in {"lesion", "blackhole"}:
                continue
            for x in range(shape[0]):
                for y in range(shape[1]):
                    for z in range(shape[2]):
                        if labels[x, y, z] != class_index:
                            continue
                        class_mem = weights[x, y, z, class_index]
                        lesion_mem = lesion_map[x, y, z]
                        black_mem = black_map[x, y, z]
                        winner = max(class_mem, lesion_mem, black_mem)
                        if winner == class_mem or winner != lesion_mem:
                            continue
                        sample = values[x, y, z, channel]
                        if field_values is not None:
                            sample = np.float32(field_values[x, y, z, channel] * sample)
                        factor = np.float32(1.0)
                        optimized = kind in {"optimized_lesion", "optimized_blackhole"}
                        apply_power = optimized or (
                            channel == flair_channel
                            and (fields is None or use_lesion_weight)
                            and winner == lesion_mem
                        )
                        if apply_power:
                            if optimized:
                                assert optimized_exponents is not None
                                distance = np.float32(
                                    np.float32(winner * winner) * sample
                                )
                                exponent = float(optimized_exponents[channel])
                            else:
                                assert exponents is not None
                                difference = np.float32(sample - previous[channel])
                                distance = np.float32(difference * difference)
                                exponent = float(exponents[channel])
                            if distance > np.float32(1.0e-10):
                                powered = toads_powertable_lookup(
                                    distance,
                                    exponent,
                                    cached_exponent=exponent,
                                )
                                factor = np.float32(
                                    np.float64(powered) / np.float64(distance)
                                )
                        weighted = np.float32(factor * lesion_mem)
                        weighted = np.float32(weighted * lesion_mem)
                        numerator = np.float32(numerator + np.float32(weighted * sample))
                        denominator = np.float32(denominator + weighted)
        result[channel] = (
            np.float32(numerator / denominator)
            if denominator > np.float32(0.0)
            else np.float32(0.0)
        )
    return result
