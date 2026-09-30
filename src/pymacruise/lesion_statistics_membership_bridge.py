"""Compose one fixed-input LesionTOADS statistics-to-membership update.

The selected wrapper refreshes statistics after thinning, computes data
energies and memberships, then updates its energy ratio. This module composes
those already validated stages; thinning, correction-field estimation,
topology, atlas registration, and stopping decisions remain caller-owned.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_data_energy_maps import (
    LesionToadsDataEnergyMaps,
    compute_lesion_toads_data_energy_maps,
)
from .lesion_objective_membership_pass import (
    LesionObjectiveMembershipPass,
    run_lesion_toads_objective_membership_pass,
)
from .lesion_similarity import LesionToadsSimilarity, compute_lesion_toads_similarity
from .lesion_statistics_refresh import (
    LesionToadsStatisticsRefresh,
    refresh_lesion_toads_statistics,
)


@dataclass(frozen=True)
class LesionToadsStatisticsMembershipBridge:
    """Intermediate and final results for the composed wrapper segment."""

    statistics: LesionToadsStatisticsRefresh
    similarity: LesionToadsSimilarity
    data_energies: LesionToadsDataEnergyMaps
    outmask_class_energy: NDArray[np.float32]
    outmask_data_energy: NDArray[np.float32]
    membership: LesionObjectiveMembershipPass
    energy_ratio: np.float32


def run_lesion_toads_statistics_membership_bridge(
    images: ArrayLike,
    initial_class_memberships: ArrayLike,
    segmentation: ArrayLike,
    initial_class_centroids: ArrayLike,
    initial_lesion_membership: ArrayLike,
    initial_black_hole_membership: ArrayLike,
    initial_outlier_membership: ArrayLike,
    initial_lesion_centroid: ArrayLike,
    initial_black_hole_centroid: ArrayLike,
    object_types: Sequence[str],
    *,
    fields: ArrayLike | None,
    centroid_options: Mapping[str, Any],
    channel_weight_options: Mapping[str, Any],
    structure_names: Sequence[str],
    prior_scale: float,
    initial_black_hole_similarity: ArrayLike,
    smoothing: float,
    class_smoothing_priors: ArrayLike,
    shape_priors: ArrayLike,
    class_shapes: ArrayLike,
    relation_types: Sequence[str],
    prior_coefficient: float,
    outlier_distance: float,
    previous_energy: float,
    use_lesion_weight: bool = False,
    flair_channel: int | None = None,
    t1_channel: int = -1,
    zero: float = 1.0e-10,
    infinity: float = 1.0e10,
) -> LesionToadsStatisticsMembershipBridge:
    """Refresh statistics, construct data terms, update memberships, ratio.

    The supported class order is exactly OBJECT, LESION, BLACKHOLE, OUTMASK.
    For OUTMASK, the selected per-channel ``maskDistance(image, 0)`` and
    ``outlier_distance`` terms are composed locally because the raw-energy
    adapter covers OBJECT/LESION/BLACKHOLE. For the other classes, centroids,
    deviations, and channel weights come from the just-refreshed statistics.
    Shape, spatial relations, and similarity-independent state are fixed
    inputs to this bridge.
    """
    if any(not isinstance(item, str) for item in object_types):
        raise ValueError("object_types must be strings")
    types = tuple(item.upper() for item in object_types)
    expected_types = ("OBJECT", "LESION", "BLACKHOLE", "OUTMASK")
    if types != expected_types:
        raise ValueError("bridge supports class order OBJECT, LESION, BLACKHOLE, OUTMASK")
    image_values = np.asarray(images, dtype=np.float32)
    membership_values = np.asarray(initial_class_memberships, dtype=np.float32)
    if image_values.ndim != 4 or membership_values.ndim != 4:
        raise ValueError("images and class memberships must be 4-D RAS arrays")
    if image_values.shape[:3] != membership_values.shape[:3]:
        raise ValueError("images and memberships must share the same spatial grid")
    if membership_values.shape[-1] != len(types):
        raise ValueError("membership class count must match object_types")
    statistics_types = tuple(item.lower() for item in types)
    if "outlier_distance" not in channel_weight_options:
        raise ValueError("channel_weight_options must include outlier_distance")
    if np.float32(channel_weight_options["outlier_distance"]) != np.float32(outlier_distance):
        raise ValueError("outlier_distance must match the statistics channel-weight input")

    statistics = refresh_lesion_toads_statistics(
        image_values,
        membership_values,
        segmentation,
        initial_class_centroids,
        initial_lesion_membership,
        initial_black_hole_membership,
        initial_outlier_membership,
        initial_lesion_centroid,
        initial_black_hole_centroid,
        statistics_types,
        fields=fields,
        centroid_options=centroid_options,
        channel_weight_options=channel_weight_options,
        use_lesion_weight=use_lesion_weight,
        flair_channel=flair_channel,
    )
    lesion_centroid = statistics.centroids.lesion_centroids
    black_hole_centroid = statistics.centroids.black_hole_centroids
    if lesion_centroid is None or black_hole_centroid is None:
        raise ValueError("the bridge requires LESION and BLACKHOLE centroids")
    standard_deviation = statistics.variances.class_standard_deviation
    similarity = compute_lesion_toads_similarity(
        statistics.centroids.class_centroids,
        standard_deviation,
        lesion_centroid,
        statistics.variances.lesion_standard_deviation,
        statistics.channel_weights,
        statistics_types,
        structure_names,
        prior_scale=prior_scale,
        t1_channel=t1_channel,
        black_hole_similarity=initial_black_hole_similarity,
    )

    class_energy_indices = (0, 1, 2)
    resolved_flair = flair_channel
    if resolved_flair is None:
        modalities = tuple(centroid_options.get("modalities") or ())
        resolved_flair = max(
            (index for index, modality in enumerate(modalities) if modality == "FLAIR"),
            default=-1,
        )
    lesion_exponents = centroid_options.get("lesion_power_exponents")
    lesion_exponent = (
        None
        if lesion_exponents is None or resolved_flair < 0
        else float(np.asarray(lesion_exponents)[resolved_flair])
    )
    energies = compute_lesion_toads_data_energy_maps(
        image_values,
        statistics.centroids.class_centroids[:, class_energy_indices].T,
        lesion_centroid,
        black_hole_centroid,
        statistics.channel_weights,
        ("object", "lesion", "blackhole"),
        use_field=fields is not None,
        intensity_fields=fields,
        flair_channel=resolved_flair,
        use_lesion_weight=use_lesion_weight,
        lesion_power_exponent=lesion_exponent,
        cached_lesion_power_exponent=lesion_exponent,
    )
    class_terms = np.zeros((*image_values.shape[:3], len(types)), dtype=np.float32)
    class_terms[..., :3] = energies.class_energy
    outmask_energy = np.zeros(image_values.shape[:3], dtype=np.float32)
    outmask_data = np.zeros(image_values.shape[:3], dtype=np.float32)
    distance = np.float32(outlier_distance)
    if not np.isfinite(distance):
        raise ValueError("outlier_distance must be finite")
    for x in range(1, image_values.shape[0] - 1):
        for y in range(1, image_values.shape[1] - 1):
            for z in range(1, image_values.shape[2] - 1):
                for channel in range(image_values.shape[3]):
                    value = image_values[x, y, z, channel]
                    mismatch = np.float32(0.0 if value == np.float32(0.0) else 1.0)
                    outmask_energy[x, y, z] = np.float32(
                        outmask_energy[x, y, z]
                        + np.float32(statistics.channel_weights[channel] * mismatch)
                    )
                    outmask_data[x, y, z] = np.float32(
                        outmask_data[x, y, z]
                        + np.float32(statistics.channel_weights[channel] * distance)
                    )
    class_terms[..., 3] = outmask_energy

    membership = run_lesion_toads_objective_membership_pass(
        membership_values,
        initial_lesion_membership,
        initial_black_hole_membership,
        initial_outlier_membership,
        class_terms,
        energies.lesion_energy,
        energies.black_hole_energy,
        outmask_data,
        class_shapes,
        class_smoothing_priors,
        types,
        relation_types,
        smoothing=smoothing,
        shape_priors=shape_priors,
        class_similarities=similarity.class_similarity,
        lesion_similarities=similarity.lesion_similarity,
        black_hole_similarities=similarity.black_hole_similarity,
        prior_coefficient=prior_coefficient,
        zero=zero,
        infinity=infinity,
    )
    previous = np.float32(previous_energy)
    energy = np.float32(membership.objective)
    numerator = np.float32(np.float32(2.0) * np.float32(previous - energy))
    denominator = np.float32(previous + energy)
    with np.errstate(divide="ignore", invalid="ignore"):
        energy_ratio = np.float32(np.abs(np.float32(numerator / denominator)))
    return LesionToadsStatisticsMembershipBridge(
        statistics,
        similarity,
        energies,
        outmask_energy,
        outmask_data,
        membership,
        energy_ratio,
    )
