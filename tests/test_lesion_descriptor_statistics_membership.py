from __future__ import annotations

from pathlib import Path

import numpy as np

from pymacruise.lesion_atlas_priors import load_lesion_atlas_priors
from pymacruise.lesion_descriptor_iteration_slice import run_lesion_toads_descriptor_iteration_slice
from pymacruise.lesion_descriptor_statistics_membership import run_lesion_toads_descriptor_statistics_membership
from pymacruise.lesion_similarity import compute_lesion_toads_similarity
from pymacruise.lesion_statistics_refresh import refresh_lesion_toads_statistics


ATLAS = Path(__file__).parents[1] / "data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"


def _inputs():
    descriptor = load_lesion_atlas_priors(str(ATLAS), ("T1_MPRAGE", "FLAIR"))
    types = tuple(item.topology for item in descriptor.structures)
    names = tuple(item.name for item in descriptor.structures)
    classes = len(types)
    shape = (5, 5, 5)
    index = np.indices(shape)
    labels = ((index[0] * 5 + index[1] * 3 + index[2]) % classes).astype(np.int16)
    images = np.stack((0.2 + 0.03 * index[0] + 0.02 * index[1] + 0.01 * index[2],
                       0.4 + 0.01 * index[0] + 0.04 * index[1] + 0.02 * index[2]), axis=-1).astype(np.float32)
    memberships = np.full((*shape, classes), np.float32(0.2 / (classes - 1)), dtype=np.float32)
    for k in range(classes):
        memberships[..., k] = np.where(labels == k, np.float32(0.8), memberships[..., k])
    lesion = np.full(shape, 0.03, dtype=np.float32)
    black_hole = np.zeros(shape, dtype=np.float32)
    outlier = np.zeros(shape, dtype=np.float32)
    centroids = np.zeros((2, classes), dtype=np.float32)
    shape_values = np.ones((*shape, classes), dtype=np.float32)
    vector = np.ones(classes, dtype=np.float32)
    centroid_options = {
        "mode": "none", "modalities": ("T1_MPRAGE", "FLAIR"),
        "tissue_power_exponents": np.zeros((2, classes), dtype=np.float32),
        "lesion_power_exponents": np.asarray((1.0, 0.5), dtype=np.float32),
    }
    return (images, memberships, labels, centroids, lesion, black_hole, outlier,
            np.asarray((0.5, 0.6), dtype=np.float32), np.zeros(2, dtype=np.float32),
            types, names, shape_values, vector, centroid_options, {"outlier_distance": 3.25})


def _run(*, use_lesion_weight=False):
    (images, memberships, labels, centroids, lesion, black_hole, outlier,
     lesion_centroid, black_centroid, types, names, shapes, vector,
     centroid_options, weight_options) = _inputs()
    normalized_types = tuple("object" if value == "obj" else value for value in types)
    return run_lesion_toads_descriptor_statistics_membership(
        images, memberships, labels, centroids, lesion, black_hole, outlier,
        lesion_centroid, black_centroid, types, names, fields=None,
        centroid_options=centroid_options, channel_weight_options=weight_options,
        prior_scale=0.65, initial_black_hole_similarity=np.ones(len(types), dtype=np.float32),
        smoothing=0.0, smoothing_priors=vector, shape_priors=np.zeros(len(types), dtype=np.float32),
        shape_values=shapes, relation_types=("direct",) * len(types),
        prior_coefficient=0.0, outlier_distance=3.25, previous_energy=10.0,
        t1_channel=0, use_lesion_weight=use_lesion_weight,
    )


def test_descriptor_statistics_membership_orchestrates_existing_steps_without_mutation() -> None:
    (images, memberships, labels, centroids, lesion, black_hole, outlier,
     lesion_centroid, black_centroid, types, names, shapes, vector,
     centroid_options, weight_options) = _inputs()
    normalized_types = tuple("object" if value == "obj" else value for value in types)
    originals = tuple(item.copy() for item in (images, memberships, labels, lesion, black_hole, outlier))
    result = _run()
    weighted = _run(use_lesion_weight=True)
    stats = refresh_lesion_toads_statistics(
        images, memberships, labels, centroids, lesion, black_hole, outlier,
        lesion_centroid, black_centroid, normalized_types, centroid_options=centroid_options,
        channel_weight_options=weight_options,
    )
    similarity = compute_lesion_toads_similarity(
        stats.centroids.class_centroids, stats.variances.class_standard_deviation,
        stats.centroids.lesion_centroids, stats.variances.lesion_standard_deviation,
        stats.channel_weights, normalized_types, names, prior_scale=0.65, t1_channel=0,
        black_hole_similarity=np.ones(len(types), dtype=np.float32),
    )
    transition = run_lesion_toads_descriptor_iteration_slice(
        images, memberships, normalized_types, stats.centroids.class_centroids.T,
        stats.centroids.lesion_centroids, stats.centroids.black_hole_centroids,
        stats.channel_weights, shapes, vector, np.zeros(len(types), dtype=np.float32),
        similarity.class_similarity, similarity.lesion_similarity, similarity.black_hole_similarity,
        ("direct",) * len(types), outlier_distance=3.25, smoothing=0.0,
        prior_coefficient=0.0, initial_lesion_membership=lesion,
        initial_black_hole_membership=black_hole, initial_outlier_membership=outlier,
    )
    assert result.stage_order == (
        "computeCentroids", "computeVariances", "computeChannelWeights",
        "computeCentroidSimilarities", "computeActiveClassDataEnergies",
        "computeMemberships", "computeEnergyRatio",
    )
    np.testing.assert_array_equal(result.statistics.centroids.class_centroids, stats.centroids.class_centroids)
    np.testing.assert_array_equal(result.statistics.variances.class_standard_deviation, stats.variances.class_standard_deviation)
    np.testing.assert_array_equal(result.similarity.class_similarity, similarity.class_similarity)
    np.testing.assert_array_equal(result.transition.data_energies.class_energy, transition.data_energies.class_energy)
    np.testing.assert_array_equal(result.transition.membership.class_memberships, transition.membership.class_memberships)
    np.testing.assert_array_equal(result.transition.membership.denominator, transition.membership.denominator)
    assert np.isfinite(result.energy_ratio)
    assert np.isfinite(weighted.energy_ratio)
    assert not np.array_equal(
        result.transition.data_energies.lesion_energy,
        weighted.transition.data_energies.lesion_energy,
    )
    for value, before in zip((images, memberships, labels, lesion, black_hole, outlier), originals):
        np.testing.assert_array_equal(value, before)
