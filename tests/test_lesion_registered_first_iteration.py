from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_next_iteration_membership import (
    run_lesion_toads_next_gm_membership_update,
)
from pymacruise.lesion_next_iteration_statistics import (
    refresh_statistics_after_next_iteration_thinning,
)
from pymacruise.lesion_registered_first_iteration import (
    run_lesion_toads_registered_first_iteration,
)
from pymacruise.lesion_registered_initial_membership import (
    LesionRegisteredInitialMembership,
    run_lesion_toads_registered_initial_membership,
)
from pymacruise.lesion_standard13_voxel import STANDARD13_TYPES


FIXTURE = Path(__file__).parent / "data" / "synthetic" / "lesion_registered_initial_membership"


def _initial_state() -> tuple[
    LesionRegisteredInitialMembership,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    images = np.asanyarray(nib.load(FIXTURE / "images_ras.nii.gz").dataobj).astype(np.float32)
    labels = np.asanyarray(nib.load(FIXTURE / "constructor_labels_ras.nii.gz").dataobj).astype(np.int16)
    shapes = np.asanyarray(nib.load(FIXTURE / "transformed_shapes_ras.nii.gz").dataobj).astype(np.float32)
    centroids = np.stack(
        (np.linspace(0.2, 0.7, 13), np.linspace(0.25, 0.8, 13))
    ).astype(np.float32)
    lesion_centroids = np.asarray((0.57, 0.63), dtype=np.float32)
    black_centroids = np.asarray((0.34, 0.4), dtype=np.float32)
    initial = run_lesion_toads_registered_initial_membership(
        images,
        labels,
        centroids,
        lesion_centroids,
        black_centroids,
        shapes,
        labels.copy(),
        np.zeros((*labels.shape, 13), dtype=np.bool_),
        np.zeros(labels.shape, dtype=np.float32),
        iterations=1,
        first_thinning_limit=1000.0,
    )
    return initial, images, labels, shapes, centroids, np.stack((lesion_centroids, black_centroids))


def test_registered_initial_thinning_composes_into_first_loop_membership() -> None:
    initial, images, _labels, shapes, centroids, auxiliary_centroids = _initial_state()
    lesion_centroids, black_centroids = auxiliary_centroids
    shape_priors = np.ones(13, dtype=np.float32)
    smoothing_priors = np.ones(13, dtype=np.float32)
    centroid_options = {
        "mode": "none",
        "modalities": ("T1_MPRAGE", "FLAIR"),
        "lesion_power_exponents": (0.4, 1.2),
        "tissue_power_exponents": np.ones((2, 13), dtype=np.float32),
    }
    channel_weight_options = {"outlier_distance": 0.1, "norm_type": "Identity"}

    result = run_lesion_toads_registered_first_iteration(
        initial,
        images,
        centroids,
        lesion_centroids,
        black_centroids,
        shapes,
        shape_priors,
        smoothing_priors,
        centroid_options=centroid_options,
        channel_weight_options=channel_weight_options,
        outlier_energy=0.1,
        lesion_power_exponent=1.2,
    )

    previous = initial.membership
    direct_statistics = refresh_statistics_after_next_iteration_thinning(
        initial.first_thinning.thinning,
        images,
        previous.class_memberships,
        previous.lesion_membership,
        previous.black_hole_membership,
        previous.outlier_membership,
        centroids,
        lesion_centroids,
        black_centroids,
        STANDARD13_TYPES,
        centroid_options=centroid_options,
        channel_weight_options=channel_weight_options,
        use_lesion_weight=True,
        flair_channel=1,
    )
    direct_relations = initial.relations.relation_membership[
        initial.relations.relation_map
    ]
    direct_membership = run_lesion_toads_next_gm_membership_update(
        direct_statistics,
        images,
        shapes,
        shape_priors,
        smoothing_priors,
        direct_relations,
        previous.class_memberships,
        previous.lesion_membership,
        previous.black_hole_membership,
        previous.outlier_membership,
        previous.lesion_relation,
        initial.distance_maps.ventricle,
        initial.distance_maps.gray_matter,
        initial.distance_maps.brainstem,
        previous_energy=float(previous.objective),
        outlier_energy=0.1,
        lesion_power_exponent=1.2,
        spread=1.0,
        source_zero_variance_nan=True,
    )

    np.testing.assert_array_equal(
        result.statistics.statistics.centroids.class_centroids,
        direct_statistics.statistics.centroids.class_centroids,
    )
    np.testing.assert_array_equal(
        result.statistics.statistics.variances.class_standard_deviation,
        direct_statistics.statistics.variances.class_standard_deviation,
    )
    np.testing.assert_array_equal(
        result.statistics.statistics.channel_weights,
        direct_statistics.statistics.channel_weights,
    )
    np.testing.assert_array_equal(
        result.membership.membership.class_memberships,
        direct_membership.membership.class_memberships,
    )
    assert result.membership.membership.objective == direct_membership.membership.objective
    assert (
        result.membership.energy_ratio.energy_ratio
        == direct_membership.energy_ratio.energy_ratio
    )
    assert result.stage_order == (
        "propagateThinning",
        "computeCentroids",
        "computeVariances",
        "computeChannelWeights",
        "computeMemberships:GMflag=true",
        "energyRatio",
    )
    assert result.membership.membership.class_memberships.shape == (*images.shape[:3], 13)
