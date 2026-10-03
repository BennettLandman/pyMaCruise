from __future__ import annotations

import numpy as np

from pymacruise.toads_initial_relation_state import toads_initialize_relation_state
from pymacruise.toads_initial_statistics_membership import compute_toads_initial_statistics_membership
from pymacruise.toads_initial_statistics_thinning import (
    TOADS_INITIAL_STATISTICS_THINNING_ORDER,
    compute_toads_initial_statistics_thinning,
)
from pymacruise.toads_thinning import toads_propagate_thinning
from test_toads_default_atlas_dispatch import _fixture


def test_exact_initial_membership_flows_into_relations_and_first_thinning() -> None:
    data, parameters, _cases, _reference = _fixture()
    shape = data["images"].shape[:3]
    x, y, z = np.indices(shape)
    segmentation = ((x + 2 * y + 3 * z) % 11).astype(np.int8)
    kwargs = {
        "image_values": data["images"],
        "initial_memberships": data["initial_memberships"],
        "centroids": parameters["centroids"],
        "segmentation": segmentation,
        "class_types": parameters["class_types"],
        "class_variance_weights": parameters["class_variance_weights"],
        "exponent_by_channel_class": parameters["exponent_by_channel_class"],
        "modality_weights": parameters["modality_weights"],
        "transformed_shapes": data["transformed_shapes"],
        "smoothing_priors": parameters["smoothing_priors"],
        "atlas_shape_priors": parameters["atlas_shape_priors"],
        "smoothing": 0.043,
        "prior_factor": parameters["prior_factor"],
        "prior_coefficient": 0.16,
        "relation_limit": 2,
        "volume_ratio": 0.5,
        "thinning_limit": 3.0,
    }
    snapshots = {
        name: value.copy()
        for name, value in {**data, "segmentation": segmentation}.items()
    }
    result = compute_toads_initial_statistics_thinning(**kwargs)

    assert result.stage_order == TOADS_INITIAL_STATISTICS_THINNING_ORDER
    expected_relation = toads_initialize_relation_state(
        segmentation, data["initial_memberships"], limit=2
    )
    np.testing.assert_array_equal(
        result.relation_state.relation.class_adjacency,
        expected_relation.relation.class_adjacency,
    )
    np.testing.assert_array_equal(
        result.relation_state.closest_labels,
        expected_relation.closest_labels,
    )
    expected_membership = compute_toads_initial_statistics_membership(
        **{key: value for key, value in kwargs.items() if key not in {"relation_limit", "volume_ratio", "thinning_limit"}},
        closest_labels=expected_relation.closest_labels,
    )
    np.testing.assert_array_equal(
        result.statistics_membership.membership_pass.memberships,
        expected_membership.membership_pass.memberships,
    )
    expected_thinning = toads_propagate_thinning(
        segmentation,
        expected_membership.membership_pass.memberships,
        expected_relation.relation.class_adjacency,
        data["transformed_shapes"],
        volume_ratio=0.5,
        limit=3.0,
    )
    np.testing.assert_array_equal(result.thinning.segmentation, expected_thinning.segmentation)
    np.testing.assert_array_equal(result.thinning.skeleton, expected_thinning.skeleton)
    assert result.thinning.processed_items == expected_thinning.processed_items
    for name, before in snapshots.items():
        np.testing.assert_array_equal(data[name] if name != "segmentation" else segmentation, before)
