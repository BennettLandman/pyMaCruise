from __future__ import annotations

import numpy as np
import pytest

from pymacruise.toads_default_atlas_membership import toads_default_atlas_membership_pass
from pymacruise.toads_initial_statistics_membership import (
    TOADS_INITIAL_STATISTICS_MEMBERSHIP_ORDER,
    compute_toads_initial_statistics_membership,
)
from test_toads_default_atlas_dispatch import _fixture


def test_source_ordered_initial_statistics_feed_exact_membership_pass() -> None:
    data, parameters, _cases, _reference = _fixture()
    shape = data["images"].shape[:3]
    # Every class is represented in the synthetic assigned-label state.
    x, y, z = np.indices(shape)
    segmentation = ((x + 2 * y + 3 * z) % 11).astype(np.int8)

    result = compute_toads_initial_statistics_membership(
        image_values=data["images"],
        initial_memberships=data["initial_memberships"],
        centroids=parameters["centroids"],
        segmentation=segmentation,
        class_types=parameters["class_types"],
        class_variance_weights=parameters["class_variance_weights"],
        exponent_by_channel_class=parameters["exponent_by_channel_class"],
        modality_weights=parameters["modality_weights"],
        transformed_shapes=data["transformed_shapes"],
        closest_labels=data["closest_labels"],
        smoothing_priors=parameters["smoothing_priors"],
        atlas_shape_priors=parameters["atlas_shape_priors"],
        smoothing=0.043,
        prior_factor=parameters["prior_factor"],
        prior_coefficient=0.16,
    )

    assert result.stage_order == TOADS_INITIAL_STATISTICS_MEMBERSHIP_ORDER
    assert result.class_variances.standard_deviations.shape == parameters["centroids"].shape
    assert result.image_variance.image_variance_weights.shape == (2,)
    assert result.centroid_similarity.similarity.shape == (11, 11)
    assert np.isfinite(result.class_variances.standard_deviations).all()
    assert np.isfinite(result.image_variance.image_variance_weights).all()
    assert np.isfinite(result.centroid_similarity.similarity).all()

    direct = toads_default_atlas_membership_pass(
        image_values=data["images"],
        initial_memberships=data["initial_memberships"],
        centroids=parameters["centroids"],
        exponent_by_channel_class=parameters["exponent_by_channel_class"],
        modality_weights=parameters["modality_weights"],
        class_variance_weights=parameters["class_variance_weights"],
        image_variances=result.image_variance.image_variance_weights,
        transformed_shapes=data["transformed_shapes"],
        closest_labels=data["closest_labels"],
        smoothing_priors=parameters["smoothing_priors"],
        centroid_similarity=result.centroid_similarity.similarity,
        atlas_shape_priors=parameters["atlas_shape_priors"],
        smoothing=0.043,
        prior_factor=parameters["prior_factor"],
        prior_coefficient=0.16,
    )
    np.testing.assert_array_equal(result.membership_pass.class_data_energy, direct.class_data_energy)
    np.testing.assert_array_equal(result.membership_pass.memberships, direct.memberships)
    np.testing.assert_array_equal(result.membership_pass.denominator, direct.denominator)
    np.testing.assert_array_equal(result.membership_pass.objective, direct.objective)


def test_initial_statistics_membership_rejects_nondefault_class_order() -> None:
    with pytest.raises(ValueError, match="default 11-class"):
        compute_toads_initial_statistics_membership(
            image_values=np.ones((3, 3, 3, 1), dtype=np.float32),
            initial_memberships=np.ones((3, 3, 3, 2), dtype=np.float32),
            centroids=np.ones((1, 2), dtype=np.float32),
            segmentation=np.zeros((3, 3, 3), dtype=np.int8),
            class_types=("MASK", "OPTIMIZED_GM"),
            class_variance_weights=np.ones((1, 2), dtype=np.float32),
            exponent_by_channel_class=np.ones((1, 2), dtype=np.float32),
            modality_weights=np.ones(1, dtype=np.float32),
            transformed_shapes=np.ones((3, 3, 3, 2), dtype=np.float32),
            closest_labels=np.zeros((3, 3, 3, 4), dtype=np.int16),
            smoothing_priors=np.ones(2, dtype=np.float32),
            atlas_shape_priors=np.ones(2, dtype=np.float32),
            smoothing=0.0,
            prior_factor=0.0,
            prior_coefficient=0.0,
        )
