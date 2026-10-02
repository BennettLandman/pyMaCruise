from __future__ import annotations

import numpy as np

from pymacruise.toads_default_atlas_membership import toads_default_atlas_membership_pass
from pymacruise.toads_default_exact_centroids import toads_default_exact_centroid_update
from pymacruise.toads_default_exact_variances import toads_default_exact_class_variances
from pymacruise.toads_math import toads_optimized_centroid_similarity, toads_optimized_image_variance
from pymacruise.toads_initial_relation_state import toads_initialize_relation_state
from pymacruise.toads_multichannel_carried_degree_iteration import (
    continue_toads_multichannel_through_carried_degree_iteration,
)
from test_toads_post_thinning_exact_field_ramp import _java_fit, _load_inputs, _next_thinning_with_fixture


def _case():
    image, memberships, labels, field, centers = _load_inputs()
    next_thinning = _next_thinning_with_fixture(current_degree=2)
    images = np.stack((image, image * np.float32(1.07) + np.float32(0.03)), axis=-1)
    fields = np.stack((field, field * np.float32(0.97) + np.float32(0.04)), axis=-1).astype(np.float32)
    centroids = np.stack((centers, centers * np.float32(1.03) + np.float32(0.01)))
    categories = np.asarray([[0.8, 1.0, 1.0], [0.9, 0.95, 1.05]], dtype=np.float32)
    exponents = np.empty((2, 11), dtype=np.float32)
    exponents[:, 0] = 1.0
    exponents[:, 1:3] = categories[:, 0, None]
    exponents[:, 3:8] = categories[:, 1, None]
    exponents[:, 8:11] = categories[:, 2, None]
    memberships = np.asarray(
        next_thinning.branch.handoff.membership_pass.memberships, dtype=np.float32
    )
    shapes = np.zeros_like(memberships)
    closest = toads_initialize_relation_state(labels, memberships, 1).closest_labels
    settings = dict(
        previous_degree=2,
        polynomial_degree=3,
        low_energy_ratio=0.1,
        category_exponents=categories,
        exponent_by_channel_class=exponents,
        modality_weights=np.asarray([0.73, 1.15], dtype=np.float32),
        class_variance_weights=np.asarray(
            [
                [0.81, 0.93, 1.17, 0.76, 0.88, 1.02, 1.14, 0.69, 0.96, 1.21, 0.73],
                [1.11, 0.72, 0.95, 1.23, 0.87, 0.74, 1.06, 0.91, 1.18, 0.67, 1.02],
            ],
            dtype=np.float32,
        ),
        transformed_shapes=shapes,
        closest_labels=closest,
        smoothing_priors=np.asarray(
            [0.65, 0.71, 0.84, 0.91, 1.02, 1.13, 0.76, 0.88, 1.05, 1.17, 0.82],
            dtype=np.float32,
        ),
        atlas_shape_priors=np.asarray(
            [0.10, 0.08, 0.09, 0.11, 0.12, 0.09, 0.10, 0.11, 0.07, 0.06, 0.07],
            dtype=np.float32,
        ),
        smoothing=0.043,
        prior_factor=0.64,
        prior_coefficient=0.16,
        prior_scale=0.05,
    )
    return next_thinning, images, fields, centroids, labels, memberships, settings


def test_two_channel_carried_degree_refits_each_field_and_matches_stage_composition(tmp_path):
    next_thinning, images, fields, centers, labels, memberships, settings = _case()
    initial_fields = fields.copy()
    result = continue_toads_multichannel_through_carried_degree_iteration(
        next_thinning,
        images,
        centers,
        fields,
        **settings,
    )

    assert result.previous_degree == 2 and result.degree == 3
    assert result.field_refit is not None
    assert result.field_refit.channel_order == (0, 1)
    assert result.stage_order[-6:] == (
        "next_iteration.computeExact3DPolynomialField.degree3",
        "next_iteration.computeExactCentroids",
        "next_iteration.computeExactVariances",
        "next_iteration.computeImageVariance",
        "next_iteration.computeCentroidSimilarity",
        "next_iteration.computeExactMemberships",
    )
    for channel in range(2):
        channel_dir = tmp_path / f"channel{channel}"
        channel_dir.mkdir()
        java = _java_fit(
            channel_dir,
            3,
            labels,
            images[..., channel],
            memberships,
            initial_fields[..., channel],
            centers[channel],
        )
        expected = np.asarray(java["field_bits"], dtype=np.int32).view(np.float32).reshape(labels.shape)
        np.testing.assert_allclose(result.field_refit.fields[..., channel], expected, rtol=0, atol=2e-6)
    np.testing.assert_array_equal(fields, result.field_refit.fields)

    expected_centers = toads_default_exact_centroid_update(
        images, memberships, centers, settings["category_exponents"],
        intensity_field=result.field_refit.fields,
    )
    expected_deviations = toads_default_exact_class_variances(
        images, memberships, expected_centers, labels,
        intensity_field=result.field_refit.fields,
    )
    expected_variance, expected_std = toads_optimized_image_variance(
        images, memberships[..., 1:], expected_centers[:, 1:],
        intensity_field=result.field_refit.fields,
    )
    expected_similarity = toads_optimized_centroid_similarity(
        expected_centers, expected_deviations, expected_variance, settings["prior_scale"]
    )
    expected_membership = toads_default_atlas_membership_pass(
        image_values=images,
        initial_memberships=memberships,
        centroids=expected_centers,
        exponent_by_channel_class=settings["exponent_by_channel_class"],
        modality_weights=settings["modality_weights"],
        class_variance_weights=settings["class_variance_weights"],
        image_variances=expected_variance,
        transformed_shapes=settings["transformed_shapes"],
        closest_labels=settings["closest_labels"],
        smoothing_priors=settings["smoothing_priors"],
        centroid_similarity=expected_similarity,
        atlas_shape_priors=settings["atlas_shape_priors"],
        smoothing=settings["smoothing"],
        prior_factor=settings["prior_factor"],
        prior_coefficient=settings["prior_coefficient"],
        intensity_field=result.field_refit.fields,
    )
    np.testing.assert_array_equal(result.centroids, expected_centers)
    np.testing.assert_array_equal(result.standard_deviations, expected_deviations)
    np.testing.assert_array_equal(result.image_variance_weights, expected_variance)
    np.testing.assert_array_equal(result.image_standard_deviation, expected_std)
    np.testing.assert_array_equal(result.centroid_similarity, expected_similarity)
    np.testing.assert_array_equal(result.membership_pass.memberships, expected_membership.memberships)
    np.testing.assert_array_equal(result.membership_pass.denominator, expected_membership.denominator)
    assert result.membership_pass.objective == expected_membership.objective


def test_multichannel_degree_zero_keeps_attached_fields_and_updates_statistics():
    next_thinning, images, fields, centers, _labels, _memberships, settings = _case()
    settings["previous_degree"] = 0
    settings["low_energy_ratio"] = 0.0
    original_fields = fields.copy()
    result = continue_toads_multichannel_through_carried_degree_iteration(
        next_thinning,
        images,
        centers,
        fields,
        **settings,
    )
    assert result.degree == 0
    assert result.field_refit is None
    assert result.stage_order[-6] == "next_iteration.computeExact3DPolynomialField.degree0.noop"
    np.testing.assert_array_equal(fields, original_fields)
