from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_alignment_startup import register_toads_from_atlas_assets
from pymacruise.toads_default_atlas_membership import toads_default_atlas_membership_pass
from pymacruise.toads_default_exact_centroids import toads_default_exact_centroid_update
from pymacruise.toads_default_exact_variances import toads_default_exact_class_variances
from pymacruise.toads_math import toads_optimized_centroid_similarity, toads_optimized_image_variance
from pymacruise.toads_registered_first_exact_iteration import (
    continue_registered_unity_startup_through_first_exact_iteration,
)
from pymacruise.toads_registered_first_iteration_cycle import (
    continue_registered_first_exact_iteration_through_next_thinning,
)
from pymacruise.toads_registered_membership_branch_cycle import (
    continue_registered_membership_through_branch,
)
from pymacruise.toads_registered_following_exact_iteration import (
    continue_registered_cycle_through_following_exact_iteration,
)
from pymacruise.toads_inhomogeneity_transition import toads_next_inhomogeneity_degree
from pymacruise.toads_registered_second_exact_iteration import (
    continue_registered_first_cycle_through_second_exact_iteration,
)
from pymacruise.toads_registered_unity_field_startup import (
    compute_registered_toads_unity_field_startup,
)
from pymacruise.toads_shape_volume import transform_toads_shape_priors
from pymacruise.toads_control import toads_energy_ratio
from pymacruise.toads_relations import toads_refresh_relation_state
from pymacruise.toads_thinning import toads_propagate_growing, toads_propagate_thinning


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "toads_registered_first_iteration_cycle"
ATLAS_DESCRIPTOR = ROOT / "data" / "synthetic" / "toads_registered_field_refit" / "atlas-11.txt"


def test_registered_first_exact_loop_reaches_iteration_two_memberships() -> None:
    run = register_toads_from_atlas_assets(
        nib.load(str(FIXTURE / "t1.nii.gz")),
        nib.load(str(FIXTURE / "brain-mask.nii.gz")),
        ATLAS_DESCRIPTOR,
        levels=1,
        max_iterations=1,
        refinement_iterations=1,
        padding=2,
    )
    settings = dict(
        exponent_by_channel_class=np.asarray([[1.0, *([2.0] * 10)]], dtype=np.float32),
        modality_weights=np.ones(1, dtype=np.float32),
        smoothing_priors=np.ones(11, dtype=np.float32),
        atlas_shape_priors=np.ones(11, dtype=np.float32),
        smoothing=np.float32(0.2 / 6.0 / 11.0),
        prior_factor=np.float32(0.2),
        prior_coefficient=np.float32(2.0 / 11.0),
        prior_scale=np.float32(0.05),
    )
    startup = compute_registered_toads_unity_field_startup(run, **settings)
    exact = continue_registered_unity_startup_through_first_exact_iteration(
        run,
        startup,
        category_exponents=np.asarray([[2.0, 2.0, 2.0]], dtype=np.float32),
        class_variance_weights=run.initial_variances.standard_deviations,
        **settings,
    )

    result = continue_registered_first_exact_iteration_through_next_thinning(
        run,
        exact,
        maximum_iterations=2,
        convergence_threshold=0.04,
        relation_limit=1,
        do_growing=False,
        subsample=2,
        max_registration_iterations=1,
        registration_energy_tolerance=0.0,
        minimum_lambda=1e-3,
    )

    assert result.exact_iteration is exact
    assert result.registration_control.handoff.membership_pass is exact.membership_pass
    assert result.registration_control.handoff.registration is not None
    expected_ratio = toads_energy_ratio(
        startup.statistics_membership.membership_pass.objective,
        exact.membership_pass.objective,
    )
    assert result.registration_control.control.energy_ratio == expected_ratio
    assert result.registered_iteration.branch.control.energy_ratio == expected_ratio
    assert not result.registered_iteration.branch.control.stop
    assert result.registered_iteration.branch.growing is not None
    assert result.registered_iteration.branch.relations is not None
    assert result.registered_iteration.branch.classification is None
    assert result.next_thinning.thinning.segmentation.shape == exact.membership_pass.memberships.shape[:3]
    assert result.next_transformed_shape_priors.shape == exact.membership_pass.memberships.shape
    assert result.stage_order[-11:] == (
        "iteration1.computeExactMemberships",
        "energyRatio",
        "atlas.registerShapes() -> registerAllShapes()",
        "convergenceCheck",
        "iterationLimitCheck",
        "propagateGrowing",
        "computeRelations",
        "computeBoundaryDistances",
        "computeClosestLabels",
        "iteration2.initThinningLabels",
        "iteration2.propagateThinning",
    )

    handoff = result.registration_control.handoff
    expected_shapes = transform_toads_shape_priors(
        np.stack(
            tuple(run.atlas_data.shapes[item.name] for item in run.atlas.structures),
            axis=-1,
        ),
        image_dimensions_xyz=exact.membership_pass.memberships.shape[:3],
        image_resolutions_mm=run.alignment.setup.signed_image_resolutions_mm,
        template_resolutions_mm=run.alignment.setup.template_resolutions_mm,
        transform=tuple(float(value) for value in handoff.transform),
        class_bounds=handoff.shape_bounds,
    )
    np.testing.assert_array_equal(result.next_transformed_shape_priors, expected_shapes)
    direct_growing = toads_propagate_growing(
        startup.thinning.segmentation,
        startup.thinning.skeleton,
        exact.membership_pass.memberships,
        np.zeros(startup.thinning.segmentation.shape, dtype=np.float32),
        startup.thinning.class_volumes,
        initial_available=startup.thinning.available,
    )
    np.testing.assert_array_equal(
        result.registered_iteration.branch.growing.segmentation,
        direct_growing.segmentation,
    )
    np.testing.assert_array_equal(
        result.registered_iteration.branch.growing.ordering,
        direct_growing.ordering,
    )
    direct_relations = toads_refresh_relation_state(
        direct_growing.segmentation,
        exact.membership_pass.memberships,
        run.relation_state.boundary_distances,
        1,
    )
    np.testing.assert_array_equal(
        result.registered_iteration.branch.relations.relation.relation_map,
        direct_relations.relation.relation_map,
    )
    np.testing.assert_array_equal(
        result.registered_iteration.branch.relations.closest_labels,
        direct_relations.closest_labels,
    )
    expected_thinning = toads_propagate_thinning(
        direct_growing.segmentation,
        exact.membership_pass.memberships,
        direct_relations.relation.class_adjacency,
        expected_shapes,
        volume_ratio=0.5,
        limit=3.0,
    )
    np.testing.assert_array_equal(result.next_thinning.thinning.segmentation, expected_thinning.segmentation)
    np.testing.assert_array_equal(result.next_thinning.thinning.skeleton, expected_thinning.skeleton)
    assert result.next_thinning.thinning.processed_items == expected_thinning.processed_items

    second_exponents = settings["exponent_by_channel_class"]
    categories = np.asarray([[2.0, 2.0, 2.0]], dtype=np.float32)
    second = continue_registered_first_cycle_through_second_exact_iteration(
        run,
        result,
        category_exponents=categories,
        exponent_by_channel_class=second_exponents,
        modality_weights=settings["modality_weights"],
        class_variance_weights=run.initial_variances.standard_deviations,
        smoothing_priors=settings["smoothing_priors"],
        atlas_shape_priors=settings["atlas_shape_priors"],
        smoothing=settings["smoothing"],
        prior_factor=settings["prior_factor"],
        prior_coefficient=settings["prior_coefficient"],
        prior_scale=settings["prior_scale"],
        maximum_distance=0.04,
    )
    assert second.next_thinning is result.next_thinning
    assert second.previous_degree == second.degree == 0
    assert second.field_refit is None
    assert second.energy_ratio == result.registration_control.control.energy_ratio
    assert second.stage_order[-6:] == (
        "next_iteration.computeExact3DPolynomialField.degree0.noop",
        "next_iteration.computeExactCentroids",
        "next_iteration.computeExactVariances",
        "next_iteration.computeImageVariance",
        "next_iteration.computeCentroidSimilarity",
        "next_iteration.computeExactMemberships",
    )

    second_images = np.moveaxis(run.template_refresh.recropped_images, 0, -1)
    second_memberships = result.registered_iteration.branch.handoff.membership_pass.memberships
    second_labels = result.next_thinning.thinning.segmentation
    second_centroids = toads_default_exact_centroid_update(
        second_images,
        second_memberships,
        exact.centroids,
        categories,
        intensity_field=exact.fields,
        smoothness=0.1,
    )
    second_deviations = toads_default_exact_class_variances(
        second_images,
        second_memberships,
        second_centroids,
        second_labels,
        intensity_field=exact.fields,
    )
    second_variance, second_std = toads_optimized_image_variance(
        second_images,
        second_memberships[..., 1:],
        second_centroids[:, 1:],
        intensity_field=exact.fields,
    )
    second_similarity = toads_optimized_centroid_similarity(
        second_centroids, second_deviations, second_variance, settings["prior_scale"]
    )
    second_membership_pass = toads_default_atlas_membership_pass(
        image_values=second_images,
        initial_memberships=second_memberships,
        centroids=second_centroids,
        exponent_by_channel_class=second_exponents,
        modality_weights=settings["modality_weights"],
        class_variance_weights=run.initial_variances.standard_deviations,
        image_variances=second_variance,
        transformed_shapes=result.next_transformed_shape_priors,
        closest_labels=result.registered_iteration.branch.relations.closest_labels,
        smoothing_priors=settings["smoothing_priors"],
        centroid_similarity=second_similarity,
        atlas_shape_priors=settings["atlas_shape_priors"],
        smoothing=settings["smoothing"],
        prior_factor=settings["prior_factor"],
        prior_coefficient=settings["prior_coefficient"],
        intensity_field=exact.fields,
    )
    np.testing.assert_array_equal(second.centroids, second_centroids)
    np.testing.assert_array_equal(second.standard_deviations, second_deviations)
    np.testing.assert_array_equal(second.image_variance_weights, second_variance)
    np.testing.assert_array_equal(second.image_standard_deviation, second_std)
    np.testing.assert_array_equal(second.centroid_similarity, second_similarity)
    np.testing.assert_array_equal(second.membership_pass.memberships, second_membership_pass.memberships)
    np.testing.assert_array_equal(second.membership_pass.denominator, second_membership_pass.denominator)
    assert second.membership_pass.objective == second_membership_pass.objective

    correction_disabled = continue_registered_first_cycle_through_second_exact_iteration(
        run,
        result,
        category_exponents=categories,
        exponent_by_channel_class=second_exponents,
        modality_weights=settings["modality_weights"],
        class_variance_weights=run.initial_variances.standard_deviations,
        smoothing_priors=settings["smoothing_priors"],
        atlas_shape_priors=settings["atlas_shape_priors"],
        smoothing=settings["smoothing"],
        prior_factor=settings["prior_factor"],
        prior_coefficient=settings["prior_coefficient"],
        prior_scale=settings["prior_scale"],
        maximum_distance=0.04,
        correct_inhomogeneity=False,
    )
    assert correction_disabled.degree == 0
    assert correction_disabled.field_refit is None
    assert correction_disabled.stage_order[-6] == (
        "next_iteration.skipCorrectionField.inhomogeneityDisabled"
    )
    np.testing.assert_array_equal(exact.fields, np.ones_like(exact.fields))
    np.testing.assert_array_equal(
        correction_disabled.membership_pass.memberships,
        second.membership_pass.memberships,
    )

    prior_branch = result.registered_iteration.branch
    prior_growing = prior_branch.growing
    prior_relations = prior_branch.relations
    assert prior_growing is not None and prior_relations is not None
    second_cycle = continue_registered_membership_through_branch(
        run,
        second.membership_pass,
        previous_energy=exact.membership_pass.objective,
        iteration=2,
        maximum_iterations=3,
        convergence_threshold=0.04,
        segmentation=result.next_thinning.thinning.segmentation,
        skeleton=result.next_thinning.thinning.skeleton,
        ordering=prior_growing.ordering,
        class_volumes=result.next_thinning.thinning.class_volumes,
        boundary_distances=prior_relations.boundary_distances,
        template_labels=np.asarray(
            [item.label for item in run.atlas.structures], dtype=np.int8
        ),
        previous_classification=run.segmentation_constructor.classification_labels,
        relation_limit=1,
        do_growing=False,
        initial_transform=result.registration_control.handoff.transform,
        initial_available=result.next_thinning.thinning.available,
        subsample=2,
        max_registration_iterations=1,
        registration_energy_tolerance=0.0,
        minimum_lambda=1e-3,
    )
    assert second_cycle.membership_pass is second.membership_pass
    assert second_cycle.registration_control.handoff.registration is not None
    assert second_cycle.registration_control.control.energy_ratio == toads_energy_ratio(
        exact.membership_pass.objective,
        second.membership_pass.objective,
    )
    if second_cycle.registered_iteration.branch.control.stop:
        assert second_cycle.next_thinning is None
        assert second_cycle.registered_iteration.branch.classification is not None
    else:
        assert second_cycle.next_thinning is not None
        assert second_cycle.registered_iteration.branch.growing is not None
        assert second_cycle.registered_iteration.branch.relations is not None
        direct_next = toads_propagate_thinning(
            second_cycle.registered_iteration.branch.growing.segmentation,
            second.membership_pass.memberships,
            second_cycle.registered_iteration.branch.relations.relation.class_adjacency,
            second_cycle.next_transformed_shape_priors,
            volume_ratio=0.5,
            limit=3.0,
        )
        np.testing.assert_array_equal(
            second_cycle.next_thinning.thinning.segmentation,
            direct_next.segmentation,
        )
    assert second_cycle.next_thinning is not None, "synthetic loop should continue into iteration three"
    third = continue_registered_cycle_through_following_exact_iteration(
        run,
        second_cycle,
        previous_degree=second.degree,
        initial_centroids=second.centroids,
        fields=exact.fields,
        category_exponents=categories,
        exponent_by_channel_class=second_exponents,
        modality_weights=settings["modality_weights"],
        class_variance_weights=run.initial_variances.standard_deviations,
        smoothing_priors=settings["smoothing_priors"],
        atlas_shape_priors=settings["atlas_shape_priors"],
        smoothing=settings["smoothing"],
        prior_factor=settings["prior_factor"],
        prior_coefficient=settings["prior_coefficient"],
        prior_scale=settings["prior_scale"],
        maximum_distance=0.04,
    )
    expected_degree = toads_next_inhomogeneity_degree(
        second.degree,
        3,
        second_cycle.registration_control.control.energy_ratio,
        np.sqrt(np.float32(0.04)),
    )
    assert third.degree == expected_degree
    assert third.energy_ratio == second_cycle.registration_control.control.energy_ratio
    assert third.next_thinning is second_cycle.next_thinning
    assert (third.field_refit is None) == (third.degree == 0)
    np.testing.assert_array_equal(exact.fields, third.field_refit.fields if third.field_refit else exact.fields)

    third_images = np.moveaxis(run.template_refresh.recropped_images, 0, -1)
    third_memberships = second_cycle.membership_pass.memberships
    third_labels = second_cycle.next_thinning.thinning.segmentation
    third_centroids = toads_default_exact_centroid_update(
        third_images,
        third_memberships,
        second.centroids,
        categories,
        intensity_field=exact.fields,
        smoothness=0.1,
    )
    third_deviations = toads_default_exact_class_variances(
        third_images,
        third_memberships,
        third_centroids,
        third_labels,
        intensity_field=exact.fields,
    )
    third_variance, third_std = toads_optimized_image_variance(
        third_images,
        third_memberships[..., 1:],
        third_centroids[:, 1:],
        intensity_field=exact.fields,
    )
    third_similarity = toads_optimized_centroid_similarity(
        third_centroids, third_deviations, third_variance, settings["prior_scale"]
    )
    third_membership_pass = toads_default_atlas_membership_pass(
        image_values=third_images,
        initial_memberships=third_memberships,
        centroids=third_centroids,
        exponent_by_channel_class=second_exponents,
        modality_weights=settings["modality_weights"],
        class_variance_weights=run.initial_variances.standard_deviations,
        image_variances=third_variance,
        transformed_shapes=second_cycle.next_transformed_shape_priors,
        closest_labels=second_cycle.registered_iteration.branch.relations.closest_labels,
        smoothing_priors=settings["smoothing_priors"],
        centroid_similarity=third_similarity,
        atlas_shape_priors=settings["atlas_shape_priors"],
        smoothing=settings["smoothing"],
        prior_factor=settings["prior_factor"],
        prior_coefficient=settings["prior_coefficient"],
        intensity_field=exact.fields,
    )
    np.testing.assert_array_equal(third.centroids, third_centroids)
    np.testing.assert_array_equal(third.standard_deviations, third_deviations)
    np.testing.assert_array_equal(third.image_variance_weights, third_variance)
    np.testing.assert_array_equal(third.image_standard_deviation, third_std)
    np.testing.assert_array_equal(third.centroid_similarity, third_similarity)
    np.testing.assert_array_equal(third.membership_pass.memberships, third_membership_pass.memberships)
    np.testing.assert_array_equal(third.membership_pass.denominator, third_membership_pass.denominator)
    assert third.membership_pass.objective == third_membership_pass.objective

    second_branch = second_cycle.registered_iteration.branch
    second_growing = second_branch.growing
    second_relations = second_branch.relations
    assert second_growing is not None and second_relations is not None
    third_cycle = continue_registered_membership_through_branch(
        run,
        third.membership_pass,
        previous_energy=second.membership_pass.objective,
        iteration=3,
        maximum_iterations=4,
        convergence_threshold=0.04,
        segmentation=second_cycle.next_thinning.thinning.segmentation,
        skeleton=second_cycle.next_thinning.thinning.skeleton,
        ordering=second_growing.ordering,
        class_volumes=second_cycle.next_thinning.thinning.class_volumes,
        boundary_distances=second_relations.boundary_distances,
        template_labels=np.asarray(
            [item.label for item in run.atlas.structures], dtype=np.int8
        ),
        previous_classification=run.segmentation_constructor.classification_labels,
        relation_limit=1,
        do_growing=False,
        initial_transform=second_cycle.registration_control.handoff.transform,
        initial_available=second_cycle.next_thinning.thinning.available,
        subsample=2,
        max_registration_iterations=1,
        registration_energy_tolerance=0.0,
        minimum_lambda=1e-3,
    )
    assert third_cycle.membership_pass is third.membership_pass
    assert third_cycle.registration_control.handoff.registration is not None
    assert third_cycle.registration_control.control.energy_ratio == toads_energy_ratio(
        second.membership_pass.objective,
        third.membership_pass.objective,
    )
    if third_cycle.registered_iteration.branch.control.stop:
        assert third_cycle.next_thinning is None
        assert third_cycle.registered_iteration.branch.classification is not None
    else:
        assert third_cycle.next_thinning is not None
        assert third_cycle.registered_iteration.branch.growing is not None
        assert third_cycle.registered_iteration.branch.relations is not None
        direct_third_next = toads_propagate_thinning(
            third_cycle.registered_iteration.branch.growing.segmentation,
            third.membership_pass.memberships,
            third_cycle.registered_iteration.branch.relations.relation.class_adjacency,
            third_cycle.next_transformed_shape_priors,
            volume_ratio=0.5,
            limit=3.0,
        )
        np.testing.assert_array_equal(
            third_cycle.next_thinning.thinning.segmentation,
            direct_third_next.segmentation,
        )
