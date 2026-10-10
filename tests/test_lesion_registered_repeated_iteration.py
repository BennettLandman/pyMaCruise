from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_descriptor_continuing_iteration import (
    run_lesion_toads_descriptor_continuing_iteration,
)
from pymacruise.lesion_descriptor_loop import (
    LesionToadsDescriptorLoopOptions,
    run_lesion_toads_descriptor_loop,
)
from pymacruise.lesion_descriptor_registered_iteration import (
    register_and_branch_lesion_toads_descriptor_iteration,
)
from pymacruise.lesion_registered_first_iteration import (
    run_lesion_toads_registered_first_iteration,
)
from pymacruise.lesion_registered_first_iteration_branch import (
    run_lesion_toads_registered_first_iteration_branch,
)
from pymacruise.lesion_registered_initial_membership import (
    run_lesion_toads_registered_initial_membership,
)
from pymacruise.lesion_registered_shapes import transform_registered_lesion_shape_priors
from pymacruise.lesion_standard13_voxel import STANDARD13_TYPES


FIXTURE = Path(__file__).parent / "data" / "synthetic" / "lesion_registered_initial_membership"


def test_registered_loop_driver_matches_two_manual_updates() -> None:
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
    shape_priors = np.ones(13, dtype=np.float32)
    smoothing_priors = np.ones(13, dtype=np.float32)
    centroid_options = {
        "mode": "none",
        "modalities": ("T1_MPRAGE", "FLAIR"),
        "lesion_power_exponents": (0.4, 1.2),
        "tissue_power_exponents": np.ones((2, 13), dtype=np.float32),
    }
    channel_weight_options = {"outlier_distance": 0.1, "norm_type": "Identity"}
    first = run_lesion_toads_registered_first_iteration(
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
    first_branch = run_lesion_toads_registered_first_iteration_branch(
        initial,
        first,
        tuple(shapes[..., class_index] for class_index in range(13)),
        np.zeros(6, dtype=np.float32),
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        object_types=STANDARD13_TYPES,
        classification=labels,
        template_labels=np.arange(13, dtype=np.int16),
        lesion_class=10,
        iteration=1,
        n_iter_max=99,
        max_distance=0.0,
        algorithm_mode="diffeomorphism",
        relation_limit=3,
        align_shapes=True,
        init_align_iter=1,
    )
    prior_branch = first_branch.branch.branch
    assert prior_branch.stop is False
    assert prior_branch.growing is not None
    assert prior_branch.relations is not None

    atlas_shape_priors = tuple(shapes[..., class_index] for class_index in range(13))
    transformed_shapes = transform_registered_lesion_shape_priors(
        first_branch.registration.registration,
        atlas_shape_priors,
        image_dimensions_xyz=tuple(images.shape[:3]),
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
    )
    # The constructor distance maps remain the source's carried maps; the
    # registered shape priors above are refreshed at the new transform.
    relation_table = prior_branch.relations.relation_membership[
        prior_branch.relations.relation_map
    ]
    second = run_lesion_toads_descriptor_continuing_iteration(
        prior_branch,
        images,
        transformed_shapes,
        object_types=STANDARD13_TYPES,
        lesion_class=10,
        volume_ratio=0.5,
        thinning_limit=1000.0,
        initial_class_centroids=first.statistics.statistics.centroids.class_centroids,
        initial_lesion_centroids=first.statistics.statistics.centroids.lesion_centroids,
        initial_black_hole_centroids=first.statistics.statistics.centroids.black_hole_centroids,
        shape_priors=shape_priors,
        smoothing_priors=smoothing_priors,
        direct_relations=relation_table,
        ventricle_distances=initial.distance_maps.ventricle,
        gray_matter_distances=initial.distance_maps.gray_matter,
        brainstem_distances=initial.distance_maps.brainstem,
        centroid_options=centroid_options,
        channel_weight_options=channel_weight_options,
        outlier_energy=0.1,
        lesion_power_exponent=1.2,
        source_zero_variance_nan=False,
    )

    assert second.stage_order == (
        "propagateThinning",
        "computeCentroids",
        "computeVariances",
        "computeChannelWeights",
        "computeMemberships:GMflag=true",
        "energyRatio",
    )
    assert second.thinning.run.segmentation.shape == images.shape[:3]
    assert transformed_shapes.shape == (*images.shape[:3], 13)
    assert second.membership.membership.class_memberships.shape == (*images.shape[:3], 13)
    assert np.isfinite(second.membership.membership.class_memberships).all()
    assert np.isfinite(second.membership.energy_ratio.energy_ratio)
    assert second.thinning is not initial.first_thinning.thinning

    second_branch = register_and_branch_lesion_toads_descriptor_iteration(
        second,
        atlas_shape_priors,
        first_branch.registration.registration.transform,
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        object_types=STANDARD13_TYPES,
        classification=labels,
        template_labels=np.arange(13, dtype=np.int16),
        lesion_class=10,
        iteration=2,
        n_iter_max=99,
        max_distance=0.0,
        algorithm_mode="diffeomorphism",
        relation_limit=3,
        align_shapes=True,
        init_align_iter=1,
    )
    assert second_branch.control.branch.stop is False
    assert second_branch.control.branch.growing is not None
    assert second_branch.control.branch.relations is not None
    assert second_branch.stage_order.index("registerShapes") < second_branch.stage_order.index(
        "convergence_check"
    )

    third_shapes = transform_registered_lesion_shape_priors(
        second_branch.registration.registration,
        atlas_shape_priors,
        image_dimensions_xyz=tuple(images.shape[:3]),
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
    )
    second_prior = second_branch.control.branch
    second_relation_table = second_prior.relations.relation_membership[
        second_prior.relations.relation_map
    ]
    third = run_lesion_toads_descriptor_continuing_iteration(
        second_prior,
        images,
        third_shapes,
        object_types=STANDARD13_TYPES,
        lesion_class=10,
        volume_ratio=0.5,
        thinning_limit=1000.0,
        initial_class_centroids=second.statistics.statistics.centroids.class_centroids,
        initial_lesion_centroids=second.statistics.statistics.centroids.lesion_centroids,
        initial_black_hole_centroids=second.statistics.statistics.centroids.black_hole_centroids,
        shape_priors=shape_priors,
        smoothing_priors=smoothing_priors,
        direct_relations=second_relation_table,
        ventricle_distances=initial.distance_maps.ventricle,
        gray_matter_distances=initial.distance_maps.gray_matter,
        brainstem_distances=initial.distance_maps.brainstem,
        centroid_options=centroid_options,
        channel_weight_options=channel_weight_options,
        outlier_energy=0.1,
        lesion_power_exponent=1.2,
        source_zero_variance_nan=False,
    )
    assert third.stage_order == second.stage_order
    assert third.thinning.run.segmentation.shape == images.shape[:3]
    assert third.membership.membership.class_memberships.shape == (*images.shape[:3], 13)
    assert np.isfinite(third.membership.membership.class_memberships).all()
    assert np.isfinite(third.membership.energy_ratio.energy_ratio)

    third_branch = register_and_branch_lesion_toads_descriptor_iteration(
        third,
        atlas_shape_priors,
        second_branch.registration.registration.transform,
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        object_types=STANDARD13_TYPES,
        classification=labels,
        template_labels=np.arange(13, dtype=np.int16),
        lesion_class=10,
        iteration=3,
        n_iter_max=3,
        max_distance=0.0,
        algorithm_mode="diffeomorphism",
        relation_limit=3,
        align_shapes=True,
        init_align_iter=1,
    )
    assert third_branch.control.branch.stop is True
    assert third_branch.control.branch.iteration_limit_stop is True
    assert third_branch.control.branch.final_classification is not None
    assert third_branch.control.branch.final_classification.classification.shape == images.shape[:3]
    assert third_branch.stage_order.index("registerShapes") < third_branch.stage_order.index(
        "convergence_check"
    )

    loop_options = LesionToadsDescriptorLoopOptions(
            images=images,
            atlas_shape_priors=atlas_shape_priors,
            shape_priors=shape_priors,
            smoothing_priors=smoothing_priors,
            object_types=STANDARD13_TYPES,
            lesion_class=10,
            volume_ratio=0.5,
            thinning_limit=1000.0,
            image_resolutions_mm=(1.0, 1.0, 1.0),
            template_resolutions_mm=(1.0, 1.0, 1.0),
            classification=labels,
            template_labels=np.arange(13, dtype=np.int16),
            n_iter_max=3,
            max_distance=0.0,
            algorithm_mode="diffeomorphism",
            relation_limit=3,
            initial_class_centroids=first.statistics.statistics.centroids.class_centroids,
            initial_lesion_centroids=first.statistics.statistics.centroids.lesion_centroids,
            initial_black_hole_centroids=first.statistics.statistics.centroids.black_hole_centroids,
            ventricle_distances=initial.distance_maps.ventricle,
            gray_matter_distances=initial.distance_maps.gray_matter,
            brainstem_distances=initial.distance_maps.brainstem,
            centroid_options=centroid_options,
            channel_weight_options=channel_weight_options,
            outlier_energy=0.1,
            lesion_power_exponent=1.2,
            source_zero_variance_nan=False,
            align_shapes=True,
            init_align_iter=1,
        )
    loop_progress = []
    loop = run_lesion_toads_descriptor_loop(
        prior_branch,
        first_branch.registration.registration,
        starting_iteration=1,
        options=loop_options,
        progress_callback=loop_progress.append,
    )
    assert tuple(item.iteration for item in loop.passes) == (2, 3)
    assert loop.stop_reason == "iteration_limit"
    assert [item["iteration"] for item in loop_progress] == [2, 3]
    assert all(np.isfinite(item["energy_ratio"]) for item in loop_progress)
    assert [item["stop"] for item in loop_progress] == [False, True]
    assert all(item["stage_order"] for item in loop_progress)
    assert all(item["elapsed_seconds"] >= 0 for item in loop_progress)
    assert all(item.transition is None and item.transformed_shape_values is None for item in loop.passes)
    assert np.isfinite([item.energy_ratio for item in loop.passes]).all()
    assert np.array_equal(
        loop.passes[0].registration_transform,
        second_branch.registration.registration.transform,
    )
    assert np.array_equal(
        loop.passes[1].registration_transform,
        third_branch.registration.registration.transform,
    )

    detailed_loop = run_lesion_toads_descriptor_loop(
        prior_branch,
        first_branch.registration.registration,
        starting_iteration=1,
        options=replace(loop_options, retain_iteration_details=True),
    )
    assert all(
        item.transition is not None and item.transformed_shape_values is not None
        for item in detailed_loop.passes
    )
    assert np.array_equal(
        detailed_loop.passes[0].transition.membership.membership.class_memberships,
        second.membership.membership.class_memberships,
    )
    assert np.array_equal(
        detailed_loop.passes[1].transition.membership.membership.class_memberships,
        third.membership.membership.class_memberships,
    )
    assert np.array_equal(
        loop.final_branch.final_classification.classification,
        third_branch.control.branch.final_classification.classification,
    )

    # The correction field is supplied in the same channel-last form used by
    # the persisted NIfTI products. It must affect both refreshed statistics
    # and the 13-class membership data terms in the continuing-loop path.
    corrected_fields = np.empty_like(images, dtype=np.float32)
    corrected_fields[..., 0] = np.float32(1.1)
    corrected_fields[..., 1] = np.float32(0.9)
    fielded_loop = run_lesion_toads_descriptor_loop(
        prior_branch,
        first_branch.registration.registration,
        starting_iteration=1,
        options=replace(
            loop_options,
            fields=corrected_fields,
            retain_iteration_details=True,
        ),
        max_passes=1,
    )
    assert fielded_loop.stop_reason == "checkpoint_boundary"
    assert fielded_loop.passes[0].transition is not None
    assert fielded_loop.passes[0].transition.membership.membership.objective.view(np.uint32) != (
        detailed_loop.passes[0].transition.membership.membership.objective.view(np.uint32)
    )

    initial_fields = np.ones_like(images, dtype=np.float32)
    refitted_loop = run_lesion_toads_descriptor_loop(
        prior_branch,
        first_branch.registration.registration,
        starting_iteration=1,
        options=replace(
            loop_options,
            fields=initial_fields,
            correction_polynomial_degree=1,
            initial_correction_degree=1,
            retain_iteration_details=True,
        ),
        max_passes=1,
    )
    refitted_transition = refitted_loop.passes[0].transition
    assert refitted_transition is not None
    assert refitted_transition.inhomogeneity is not None
    assert refitted_transition.inhomogeneity.degree == 1
    assert refitted_transition.stage_order.index("degreeGate") < refitted_transition.stage_order.index(
        "computeCentroids"
    )
    assert refitted_loop.final_correction_degree == 1
    assert refitted_loop.final_fields is not None
    assert np.isfinite(refitted_loop.final_fields).all()
    np.testing.assert_array_equal(initial_fields, np.ones_like(initial_fields))
    np.testing.assert_array_equal(
        refitted_transition.inhomogeneity.fields[0],
        refitted_loop.final_fields[..., 0],
    )

    converged = run_lesion_toads_descriptor_loop(
        prior_branch,
        first_branch.registration.registration,
        starting_iteration=1,
        options=replace(loop_options, max_distance=1.0e6),
    )
    assert len(converged.passes) == 1
    assert converged.stop_reason == "convergence"
    assert converged.final_branch.convergence_stop is True
    assert converged.final_branch.iteration_limit_stop is False
