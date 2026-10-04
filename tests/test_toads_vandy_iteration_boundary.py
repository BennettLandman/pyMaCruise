from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_alignment_startup import register_toads_from_atlas_assets
from pymacruise.toads_atlas import (
    load_toads_atlas,
    load_toads_intensity_variance_factors,
    load_toads_optimized_factors,
)
from pymacruise.toads_default_exact_variances import toads_default_exact_class_variances
from pymacruise.toads_registration import exact_registration_measure_coefficients
from pymacruise.toads_registration_sampling import (
    vectorized_toads_registration_coefficients,
    vectorized_toads_registration_energy,
)
from pymacruise.toads_shape_volume import transform_toads_shape_priors
from pymacruise.toads_atlas import toads_atlas_shapes_in_structure_order
from pymacruise.toads_registered_loop import (
    ToadsRegisteredLoopConfig,
    run_registered_toads_loop,
)


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data/synthetic/toads_downsampled_vandy_atlas"
REFERENCE = FIXTURE / "two_iteration_java_reference.npz"
PROVENANCE = json.loads((FIXTURE / "two_iteration_java_reference.json").read_text())
ATLAS = ROOT.parent / "data/atlases/cruise-atlas-10obj-toads2008.txt"


def _registered_two_iteration_run():
    atlas = load_toads_atlas(ATLAS)
    t1 = nib.load(str(FIXTURE / "t1.nii.gz"))
    mask = nib.load(str(FIXTURE / "brain-mask.nii.gz"))
    alignment = register_toads_from_atlas_assets(
        t1,
        mask,
        atlas,
        levels=1,
        max_iterations=1,
        refinement_iterations=1,
        modality="T1_MPRAGE",
        padding=10,
        relation_limit=3,
        energy_tolerance=1e-6,
        minimum_lambda=1e-3,
    )
    category_exponents = load_toads_optimized_factors(atlas.descriptor, ["T1_MPRAGE"])
    class_indices = {"optimized_csf": 0, "optimized_gm": 1, "optimized_wm": 2}
    exponent_by_class = np.empty((1, len(atlas.structures)), dtype=np.float32)
    for index, structure in enumerate(atlas.structures):
        topology = structure.topology.lower()
        exponent_by_class[0, index] = (
            np.float32(1.0)
            if topology == "mask"
            else category_exponents[0, class_indices[topology]]
        )
    config = ToadsRegisteredLoopConfig(
        category_exponents=category_exponents,
        exponent_by_channel_class=exponent_by_class,
        modality_weights=np.ones(1, dtype=np.float32),
        class_variance_weights=load_toads_intensity_variance_factors(
            atlas.descriptor, ["T1_MPRAGE"]
        ),
        smoothing_priors=np.ones(len(atlas.structures), dtype=np.float32),
        atlas_shape_priors=np.ones(len(atlas.structures), dtype=np.float32),
        smoothing=0.2 / 6.0 / len(atlas.structures),
        prior_factor=0.5,
        prior_coefficient=2.0 / len(atlas.structures),
        maximum_iterations=2,
        convergence_threshold=1e-4,
        relation_limit=3,
        do_growing=True,
        maximum_polynomial_degree=3,
        maximum_distance=1e-4,
        max_registration_iterations=1,
        retain_iteration_history=True,
    )
    return atlas, run_registered_toads_loop(alignment, config)


def test_vandy_iteration_memberships_and_next_thinning_match_java_with_known_drift() -> None:
    """Compare two loop-boundary states to a pinned selected-source Java run.

    The second-update envelope is diagnostic, not a membership-equivalence
    tolerance: iteration-2 membership drift is still a source-alignment issue.
    """
    assert PROVENANCE["source_scope"] == "selected VandyCRUISE two-iteration run on committed synthetic input"
    assert PROVENANCE["source_sha256"]["AlgorithmMultipleToads.java"] == (
        "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    )
    atlas, result = _registered_two_iteration_run()
    assert len(atlas.structures) == 11
    assert result.iteration_count == 2

    with np.load(REFERENCE, allow_pickle=False) as java:
        first = result.exact_iterations[0]
        second = result.exact_iterations[1]
        first_memberships = np.asarray(first.membership_pass.memberships, dtype=np.float32)
        java_first = java["iteration1_memberships"]
        np.testing.assert_allclose(first_memberships, java_first, rtol=0, atol=4e-6)
        assert abs(
            float(first.membership_pass.objective)
            - PROVENANCE["java_trace"]["iteration1_objective"]
        ) < 2e-6

        first_cycle = result.registered_cycles[0]
        java_registration = PROVENANCE["java_trace"]["iteration1_registration"]
        initial_transform = np.asarray(
            result.alignment_run.post_refresh_registration.final_transform,
            dtype=np.float32,
        )
        np.testing.assert_allclose(
            initial_transform,
            java["iteration1_java_pre_registration_transform"],
            rtol=0,
            atol=3e-8,
        )
        class_bounds_low = np.zeros((11, 3), dtype=np.int32)
        class_bounds_high = np.broadcast_to(
            np.asarray(first_memberships.shape[:3], dtype=np.int32), (11, 3)
        ).copy()
        registered_classes = np.ones(11, dtype=bool)
        registered_classes[0] = False
        registration_kwargs = {
            "image_resolutions_mm": result.alignment_run.alignment.setup.signed_image_resolutions_mm,
            "template_resolutions_mm": result.alignment_run.alignment.setup.template_resolutions_mm,
            "transform": tuple(float(value) for value in initial_transform),
            "class_bounds_low_xyz": class_bounds_low,
            "class_bounds_high_xyz": class_bounds_high,
            "registered_classes": registered_classes,
            "offset": 0,
            "subsample": 2,
            "precomputed_affine": True,
        }
        atlas_shapes = toads_atlas_shapes_in_structure_order(
            atlas, result.alignment_run.atlas_data
        )
        java_transform_priors = transform_toads_shape_priors(
            atlas_shapes,
            image_dimensions_xyz=first_memberships.shape[:3],
            image_resolutions_mm=result.alignment_run.alignment.setup.signed_image_resolutions_mm,
            template_resolutions_mm=result.alignment_run.alignment.setup.template_resolutions_mm,
            transform=java_registration["pre_update_transform"],
            class_bounds=result.alignment_run.post_refresh_registration.final_shape_bounds,
        )
        np.testing.assert_array_equal(
            java_transform_priors, java["iteration1_registered_shape_priors"]
        )
        python_coefficients = vectorized_toads_registration_coefficients(
            first_memberships,
            atlas_shapes,
            **registration_kwargs,
        )
        python_gradient, python_hessian = exact_registration_measure_coefficients(
            python_coefficients.gradient,
            python_coefficients.hessian,
            python_coefficients.norm_gradient,
            python_coefficients.norm_hessian,
            python_coefficients.norm,
            source_faithful=True,
        )
        assert np.max(
            np.abs(python_gradient - java["iteration1_java_registration_gradient"])
        ) < 4e-9
        assert np.max(
            np.abs(python_hessian - java["iteration1_java_registration_hessian"])
        ) < 3e-6
        python_registration_energy = vectorized_toads_registration_energy(
            first_memberships,
            atlas_shapes,
            **registration_kwargs,
        )
        java_initial_transform = java["iteration1_java_pre_registration_transform"]
        shared_java_transform_kwargs = {
            **registration_kwargs,
            "transform": tuple(float(value) for value in java_initial_transform),
        }
        java_source_energy = np.float32(java["iteration1_java_registration_energy"][0])
        java_membership_energy = vectorized_toads_registration_energy(
            java["iteration1_memberships"],
            atlas_shapes,
            **shared_java_transform_kwargs,
        )
        python_membership_java_transform_energy = vectorized_toads_registration_energy(
            first_memberships,
            atlas_shapes,
            **shared_java_transform_kwargs,
        )
        assert java_membership_energy == java_source_energy
        assert python_membership_java_transform_energy == java_source_energy
        assert abs(float(python_registration_energy) - float(java_source_energy)) < 2e-7

        java_transform = java["iteration1_registered_transform"]
        python_transform = np.asarray(
            first_cycle.registration_control.handoff.transform, dtype=np.float32
        )
        transform_delta = python_transform - java_transform
        expected_transform_delta = np.asarray(
            [2e-4, 2e-4, 2e-4, 2e-4, 2e-4, 2e-4],
            dtype=np.float32,
        )
        np.testing.assert_allclose(transform_delta, expected_transform_delta, rtol=0, atol=5e-7)

        java_grown = java["iteration1_grown_segmentation"]
        python_grown = first_cycle.registered_iteration.branch.growing.segmentation
        assert np.count_nonzero(python_grown != java_grown) == 55
        java_closest = java["iteration1_closest_labels"]
        python_closest = first_cycle.registered_iteration.branch.relations.closest_labels
        assert np.count_nonzero(python_closest != java_closest) == 843

        next_thinning = result.registered_cycles[0].next_thinning.thinning.segmentation
        java_next_thinning = java["iteration2_thinning_segmentation"]
        assert next_thinning.shape == java_next_thinning.shape == (53, 60, 51)
        labels = np.asarray([item.label for item in atlas.structures], dtype=np.int16)
        dice_by_label = {}
        for index, label in enumerate(labels):
            python_class = next_thinning == index
            java_class = java_next_thinning == index
            denominator = int(np.count_nonzero(python_class) + np.count_nonzero(java_class))
            dice_by_label[int(label)] = (
                1.0
                if denominator == 0
                else 2.0 * np.count_nonzero(python_class & java_class) / denominator
            )
        changed = int(np.count_nonzero(next_thinning != java_next_thinning))
        assert changed == PROVENANCE["regression_envelope"]["iteration2_thinning_different_voxels"] == 42
        assert min(dice_by_label.values()) >= 0.995

        java_centroids = java["iteration2_centroids"]
        np.testing.assert_allclose(second.centroids, java_centroids, rtol=0, atol=1e-6)
        java_deviations = java["iteration2_standard_deviations"]
        deviation_delta = np.abs(second.standard_deviations - java_deviations)
        assert np.count_nonzero(deviation_delta > 1e-6) == 5
        assert np.max(deviation_delta) <= 1.1e-4

        # Reuse the Python source-matched variance kernel on identical Java
        # state to distinguish an input-state split from a kernel mismatch.
        image_values = np.moveaxis(
            np.asarray(result.alignment_run.template_refresh.recropped_images, dtype=np.float32),
            0,
            -1,
        )
        shared_state_deviations = toads_default_exact_class_variances(
            image_values,
            java["iteration1_memberships"],
            java_centroids,
            java_next_thinning,
            intensity_field=result.startup.fields,
        )
        np.testing.assert_allclose(
            shared_state_deviations, java_deviations, rtol=0, atol=2e-7
        )

        second_memberships = np.asarray(second.membership_pass.memberships, dtype=np.float32)
        java_second = java["iteration2_memberships"]
        difference = np.abs(second_memberships - java_second)
        mean_absolute = float(np.mean(difference))
        maximum_absolute = float(np.max(difference))
        count_over_1e5 = int(np.count_nonzero(difference > 1e-5))
        assert mean_absolute <= PROVENANCE["regression_envelope"]["iteration2_membership_mae_max"]
        assert maximum_absolute <= PROVENANCE["regression_envelope"]["iteration2_membership_max_error_max"]
        assert count_over_1e5 <= PROVENANCE["regression_envelope"]["iteration2_membership_count_over_1e5_max"]
        assert abs(
            float(second.membership_pass.objective)
            - PROVENANCE["java_trace"]["iteration2_objective"]
        ) <= PROVENANCE["regression_envelope"]["iteration2_objective_absolute_difference_max"]
