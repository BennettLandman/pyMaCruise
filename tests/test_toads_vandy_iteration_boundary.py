from __future__ import annotations

import hashlib
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
    register_toads_atlas_rigid,
    vectorized_toads_registration_coefficients,
    vectorized_toads_registration_energy,
)
from pymacruise.toads_shape_volume import transform_toads_shape_priors
from pymacruise.toads_atlas import toads_atlas_shapes_in_structure_order
from pymacruise.toads_registered_loop import (
    ToadsRegisteredLoopConfig,
    run_registered_toads_loop,
    toads_effective_prior_coefficient,
    toads_effective_smoothing,
)


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data/synthetic/toads_downsampled_vandy_atlas"
REFERENCE = FIXTURE / "two_iteration_java_reference.npz"
PROVENANCE = json.loads((FIXTURE / "two_iteration_java_reference.json").read_text())
TWELVE_ITERATION_PROVENANCE = json.loads(
    (FIXTURE / "twelve_iteration_java_reference.json").read_text()
)
TWELVE_ITERATION_MEMBERSHIPS = FIXTURE / "twelve_iteration_java_memberships.npz"
POST_REFRESH_REFERENCE = json.loads(
    (FIXTURE / "post_refresh_registration_java_reference.json").read_text()
)
POST_REFRESH_HARNESS = (
    ROOT / "java_reference" / "ToadsVandyPostRefreshRegistrationHarness.java"
)
ATLAS = ROOT.parent / "data/atlases/cruise-atlas-10obj-toads2008.txt"


def _registered_two_iteration_run(maximum_iterations: int = 2):
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
        smoothing=toads_effective_smoothing(0.2, len(atlas.structures)),
        prior_factor=0.5,
        prior_coefficient=toads_effective_prior_coefficient(2.0, len(atlas.structures)),
        maximum_iterations=maximum_iterations,
        convergence_threshold=1e-4,
        relation_limit=3,
        do_growing=True,
        maximum_polynomial_degree=3,
        maximum_distance=1e-4,
        max_registration_iterations=1,
        retain_iteration_history=True,
    )
    return atlas, run_registered_toads_loop(alignment, config)


def test_vandy_iteration_growth_and_next_thinning_match_java_reference() -> None:
    """Compare two loop-boundary states to a pinned selected-source Java run.

    Compare the selected two-iteration synthetic Vandy reference, with the
    whole Java composite and natural convergence kept as separate claims.
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
        np.testing.assert_array_equal(first_memberships, java_first)
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
        java_initial_transform = java["iteration1_java_pre_registration_transform"]
        registration_kwargs = {
            "image_resolutions_mm": result.alignment_run.alignment.setup.signed_image_resolutions_mm,
            "template_resolutions_mm": result.alignment_run.alignment.setup.template_resolutions_mm,
            "transform": tuple(float(value) for value in java_initial_transform),
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
        assert hashlib.sha256(POST_REFRESH_HARNESS.read_bytes()).hexdigest() == (
            POST_REFRESH_REFERENCE["harness_sha256"]
        )
        java_post_refresh_start = result.alignment_run.template_refresh.recropped_transform
        np.testing.assert_array_equal(
            java_post_refresh_start.view(np.uint32),
            np.asarray(
                POST_REFRESH_REFERENCE["initial_transform_bits"], dtype=np.int32
            ).view(np.uint32),
        )
        source_coefficients = vectorized_toads_registration_coefficients(
            java_first,
            atlas_shapes,
            **{
                **registration_kwargs,
                "transform": tuple(float(value) for value in java_post_refresh_start),
            },
        )
        source_gradient, source_hessian = exact_registration_measure_coefficients(
            source_coefficients.gradient,
            source_coefficients.hessian,
            source_coefficients.norm_gradient,
            source_coefficients.norm_hessian,
            source_coefficients.norm,
            source_faithful=True,
        )
        np.testing.assert_array_equal(
            source_gradient.view(np.uint32),
            np.asarray(POST_REFRESH_REFERENCE["gradient_bits"], dtype=np.int32).view(np.uint32),
        )
        np.testing.assert_array_equal(
            source_hessian.view(np.uint32),
            np.asarray(POST_REFRESH_REFERENCE["hessian_bits"], dtype=np.int32).view(np.uint32),
        )
        source_energy = vectorized_toads_registration_energy(
            java_first,
            atlas_shapes,
            **{
                **registration_kwargs,
                "transform": tuple(float(value) for value in java_post_refresh_start),
            },
        )
        expected_source_energy = np.asarray(
            [POST_REFRESH_REFERENCE["energy_bits"]], dtype=np.int32
        ).view(np.float32)[0]
        assert abs(
            int(np.asarray(source_energy).view(np.uint32))
            - int(np.asarray(expected_source_energy).view(np.uint32))
        ) <= 1
        source_registration = register_toads_atlas_rigid(
            java_first,
            atlas_shapes,
            image_resolutions_mm=result.alignment_run.alignment.setup.signed_image_resolutions_mm,
            template_resolutions_mm=result.alignment_run.alignment.setup.template_resolutions_mm,
            initial_transform=tuple(float(value) for value in java_post_refresh_start),
            class_bounds_low_xyz=class_bounds_low,
            class_bounds_high_xyz=class_bounds_high,
            registered_classes=registered_classes,
            offset=0,
            subsample=2,
            max_iterations=1,
            energy_tolerance=1e-6,
            minimum_lambda=1e-3,
            precomputed_affine=True,
        )
        np.testing.assert_array_equal(
            source_registration.transform.view(np.uint32),
            np.asarray(POST_REFRESH_REFERENCE["final_transform_bits"], dtype=np.int32).view(np.uint32),
        )
        python_transform_bits = result.alignment_run.post_refresh_registration.final_transform.view(np.uint32).astype(np.int64)
        java_transform_bits = np.asarray(
            POST_REFRESH_REFERENCE["final_transform_bits"], dtype=np.int32
        ).view(np.uint32).astype(np.int64)
        assert int(np.abs(python_transform_bits - java_transform_bits).max()) <= 1
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
        # The selected Vandy atlas constructs its affine mapper with
        # precompute=true. At a shared transform, source mapping plus
        # coefficient accumulation should therefore stay within one float32
        # ULP; a center-relative mapper can be many ULPs away on this grid.
        for actual, expected in (
            (python_gradient, java["iteration1_java_registration_gradient"]),
            (python_hessian, java["iteration1_java_registration_hessian"]),
        ):
            ulp_distance = np.abs(
                actual.view(np.int32).astype(np.int64)
                - np.asarray(expected, dtype=np.float32).view(np.int32).astype(np.int64)
            )
            assert int(ulp_distance.max()) <= 1
        python_registration_energy = vectorized_toads_registration_energy(
            first_memberships,
            atlas_shapes,
            **registration_kwargs,
        )
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
        np.testing.assert_array_equal(python_transform, java_transform)

        java_grown = java["iteration1_grown_segmentation"]
        python_grown = first_cycle.registered_iteration.branch.growing.segmentation
        assert np.count_nonzero(python_grown != java_grown) == 0
        java_closest = java["iteration1_closest_labels"]
        python_closest = first_cycle.registered_iteration.branch.relations.closest_labels
        assert np.count_nonzero(python_closest != java_closest) == 0

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
        assert changed == PROVENANCE["python_observed_after_ordering_carry_fix"]["iteration2_thinning_different_voxels"] == 0
        assert min(dice_by_label.values()) == 1.0

        java_centroids = java["iteration2_centroids"]
        np.testing.assert_allclose(second.centroids, java_centroids, rtol=0, atol=1e-6)
        java_deviations = java["iteration2_standard_deviations"]
        deviation_delta = np.abs(second.standard_deviations - java_deviations)
        assert np.count_nonzero(deviation_delta > 1e-6) == 0
        assert np.max(deviation_delta) <= 2e-7

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
        corrected_python = PROVENANCE["python_observed_after_ordering_carry_fix"]
        assert mean_absolute == corrected_python["iteration2_membership_mae"]
        assert maximum_absolute == corrected_python["iteration2_membership_max_error"]
        assert count_over_1e5 == corrected_python["iteration2_membership_count_over_1e5"]
        objective_difference = abs(
            float(second.membership_pass.objective)
            - PROVENANCE["java_trace"]["iteration2_objective"]
        )
        assert objective_difference == corrected_python["iteration2_objective_absolute_difference"]
        assert np.float32(second.membership_pass.objective) == np.float32(PROVENANCE["java_trace"]["iteration2_objective"])


def test_vandy_source_faithful_coefficients_track_java_and_meet_label_dice_gate() -> None:
    """Check corrected source coefficients and per-label output acceptance.

    Python now reaches its configured energy threshold before the selected
    Java iteration cap. Compare final labels by the agreed per-label Dice
    metric; do not require matching iteration counts or hashes.
    """
    reference = TWELVE_ITERATION_PROVENANCE
    assert reference["source_sha256"]["AlgorithmMultipleToads.java"] == (
        PROVENANCE["source_sha256"]["AlgorithmMultipleToads.java"]
    )
    assert reference["input_sha256"] == PROVENANCE["input_sha256"]
    assert hashlib.sha256(TWELVE_ITERATION_MEMBERSHIPS.read_bytes()).hexdigest() == (
        reference["java_membership_archive"]["sha256"]
    )

    atlas, result = _registered_two_iteration_run(maximum_iterations=12)
    assert len(atlas.structures) == 11
    python_run = reference["python_run_after_float32_coefficient_fix"]
    assert reference["iterations"] == 12
    assert result.iteration_count == python_run["iterations"] == 11
    assert result.iteration_limit_reached is False
    assert result.energy_threshold_reached is True
    assert float(result.terminal_energy_ratio) == python_run["terminal_energy_ratio"]

    with np.load(TWELVE_ITERATION_MEMBERSHIPS, allow_pickle=False) as java:
        for iteration, exact in enumerate(result.exact_iterations, start=1):
            python_memberships = np.asarray(
                exact.membership_pass.memberships, dtype=np.float32
            )
            java_memberships = np.asarray(java[f"iteration{iteration}"], dtype=np.float32)
            assert python_memberships.shape == java_memberships.shape == (53, 60, 51, 11)
            difference = np.abs(python_memberships - java_memberships)
            diagnostic = reference["python_java_membership_diagnostics"][iteration - 1]
            assert float(difference.mean()) == diagnostic["mean_abs"]
            assert float(difference.max()) == diagnostic["max_abs"]
            assert int(np.count_nonzero(difference > 1e-6)) == diagnostic["count_gt_1e6"]
            assert int(np.count_nonzero(difference > 1e-5)) == diagnostic["count_gt_1e5"]
            python_argmax = np.argmax(python_memberships, axis=-1)
            java_argmax = np.argmax(java_memberships, axis=-1)
            assert float(np.mean(python_argmax == java_argmax)) == diagnostic[
                "argmax_agreement"
            ]

    labels = np.asarray(result.native_classification, dtype=np.uint8)
    assert labels.shape == (36, 45, 34)
    labels_source_order = np.asfortranarray(labels)
    actual_hash = hashlib.sha256(labels_source_order.tobytes(order="F")).hexdigest()
    assert actual_hash == python_run["python_classification_sha256"]

    java_path = FIXTURE / reference["java_final_classification_nifti"]
    assert hashlib.sha256(java_path.read_bytes()).hexdigest() == (
        reference["java_final_classification_nifti_sha256"]
    )
    java_labels = np.asanyarray(nib.load(str(java_path)).dataobj, dtype=np.uint8)
    java_source_hash = hashlib.sha256(
        np.asfortranarray(java_labels).tobytes(order="F")
    ).hexdigest()
    assert java_source_hash == reference["java_final_classification_sha256"]
    dice_by_label = {}
    for structure in atlas.structures:
        label = int(structure.label)
        python_class = labels == label
        java_class = java_labels == label
        denominator = int(np.count_nonzero(python_class) + np.count_nonzero(java_class))
        dice_by_label[label] = (
            1.0
            if denominator == 0
            else 2.0 * np.count_nonzero(python_class & java_class) / denominator
        )
    java_python = reference["java_python_final_label_dice_after_float32_coefficient_fix"]
    assert dice_by_label == {
        int(key): value for key, value in java_python["per_label"].items()
    }
    assert min(dice_by_label.values()) == java_python["minimum_dice"]
    assert min(dice_by_label.values()) >= 0.97
    assert int(np.count_nonzero(labels != java_labels)) == java_python["changed_voxels"] == 4
