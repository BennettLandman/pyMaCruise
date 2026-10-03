from __future__ import annotations

import json
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_alignment_startup import register_toads_from_atlas_assets
from pymacruise.toads_atlas import load_toads_intensity_variance_factors
from pymacruise.toads_default_atlas_membership import toads_default_atlas_membership_pass
from pymacruise.toads_default_exact_centroids import toads_default_exact_centroid_update
from pymacruise.toads_default_exact_variances import toads_default_exact_class_variances
from pymacruise.toads_math import toads_optimized_centroid_similarity, toads_optimized_image_variance
from pymacruise.toads_registered_field_refit import (
    compute_registered_toads_corrected_initial_membership,
    fit_registered_toads_correction_fields,
)
from pymacruise.toads_registered_initial_membership import compute_registered_toads_initial_membership
from pymacruise.toads_registered_unity_field_startup import compute_registered_toads_unity_field_startup
from pymacruise.toads_registered_first_exact_iteration import (
    continue_registered_unity_startup_through_first_exact_iteration,
)
from pymacruise.toads_thinning import toads_propagate_thinning
ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "toads_registered_field_refit"
HARNESS = ROOT / "java_reference" / "ToadsExactPolynomialFieldReferenceHarness.java"
UPSTREAM = ROOT / "java_reference" / "upstream"


def _field_input_bytes(labels, image, memberships, field, centers) -> bytes:
    return b"".join((
        np.asarray(labels, dtype=np.int8).tobytes(order="C"),
        np.asarray(image, dtype=">f4").tobytes(order="C"),
        np.asarray(memberships, dtype=">f4").tobytes(order="C"),
        np.asarray(field, dtype=">f4").tobytes(order="C"),
        np.asarray(centers, dtype=">f4").tobytes(order="C"),
    ))


def test_registered_startup_field_fit_matches_degree_three_java_reference(tmp_path: Path) -> None:
    descriptor = FIXTURE / "atlas-11.txt"
    t1 = nib.load(str(FIXTURE / "t1.nii.gz"))
    mask = nib.load(str(FIXTURE / "brain-mask.nii.gz"))
    assert nib.aff2axcodes(t1.affine) == nib.aff2axcodes(mask.affine) == ("R", "A", "S")
    run = register_toads_from_atlas_assets(
        t1, mask, descriptor, levels=1, max_iterations=1,
        refinement_iterations=1, padding=2,
    )
    settings = dict(
        exponent_by_channel_class=np.full((1, 11), 2.0, dtype=np.float32),
        modality_weights=np.ones(1, dtype=np.float32),
        smoothing_priors=np.ones(11, dtype=np.float32),
        atlas_shape_priors=np.ones(11, dtype=np.float32),
        smoothing=np.float32(0.2 / 6.0 / 11.0),
        prior_factor=np.float32(0.2),
        prior_coefficient=np.float32(2.0 / 11.0),
        prior_scale=np.float32(0.05),
    )
    result = compute_registered_toads_corrected_initial_membership(run, **settings)
    field_refit = result.field_refit
    corrected = result.statistics_membership
    images = np.moveaxis(run.template_refresh.recropped_images, 0, -1)
    memberships = run.segmentation_initialization.memberships
    segmentation = run.segmentation_constructor.segmentation
    centroids = run.segmentation_initialization.centroids

    build = tmp_path / "classes"
    build.mkdir()
    java_sources = sorted((UPSTREAM / "Jama").glob("*.java")) + sorted((UPSTREAM / "Jama" / "util").glob("*.java"))
    java_sources += [UPSTREAM / "gov/nih/mipav/model/structures/jama/JamaMatrix.java", UPSTREAM / "edu/jhmi/rad/medic/utilities/PowerTable.java", HARNESS]
    subprocess.run(["javac", "-d", str(build), *map(str, java_sources)], check=True, capture_output=True, text=True)
    input_path = tmp_path / "registered_exact_field.bin"
    input_path.write_bytes(_field_input_bytes(
        segmentation, images[..., 0], memberships, np.ones(segmentation.shape, dtype=np.float32), centroids[0],
    ))
    java = json.loads(subprocess.run(
        ["java", "-cp", str(build), "ToadsExactPolynomialFieldReferenceHarness",
         str(input_path), *(str(n) for n in segmentation.shape), "11", "3"],
        check=True, capture_output=True, text=True,
    ).stdout)

    assert field_refit.degree == 3
    assert field_refit.channel_order == (0,)
    assert result.stage_order == (
        "registered_startup.computeExact3DPolynomialField.degree3",
        "registered_startup.computeVariances",
        "registered_startup.computeImageVariance",
        "registered_startup.computeCentroidSimilarity",
        "registered_startup.computeExactMemberships",
        "registered_startup.initThinningLabels",
        "registered_startup.propagateThinning",
    )
    expected_field = np.asarray(java["field_bits"], dtype=np.int32).view(np.float32).reshape(segmentation.shape)
    np.testing.assert_array_equal(field_refit.fits[0].normal_matrix.ravel(), np.asarray(java["normal"], dtype=np.float64))
    np.testing.assert_array_equal(field_refit.fits[0].right_hand_side, np.asarray(java["rhs"], dtype=np.float64))
    np.testing.assert_allclose(field_refit.fields[..., 0], expected_field, rtol=0, atol=2e-6)
    np.testing.assert_array_equal(field_refit.fields[..., 0][segmentation == 0], 1.0)

    direct = compute_registered_toads_initial_membership(
        run, **settings, intensity_field=field_refit.fields,
    )
    uncorrected = compute_registered_toads_initial_membership(run, **settings)
    np.testing.assert_array_equal(corrected.membership_pass.memberships, direct.membership_pass.memberships)
    np.testing.assert_array_equal(corrected.membership_pass.class_data_energy, direct.membership_pass.class_data_energy)
    assert not np.array_equal(
        corrected.membership_pass.memberships,
        uncorrected.membership_pass.memberships,
    )
    expected_thinning = toads_propagate_thinning(
        segmentation,
        corrected.membership_pass.memberships,
        run.relation_state.relation.class_adjacency,
        run.transformed_shape_priors,
        volume_ratio=0.5,
        limit=3.0,
    )
    np.testing.assert_array_equal(result.thinning.segmentation, expected_thinning.segmentation)
    np.testing.assert_array_equal(result.thinning.skeleton, expected_thinning.skeleton)
    assert result.thinning.processed_items == expected_thinning.processed_items


def test_registered_wrapper_startup_keeps_unity_field_and_degree_zero(tmp_path: Path) -> None:
    descriptor = FIXTURE / "atlas-11.txt"
    t1 = nib.load(str(FIXTURE / "t1.nii.gz"))
    mask = nib.load(str(FIXTURE / "brain-mask.nii.gz"))
    run = register_toads_from_atlas_assets(
        t1, mask, descriptor, levels=1, max_iterations=1,
        refinement_iterations=1, padding=2,
    )
    settings = dict(
        exponent_by_channel_class=np.full((1, 11), 2.0, dtype=np.float32),
        modality_weights=np.ones(1, dtype=np.float32),
        smoothing_priors=np.ones(11, dtype=np.float32),
        atlas_shape_priors=np.ones(11, dtype=np.float32),
        smoothing=np.float32(0.2 / 6.0 / 11.0),
        prior_factor=np.float32(0.2),
        prior_coefficient=np.float32(2.0 / 11.0),
        prior_scale=np.float32(0.05),
    )
    result = compute_registered_toads_unity_field_startup(run, **settings)

    image_shape = run.template_refresh.recropped_images.shape[1:]
    assert result.degree_state == 0
    assert result.energy_ratio == np.float32(1.0)
    assert result.fields.shape == image_shape + (1,)
    np.testing.assert_array_equal(result.fields, np.ones(result.fields.shape, dtype=np.float32))
    assert result.stage_order == (
        "registered_startup.initializeCorrectionFieldsToUnity",
        "registered_startup.computeVariances",
        "registered_startup.computeImageVariance",
        "registered_startup.computeCentroidSimilarity",
        "registered_startup.computeExactMemberships",
        "iteration1.propagateThinning",
    )

    direct_membership = compute_registered_toads_initial_membership(
        run, **settings, intensity_field=result.fields,
    )
    np.testing.assert_array_equal(
        result.statistics_membership.membership_pass.memberships,
        direct_membership.membership_pass.memberships,
    )
    direct_thinning = toads_propagate_thinning(
        run.segmentation_constructor.segmentation,
        direct_membership.membership_pass.memberships,
        run.relation_state.relation.class_adjacency,
        run.transformed_shape_priors,
        volume_ratio=0.5,
        limit=3.0,
    )
    np.testing.assert_array_equal(result.thinning.segmentation, direct_thinning.segmentation)
    np.testing.assert_array_equal(result.thinning.skeleton, direct_thinning.skeleton)


def test_registered_unity_startup_continues_through_first_degree_zero_exact_loop() -> None:
    descriptor = FIXTURE / "atlas-11.txt"
    t1 = nib.load(str(FIXTURE / "t1.nii.gz"))
    mask = nib.load(str(FIXTURE / "brain-mask.nii.gz"))
    run = register_toads_from_atlas_assets(
        t1, mask, descriptor, levels=1, max_iterations=1,
        refinement_iterations=1, padding=2,
    )
    settings = dict(
        exponent_by_channel_class=np.full((1, 11), 2.0, dtype=np.float32),
        modality_weights=np.ones(1, dtype=np.float32),
        smoothing_priors=np.ones(11, dtype=np.float32),
        atlas_shape_priors=np.ones(11, dtype=np.float32),
        smoothing=np.float32(0.2 / 6.0 / 11.0),
        prior_factor=np.float32(0.2),
        prior_coefficient=np.float32(2.0 / 11.0),
        prior_scale=np.float32(0.05),
    )
    startup = compute_registered_toads_unity_field_startup(run, **settings)
    result = continue_registered_unity_startup_through_first_exact_iteration(
        run,
        startup,
        category_exponents=np.asarray([[1.0, 1.0, 1.0]], dtype=np.float32),
            class_variance_weights=load_toads_intensity_variance_factors(FIXTURE / "atlas-11.txt", [run.modality]),
        **settings,
    )

    assert result.degree == 0
    assert result.fields is startup.fields
    np.testing.assert_array_equal(result.fields, np.ones_like(result.fields))
    assert result.stage_order[-6:] == (
        "iteration1.computeCorrectionField.degree0.noop",
        "iteration1.computeExactCentroids",
        "iteration1.computeExactVariances",
        "iteration1.computeImageVariance",
        "iteration1.computeCentroidSimilarity",
        "iteration1.computeExactMemberships",
    )
    assert result.membership_pass.memberships.shape == startup.thinning.segmentation.shape + (11,)
    assert not np.array_equal(
        result.membership_pass.memberships,
        startup.statistics_membership.membership_pass.memberships,
    )
    images = np.moveaxis(run.template_refresh.recropped_images, 0, -1)
    memberships = startup.statistics_membership.membership_pass.memberships
    fields = startup.fields
    centers = toads_default_exact_centroid_update(
        images,
        memberships,
        run.segmentation_initialization.centroids,
        np.asarray([[1.0, 1.0, 1.0]], dtype=np.float32),
        intensity_field=fields,
    )
    deviations = toads_default_exact_class_variances(
        images,
        memberships,
        centers,
        startup.thinning.segmentation,
        intensity_field=fields,
    )
    image_variance, image_std = toads_optimized_image_variance(
        images,
        memberships[..., 1:],
        centers[:, 1:],
        intensity_field=fields,
    )
    similarity = toads_optimized_centroid_similarity(
        centers, deviations, image_variance, settings["prior_scale"]
    )
    expected_membership = toads_default_atlas_membership_pass(
        image_values=images,
        initial_memberships=memberships,
        centroids=centers,
        exponent_by_channel_class=settings["exponent_by_channel_class"],
        modality_weights=settings["modality_weights"],
            class_variance_weights=load_toads_intensity_variance_factors(FIXTURE / "atlas-11.txt", [run.modality]),
        image_variances=image_variance,
        transformed_shapes=run.transformed_shape_priors,
        closest_labels=run.relation_state.closest_labels,
        smoothing_priors=settings["smoothing_priors"],
        centroid_similarity=similarity,
        atlas_shape_priors=settings["atlas_shape_priors"],
        smoothing=settings["smoothing"],
        prior_factor=settings["prior_factor"],
        prior_coefficient=settings["prior_coefficient"],
        intensity_field=fields,
    )
    np.testing.assert_array_equal(result.centroids, centers)
    np.testing.assert_array_equal(result.class_standard_deviations, deviations)
    np.testing.assert_array_equal(result.image_variance_weights, image_variance)
    np.testing.assert_array_equal(result.image_standard_deviation, image_std)
    np.testing.assert_array_equal(result.centroid_similarity, similarity)
    np.testing.assert_array_equal(result.membership_pass.memberships, expected_membership.memberships)
    np.testing.assert_array_equal(result.membership_pass.denominator, expected_membership.denominator)
    assert result.membership_pass.objective == expected_membership.objective


@pytest.mark.parametrize("degree", [0, 4, True])
def test_registered_field_refit_rejects_unsupported_degree(degree) -> None:
    with pytest.raises(ValueError, match="degree must be"):
        fit_registered_toads_correction_fields(object(), degree=degree)
