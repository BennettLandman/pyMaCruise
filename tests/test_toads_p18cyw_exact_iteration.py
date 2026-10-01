from __future__ import annotations

import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.toads_default_membership_registration_control import (
    continue_default_membership_through_registration_control,
)
from pymacruise.toads_default_membership_registered_iteration import (
    continue_default_membership_registered_iteration,
)
from pymacruise.toads_p18cyv_next_iteration import (
    continue_p18cyv_iteration_through_next_thinning,
)
from pymacruise.toads_p18cyw_exact_iteration import (
    continue_p18cyv_thinning_through_exact_iteration_degree1,
)
from test_toads_post_thinning_exact_iteration import (
    FIELD_HARNESS,
    ITERATION_HARNESS,
    _continuation_state,
    _field_input_bytes,
    _iteration_input_bytes,
    _bits,
)
from test_toads_p18cyt_registered_iteration import _p18cyt_membership

ROOT = Path(__file__).parent
CASE = json.loads(
    (ROOT / "data/synthetic/toads_p18cyw_exact_iteration/case.json").read_text()
)
RAS_FIXTURE = ROOT / "data/synthetic/toads_post_registration_branch/membership_maps.nii.gz"
FIELD_SOURCE_SHA256 = "6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8"
WRAPPER_SOURCE_SHA256 = "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
ITERATION_SOURCE_SHA256 = "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"


def _p18cyv_state():
    next_thinning, image, centroids, field, closest = _continuation_state()
    p18cyt = _p18cyt_membership(next_thinning, image, centroids, field, closest)
    current_membership = p18cyt.membership_pass
    nx, ny, nz = image.shape
    gx, gy, gz = np.meshgrid(
        np.linspace(-1, 1, nx, dtype=np.float32),
        np.linspace(-1, 1, ny, dtype=np.float32),
        np.linspace(-1, 1, nz, dtype=np.float32),
        indexing="ij",
    )
    atlas = np.empty((*image.shape, 11), dtype=np.float32)
    for k in range(11):
        radius = (
            (gx - np.float32(-0.45 + 0.09 * k)) ** 2
            + 1.3 * (gy - np.float32(0.3 * np.sin(k))) ** 2
            + 0.8 * (gz - np.float32(0.25 * np.cos(k))) ** 2
        )
        atlas[..., k] = np.asarray(
            np.exp(-radius / np.float32(0.42)) * 0.75 + 0.05, dtype=np.float32
        )
    objective = np.float32(current_membership.objective)
    previous_energy = np.float32(objective * np.float32(1.01))
    registration = continue_default_membership_through_registration_control(
        current_membership,
        atlas,
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        initial_transform=(0, 0, 0, 0, 0, 0),
        previous_energy=previous_energy,
        iteration=3,
        maximum_iterations=5,
        convergence_threshold=1e-4,
        do_growing=False,
        align_shapes=False,
    )
    prior_branch = next_thinning.branch
    old_relations = prior_branch.relations
    assert old_relations is not None and prior_branch.growing is not None
    boundary = np.repeat(
        np.repeat(
            np.repeat(
                old_relations.boundary_distances,
                CASE["spatial_repeat_factor"],
                axis=0,
            ),
            CASE["spatial_repeat_factor"],
            axis=1,
        ),
        CASE["spatial_repeat_factor"],
        axis=2,
    )
    registered_iteration = continue_default_membership_registered_iteration(
        registration,
        previous_energy=previous_energy,
        iteration=3,
        maximum_iterations=5,
        convergence_threshold=1e-4,
        segmentation=next_thinning.thinning.segmentation,
        skeleton=next_thinning.thinning.skeleton,
        ordering=prior_branch.growing.ordering,
        class_volumes=next_thinning.thinning.class_volumes,
        boundary_distances=boundary,
        template_labels=np.asarray([0, 4, 7, 9, 11, 13, 17, 19, 21, 23, 25], dtype=np.int8),
        previous_classification=np.full(image.shape, -1, dtype=np.int8),
        relation_limit=1,
        do_growing=False,
        initial_available=next_thinning.thinning.available,
    )
    shapes = np.zeros_like(current_membership.memberships, dtype=np.float32)
    p18cyv = continue_p18cyv_iteration_through_next_thinning(registered_iteration, shapes)
    return p18cyt, p18cyv, image


def test_p18cyv_output_continues_through_java_referenced_exact_iteration(tmp_path: Path) -> None:
    p18cyt, p18cyv, image = _p18cyv_state()
    source_image = nib.load(str(RAS_FIXTURE))
    assert nib.aff2axcodes(source_image.affine) == tuple(CASE["source_fixture_axis_codes"])
    affine = source_image.affine @ np.diag([1 / CASE["spatial_repeat_factor"]] * 3 + [1])
    assert nib.aff2axcodes(affine) == ("R", "A", "S")
    assert nib.Nifti1Image(image, affine).shape == image.shape

    membership = p18cyv.next_thinning.branch.handoff.membership_pass.memberships
    labels = p18cyv.next_thinning.thinning.segmentation
    assert labels.shape == image.shape == membership.shape[:3]
    relations = p18cyv.registered_iteration.branch.relations
    assert relations is not None
    closest = relations.closest_labels
    assert closest.shape == (*image.shape, 4)

    field = np.asarray(p18cyt.field_update.fit.field, dtype=np.float32).copy()
    original_field = field.copy()
    centroids = np.asarray(p18cyt.centroids[0], dtype=np.float32).copy()
    exponents = np.empty((1, 11), dtype=np.float32)
    exponents[0, 0] = 1.0
    exponents[0, 1:3] = 0.8
    exponents[0, 3:8] = 1.0
    exponents[0, 8:11] = 1.0
    shape_values = np.zeros((*image.shape, 11), dtype=np.float32)
    modality_weights = np.asarray([0.73], dtype=np.float32)
    variance_weights = np.asarray(
        [[0.81, 0.93, 1.17, 0.76, 0.88, 1.02, 1.14, 0.69, 0.96, 1.21, 0.73]],
        dtype=np.float32,
    )
    smoothing_priors = np.asarray(
        [0.65, 0.71, 0.84, 0.91, 1.02, 1.13, 0.76, 0.88, 1.05, 1.17, 0.82],
        dtype=np.float32,
    )
    atlas_priors = np.asarray(
        [0.10, 0.08, 0.09, 0.11, 0.12, 0.09, 0.10, 0.11, 0.07, 0.06, 0.07],
        dtype=np.float32,
    )
    arguments = dict(
        category_exponents=np.asarray([0.8, 1.0, 1.0], dtype=np.float32),
        exponent_by_channel_class=exponents,
        modality_weights=modality_weights,
        class_variance_weights=variance_weights,
        transformed_shapes=shape_values,
        closest_labels=closest,
        smoothing_priors=smoothing_priors,
        atlas_shape_priors=atlas_priors,
        smoothing=0.043,
        prior_factor=0.64,
        prior_coefficient=0.16,
        prior_scale=0.05,
    )

    field_input = tmp_path / "p18cyw_field_input.bin"
    field_input.write_bytes(_field_input_bytes(labels, image, membership, original_field, centroids))
    field_reference = json.loads(subprocess.run(
        ["sh", str(FIELD_HARNESS), str(field_input), *(str(v) for v in image.shape), "11"],
        check=True, capture_output=True, text=True,
    ).stdout)
    result = continue_p18cyv_thinning_through_exact_iteration_degree1(
        p18cyv,
        image,
        centroids,
        field,
        **arguments,
    )
    exact = result.exact_iteration
    assert exact.field_update.degree == 1
    assert result.p18cyv is p18cyv
    assert exact.field_update.next_thinning is p18cyv.next_thinning
    assert result.stage_order[-7:] == (
        "next_iteration.propagateThinning",
        "next_iteration.computeExact3DPolynomialField.degree1",
        "next_iteration.computeExactCentroids",
        "next_iteration.computeExactVariances",
        "next_iteration.computeImageVariance",
        "next_iteration.computeCentroidSimilarity",
        "next_iteration.computeExactMemberships",
    )
    assert field_reference["reference_kind"] == "formula_harness_not_source_class_execution"
    assert field_reference["inhomogeneity_upstream_source_sha256"] == FIELD_SOURCE_SHA256
    assert field_reference["wrapper_upstream_source_sha256"] == WRAPPER_SOURCE_SHA256
    java_field = np.asarray(field_reference["field_bits"], dtype=np.int32).view(np.float32).reshape(image.shape)
    np.testing.assert_allclose(exact.field_update.fit.field, java_field, rtol=0, atol=2e-7)
    np.testing.assert_allclose(
        exact.field_update.fit.coefficients,
        field_reference["coefficients"],
        rtol=2e-6,
        atol=2e-8,
    )

    iteration_input = tmp_path / "p18cyw_iteration_input.bin"
    iteration_output = tmp_path / "p18cyw_iteration_output.bin"
    power_output = tmp_path / "p18cyw_power_output.bin"
    iteration_input.write_bytes(_iteration_input_bytes(
        labels, image, membership, exact.field_update.fit.field, centroids, closest,
    ))
    reference = json.loads(subprocess.run(
        [
            "sh", str(ITERATION_HARNESS), str(iteration_input), *(str(v) for v in image.shape),
            str(iteration_output), str(power_output),
        ],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert reference["source_sha256"] == ITERATION_SOURCE_SHA256
    assert reference["reference_kind"] == "source_formula_harness"
    assert reference["source_order"] == [
        "computeCorrectionField.degree1", "computeExactCentroids",
        "computeExactVariances", "computeImageVariance", "centroid similarity",
        "computeExactMemberships",
    ]
    assert reference["centroid_bits"] == _bits(exact.centroids)
    assert reference["standard_deviation_bits"] == _bits(exact.standard_deviations)
    assert reference["similarity_bits"] == _bits(exact.centroid_similarity)
    assert reference["image_variance_bits"] == int(exact.image_variance_weights.view(np.int32)[0])
    assert reference["image_std_bits"] == int(exact.image_standard_deviation.view(np.int32)[0])

    reference_arrays = np.fromfile(iteration_output, dtype=">f4").astype(np.float32)
    energy_count = exact.membership_pass.class_data_energy.size
    membership_count = exact.membership_pass.memberships.size
    denominator_count = exact.membership_pass.denominator.size
    java_energy = reference_arrays[:energy_count].reshape(exact.membership_pass.class_data_energy.shape)
    java_membership = reference_arrays[energy_count:energy_count + membership_count].reshape(exact.membership_pass.memberships.shape)
    java_denominator = reference_arrays[energy_count + membership_count:energy_count + membership_count + denominator_count].reshape(exact.membership_pass.denominator.shape)
    np.testing.assert_array_equal(_bits(java_energy), _bits(exact.membership_pass.class_data_energy))
    np.testing.assert_array_equal(_bits(java_membership), _bits(exact.membership_pass.memberships))
    np.testing.assert_array_equal(_bits(java_denominator), _bits(exact.membership_pass.denominator))
    assert reference["objective_bits"] == int(_bits([exact.membership_pass.objective])[0])
