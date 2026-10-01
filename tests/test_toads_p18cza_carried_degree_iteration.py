from __future__ import annotations

import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.toads_inhomogeneity_transition import toads_next_inhomogeneity_degree
from pymacruise.toads_p18cza_carried_degree_iteration import (
    continue_p18cyz_thinning_through_carried_degree1_iteration,
)
from pymacruise.toads_p18cyx_next_loop_thinning import (
    continue_p18cyx_branch_through_next_thinning,
)
from test_toads_p18cyx_registered_cycle import _p18cyw_state
from test_toads_p18cyz_next_loop_thinning import _make_p18cyx
from test_toads_post_thinning_exact_iteration import (
    FIELD_HARNESS,
    ITERATION_HARNESS,
    _field_input_bytes,
    _iteration_input_bytes,
    _bits,
)

ROOT = Path(__file__).parent
CASE = json.loads(
    (ROOT / "data/synthetic/toads_p18cza_carried_degree_iteration/case.json").read_text()
)
RAS_FIXTURE = ROOT / "data" / "synthetic" / CASE["source_fixture"]
FIELD_SOURCE_SHA256 = "6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8"
ITERATION_SOURCE_SHA256 = "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"


def test_p18cyz_thinning_composes_through_carried_degree_one_java_references(tmp_path: Path) -> None:
    previous, image = _p18cyw_state()
    cycle = _make_p18cyx(previous, image)
    p18cyz = cycle.registered_iteration.branch
    shapes = np.zeros_like(p18cyz.handoff.membership_pass.memberships, dtype=np.float32)
    source_args = dict(
        category_exponents=np.asarray([0.8, 1.0, 1.0], dtype=np.float32),
        exponent_by_channel_class=np.asarray(
            [[1.0, 0.8, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]],
            dtype=np.float32,
        ),
        modality_weights=np.asarray([0.73], dtype=np.float32),
        class_variance_weights=np.asarray(
            [[0.81, 0.93, 1.17, 0.76, 0.88, 1.02, 1.14, 0.69, 0.96, 1.21, 0.73]],
            dtype=np.float32,
        ),
        transformed_shapes=shapes,
        closest_labels=p18cyz.relations.closest_labels,
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
    source_image = nib.load(str(RAS_FIXTURE))
    assert nib.aff2axcodes(source_image.affine) == tuple(CASE["source_fixture_axis_codes"])
    ras_affine = source_image.affine @ np.diag(
        [1 / CASE["spatial_repeat_factor"]] * 3 + [1]
    )
    assert nib.aff2axcodes(ras_affine) == ("R", "A", "S")

    field = np.asarray(previous.exact_iteration.field_update.fit.field, dtype=np.float32).copy()
    original_field = field.copy()
    centroids = np.asarray(previous.exact_iteration.centroids[0], dtype=np.float32).copy()
    p18cyz_result = continue_p18cyx_branch_through_next_thinning(cycle, shapes)
    ratio = np.float32(p18cyz_result.next_thinning.branch.control.energy_ratio)
    scheduled_degree = toads_next_inhomogeneity_degree(
        CASE["previous_polynomial_degree"],
        CASE["maximum_polynomial_degree"],
        ratio,
        CASE["low_energy_ratio"],
    )
    assert float(ratio) == CASE["observed_energy_ratio"]
    assert scheduled_degree == 1

    result = continue_p18cyz_thinning_through_carried_degree1_iteration(
        p18cyz_result.next_thinning,
        image,
        centroids,
        field,
        previous_degree=CASE["previous_polynomial_degree"],
        polynomial_degree=CASE["maximum_polynomial_degree"],
        low_energy_ratio=CASE["low_energy_ratio"],
        **source_args,
    )
    exact = result.exact_iteration
    labels = p18cyz_result.next_thinning.thinning.segmentation
    memberships = p18cyz_result.next_thinning.branch.handoff.membership_pass.memberships
    closest = p18cyz.relations.closest_labels
    assert result.previous_degree == result.degree == 1
    assert result.energy_ratio == ratio
    assert exact.field_update.next_thinning is p18cyz_result.next_thinning
    assert exact.stage_order[-7:] == (
        "next_iteration.propagateThinning",
        "next_iteration.computeExact3DPolynomialField.degree1",
        "next_iteration.computeExactCentroids",
        "next_iteration.computeExactVariances",
        "next_iteration.computeImageVariance",
        "next_iteration.computeCentroidSimilarity",
        "next_iteration.computeExactMemberships",
    )

    field_input = tmp_path / "p18cza_field_input.bin"
    field_input.write_bytes(_field_input_bytes(labels, image, memberships, original_field, centroids))
    field_reference = json.loads(subprocess.run(
        ["sh", str(FIELD_HARNESS), str(field_input), *(str(v) for v in image.shape), "11"],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert field_reference["reference_kind"] == "formula_harness_not_source_class_execution"
    assert field_reference["inhomogeneity_upstream_source_sha256"] == FIELD_SOURCE_SHA256
    java_field = np.asarray(field_reference["field_bits"], dtype=np.int32).view(np.float32).reshape(image.shape)
    np.testing.assert_allclose(exact.field_update.fit.field, java_field, rtol=0, atol=2e-7)
    np.testing.assert_allclose(exact.field_update.fit.coefficients, field_reference["coefficients"], rtol=2e-6, atol=2e-8)

    iteration_input = tmp_path / "p18cza_iteration_input.bin"
    iteration_output = tmp_path / "p18cza_iteration_output.bin"
    power_output = tmp_path / "p18cza_power_output.bin"
    iteration_input.write_bytes(_iteration_input_bytes(
        labels, image, memberships, exact.field_update.fit.field, centroids, closest,
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

    arrays = np.fromfile(iteration_output, dtype=">f4").astype(np.float32)
    mem_pass = exact.membership_pass
    energy_count = mem_pass.class_data_energy.size
    membership_count = mem_pass.memberships.size
    denominator_count = mem_pass.denominator.size
    java_energy = arrays[:energy_count].reshape(mem_pass.class_data_energy.shape)
    java_membership = arrays[energy_count:energy_count + membership_count].reshape(mem_pass.memberships.shape)
    java_denominator = arrays[energy_count + membership_count:energy_count + membership_count + denominator_count].reshape(mem_pass.denominator.shape)
    np.testing.assert_array_equal(_bits(java_energy), _bits(mem_pass.class_data_energy))
    np.testing.assert_array_equal(_bits(java_membership), _bits(mem_pass.memberships))
    np.testing.assert_array_equal(_bits(java_denominator), _bits(mem_pass.denominator))
    assert reference["objective_bits"] == _bits([mem_pass.objective])[0]
