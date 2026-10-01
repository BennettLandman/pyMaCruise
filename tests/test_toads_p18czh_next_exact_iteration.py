from __future__ import annotations

import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.toads_p18czf_registered_branch import (
    continue_p18cze_membership_through_registered_branch,
)
from pymacruise.toads_p18czg_next_loop_thinning import (
    continue_p18czf_branch_through_next_loop_thinning,
)
from pymacruise.toads_p18czh_next_exact_iteration import (
    continue_p18czg_thinning_through_next_exact_iteration_degree1,
)
from pymacruise.toads_p18cze_next_exact_iteration import (
    continue_p18czd_thinning_through_next_exact_iteration_degree1,
)
from test_toads_p18czb_registered_handoff import _atlas
from test_toads_p18czd_next_loop_thinning import _float_bits
from test_toads_p18cze_next_exact_iteration import _parameters, _p18czd_state
from test_toads_post_thinning_exact_iteration import (
    FIELD_HARNESS,
    ITERATION_HARNESS,
    _field_input_bytes,
    _iteration_input_bytes,
    _bits,
)

ROOT = Path(__file__).parent
CASE = json.loads(
    (ROOT / "data/synthetic/toads_p18czh_next_exact_iteration/case.json").read_text()
)
RAS_FIXTURE = ROOT / "data" / "synthetic" / CASE["source_fixture"]
FIELD_SOURCE_SHA256 = "6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8"
ITERATION_SOURCE_SHA256 = "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"


def _p18czg_state():
    image, _, p18czd = _p18czd_state()
    memberships = p18czd.next_thinning.branch.handoff.membership_pass.memberships
    exact = continue_p18czd_thinning_through_next_exact_iteration_degree1(
        p18czd,
        image,
        polynomial_degree=CASE["maximum_polynomial_degree"],
        low_energy_ratio=CASE["low_energy_ratio"],
        **_parameters(memberships.shape),
    )
    p18czf = continue_p18cze_membership_through_registered_branch(
        exact,
        _atlas(image),
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        iteration=6,
        maximum_iterations=8,
        convergence_threshold=0.0,
        do_growing=False,
        template_labels=np.asarray([0, 4, 7, 9, 11, 13, 17, 19, 21, 23, 25], dtype=np.int8),
        previous_classification=np.full(image.shape, -1, dtype=np.int8),
        relation_limit=1,
        subsample=2,
        max_registration_iterations=5,
        registration_energy_tolerance=0.0,
        minimum_lambda=1e-3,
    )
    shapes = np.zeros_like(
        p18czf.registered_iteration.branch.handoff.membership_pass.memberships,
        dtype=np.float32,
    )
    p18czg = continue_p18czf_branch_through_next_loop_thinning(p18czf, shapes)
    return image, p18czf, p18czg


def test_p18czg_thinning_composes_through_next_degree_one_exact_iteration(tmp_path: Path) -> None:
    image, p18czf, p18czg = _p18czg_state()
    source = nib.load(str(RAS_FIXTURE))
    assert nib.aff2axcodes(source.affine) == tuple(CASE["source_fixture_axis_codes"])
    affine = source.affine @ np.diag([1 / CASE["spatial_repeat_factor"]] * 3 + [1])
    assert nib.aff2axcodes(affine) == ("R", "A", "S")

    previous_exact = p18czf.exact_iteration.exact_iteration
    original_field = np.asarray(
        previous_exact.exact_iteration.field_update.fit.field, dtype=np.float32
    ).copy()
    initial_centroids = np.asarray(
        previous_exact.exact_iteration.centroids[0], dtype=np.float32
    ).copy()
    branch = p18czg.next_thinning.branch
    memberships = branch.handoff.membership_pass.memberships
    result = continue_p18czg_thinning_through_next_exact_iteration_degree1(
        p18czg,
        image,
        polynomial_degree=CASE["maximum_polynomial_degree"],
        low_energy_ratio=CASE["low_energy_ratio"],
        **_parameters(memberships.shape),
    )
    exact = result.exact_iteration
    assert exact.previous_degree == CASE["previous_degree"]
    assert exact.degree == CASE["previous_degree"]
    assert exact.energy_ratio == branch.control.energy_ratio
    assert exact.exact_iteration.field_update.next_thinning is p18czg.next_thinning
    assert exact.stage_order[-7:] == (
        "next_iteration.propagateThinning",
        "next_iteration.computeExact3DPolynomialField.degree1",
        "next_iteration.computeExactCentroids",
        "next_iteration.computeExactVariances",
        "next_iteration.computeImageVariance",
        "next_iteration.computeCentroidSimilarity",
        "next_iteration.computeExactMemberships",
    )
    assert nib.Nifti1Image(exact.exact_iteration.field_update.fit.field, affine).shape == image.shape

    labels = p18czg.next_thinning.thinning.segmentation
    field_input = tmp_path / "p18czh_field_input.bin"
    field_input.write_bytes(
        _field_input_bytes(labels, image, memberships, original_field, initial_centroids)
    )
    field_reference = json.loads(subprocess.run(
        ["sh", str(FIELD_HARNESS), str(field_input), *(str(v) for v in image.shape), "11"],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert field_reference["reference_kind"] == "formula_harness_not_source_class_execution"
    assert field_reference["inhomogeneity_upstream_source_sha256"] == FIELD_SOURCE_SHA256
    java_field = np.asarray(field_reference["field_bits"], dtype=np.int32).view(np.float32).reshape(image.shape)
    np.testing.assert_allclose(exact.exact_iteration.field_update.fit.field, java_field, rtol=0, atol=2e-7)
    np.testing.assert_allclose(
        exact.exact_iteration.field_update.fit.coefficients,
        field_reference["coefficients"],
        rtol=2e-6,
        atol=2e-8,
    )

    iteration_input = tmp_path / "p18czh_iteration_input.bin"
    iteration_output = tmp_path / "p18czh_iteration_output.bin"
    power_output = tmp_path / "p18czh_power_output.bin"
    iteration_input.write_bytes(_iteration_input_bytes(
        labels,
        image,
        memberships,
        exact.exact_iteration.field_update.fit.field,
        initial_centroids,
        branch.relations.closest_labels,
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
    computed = exact.exact_iteration
    assert reference["centroid_bits"] == _bits(computed.centroids)
    assert reference["standard_deviation_bits"] == _bits(computed.standard_deviations)
    assert reference["similarity_bits"] == _bits(computed.centroid_similarity)
    assert reference["image_variance_bits"] == int(computed.image_variance_weights.view(np.int32)[0])
    assert reference["image_std_bits"] == int(computed.image_standard_deviation.view(np.int32)[0])
    arrays = np.fromfile(iteration_output, dtype=">f4").astype(np.float32)
    membership_pass = computed.membership_pass
    energy_count = membership_pass.class_data_energy.size
    membership_count = membership_pass.memberships.size
    denominator_count = membership_pass.denominator.size
    java_energy = arrays[:energy_count].reshape(membership_pass.class_data_energy.shape)
    java_membership = arrays[energy_count:energy_count + membership_count].reshape(membership_pass.memberships.shape)
    java_denominator = arrays[energy_count + membership_count:energy_count + membership_count + denominator_count].reshape(membership_pass.denominator.shape)
    np.testing.assert_array_equal(_bits(java_energy), _bits(membership_pass.class_data_energy))
    np.testing.assert_array_equal(_bits(java_membership), _bits(membership_pass.memberships))
    np.testing.assert_array_equal(_bits(java_denominator), _bits(membership_pass.denominator))
    assert reference["objective_bits"] == _bits([membership_pass.objective])[0]
