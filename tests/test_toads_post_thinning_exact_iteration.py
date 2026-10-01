from __future__ import annotations

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import re
import subprocess

import numpy as np
import pytest

from pymacruise.toads_post_registration_next_thinning import toads_post_registration_next_thinning
from pymacruise.toads_post_thinning_exact_iteration import toads_post_thinning_exact_iteration_degree1
from pymacruise.toads_default_membership_registration_control import (
    continue_default_membership_through_registration_control,
)
from pymacruise.toads_math import toads_powertable_lookup
from pymacruise.toads_registration_sampling import (
    vectorized_toads_registration_coefficients,
    vectorized_toads_registration_energy,
)
from pymacruise.toads_registration import exact_registration_measure_coefficients
from test_toads_post_registration_next_thinning import _registered_branch


ROOT = Path(__file__).parent
FIELD_HARNESS = ROOT.parent / "scripts/run_java_toads_post_thinning_exact_field_reference.sh"
ITERATION_HARNESS = ROOT.parent / "scripts/run_java_toads_post_thinning_exact_iteration_reference.sh"
SOURCE_SHA256 = "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
POWER_TABLE_SHA256 = "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"
POWER_TABLE_COPY = ROOT / "java_reference/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java"


def _field_input_bytes(labels, image, memberships, field, centroids) -> bytes:
    return b"".join((
        np.asarray(labels, dtype=np.int8).tobytes(order="C"),
        np.asarray(image, dtype=">f4").tobytes(order="C"),
        np.asarray(memberships, dtype=">f4").tobytes(order="C"),
        np.asarray(field, dtype=">f4").tobytes(order="C"),
        np.asarray(centroids, dtype=">f4").tobytes(order="C"),
    ))


def _iteration_input_bytes(labels, image, memberships, field, centroids, closest) -> bytes:
    return b"".join((
        np.asarray(labels, dtype=np.int8).tobytes(order="C"),
        np.asarray(image, dtype=">f4").tobytes(order="C"),
        np.asarray(memberships, dtype=">f4").tobytes(order="C"),
        np.asarray(field, dtype=">f4").tobytes(order="C"),
        np.asarray(centroids, dtype=">f4").tobytes(order="C"),
        np.asarray(closest, dtype=np.int8).tobytes(order="C"),
    ))


def _bits(values) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).astype(int).tolist()


def _float_from_bits(bits: int) -> np.float32:
    return np.asarray([bits], dtype=np.uint32).view(np.float32)[0]


def _continuation_state():
    current = _registered_branch()
    energy = current.handoff.membership_pass.objective
    current_energy = np.float32(energy * np.float32(1.005))
    branch = _registered_branch(
        previous_energy=current_energy,
        convergence_threshold=np.float32(0.0001),
    )
    expand = lambda value: np.repeat(np.repeat(np.repeat(value, 3, axis=0), 3, axis=1), 3, axis=2)
    grown = replace(
        branch.growing,
        segmentation=expand(branch.growing.segmentation),
        skeleton=expand(branch.growing.skeleton),
        available=expand(branch.growing.available),
        ordering=expand(branch.growing.ordering),
        class_volumes=np.asarray(branch.growing.class_volumes * np.float32(27.0), dtype=np.float32),
    )
    membership = branch.handoff.membership_pass
    membership = replace(
        membership,
        class_data_energy=expand(membership.class_data_energy),
        memberships=expand(membership.memberships),
        denominator=expand(membership.denominator),
    )
    branch = replace(branch, handoff=replace(branch.handoff, membership_pass=membership), growing=grown)
    next_thinning = toads_post_registration_next_thinning(
        branch, np.zeros_like(membership.memberships, dtype=np.float32)
    )
    nx, ny, nz = next_thinning.thinning.segmentation.shape
    x, y, z = np.meshgrid(
        np.linspace(0.0, 1.0, nx, dtype=np.float32),
        np.linspace(0.0, 1.0, ny, dtype=np.float32),
        np.linspace(0.0, 1.0, nz, dtype=np.float32),
        indexing="ij",
    )
    image = np.asarray(0.35 + 0.31 * x + 0.19 * y + 0.13 * z + 0.07 * x * y, dtype=np.float32)
    centroids = np.asarray([0.0, 0.20, 0.31, 0.39, 0.46, 0.52, 0.58, 0.63, 0.69, 0.75, 0.82], dtype=np.float32)
    field = np.ones((nx, ny, nz), dtype=np.float32)
    closest = np.empty((nx, ny, nz, 4), dtype=np.int8)
    for i in range(4):
        closest[..., i] = (x * 10 + y * 4 + z * 3 + i) % 11
    return next_thinning, image, centroids, field, closest


def test_selected_single_channel_exact_refit_statistics_membership_sequence(tmp_path: Path) -> None:
    next_thinning, image, centroids, field, closest = _continuation_state()
    shape = np.zeros((*image.shape, 11), dtype=np.float32)
    exponent_map = np.empty((1, 11), dtype=np.float32)
    exponent_map[0, 0] = 1.0
    exponent_map[0, 1:3] = 0.8
    exponent_map[0, 3:] = 1.0
    memberships = next_thinning.branch.handoff.membership_pass.memberships
    field_input = tmp_path / "field_input.bin"
    field_input.write_bytes(_field_input_bytes(
        next_thinning.thinning.segmentation, image, memberships, field, centroids,
    ))
    field_result = json.loads(subprocess.run(
        ["sh", str(FIELD_HARNESS), str(field_input), *(str(value) for value in image.shape), "11"],
        check=True, capture_output=True, text=True,
    ).stdout)
    java_field = np.asarray(field_result["field_bits"], dtype=np.int32).view(np.float32).reshape(image.shape)
    result = toads_post_thinning_exact_iteration_degree1(
        next_thinning,
        image,
        centroids,
        field,
        category_exponents=np.asarray([0.8, 1.0, 1.0], dtype=np.float32),
        exponent_by_channel_class=exponent_map,
        modality_weights=np.asarray([0.73], dtype=np.float32),
        class_variance_weights=np.asarray([[0.81, 0.93, 1.17, 0.76, 0.88, 1.02, 1.14, 0.69, 0.96, 1.21, 0.73]], dtype=np.float32),
        transformed_shapes=shape,
        closest_labels=closest,
        smoothing_priors=np.asarray([0.65, 0.71, 0.84, 0.91, 1.02, 1.13, 0.76, 0.88, 1.05, 1.17, 0.82], dtype=np.float32),
        atlas_shape_priors=np.asarray([0.10, 0.08, 0.09, 0.11, 0.12, 0.09, 0.10, 0.11, 0.07, 0.06, 0.07], dtype=np.float32),
        smoothing=0.043,
        prior_factor=0.64,
        prior_coefficient=0.16,
        prior_scale=0.05,
    )
    field_delta = float(np.max(np.abs(java_field - result.field_update.fit.field)))
    assert field_delta <= 2e-7, f"max exact-field difference was {field_delta:.9g}"
    iteration_input = tmp_path / "iteration_input.bin"
    iteration_output = tmp_path / "iteration_output.bin"
    iteration_input.write_bytes(_iteration_input_bytes(
        next_thinning.thinning.segmentation, image, memberships, result.field_update.fit.field, centroids, closest,
    ))
    power_output = tmp_path / "power_lookup_output.bin"
    java = json.loads(subprocess.run(
        ["sh", str(ITERATION_HARNESS), str(iteration_input), *(str(value) for value in image.shape), str(iteration_output), str(power_output)],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert result.stage_order[-5:] == (
        "next_iteration.computeExactCentroids",
        "next_iteration.computeExactVariances",
        "next_iteration.computeImageVariance",
        "next_iteration.computeCentroidSimilarity",
        "next_iteration.computeExactMemberships",
    )
    assert result.field_update.degree == 1
    assert result.centroids.shape == result.standard_deviations.shape == (1, 11)
    assert result.centroid_similarity.shape == (11, 11)
    assert result.membership_pass.memberships.shape == (*image.shape, 11)
    assert np.isfinite(result.centroids).all()
    assert np.isfinite(result.standard_deviations).all()
    assert np.isfinite(result.membership_pass.memberships).all()
    assert field_result["inhomogeneity_upstream_source_sha256"] == "6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8"
    assert field_result["wrapper_upstream_source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert java["source_sha256"] == SOURCE_SHA256
    assert java["reference_kind"] == "source_formula_harness"
    assert java["power_table_source_sha256"] == "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678"
    assert java["source_order"] == [
        "computeCorrectionField.degree1", "computeExactCentroids",
        "computeExactVariances", "computeImageVariance", "centroid similarity",
        "computeExactMemberships",
    ]
    assert java["shape"] == list(image.shape)
    assert java["centroid_bits"] == _bits(result.centroids)
    assert java["standard_deviation_bits"] == _bits(result.standard_deviations)
    assert java["similarity_bits"] == _bits(result.centroid_similarity)
    assert java["image_variance_bits"] == int(result.image_variance_weights.view(np.int32)[0])
    assert java["image_std_bits"] == int(result.image_standard_deviation.view(np.int32)[0])
    java_arrays = np.fromfile(iteration_output, dtype=">f4").astype(np.float32)
    energy_count = result.membership_pass.class_data_energy.size
    membership_count = result.membership_pass.memberships.size
    denominator_count = result.membership_pass.denominator.size
    java_energy = java_arrays[:energy_count].reshape(result.membership_pass.class_data_energy.shape)
    java_memberships = java_arrays[energy_count:energy_count + membership_count].reshape(result.membership_pass.memberships.shape)
    java_denominator = java_arrays[energy_count + membership_count:energy_count + membership_count + denominator_count].reshape(result.membership_pass.denominator.shape)
    java_power = np.fromfile(power_output, dtype=">f8").reshape(memberships.shape)
    python_power = np.zeros(memberships.shape, dtype=np.float64)
    corrected = np.asarray(result.field_update.fit.field * image, dtype=np.float32)
    for class_index in range(1, 11):
        exponent = float(exponent_map[0, class_index])
        residual = np.asarray(corrected[1:-1, 1:-1, 1:-1] - result.centroids[0, class_index], dtype=np.float32)
        distance = np.asarray(residual * residual, dtype=np.float32)
        python_power[1:-1, 1:-1, 1:-1, class_index] = toads_powertable_lookup(
            distance, exponent, cached_exponent=exponent,
        )
    power_delta = float(np.max(np.abs(java_power - python_power)))
    energy_delta = float(np.max(np.abs(java_energy - result.membership_pass.class_data_energy)))
    denominator_delta = float(np.max(np.abs(java_denominator - result.membership_pass.denominator)))
    membership_delta = float(np.max(np.abs(java_memberships - result.membership_pass.memberships)))
    print(
        "P18 exact continuation observed max absolute differences: "
        f"field={field_delta:.9g}, PowerTable lookup={power_delta:.17g}, "
        f"same-input class energy={energy_delta:.9g}, denominator={denominator_delta:.9g}, "
        f"memberships={membership_delta:.9g}"
    )
    # Math.pow and NumPy power differ in double precision by a measured
    # 2.78e-17 maximum here; membership-stage float32 arrays happen to round
    # identically for this fixture. That does not imply bitwise identity for
    # other inputs or runtimes.
    np.testing.assert_array_equal(_bits(java_energy), _bits(result.membership_pass.class_data_energy))
    np.testing.assert_array_equal(_bits(java_denominator), _bits(result.membership_pass.denominator))
    np.testing.assert_array_equal(_bits(java_memberships), _bits(result.membership_pass.memberships))
    assert java["objective_bits"] == int(np.asarray([result.membership_pass.objective], dtype=np.float32).view(np.int32)[0])
    harness = ROOT / "java_reference/ToadsPostThinningExactIterationReferenceHarness.java"
    assert hashlib.sha256(harness.read_bytes()).hexdigest() == "e005b859f9924dfbe3f91a42a0faf60e5342f0b7bd2f507db4397561b96f886c"
    assert hashlib.sha256(POWER_TABLE_COPY.read_bytes()).hexdigest() == POWER_TABLE_SHA256
    assert power_delta <= 1e-14, f"max PowerTable-vs-NumPy lookup difference was {power_delta:.17g}"
    assert energy_delta == 0.0, f"max same-input membership class-energy difference was {energy_delta:.9g}"


def test_p18cyt_registration_methods_and_source_inspected_control_order(tmp_path: Path) -> None:
    """Check Java registration methods and the inspected one-channel wrapper seam."""
    next_thinning, image, centroids, field, closest = _continuation_state()
    shape = np.zeros((*image.shape, 11), dtype=np.float32)
    exponent_map = np.empty((1, 11), dtype=np.float32)
    exponent_map[0, 0] = 1.0
    exponent_map[0, 1:3] = 0.8
    exponent_map[0, 3:] = 1.0
    memberships = next_thinning.branch.handoff.membership_pass.memberships
    result = toads_post_thinning_exact_iteration_degree1(
        next_thinning,
        image,
        centroids,
        field,
        category_exponents=np.asarray([0.8, 1.0, 1.0], dtype=np.float32),
        exponent_by_channel_class=exponent_map,
        modality_weights=np.asarray([0.73], dtype=np.float32),
        class_variance_weights=np.asarray(
            [[0.81, 0.93, 1.17, 0.76, 0.88, 1.02, 1.14, 0.69, 0.96, 1.21, 0.73]],
            dtype=np.float32,
        ),
        transformed_shapes=shape,
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
    membership_pass = result.membership_pass
    # This is the actual P18cyt one-channel output: one intensity input and
    # eleven statistical membership maps. Give each atlas class a smooth,
    # nonconstant prior so rigid registration exercises spatial derivatives.
    nx, ny, nz = membership_pass.memberships.shape[:3]
    gx, gy, gz = np.meshgrid(
        np.linspace(-1.0, 1.0, nx, dtype=np.float32),
        np.linspace(-1.0, 1.0, ny, dtype=np.float32),
        np.linspace(-1.0, 1.0, nz, dtype=np.float32),
        indexing="ij",
    )
    atlas = np.empty((nx, ny, nz, 11), dtype=np.float32)
    for class_index in range(11):
        cx = np.float32(-0.45 + 0.09 * class_index)
        cy = np.float32(0.3 * np.sin(class_index))
        cz = np.float32(0.25 * np.cos(class_index))
        radius = np.float32((gx - cx) ** 2 + 1.3 * (gy - cy) ** 2 + 0.8 * (gz - cz) ** 2)
        atlas[..., class_index] = np.asarray(
            np.exp(-radius / np.float32(0.42)) * np.float32(0.75) + np.float32(0.05),
            dtype=np.float32,
        )

    objective = np.float32(membership_pass.objective)
    previous_energy = np.float32(objective * np.float32(1.01))
    ratio_harness = ROOT / "java_reference/ToadsP18cyqEnergyControlReferenceHarness.java"
    assert hashlib.sha256(ratio_harness.read_bytes()).hexdigest() == (
        "d831bae542aaa7ae79d9c741d93aa7fcc3265e9e32bc5dffc88abb30d12dfcdc"
    )
    subprocess.run(
        ["javac", "-d", str(tmp_path), str(ratio_harness)],
        check=True,
        capture_output=True,
        text=True,
    )

    def float_bits(value: float) -> str:
        return str(int(np.asarray([value], dtype=np.float32).view(np.int32)[0]))

    java_control = subprocess.run(
        [
            "java", "-cp", str(tmp_path), "ToadsP18cyqEnergyControlReferenceHarness",
            float_bits(previous_energy), float_bits(objective), float_bits(1e-4), "3", "5", "false",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    java_control_reference = json.loads(java_control.stdout)
    assert java_control_reference["call_order"] == [
        "computeExactMemberships", "energyRatio", "registerShapes",
        "convergenceCheck", "iterationLimitCheck",
    ]
    common = dict(
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        initial_transform=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        previous_energy=previous_energy,
        iteration=3,
        maximum_iterations=5,
        convergence_threshold=1e-4,
        do_growing=False,
        offset=0,
        subsample=2,
        max_registration_iterations=5,
        registration_energy_tolerance=0.0,
        minimum_lambda=1e-3,
    )
    aligned = continue_default_membership_through_registration_control(
        membership_pass, atlas, **common,
    )
    ratio = np.float32(
        np.abs(
            np.float32(
                np.float32(np.float32(2.0) * np.float32(previous_energy - objective))
                / np.float32(previous_energy + objective)
            )
        )
    )
    assert aligned.control.energy_ratio.view(np.int32) == ratio.view(np.int32)
    assert java_control_reference["ratio_bits"] == int(ratio.view(np.int32))
    assert java_control_reference["stop"] is aligned.control.stop
    assert java_control_reference["run_growing"] is aligned.control.run_growing
    assert aligned.stage_order == (
        "computeExactMemberships",
        "energyRatio",
        "atlas.registerShapes() -> registerAllShapes()",
        "convergenceCheck",
        "iterationLimitCheck",
    )
    assert aligned.handoff.membership_pass is membership_pass
    assert aligned.handoff.registration is not None
    assert not aligned.control.stop

    unaligned = continue_default_membership_through_registration_control(
        membership_pass, atlas, align_shapes=False, **common,
    )
    assert unaligned.handoff.registration is None
    assert unaligned.control.energy_ratio.view(np.int32) == ratio.view(np.int32)
    assert unaligned.stage_order == (
        "computeExactMemberships", "energyRatio", "convergenceCheck", "iterationLimitCheck",
    )
    limited = continue_default_membership_through_registration_control(
        membership_pass,
        atlas,
        align_shapes=False,
        **{**common, "iteration": 5},
    )
    assert limited.control.stop and not limited.control.run_growing

    mem_path, shape_path = tmp_path / "p18cyt_memberships.bin", tmp_path / "p18cyt_atlas.bin"
    mem_path.write_bytes(np.asarray(membership_pass.memberships, dtype="<f4").tobytes(order="C"))
    shape_path.write_bytes(np.asarray(atlas, dtype="<f4").tobytes(order="C"))
    harness = ROOT / "java_reference/generate_toads_p18cyp_actual_registration_reference.sh"
    completed = subprocess.run(
        ["sh", str(harness), str(mem_path), str(shape_path), str(nx), str(ny), str(nz)],
        check=True,
        capture_output=True,
        text=True,
    )
    java = json.loads(completed.stdout.splitlines()[-1])
    assert java["source"] == (
        "VandyCRUISE MultipleToadDeformableAtlas.updateShapeRegistration(rigid) + "
        "registerShapes() + private registerAllShapes()"
    )
    registration = aligned.handoff.registration
    assert registration is not None
    # The harness runs the same registration source with these exact P18cyt
    # memberships and atlas values, exposing per-iteration derivative and
    # line-search traces on stdout before its terminal transform/energy.
    trace = completed.stdout.splitlines()
    number = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?"
    java_hessians = [
        np.asarray([float(v) for v in re.findall(number, line.split("vector:", 1)[1])], dtype=np.float32)
        for line in trace if line.startswith("H: vector:")
    ]
    java_gradients = [
        np.asarray([float(v) for v in re.findall(number, line.split("vector:", 1)[1])], dtype=np.float32)
        for line in trace if line.startswith("G: vector:")
    ]
    java_trials = [
        tuple(float(value) for value in match.groups())
        for line in trace if (match := re.search(r"a: ([^ (]+)\(([^)]+)\)->\(([^)]+)\)", line))
    ]
    assert len(java_hessians) == len(java_gradients) == len(java_trials) == len(registration.iterations) == 5
    low = np.zeros((11, 3), dtype=np.int32)
    high = np.broadcast_to(np.asarray((nx, ny, nz), dtype=np.int32), (11, 3)).copy()
    registered_classes = np.asarray([False] + [True] * 10)
    current_transform = np.zeros(6, dtype=np.float32)
    for index, (iteration_result, java_hessian, java_gradient, java_trial) in enumerate(zip(
        registration.iterations, java_hessians, java_gradients, java_trials, strict=True,
    )):
        coeff = vectorized_toads_registration_coefficients(
            membership_pass.memberships,
            atlas,
            image_resolutions_mm=(1.0, 1.0, 1.0),
            template_resolutions_mm=(1.0, 1.0, 1.0),
            transform=tuple(float(v) for v in current_transform),
            class_bounds_low_xyz=low,
            class_bounds_high_xyz=high,
            registered_classes=registered_classes,
            offset=0,
            subsample=2,
        )
        gradient, hessian = exact_registration_measure_coefficients(
            coeff.gradient, coeff.hessian, coeff.norm_gradient, coeff.norm_hessian, coeff.norm,
            source_faithful=True,
        )
        np.testing.assert_allclose(java_hessian, hessian, rtol=3e-5, atol=3e-6)
        np.testing.assert_allclose(java_gradient, gradient, rtol=3e-5, atol=3e-6)
        energy_before = vectorized_toads_registration_energy(
            membership_pass.memberships,
            atlas,
            image_resolutions_mm=(1.0, 1.0, 1.0),
            template_resolutions_mm=(1.0, 1.0, 1.0),
            transform=tuple(float(v) for v in current_transform),
            class_bounds_low_xyz=low,
            class_bounds_high_xyz=high,
            registered_classes=registered_classes,
            offset=0,
            subsample=2,
        )
        assert java_trial[0] == pytest.approx(float(np.float32(1.5**index)), rel=2e-6)
        assert java_trial[1] == pytest.approx(float(energy_before), rel=3e-5, abs=3e-6)
        np.testing.assert_allclose(java_trial[2], iteration_result.energy, rtol=3e-5, atol=3e-6)
        current_transform = iteration_result.transform.copy()
    java_transform = np.asarray(java["transform_bits"], dtype=np.int32).view(np.float32)
    java_energy = np.asarray([java["energy_bits"]], dtype=np.int32).view(np.float32)
    np.testing.assert_allclose(java_transform, registration.transform, rtol=3e-5, atol=3e-6)
    np.testing.assert_allclose(java_energy, [registration.energy], rtol=3e-5, atol=3e-6)
