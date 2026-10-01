from __future__ import annotations

import json
import hashlib
from dataclasses import replace
from pathlib import Path
import re
import subprocess

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_initial_registration_control import register_initial_memberships_and_control
from pymacruise.toads_initial_thinning_composition import toads_initial_fielded_membership_to_thinning
from pymacruise.toads_initial_thinning_statistics_membership import continue_initial_thinning_through_memberships
from pymacruise.toads_default_registered_iteration import continue_registered_toads_iteration
from pymacruise.toads_control import toads_iteration_control
from pymacruise.toads_registration import exact_registration_measure_coefficients
from pymacruise.toads_registration_sampling import (
    vectorized_toads_registration_coefficients,
    vectorized_toads_registration_energy,
)


ROOT = Path(__file__).parent / "data/synthetic/toads_initial_thinning_startup"
JAVA_P18CYO_REFERENCE = ROOT / "continuation_java_reference.json"


def _load(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / name))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    return np.asarray(image.dataobj)


def _p18cyo():
    spec = json.loads((ROOT / "spec.json").read_text())
    p = spec["parameters"]
    image = _load("images_ras.nii.gz").astype(np.float32)
    segmentation = _load("segmentation_ras.nii.gz").astype(np.uint8)
    memberships = _load("initial_memberships_ras.nii.gz").astype(np.float32)
    fields = _load("fields_ras.nii.gz").astype(np.float32)
    bias = _load("psf_bias_ras.nii.gz").astype(np.float32)
    transformed = _load("transformed_shapes_ras.nii.gz").astype(np.float32)
    centers = np.asarray(p["centroids"], dtype=np.float32)
    exponents = np.asarray(p["exponent_by_channel_class"], dtype=np.float32)
    startup = toads_initial_fielded_membership_to_thinning(
        image_values=image, segmentation=segmentation, initial_memberships=memberships,
        centroids=centers, exponent_by_channel_class=exponents,
        modality_weights=p["modality_weights"], class_variance_weights=p["class_variance_weights"],
        image_variances=p["image_variances"], transformed_shapes=transformed,
        smoothing_priors=p["smoothing_priors"], centroid_similarity=p["centroid_similarity"],
        atlas_shape_priors=p["atlas_shape_priors"], smoothing=spec["membership_smoothing"],
        prior_factor=spec["prior_factor"], prior_coefficient=spec["prior_coefficient"],
        channel_fields=[fields[..., 0], fields[..., 1]], limit=spec["relation_limit"],
        volume_ratio=spec["thinning_volume_ratio"], thinning_limit=spec["thinning_limit"], psf_bias=bias,
    )
    category_exponents = np.asarray([exponents[:, 1], exponents[:, 3], exponents[:, 8]], dtype=np.float32).T
    continued = continue_initial_thinning_through_memberships(
        startup, image_values=image, initial_centroids=centers,
        category_exponents=category_exponents, transformed_shapes=transformed,
        image_variance_weights=p["image_variances"], smoothing_priors=p["smoothing_priors"],
        atlas_shape_priors=p["atlas_shape_priors"], exponent_by_channel_class=exponents,
        modality_weights=p["modality_weights"], prior_factor=spec["prior_factor"],
        prior_coefficient=spec["prior_coefficient"], smoothing=spec["membership_smoothing"],
        prior_scale=0.05, intensity_field=fields, psf_bias=bias,
    )
    return continued, transformed


def test_p18cyo_membership_and_objective_feed_registration_then_control(tmp_path: Path) -> None:
    continuation, atlas = _p18cyo()
    source_reference = json.loads(JAVA_P18CYO_REFERENCE.read_text())
    membership_digest = hashlib.sha256(
        np.asarray(continuation.membership_pass.memberships, dtype=">f4").tobytes(order="C")
    ).hexdigest()
    assert membership_digest == source_reference["updated_membership_sha256"]
    assert int(np.asarray(continuation.membership_pass.objective, dtype=np.float32).view(np.int32)) == source_reference["objective_bits"]
    membership_before = continuation.membership_pass.memberships.copy()
    objective = np.float32(continuation.membership_pass.objective)
    result = register_initial_memberships_and_control(
        continuation, atlas, image_resolutions_mm=(1, 1, 1), template_resolutions_mm=(1, 1, 1),
        initial_transform=(0, 0, 0, 0, 0, 0), previous_energy=np.float32(objective * 1.01),
        iteration=3, maximum_iterations=5, convergence_threshold=1e-4, do_growing=False,
    )
    assert result.handoff.membership_pass is continuation.membership_pass
    np.testing.assert_array_equal(result.handoff.membership_pass.memberships, membership_before)
    assert result.handoff.registration is not None
    assert len(result.handoff.registration.iterations) in range(1, 6)
    assert result.handoff.transform.shape == (6,)
    assert result.handoff.rotation.shape == (3, 3)
    assert len(result.handoff.shape_bounds) == 11
    assert result.control.energy_ratio == np.float32(
        np.abs(np.float32(np.float32(2) * np.float32(np.float32(np.float32(objective * 1.01) - objective)) /
                          np.float32(np.float32(objective * 1.01) + objective)))
    )
    assert result.stage_order[-4:] == (
        "atlas.registerShapes() -> registerAllShapes()", "energyRatio", "convergenceCheck", "iterationLimitCheck",
    )

    # Run the actual supplied VandyCRUISE atlas implementation headlessly on
    # these exact arrays. Stubs cover only MIPAV token parsing and progress UI.
    mem_path, shape_path = tmp_path / "memberships.bin", tmp_path / "shapes.bin"
    mem_path.write_bytes(np.asarray(membership_before, dtype="<f4").tobytes(order="C"))
    shape_path.write_bytes(np.asarray(atlas, dtype="<f4").tobytes(order="C"))
    harness = Path(__file__).parent / "java_reference/generate_toads_p18cyp_actual_registration_reference.sh"
    completed = subprocess.run(
        ["sh", str(harness), str(mem_path), str(shape_path), *map(str, membership_before.shape[:3])],
        check=True, capture_output=True, text=True,
    )
    java = json.loads(completed.stdout.splitlines()[-1])
    assert java["source"] == (
        "VandyCRUISE MultipleToadDeformableAtlas.updateShapeRegistration(rigid) + "
        "registerShapes() + private registerAllShapes()"
    )
    python_registration = result.handoff.registration
    assert python_registration is not None
    trace = completed.stdout.splitlines()
    java_hessians = [
        np.asarray([float(v) for v in re.findall(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?", line.split("vector:", 1)[1])], dtype=np.float32)
        for line in trace if line.startswith("H: vector:")
    ]
    java_gradients = [
        np.asarray([float(v) for v in re.findall(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?", line.split("vector:", 1)[1])], dtype=np.float32)
        for line in trace if line.startswith("G: vector:")
    ]
    java_trials = [
        tuple(float(value) for value in match.groups())
        for line in trace if (match := re.search(
            r"a: ([^ (]+)\(([^)]+)\)->\(([^)]+)\)", line
        ))
    ]
    assert len(java_hessians) == len(java_gradients) == len(java_trials) == len(python_registration.iterations) == 5
    dims = membership_before.shape[:3]
    low = np.zeros((11, 3), dtype=np.int32)
    high = np.broadcast_to(np.asarray(dims, dtype=np.int32), (11, 3)).copy()
    registered = np.ones(11, dtype=bool)
    registered[0] = False
    python_transform = np.zeros(6, dtype=np.float32)
    for index, (iteration_result, java_hessian, java_gradient, java_trial) in enumerate(zip(
        python_registration.iterations, java_hessians, java_gradients, java_trials, strict=True,
    )):
        coefficients = vectorized_toads_registration_coefficients(
            membership_before, atlas, image_resolutions_mm=(1, 1, 1),
            template_resolutions_mm=(1, 1, 1), transform=tuple(float(v) for v in python_transform),
            class_bounds_low_xyz=low, class_bounds_high_xyz=high,
            registered_classes=registered, offset=0, subsample=2,
        )
        gradient, hessian = exact_registration_measure_coefficients(
            coefficients.gradient, coefficients.hessian, coefficients.norm_gradient,
            coefficients.norm_hessian, coefficients.norm, source_faithful=True,
        )
        np.testing.assert_allclose(java_hessian, hessian, rtol=3e-5, atol=3e-6)
        np.testing.assert_allclose(java_gradient, gradient, rtol=3e-5, atol=3e-6)
        energy_before = vectorized_toads_registration_energy(
            membership_before, atlas, image_resolutions_mm=(1, 1, 1),
            template_resolutions_mm=(1, 1, 1), transform=tuple(float(v) for v in python_transform),
            class_bounds_low_xyz=low, class_bounds_high_xyz=high,
            registered_classes=registered, offset=0, subsample=2,
        )
        expected_lambda = np.float32(1.5**index)
        assert java_trial[0] == pytest.approx(float(expected_lambda), rel=2e-6)
        assert java_trial[1] == pytest.approx(float(energy_before), rel=3e-5, abs=3e-6)
        source_trial = np.asarray([
            np.float32(python_transform[i] + np.float32(np.float32(expected_lambda / hessian[i]) * gradient[i]))
            for i in range(6)
        ], dtype=np.float32)
        np.testing.assert_allclose(source_trial, iteration_result.transform, rtol=3e-5, atol=3e-6)
        assert java_trial[2] == pytest.approx(float(iteration_result.energy), rel=3e-5, abs=3e-6)
        python_transform = iteration_result.transform.copy()
    java_transform = np.asarray(java["transform_bits"], dtype=np.int32).view(np.float32)
    np.testing.assert_allclose(java_transform, python_registration.transform, rtol=3e-5, atol=3e-6)
    java_energy = np.asarray([java["energy_bits"]], dtype=np.int32).view(np.float32)
    np.testing.assert_allclose(java_energy, [python_registration.energy], rtol=3e-5, atol=3e-6)
    java_bounds = np.asarray(java["bounds"], dtype=np.int32)
    python_bounds = np.asarray(
        [[*bound.lower_xyz, *bound.upper_xyz] for bound in result.handoff.shape_bounds],
        dtype=np.int32,
    )
    np.testing.assert_array_equal(java_bounds, python_bounds)

    # Continue that exact registered membership state through the selected
    # non-stopping branch. Growth/relation kernels have their own Java-backed
    # references; this assertion pins their source order after registration.
    startup_thinning = continuation.startup_and_thinning
    thinning = startup_thinning.thinning
    relation_state = startup_thinning.startup.relation_state
    ordering = np.zeros(thinning.segmentation.shape, dtype=np.float32)
    for x in range(ordering.shape[0]):
        for y in range(ordering.shape[1]):
            for z in range(ordering.shape[2]):
                ordering[x, y, z] = np.float32(-((13 * x + 7 * y + 3 * z) % 19) / 10.0)
    branch_inputs = (
        thinning.segmentation.copy(),
        thinning.skeleton.copy(),
        ordering.copy(),
        thinning.class_volumes.copy(),
        relation_state.boundary_distances.copy(),
        result.handoff.membership_pass.memberships.copy(),
    )
    transition = continue_registered_toads_iteration(
        result,
        np.float32(objective * np.float32(1.01)),
        3,
        5,
        1e-4,
        thinning.segmentation,
        thinning.skeleton,
        ordering,
        thinning.class_volumes,
        do_growing=False,
        boundary_distances=relation_state.boundary_distances,
        relation_limit=json.loads((ROOT / "spec.json").read_text())["relation_limit"],
        available=thinning.available,
    )
    assert transition.registration_control is result
    assert transition.branch.control == result.control
    assert not transition.branch.control.stop
    assert transition.branch.growing.events
    assert transition.stage_order[-8:] == (
        "atlas.registerShapes() -> registerAllShapes()",
        "energyRatio",
        "convergenceCheck",
        "iterationLimitCheck",
        "propagateGrowing",
        "computeRelations",
        "computeBoundaryDistances",
        "computeClosestLabels",
    )
    for current, before in zip(
        (
            thinning.segmentation,
            thinning.skeleton,
            ordering,
            thinning.class_volumes,
            relation_state.boundary_distances,
            result.handoff.membership_pass.memberships,
        ),
        branch_inputs,
        strict=True,
    ):
        np.testing.assert_array_equal(current, before)

    # Pin the other wrapper branch on the same registered membership state:
    # at the iteration limit, disabled final growth goes straight to label
    # transfer, with no relation refresh.
    stopping_control = toads_iteration_control(
        np.float32(objective * np.float32(1.01)),
        result.handoff.membership_pass.objective,
        5,
        5,
        1e-4,
        do_growing=False,
    )
    stopping_registration = replace(result, control=stopping_control)
    template_labels = np.asarray([0, 4, 7, 9, 11, 13, 17, 19, 21, 23, 25], dtype=np.int8)
    previous_classification = np.full(thinning.segmentation.shape, -1, dtype=np.int8)
    stopped = continue_registered_toads_iteration(
        stopping_registration,
        np.float32(objective * np.float32(1.01)),
        5,
        5,
        1e-4,
        thinning.segmentation,
        thinning.skeleton,
        ordering,
        thinning.class_volumes,
        do_growing=False,
        template_labels=template_labels,
        previous_classification=previous_classification,
        available=thinning.available,
    )
    assert stopped.branch.control.stop
    assert stopped.branch.growing is None
    np.testing.assert_array_equal(
        stopped.branch.classification,
        template_labels[thinning.segmentation],
    )
    assert stopped.stage_order[-1] == "updateClassificationFromLabels"
