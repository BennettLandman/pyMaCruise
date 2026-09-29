from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.toads_registration import exact_registration_measure_coefficients
from pymacruise.toads_registration import exact_registration_gradient_descent_step
from pymacruise.toads_registration_sampling import (
    register_toads_atlas_rigid,
    vectorized_toads_registration_coefficients,
    vectorized_toads_registration_energy,
)


REFERENCE = Path(__file__).parent / "data/synthetic/toads_rigid_registration_update/java_reference.json"


def _bits(values) -> list[int]:
    return [int(value) for value in np.asarray(values, dtype=np.float32).ravel().view(np.int32)]


def _from_bits(values) -> np.ndarray:
    return np.asarray(values, dtype=np.int32).view(np.float32)


def _assert_float32_close(actual, expected_bits, *, atol: float = 2e-7) -> None:
    np.testing.assert_allclose(
        np.asarray(actual, dtype=np.float32),
        _from_bits(expected_bits),
        rtol=2e-6,
        atol=atol,
    )


def _inputs():
    atlas_shape = (5, 5, 5, 3)
    x, y, z, k = np.indices(atlas_shape, dtype=np.float32)
    atlas = (0.18 + 0.035 * x + 0.025 * y + 0.02 * z + 0.03 * k + 0.004 * x * y).astype(np.float32)
    dimensions = (4, 4, 4)
    ix, iy, iz = np.indices(dimensions, dtype=np.float32)
    memberships = np.empty(dimensions + (3,), dtype=np.float32)
    memberships[..., 0] = 0.1
    memberships[..., 1] = 0.25 + 0.04 * ix + 0.02 * iy
    memberships[..., 2] = 0.2 + 0.015 * iz + 0.01 * ix
    low = np.zeros((3, 3), dtype=np.int64)
    high = np.full((3, 3), 4, dtype=np.int64)
    kwargs = {
        "image_resolutions_mm": (1.0, 1.0, 1.0),
        "template_resolutions_mm": (1.0, 1.0, 1.0),
        "class_bounds_low_xyz": low,
        "class_bounds_high_xyz": high,
        "max_iterations": 1,
        "energy_tolerance": 0.0,
    }
    return memberships, atlas, kwargs


def test_full_rigid_registration_update_matches_java_source_harness() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source"] == "VandyCRUISE MultipleToadDeformableAtlas.computeExactRegistrationCoefficients + registerGradientDescent"
    assert reference["source_sha256"] == "e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2"
    assert reference["parametric_transform_sha256"] == "22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac"
    assert reference["rotation_matrix_sha256"] == "144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0"
    assert reference["image_functions_sha256"] == "be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892"
    assert set(reference["cases"]) == {"accepted", "backtracked"}
    memberships, atlas, kwargs = _inputs()
    for name, expected in reference["cases"].items():
        initial = _from_bits(expected["initial_transform"])
        coefficients = vectorized_toads_registration_coefficients(
            memberships,
            atlas,
            image_resolutions_mm=kwargs["image_resolutions_mm"],
            template_resolutions_mm=kwargs["template_resolutions_mm"],
            transform=tuple(initial),
            class_bounds_low_xyz=kwargs["class_bounds_low_xyz"],
            class_bounds_high_xyz=kwargs["class_bounds_high_xyz"],
        )
        gradient, hessian = exact_registration_measure_coefficients(
            coefficients.gradient,
            coefficients.hessian,
            coefficients.norm_gradient,
            coefficients.norm_hessian,
            coefficients.norm,
            source_faithful=True,
        )
        _assert_float32_close([coefficients.cost, coefficients.norm], expected["cost_norm"])
        _assert_float32_close(coefficients.gradient, expected["raw_gradient"])
        _assert_float32_close(coefficients.hessian, expected["raw_hessian"])
        _assert_float32_close(coefficients.norm_gradient, expected["norm_gradient"])
        _assert_float32_close(coefficients.norm_hessian, expected["norm_hessian"])
        _assert_float32_close(gradient, expected["gradient"])
        _assert_float32_close(hessian, expected["hessian"])

        initial_energy = vectorized_toads_registration_energy(
            memberships,
            atlas,
            image_resolutions_mm=kwargs["image_resolutions_mm"],
            template_resolutions_mm=kwargs["template_resolutions_mm"],
            transform=tuple(initial),
            class_bounds_low_xyz=kwargs["class_bounds_low_xyz"],
            class_bounds_high_xyz=kwargs["class_bounds_high_xyz"],
        )
        _assert_float32_close([initial_energy], [expected["initial_energy"]])
        python_trials = []

        def record_trial(candidate):
            energy = vectorized_toads_registration_energy(
                memberships,
                atlas,
                image_resolutions_mm=kwargs["image_resolutions_mm"],
                template_resolutions_mm=kwargs["template_resolutions_mm"],
                transform=tuple(candidate),
                class_bounds_low_xyz=kwargs["class_bounds_low_xyz"],
                class_bounds_high_xyz=kwargs["class_bounds_high_xyz"],
            )
            python_trials.append((candidate.copy(), energy))
            return energy

        line_search = exact_registration_gradient_descent_step(
            initial,
            gradient,
            hessian,
            initial_energy,
            record_trial,
        )
        assert len(python_trials) == len(expected["trials"])
        expected_lambda = np.float32(1.0)
        for lambda_bits, trial_bits, energy_bits, _improves in expected["trials"]:
            _assert_float32_close([expected_lambda], [lambda_bits], atol=1e-7)
            actual_trial, actual_trial_energy = python_trials.pop(0)
            _assert_float32_close(actual_trial, trial_bits)
            trial = _from_bits(trial_bits)
            trial_energy = vectorized_toads_registration_energy(
                memberships,
                atlas,
                image_resolutions_mm=kwargs["image_resolutions_mm"],
                template_resolutions_mm=kwargs["template_resolutions_mm"],
                transform=tuple(trial),
                class_bounds_low_xyz=kwargs["class_bounds_low_xyz"],
                class_bounds_high_xyz=kwargs["class_bounds_high_xyz"],
            )
            _assert_float32_close([trial_energy], [energy_bits])
            _assert_float32_close([actual_trial_energy], [energy_bits])
            assert bool(actual_trial_energy > initial_energy) == bool(_improves)
            expected_lambda = np.float32(
                expected_lambda * np.float32(1.5)
                if _improves
                else expected_lambda / np.float32(1.5)
            )
        assert not python_trials
        _assert_float32_close([line_search.lambda_value], [_bits([expected_lambda])[0]])

        result = register_toads_atlas_rigid(
            memberships,
            atlas,
            initial_transform=tuple(initial),
            **kwargs,
        )
        iteration = result.iterations[0]
        _assert_float32_close(result.transform, expected["transform"])
        _assert_float32_close([result.energy], [expected["energy"]])
        _assert_float32_close([iteration.lambda_value], [expected["lambda"]])
        _assert_float32_close([iteration.relative_change], [expected["relative_change"]])
        assert iteration.accepted is expected["accepted"]
        assert iteration.rejected_trials == expected["rejected"]

    assert reference["cases"]["accepted"]["accepted"] is True
    assert reference["cases"]["accepted"]["rejected"] == 0
    assert reference["cases"]["backtracked"]["accepted"] is False
    assert reference["cases"]["backtracked"]["rejected"] == 11
