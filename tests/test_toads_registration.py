from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_registration import (
    exact_membership_registration_energy,
    exact_membership_registration_coefficients,
    exact_registration_measure_coefficients,
    exact_registration_coefficient_contributions,
    exact_registration_gradient_descent_step,
    rigid_image_to_template_jacobian,
    synthetic_registration_objective_inputs,
    trilinear_prior_derivatives,
)


ROOT = Path(__file__).parent
REFERENCE = ROOT / "data" / "synthetic" / "toads_registration_energy" / "java_reference.json"
COEFFICIENT_REFERENCE = (
    ROOT / "data" / "synthetic" / "toads_registration_energy" / "java_coefficient_reference.json"
)
HARNESS_HASH = "3c2e8bc403a1e04f201cbc29f6aa468175b202c4d2013bc1aefa4429d7bf8b41"
COEFFICIENT_HARNESS_HASH = "cc62e5e25bd3297a0d03a981c951cea0f7bb137810fe0791b2952114c39af1a5"
DERIVATIVE_HARNESS_HASH = "127099ca604f07537992594a9dbe871ce0efe0fdf3a3106277ae91c72f85a5b0"
RIGID_JACOBIAN_HARNESS_HASH = "e04a9319eae051362547f0461100512a0bcf58e0d7f70763c300238f939bbcc4"
COEFFICIENT_ACCUMULATION_HARNESS_HASH = "20ac7238c1eaac7c6a9b72a2a974cca58582c0eaf2ab8fd7dead86755d7aa2d9"
LINE_SEARCH_HARNESS_HASH = "bbe218ebe24b429b032ff023fb53bf0e65da64997d8567c54acd69a8503c3253"
FULL_COEFFICIENT_HARNESS_HASH = "8f0d04367bc7e70a2ca9f9b48ad6331903f256a8c30f811e528530bb4d889246"


def test_exact_membership_energy_matches_source_formula_harness() -> None:
    expected = json.loads(REFERENCE.read_text())
    memberships, priors = synthetic_registration_objective_inputs()
    cases = {
        "all_registered": (None, 0, 1),
        "selected_classes_subsampled": ([False, True, False, True], 1, 2),
    }
    for name, (registered, offset, subsample) in cases.items():
        energy = exact_membership_registration_energy(
            memberships,
            priors,
            registered_classes=registered,
            offset=offset,
            subsample=subsample,
        )
        assert energy == np.float32(expected[name]["energy"]), name


def test_exact_membership_energy_preserves_source_empty_norm_nan() -> None:
    memberships = np.ones((3, 3, 3, 2), dtype=np.float32)
    priors = np.zeros_like(memberships)
    assert np.isnan(exact_membership_registration_energy(memberships, priors))


def test_registration_reference_sources_have_recorded_hashes() -> None:
    for name, expected in (
        ("ToadsRegistrationEnergyReferenceHarness.java", HARNESS_HASH),
        ("ToadsRegistrationCoefficientReferenceHarness.java", COEFFICIENT_HARNESS_HASH),
        ("ToadsPriorDerivativeReferenceHarness.java", DERIVATIVE_HARNESS_HASH),
        ("ToadsRigidJacobianReferenceHarness.java", RIGID_JACOBIAN_HARNESS_HASH),
        ("ToadsRegistrationCoefficientAccumulationHarness.java", COEFFICIENT_ACCUMULATION_HARNESS_HASH),
        ("ToadsRegistrationLineSearchHarness.java", LINE_SEARCH_HARNESS_HASH),
        ("ToadsRegistrationFullCoefficientsHarness.java", FULL_COEFFICIENT_HARNESS_HASH),
    ):
        harness = ROOT / "java_reference" / name
        assert hashlib.sha256(harness.read_bytes()).hexdigest() == expected


def test_registration_coefficients_preserve_source_normalization_call_order() -> None:
    expected = json.loads(COEFFICIENT_REFERENCE.read_text())
    gradient, hessian = exact_registration_measure_coefficients(
        [0.75, -0.125, 3.5],
        [2.0, 0.5, 8.25],
        [12.0, 4.0, -1.0],
        [9.0, 8.0, 7.0],
        2.5,
    )
    np.testing.assert_array_equal(gradient, np.asarray(expected["gradient"], dtype=np.float32))
    np.testing.assert_array_equal(hessian, np.asarray(expected["hessian"], dtype=np.float32))


def test_registration_coefficient_correction_is_opt_in() -> None:
    gradient, hessian = exact_registration_measure_coefficients(
        [0.75, -0.125, 3.5],
        [2.0, 0.5, 8.25],
        [12.0, 4.0, -1.0],
        [9.0, 8.0, 7.0],
        2.5,
        source_faithful=False,
    )
    np.testing.assert_array_equal(gradient, np.asarray([0.3, -0.05, 1.4], dtype=np.float32))
    np.testing.assert_array_equal(hessian, np.asarray([0.8, 0.2, 3.3], dtype=np.float32))


def test_prior_derivatives_match_selected_java_formula_excerpt() -> None:
    expected = json.loads(
        (REFERENCE.parent / "java_prior_derivative_reference.json").read_text()
    )
    x, y, z = np.indices((5, 5, 5), dtype=np.int32)
    shape = (((7 * x + 3 * y + 5 * z + 2 * x * y + 4 * y * z) % 23) / np.float32(22.0)).astype(np.float32)
    points = ((0.25, 1.5, 2.75), (3.0, 3.0, 3.0), (3.25, 1.0, 1.0), (-0.1, 1.0, 1.0))
    for point, reference in zip(points, expected["derivatives"], strict=True):
        np.testing.assert_array_equal(
            trilinear_prior_derivatives(shape, point),
            np.asarray(reference, dtype=np.float32),
        )


def test_rigid_coordinate_jacobian_matches_selected_java_formula_excerpt() -> None:
    expected = json.loads(
        (REFERENCE.parent / "java_rigid_jacobian_reference.json").read_text()
    )
    actual = rigid_image_to_template_jacobian(
        (2, 1, 3),
        (0.14, 0.07, -0.05, 1.0, 2.0, 3.0),
        (2.5, 2.0, 1.5),
        (0.8, 1.2, 2.0),
        (0.9, 1.1, 1.3),
    )
    np.testing.assert_array_equal(actual, np.asarray(expected["jacobian"], dtype=np.float32))


def test_rigid_coordinate_jacobian_rejects_singular_rotation_derivative() -> None:
    with pytest.raises(ValueError, match="singular"):
        rigid_image_to_template_jacobian(
            (1, 1, 1),
            (1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (1.0, 1.0, 1.0),
            (1.0, 1.0, 1.0),
            (1.0, 1.0, 1.0),
        )


def test_registration_coefficient_contributions_match_source_formula_excerpt() -> None:
    expected = json.loads(
        (REFERENCE.parent / "java_coefficient_accumulation_reference.json").read_text()
    )
    actual = exact_registration_coefficient_contributions(
        [0.9, 0.25, 0.75, 0.6],
        [1.0, 0.5, 0.0, 0.8],
        [[0, 0, 0], [0.1, -0.2, 0.3], [-0.4, 0.5, 0.2], [0.15, 0.05, -0.1]],
        [[1, 2, 3, 4, 0, 0], [0.5, -1, 0.25, 0, 5, 0], [-0.75, 0.1, 0.3, 0, 0, 6]],
        registered_classes=[True, True, False, True],
    )
    assert actual.cost == np.float32(expected["cost"])
    assert actual.norm == np.float32(expected["norm"])
    for field in ("gradient", "hessian", "norm_gradient", "norm_hessian"):
        np.testing.assert_array_equal(getattr(actual, field), np.asarray(expected[field], dtype=np.float32))


def test_registration_gradient_descent_accept_and_backtrack_match_source() -> None:
    expected = json.loads(
        (REFERENCE.parent / "java_line_search_reference.json").read_text()
    )
    accepted = exact_registration_gradient_descent_step(
        [0.0], [1.0], [1.0], 0.0,
        lambda trial: 1.0 if trial[0] > 0 else 0.0,
    )
    rejected = exact_registration_gradient_descent_step(
        [0.0], [1.0], [1.0], 2.0,
        lambda _trial: 1.0,
    )
    for name, actual in (("accepted", accepted), ("rejected", rejected)):
        record = expected[name]
        np.testing.assert_array_equal(actual.transform, np.asarray(record[0], dtype=np.float32))
        assert actual.energy == np.float32(record[1])
        assert actual.lambda_value == np.float32(record[2])
        assert actual.relative_change == np.float32(record[3])
        assert actual.accepted is bool(record[4])
        assert actual.rejected_trials == int(record[5])


def test_full_grid_registration_coefficients_match_source_loop_order() -> None:
    x, y, z, k = np.indices((3, 2, 2, 4), dtype=np.int32)
    memberships = (((7*x + 3*y + 5*z + 2*k) % 17) / np.float32(16.0)).astype(np.float32)
    priors = (((2*x + 5*y + 3*z + 4*k) % 11) / np.float32(10.0)).astype(np.float32)
    prior_gradients = np.empty(memberships.shape + (3,), dtype=np.float32)
    for axis in range(3):
        if axis == 0:
            value = x + 2*y + 3*z + k
        elif axis == 1:
            value = 2*x + y + z + 2*k
        else:
            value = 3*x + y + 2*z + k
        prior_gradients[..., axis] = ((value % 7 - 3) / np.float32(10.0)).astype(np.float32)
    xj, yj, zj, axisj, parameter = np.indices((3, 2, 2, 3, 6), dtype=np.int32)
    jacobians = (((3*xj + 5*yj + 7*zj + 2*axisj + 4*parameter) % 13 - 6) / np.float32(10.0)).astype(np.float32)
    expected = json.loads((REFERENCE.parent / "java_full_coefficients_reference.json").read_text())
    for name, offset, subsample in (("full_grid", 0, 1), ("offset_subsample", 1, 2)):
        actual = exact_membership_registration_coefficients(
            memberships, priors, prior_gradients, jacobians,
            registered_classes=[True, True, False, True],
            offset=offset, subsample=subsample,
        )
        reference = expected[name]
        np.testing.assert_array_equal(
            np.asarray([actual.cost, actual.norm], dtype=np.float32),
            np.asarray(reference["cost_norm"], dtype=np.float32),
        )
        for field in ("gradient", "hessian", "norm_gradient", "norm_hessian"):
            np.testing.assert_array_equal(getattr(actual, field), np.asarray(reference[field], dtype=np.float32))


@pytest.mark.parametrize(
    ("memberships", "priors", "kwargs", "message"),
    [
        (np.ones((3, 3, 3)), np.ones((3, 3, 3)), {}, "4D"),
        (np.ones((3, 3, 3, 2)), np.ones((3, 3, 3, 3)), {}, "share a 4D"),
        (np.ones((3, 3, 3, 2)), np.ones((3, 3, 3, 2)), {"offset": 3}, "inside"),
        (np.ones((3, 3, 3, 2)), np.ones((3, 3, 3, 2)), {"subsample": 0}, "positive"),
        (np.ones((3, 3, 3, 2)), np.ones((3, 3, 3, 2)), {"registered_classes": [True]}, "one value per class"),
    ],
)
def test_exact_membership_energy_rejects_invalid_inputs(
    memberships: np.ndarray,
    priors: np.ndarray,
    kwargs: dict[str, object],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        exact_membership_registration_energy(memberships, priors, **kwargs)
