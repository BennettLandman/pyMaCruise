from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_registration import (
    exact_membership_registration_energy,
    exact_registration_measure_coefficients,
    synthetic_registration_objective_inputs,
)


ROOT = Path(__file__).parent
REFERENCE = ROOT / "data" / "synthetic" / "toads_registration_energy" / "java_reference.json"
COEFFICIENT_REFERENCE = (
    ROOT / "data" / "synthetic" / "toads_registration_energy" / "java_coefficient_reference.json"
)
HARNESS_HASH = "3c2e8bc403a1e04f201cbc29f6aa468175b202c4d2013bc1aefa4429d7bf8b41"
COEFFICIENT_HARNESS_HASH = "cc62e5e25bd3297a0d03a981c951cea0f7bb137810fe0791b2952114c39af1a5"


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
