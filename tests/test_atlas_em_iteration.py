from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.atlas_em_iteration import run_atlas_em_fixed_prior_iterations


ROOT = Path(__file__).parent
HARNESS = ROOT / "java_reference" / "AtlasEMIterationReferenceHarness.java"


def _run_java(tmp_path: Path) -> dict:
    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        pytest.skip("JDK is required for pinned Java method references")
    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run([javac, "-d", str(classes), str(HARNESS)], check=True, capture_output=True, text=True)
    output = subprocess.run([java, "-cp", str(classes), "AtlasEMIterationReferenceHarness"], check=True, capture_output=True, text=True)
    return json.loads(output.stdout)


def _bits(values) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def _input():
    image = np.empty((5, 5, 5, 1), np.float32)
    mem = np.empty((5, 5, 5, 3), np.float32)
    prior = np.zeros_like(mem)
    mask = np.zeros((5, 5, 5), bool)
    for x in range(5):
        for y in range(5):
            for z in range(5):
                mask[x, y, z] = x > 0 and x < 4 and y > 0 and y < 4 and z > 0 and z < 4
                image[x, y, z, 0] = np.float32(x * 11 + y * 7 + z * 3 + 1) / np.float32(80)
                for k in range(3):
                    mem[x, y, z, k] = np.float32(1 + x + y * 2 + z * 3 + k) / np.float32(23)
                    if mask[x, y, z]:
                        prior[x, y, z, k] = np.float32(1 + x * 2 + y + z * 3 + k) / np.float32(20)
    means = np.asarray([[0.2, 0.5, 0.8]], np.float32)
    covariance = np.asarray([[[0.2]], [[0.3]], [[0.4]]], np.float32)
    return image, mem, mask, means, covariance, prior


def test_fixed_prior_em_iteration_loop_matches_java_excerpt(tmp_path: Path) -> None:
    reference = _run_java(tmp_path)
    result = run_atlas_em_fixed_prior_iterations(*_input(), max_iterations=4, max_membership_change=0.01, smoothing=0.125)
    assert reference["source_sha256"] == "b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a"
    assert result.iterations == reference["iterations"]
    assert result.converged is reference["converged"]
    expected_energies = np.asarray([int(value) & 0xFFFFFFFF for value in reference["energies"]], dtype=np.uint32)
    actual_energies = _bits(result.energy_history)
    # Numba and the JVM can differ by a few float32 ULPs in the diagnostic
    # energy because their double precision log implementations round at
    # different points. The posterior memberships and source stop metric stay
    # bit-identical, so keep the energy check tight without requiring a JVM's
    # exact libm rounding from the compiled kernel.
    assert all(abs(int(actual) - int(expected)) <= 4 for actual, expected in zip(actual_energies, expected_energies))
    assert _bits(result.maximum_change_history) == [int(value) & 0xFFFFFFFF for value in reference["changes"]]
    assert _bits(result.centroids.ravel()) == [int(value) & 0xFFFFFFFF for row in reference["centroids"] for value in row]
    assert _bits(result.covariances.ravel()) == [int(value) & 0xFFFFFFFF for row in reference["covariances"] for value in row]
    assert _bits(result.memberships.ravel()) == [int(value) & 0xFFFFFFFF for value in reference["memberships"]]


def test_zero_iteration_cap_preserves_initialized_state() -> None:
    inputs = _input()
    result = run_atlas_em_fixed_prior_iterations(*inputs, max_iterations=0)
    assert result.iterations == 0
    assert not result.converged
    assert result.energy_history.size == result.maximum_change_history.size == 0
    np.testing.assert_array_equal(result.memberships, inputs[1])
    np.testing.assert_array_equal(result.centroids, inputs[3])
    np.testing.assert_array_equal(result.covariances, inputs[4])


def test_empty_class_covariance_is_floored_and_recorded() -> None:
    image, memberships, mask, means, covariance, priors = _input()
    memberships[..., 2] = 0.0
    result = run_atlas_em_fixed_prior_iterations(
        image,
        memberships,
        mask,
        means,
        covariance,
        priors,
        max_iterations=1,
        smoothing=0.125,
        covariance_floor=1.0e-6,
    )
    assert result.regularized_classes == (2,)
    assert result.covariances[2, 0, 0] == np.float32(1.0e-6)
    assert np.isfinite(result.memberships).all()
