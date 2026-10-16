from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.atlas_em_statistics import (
    compute_atlas_em_centroids,
    compute_atlas_em_covariances,
    compute_atlas_em_mixing_coefficients,
)

ROOT = Path(__file__).parent
HARNESS = ROOT / "java_reference" / "AtlasEMStatisticsReferenceHarness.java"


def _input():
    nx, ny, nz, nc, classes = 3, 2, 2, 2, 3
    image = np.empty((nx, ny, nz, nc), np.float32)
    field = np.empty_like(image)
    memberships = np.empty((nx, ny, nz, classes), np.float32)
    mask = np.empty((nx, ny, nz), bool)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                mask[x, y, z] = (x + 2 * y + z) % 3 != 0
                for c in range(nc):
                    image[x, y, z, c] = np.float32(np.float32(0.25 + (x * 17 + y * 7 + z * 3 + c * 11) % 29) / np.float32(7.0))
                    field[x, y, z, c] = np.float32(np.float32(0.8) + np.float32(x + y + z + c) / np.float32(10.0))
                for k in range(classes):
                    memberships[x, y, z, k] = np.float32(1 + x * 3 + y * 2 + z + k * 5) / np.float32(17.0)
    memberships[..., 2] = 0.0
    return image, memberships, mask, field


def _run_java(tmp_path: Path) -> dict:
    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        pytest.skip("JDK is required for pinned Java method references")
    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run([javac, "-d", str(classes), str(HARNESS)], check=True, capture_output=True, text=True)
    output = subprocess.run([java, "-cp", str(classes), "AtlasEMStatisticsReferenceHarness"], check=True, capture_output=True, text=True)
    return json.loads(output.stdout)


def _float_bits(values) -> list:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def _unsigned_tree(values):
    if isinstance(values, list):
        return [_unsigned_tree(item) for item in values]
    return int(values) & 0xFFFFFFFF


def test_fixed_input_atlas_em_statistics_match_java_source_excerpt(tmp_path: Path) -> None:
    ref = _run_java(tmp_path)
    image, mem, mask, field = _input()
    initial = np.zeros((2, 3), dtype=np.float32)
    initial[:, 2] = (0.7, 0.8)

    priors = compute_atlas_em_mixing_coefficients(mem, mask)
    equal = compute_atlas_em_mixing_coefficients(mem, mask, mode="Equal")
    means = compute_atlas_em_centroids(image, mem, mask, initial=initial)
    covariance = compute_atlas_em_covariances(image, mem, mask, means)

    assert ref["source_sha256"] == "b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a"
    assert _float_bits(priors) == ref["prior"]
    assert _float_bits(equal) == ref["equal_prior"]
    assert _float_bits(means) == ref["centroid"]
    assert _float_bits(covariance) == _unsigned_tree(ref["covariance"])
    assert np.all(covariance[2] == 0.0)  # empty class matches Java zero covariance


def test_atlas_em_statistics_validate_shapes_and_finite_inputs() -> None:
    with pytest.raises(ValueError, match="finite"):
        compute_atlas_em_mixing_coefficients(np.full((2, 2, 2, 1), np.nan), np.ones((2, 2, 2), bool))
