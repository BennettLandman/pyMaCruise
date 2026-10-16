from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.atlas_em_statistics import update_atlas_em_gaussian_memberships


ROOT = Path(__file__).parent
HARNESS = ROOT / "java_reference" / "AtlasEMGaussianMembershipReferenceHarness.java"


def _run_java(tmp_path: Path, *arguments: str) -> dict:
    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        pytest.skip("JDK is required for pinned Java method references")
    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run([javac, "-d", str(classes), str(HARNESS)], check=True, capture_output=True, text=True)
    output = subprocess.run([java, "-cp", str(classes), "AtlasEMGaussianMembershipReferenceHarness", *arguments], check=True, capture_output=True, text=True)
    return json.loads(output.stdout)


def _bits(values) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def _input():
    nx = ny = nz = 4
    classes = 3
    image = np.empty((nx, ny, nz, 1), np.float32)
    memberships = np.empty((nx, ny, nz, classes), np.float32)
    prior = np.zeros_like(memberships)
    mask = np.zeros((nx, ny, nz), bool)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                mask[x, y, z] = x > 0 and x < nx - 1 and y > 0 and y < ny - 1 and z > 0 and z < nz - 1
                image[x, y, z, 0] = np.float32(x * 7 + y * 5 + z * 3 + 1) / np.float32(20)
                for k in range(classes):
                    memberships[x, y, z, k] = np.float32(1 + x + y * 2 + z * 3 + k) / np.float32(13)
                    if mask[x, y, z]:
                        prior[x, y, z, k] = np.float32(1 + x * 2 + y + z * 3 + k) / np.float32(20)
    means = np.asarray([[0.15, 0.47, 0.83]], np.float32)
    covariance = np.asarray([[[0.12]], [[0.2]], [[0.31]]], np.float32)
    return image, memberships, mask, means, covariance, prior


def test_one_channel_atlas_mode_update_matches_java_excerpt(tmp_path: Path) -> None:
    reference = _run_java(tmp_path)
    result, metrics = update_atlas_em_gaussian_memberships(*_input()[:3], *_input()[3:], smoothing=0.625)
    assert reference["source_sha256"] == "b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a"
    assert _bits(result.ravel()) == [int(value) & 0xFFFFFFFF for value in reference["memberships"]]
    assert _bits(metrics) == [int(value) & 0xFFFFFFFF for value in reference["metrics"]]


def test_gaussian_membership_update_requires_source_padding_and_positive_covariance() -> None:
    image, memberships, mask, means, covariance, prior = _input()
    mask[0, 1, 1] = True
    with pytest.raises(ValueError, match="boundary faces"):
        update_atlas_em_gaussian_memberships(image, memberships, mask, means, covariance, prior)
    _, memberships, mask, means, covariance, prior = _input()
    covariance[0, 0, 0] = 0.0
    with pytest.raises(ValueError, match="covariance must be positive"):
        update_atlas_em_gaussian_memberships(image, memberships, mask, means, covariance, prior)


def test_source_zero_membership_energy_nan_is_preserved(tmp_path: Path) -> None:
    reference = _run_java(tmp_path, "zero")
    image, memberships, mask, means, covariance, prior = _input()
    memberships[1, 1, 1, 0] = 0.0
    prior[1, 1, 1, 0] = 0.0
    _, metrics = update_atlas_em_gaussian_memberships(image, memberships, mask, means, covariance, prior, smoothing=0.625)
    java_energy_bits = int(reference["metrics"][0]) & 0xFFFFFFFF
    assert np.isnan(metrics[0])
    assert java_energy_bits in {0x7FC00000, 0xFFC00000}
