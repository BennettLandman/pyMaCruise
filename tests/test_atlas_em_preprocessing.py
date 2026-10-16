from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.atlas_em_preprocessing import (
    initialize_atlas_em_centroids,
    initialize_atlas_em_memberships,
    normalize_atlas_em_images,
)


ROOT = Path(__file__).parent
HARNESS = ROOT / "java_reference" / "AtlasEMPreprocessReferenceHarness.java"


def _run_java(tmp_path: Path) -> dict:
    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        pytest.skip("JDK is required for pinned Java method references")
    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run([javac, "-d", str(classes), str(HARNESS)], check=True, capture_output=True, text=True)
    output = subprocess.run([java, "-cp", str(classes), "AtlasEMPreprocessReferenceHarness"], check=True, capture_output=True, text=True)
    return json.loads(output.stdout)


def _bits(array) -> list[int]:
    return np.asarray(array, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def test_atlas_em_initial_centroids_and_memberships_match_java_excerpt(tmp_path: Path) -> None:
    reference = _run_java(tmp_path)
    low = np.asarray([0.25, 0.5], np.float32)
    high = np.asarray([2.25, 4.5], np.float32)
    priors = np.asarray([[0.0, 0.3, 0.3, 0.8], [0.1, 0.2, 0.4, 0.7]], np.float32)
    centroids = initialize_atlas_em_centroids(low, high, priors)

    image = np.empty((4, 3, 3, 2), np.float32)
    for x in range(4):
        for y in range(3):
            for z in range(3):
                for c in range(2):
                    image[x, y, z, c] = np.float32(x * 13 + y * 5 + z * 3 + c * 7 + 1) / np.float32(17)
    memberships = initialize_atlas_em_memberships(image, centroids, interior=(slice(0, 3), slice(0, 2), slice(0, 2)))
    default_crop_memberships = initialize_atlas_em_memberships(image, centroids)

    assert reference["source_sha256"] == "acec23f916b3134e53381e83ff3de7ca569f785e6b866bfe135fb409c2926d2b"
    assert _bits(centroids.ravel()) == [int(v) & 0xFFFFFFFF for row in reference["centroid"] for v in row]
    assert _bits(memberships.ravel()) == [int(v) & 0xFFFFFFFF for v in reference["memberships"]]
    assert _bits(default_crop_memberships.ravel()) == _bits(memberships.ravel())
    assert np.all(memberships[3, ...] == 0.0)
    assert np.all(memberships[:, 2, ...] == 0.0)


def test_atlas_em_normalization_uses_positive_clip_and_robust_span() -> None:
    image = np.empty((13, 11, 9, 1), np.float32)
    for x in range(13):
        for y in range(11):
            for z in range(9):
                value = np.float32(0.5 + ((x + 13 * y + 143 * z) * 37 % 501) / 10.0)
                if x == 0 and y == 0 and z == 0:
                    value = -2.0
                image[x, y, z, 0] = value
    result = normalize_atlas_em_images(image)
    assert result.images.dtype == np.float32
    assert result.images[0, 0, 0, 0] == 0.0
    assert result.robust_minima[0] > 0.0
    assert result.robust_maxima[0] > result.robust_minima[0]
    np.testing.assert_array_equal(result.images[..., 0], np.maximum(image[..., 0], 0) / (result.robust_maxima[0] - result.robust_minima[0]))


def test_atlas_em_initialization_rejects_degenerate_robust_span() -> None:
    with pytest.raises(ValueError, match="zero or invalid robust span"):
        normalize_atlas_em_images(np.ones((3, 3, 3, 1), np.float32))
