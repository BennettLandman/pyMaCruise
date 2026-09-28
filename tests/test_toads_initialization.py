from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_initialization import (
    initial_exact_memberships,
    initial_intensity_centroids,
)


FIXTURES = Path(__file__).parent / "data" / "synthetic" / "toads_initialization"


def _fixture_images(centroids: np.ndarray) -> np.ndarray:
    nx, ny, nz = 6, 6, 6
    images = np.empty((2, nx, ny, nz), dtype=np.float32)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                code0 = (x * 11 + y * 7 + z * 3) % 23
                code1 = (x * 5 + y * 13 + z * 17) % 29
                images[0, x, y, z] = np.float32(
                    np.float32(0.2) + np.float32(code0) / np.float32(20.0)
                )
                images[1, x, y, z] = np.float32(
                    np.float32(0.2) + np.float32(code1) / np.float32(25.0)
                )
    images[:, 2, 2, 2] = centroids[:, 1]
    return images


def test_centroids_and_exact_memberships_match_java_harness() -> None:
    spec = json.loads((FIXTURES / "spec.json").read_text())
    reference = json.loads((FIXTURES / "java_reference.json").read_text())
    assert reference["shape"] == spec["shape_xyz"]
    assert reference["crop_border"] == spec["crop_border_voxels"]

    centroids = initial_intensity_centroids(
        spec["low"], spec["high"], spec["intensity_priors"]
    )
    expected_centroids = np.asarray(reference["centroids"], dtype=np.float32)
    np.testing.assert_array_equal(centroids, expected_centroids)

    images = _fixture_images(centroids)
    memberships = initial_exact_memberships(
        images, centroids, crop_border=spec["crop_border_voxels"]
    )
    expected_memberships = np.asarray(reference["memberships"], dtype=np.float32).reshape(
        (*spec["shape_xyz"], spec["classes"])
    )
    np.testing.assert_array_equal(memberships, expected_memberships)

    # The Java implementation leaves crop padding and the final source plane at zero.
    assert np.count_nonzero(memberships[0]) == 0
    assert np.count_nonzero(memberships[4:, :, :, :]) == 0
    assert np.count_nonzero(memberships[:, 4:, :, :]) == 0
    assert np.count_nonzero(memberships[:, :, 4:, :]) == 0
    # Two classes share channel-0 centroid; source multiplies their reciprocal
    # scores by their exact-equality group size.
    assert centroids[0, 1] == centroids[0, 2]
    assert memberships[2, 2, 2, 1] > 0.999
    # With source grouping weights, total membership can exceed one when
    # channel-0 centroids tie; do not impose a conventional FCM sum-to-one rule.
    assert memberships[2, 2, 2].sum() > 1.0


@pytest.mark.parametrize(
    "low,high,priors",
    [
        ([1.0], [1.0], [[0.2]]),
        ([2.0, 3.0], [4.0], [[0.2]]),
        ([1.0], [2.0], [[np.nan]]),
    ],
)
def test_centroid_inputs_are_validated(low, high, priors) -> None:
    with pytest.raises(ValueError):
        initial_intensity_centroids(low, high, priors)


def test_membership_inputs_and_crop_dimensions_are_validated() -> None:
    images = np.ones((1, 3, 3, 3), dtype=np.float32)
    centroids = np.ones((1, 2), dtype=np.float32)
    with pytest.raises(ValueError, match="positive robust intensity span"):
        initial_intensity_centroids([2], [2], [[0.2]])
    with pytest.raises(ValueError, match="too small"):
        initial_exact_memberships(images, centroids, crop_border=1)
    with pytest.raises(ValueError, match="finite"):
        initial_exact_memberships(np.full_like(images, np.inf), centroids, crop_border=0)
