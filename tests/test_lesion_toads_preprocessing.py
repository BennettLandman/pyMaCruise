from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_toads_preprocessing import (
    initialize_lesion_toads_intensities,
)


FIXTURES = Path(__file__).parent / "data" / "synthetic" / "lesion_toads_preprocessing"


def _load_channels() -> np.ndarray:
    t1 = np.asanyarray(nib.load(FIXTURES / "t1_mprage.nii.gz").dataobj)
    flair = np.asanyarray(nib.load(FIXTURES / "flair.nii.gz").dataobj)
    return np.stack((t1, flair)).astype(np.float32)


def test_lesion_toads_initialization_matches_java_reference() -> None:
    recipe = json.loads((FIXTURES / "spec.json").read_text())
    reference = json.loads((FIXTURES / "java_reference.json").read_text())
    result = initialize_lesion_toads_intensities(
        _load_channels(),
        recipe["lesion_priors"],
        recipe["black_hole_priors"],
        robust_fraction=recipe["robust_fraction"],
        robust_scales=recipe["robust_scales"],
    )

    assert list(result.normalized_images.shape[1:]) == reference["shape"]
    np.testing.assert_array_equal(
        result.robust_minima, np.asarray(reference["robust_minima"], dtype=np.float32)
    )
    np.testing.assert_array_equal(
        result.robust_maxima, np.asarray(reference["robust_maxima"], dtype=np.float32)
    )
    np.testing.assert_array_equal(
        result.lesion_centroid,
        np.asarray(reference["lesion_centroid"], dtype=np.float32),
    )
    np.testing.assert_array_equal(
        result.black_hole_centroid,
        np.asarray(reference["black_hole_centroid"], dtype=np.float32),
    )
    for channel in range(result.normalized_images.shape[0]):
        actual = result.normalized_images[channel].ravel(order="F")
        expected = np.asarray(reference["channels"][channel], dtype=np.float32)
        np.testing.assert_array_equal(actual, expected)


def test_lesion_toads_clamps_negative_values_before_bounds() -> None:
    images = _load_channels()
    assert np.any(images < 0)
    result = initialize_lesion_toads_intensities(images, [0.4, 1.2], [0.1, 0.6])
    assert np.all(result.normalized_images[images < 0] == np.float32(0.0))


@pytest.mark.parametrize("value", [0.0, 7.25])
def test_lesion_toads_zero_robust_span_is_explicit_deviation(value: float) -> None:
    image = np.full((1, 3, 4, 5), value, dtype=np.float32)
    with pytest.raises(ValueError, match="divide by zero"):
        initialize_lesion_toads_intensities(image, [0.4], [0.1])


def test_lesion_toads_rejects_channel_and_prior_mismatch() -> None:
    image = np.empty((2, 3, 4, 5), dtype=np.float32)
    image[0] = np.arange(60, dtype=np.float32).reshape(3, 4, 5)
    image[1] = image[0] * np.float32(2.0)
    with pytest.raises(ValueError, match="lesion priors"):
        initialize_lesion_toads_intensities(image, [0.4], [0.1, 0.6])
