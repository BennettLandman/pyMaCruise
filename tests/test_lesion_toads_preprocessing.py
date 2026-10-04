from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest
from nibabel.orientations import aff2axcodes

from pymacruise.lesion_toads_preprocessing import (
    initialize_lesion_toads_class_centroids,
    initialize_lesion_toads_centroids,
    initialize_lesion_toads_intensities,
    normalize_lesion_toads_images,
)


FIXTURES = Path(__file__).parent / "data" / "synthetic" / "lesion_toads_preprocessing"


def _load_channels() -> np.ndarray:
    t1 = np.asanyarray(nib.load(FIXTURES / "t1_mprage.nii.gz").dataobj)
    flair = np.asanyarray(nib.load(FIXTURES / "flair.nii.gz").dataobj)
    return np.stack((t1, flair)).astype(np.float32)


def test_lesion_toads_initialization_matches_java_reference() -> None:
    recipe = json.loads((FIXTURES / "spec.json").read_text())
    reference = json.loads((FIXTURES / "java_reference.json").read_text())
    assert aff2axcodes(nib.load(FIXTURES / "t1_mprage.nii.gz").affine) == (
        "R",
        "A",
        "S",
    )
    assert aff2axcodes(nib.load(FIXTURES / "flair.nii.gz").affine) == (
        "R",
        "A",
        "S",
    )
    harness = (
        Path(__file__).parent
        / "java_reference"
        / "LesionToadsPreprocessingReferenceHarness.java"
    )
    assert hashlib.sha256(harness.read_bytes()).hexdigest() == (
        "f0c7d925ff96cc0414c04714fe0335ce03ff023f3a388b1e618a349d1262ae53"
    )
    result = initialize_lesion_toads_intensities(
        _load_channels(),
        recipe["lesion_priors"],
        recipe["black_hole_priors"],
        robust_fraction=recipe["robust_fraction"],
        robust_scales=recipe["robust_scales"],
    )

    assert reference["robustMaxMin"] is False
    assert reference["call_order"] == [
        "LesionToadPreprocess.constructor",
        "AlgorithmLesionToads.normalizeImages",
        "findCroppingBoundaries",
        "cropImages",
        "initialCentroids",
        "initialLesionCentroid",
        "initialBlackHoleCentroid",
    ]
    assert reference["source_sha256"] == {
        "LesionToadPreprocess.java": "72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08",
        "AlgorithmLesionToads.java": "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb",
        "ImageFunctions.java": "be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892",
    }
    assert list(result.normalized_images.shape[1:]) == reference["shape"]
    np.testing.assert_array_equal(
        result.intensity_minima,
        np.asarray(reference["intensity_minima"], dtype=np.float32),
    )
    np.testing.assert_array_equal(
        result.intensity_maxima,
        np.asarray(reference["intensity_maxima"], dtype=np.float32),
    )
    np.testing.assert_array_equal(
        result.centroid_robust_minima,
        np.asarray(reference["centroid_robust_minima"], dtype=np.float32),
    )
    np.testing.assert_array_equal(
        result.centroid_robust_maxima,
        np.asarray(reference["centroid_robust_maxima"], dtype=np.float32),
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
        assert np.all(actual >= 0.0) and np.all(actual <= 1.0)


def test_lesion_toads_clamps_negative_values_before_bounds() -> None:
    images = _load_channels()
    assert np.any(images < 0)
    result = normalize_lesion_toads_images(images)
    assert np.all(result.normalized_images[images < 0] == np.float32(0.0))


@pytest.mark.parametrize("value", [0.0, 7.25])
def test_lesion_toads_zero_ordinary_span_is_explicit_deviation(value: float) -> None:
    image = np.full((1, 3, 4, 5), value, dtype=np.float32)
    with pytest.raises(ValueError, match="ordinary intensity span.*divide by zero"):
        initialize_lesion_toads_intensities(image, [0.4], [0.1])


def test_lesion_toads_rejects_channel_and_prior_mismatch() -> None:
    image = np.empty((2, 3, 4, 5), dtype=np.float32)
    image[0] = np.arange(60, dtype=np.float32).reshape(3, 4, 5)
    image[1] = image[0] * np.float32(2.0)
    with pytest.raises(ValueError, match="lesion priors"):
        initialize_lesion_toads_intensities(image, [0.4], [0.1, 0.6])


def test_centroid_bounds_use_normalized_cropped_inputs_separately() -> None:
    image = _load_channels()
    normalized = normalize_lesion_toads_images(image)
    cropped_normalized = normalized.normalized_images[:, 1:-1, 1:-1, 1:-1]
    centroids = initialize_lesion_toads_centroids(
        cropped_normalized, [0.4, 1.2], [0.1, 0.6]
    )
    assert centroids.robust_minima.shape == (2,)
    assert centroids.robust_maxima.shape == (2,)
    assert np.all(centroids.robust_maxima >= centroids.robust_minima)


def test_selected_wrapper_uses_robust_bounds_only_for_cropped_centroids() -> None:
    reference = json.loads((FIXTURES / "java_reference.json").read_text())
    normalized = normalize_lesion_toads_images(_load_channels())
    # This fixed asymmetric interior region stands in for the wrapper's
    # normalized, already-cropped image argument. Crop selection is not under
    # test; the Java harness uses these same voxel limits for its centroid
    # method reference.
    cropped = normalized.normalized_images[:, 1:-1, 2:-2, 1:-1]
    result = initialize_lesion_toads_centroids(cropped, [0.4, 1.2], [0.1, 0.6])
    full_grid = initialize_lesion_toads_centroids(
        normalized.normalized_images, [0.4, 1.2], [0.1, 0.6]
    )

    np.testing.assert_array_equal(
        result.robust_minima,
        np.asarray(reference["cropped_centroid_robust_minima"], dtype=np.float32),
    )
    np.testing.assert_array_equal(
        result.robust_maxima,
        np.asarray(reference["cropped_centroid_robust_maxima"], dtype=np.float32),
    )
    np.testing.assert_array_equal(
        result.lesion_centroid,
        np.asarray(reference["cropped_lesion_centroid"], dtype=np.float32),
    )
    np.testing.assert_array_equal(
        result.black_hole_centroid,
        np.asarray(reference["cropped_black_hole_centroid"], dtype=np.float32),
    )
    assert not np.array_equal(result.robust_minima, full_grid.robust_minima)
    assert not np.array_equal(result.robust_maxima, full_grid.robust_maxima)


def test_class_centroids_match_selected_java_initial_centroids() -> None:
    recipe = json.loads((FIXTURES / "spec.json").read_text())
    reference = json.loads((FIXTURES / "java_reference.json").read_text())
    normalized = normalize_lesion_toads_images(_load_channels())
    cropped = normalized.normalized_images[:, 1:-1, 2:-2, 1:-1]
    result = initialize_lesion_toads_class_centroids(
        cropped,
        np.asarray(recipe["class_intensity_priors"], dtype=np.float32),
    )

    np.testing.assert_array_equal(
        result.robust_minima,
        np.asarray(reference["cropped_centroid_robust_minima"], dtype=np.float32),
    )
    np.testing.assert_array_equal(
        result.robust_maxima,
        np.asarray(reference["cropped_centroid_robust_maxima"], dtype=np.float32),
    )
    np.testing.assert_array_equal(
        result.class_centroids.view(np.int32),
        np.asarray(reference["initial_class_centroids"], dtype=np.float32).view(np.int32),
    )
