from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.toads_math import toads_optimized_centroid_similarity


ROOT = Path(__file__).parent / "data" / "synthetic"
REFERENCE = ROOT / "toads_optimized_similarity" / "java_reference.json"
SEGMENTATION_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"


def test_optimized_centroid_similarity_matches_java_formula_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["segmentation_source_sha256"] == SEGMENTATION_SHA256
    assert reference["prior_scale"] == 0.05
    variance_reference = json.loads(
        (ROOT / "toads_optimized_variance" / "java_reference.json").read_text()
    )
    class_variance_reference = json.loads(
        (ROOT / "toads_optimized_class_variance" / "java_reference.json").read_text()
    )
    np.testing.assert_array_equal(
        reference["image_variance_weights"],
        variance_reference["without_field"]["normalized_inverse_variance"],
    )
    np.testing.assert_array_equal(
        np.asarray(reference["class_standard_deviations"], dtype=np.float32),
        np.asarray(class_variance_reference["without_field"], dtype=np.float32).reshape(3, 3),
    )
    similarity = toads_optimized_centroid_similarity(
        reference["centroids"],
        reference["class_standard_deviations"],
        reference["image_variance_weights"],
        reference["prior_scale"],
    )
    np.testing.assert_array_equal(
        similarity,
        np.asarray(reference["similarity"], dtype=np.float32).reshape(3, 3),
    )
    np.testing.assert_array_equal(np.diag(similarity), np.zeros(3, dtype=np.float32))
    np.testing.assert_array_equal(similarity, similarity.T)


def test_similarity_preserves_source_ieee_behavior_for_zero_deviations() -> None:
    centers = np.asarray([[0.0, 1.0]], dtype=np.float32)
    deviations = np.zeros((1, 2), dtype=np.float32)
    similarity = toads_optimized_centroid_similarity(
        centers, deviations, np.ones(1, dtype=np.float32)
    )
    assert similarity[0, 1] == np.float32(0.0)
    assert similarity[1, 0] == np.float32(0.0)

    identical_centers = np.zeros((1, 2), dtype=np.float32)
    nonfinite = toads_optimized_centroid_similarity(
        identical_centers, deviations, np.ones(1, dtype=np.float32)
    )
    assert np.isnan(nonfinite[0, 1])
    assert np.isnan(nonfinite[1, 0])
