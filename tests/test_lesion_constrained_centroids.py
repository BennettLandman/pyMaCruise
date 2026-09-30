from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_statistics import (
    compute_constrained_lesion_toads_centroids,
    compute_weighted_lesion_toads_centroids,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_constrained_centroids"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
WEIGHTED_REFERENCE = json.loads((ROOT / "weighted_java_reference.json").read_text())
INITIAL = np.asarray([[0.15, 0.55, 1.05, 1.65], [0.30, 0.70, 1.15, 1.80]], dtype=np.float32)


def _bits(values: np.ndarray) -> list[list[int]]:
    return values.view(np.uint32).astype(np.uint64).tolist()


def test_constrained_centroid_update_matches_selected_java() -> None:
    assert REFERENCE["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    image = nib.load(str(ROOT / "centroid_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (5, 4, 3, 9)
    labels = frames[..., 0].astype(np.int16)
    memberships = frames[..., 1:5]
    values = frames[..., 5:7]
    fields = frames[..., 7:9]
    initial_copy = INITIAL.copy()

    no_field = compute_constrained_lesion_toads_centroids(
        values, memberships, labels, INITIAL
    )
    with_field = compute_constrained_lesion_toads_centroids(
        values, memberships, labels, INITIAL, fields=fields
    )
    assert _bits(no_field) == REFERENCE["without_field_bits"]
    assert _bits(with_field) == REFERENCE["with_field_bits"]
    np.testing.assert_array_equal(INITIAL, initial_copy)
    for channel in range(2):
        assert no_field[channel, 0] == np.float32(
            INITIAL[channel, 0]
            + np.float32(
                np.float32(0.45)
                * np.float32(INITIAL[channel, 1] - INITIAL[channel, 0])
            )
        )
        assert no_field[channel, 1] == np.float32(
            INITIAL[channel, 1]
            + np.float32(
                np.float32(0.45)
                * np.float32(INITIAL[channel, 2] - INITIAL[channel, 1])
            )
        )
    for channel in range(2):
        empty_class_lower = np.float32(
            INITIAL[channel, 3]
            - np.float32(
                np.float32(0.45)
                * np.float32(INITIAL[channel, 3] - INITIAL[channel, 2])
            )
        )
        assert no_field[channel, 3] == empty_class_lower


def test_constrained_centroid_update_rejects_invalid_grids() -> None:
    image = np.zeros((3, 3, 3, 1), dtype=np.float32)
    memberships = np.ones((3, 3, 3, 2), dtype=np.float32)
    labels = np.zeros((3, 3, 3), dtype=np.int16)
    centroids = np.asarray([[0.5, 1.0]], dtype=np.float32)
    with pytest.raises(ValueError, match="centroids must have shape"):
        compute_constrained_lesion_toads_centroids(image, memberships, labels, np.zeros((2, 2)))
    with pytest.raises(ValueError, match="labels must index"):
        compute_constrained_lesion_toads_centroids(image, memberships, labels + 2, centroids)
    with pytest.raises(ValueError, match="fields must match"):
        compute_constrained_lesion_toads_centroids(
            image, memberships, labels, centroids, fields=np.zeros_like(image[..., :0])
        )


def test_weighted_centroid_update_matches_selected_java() -> None:
    assert WEIGHTED_REFERENCE["source_sha256"] == REFERENCE["source_sha256"]
    image = nib.load(str(ROOT / "centroid_inputs.nii.gz"))
    frames = np.asarray(image.dataobj, dtype=np.float32)
    labels = frames[..., 0].astype(np.int16)
    memberships = frames[..., 1:5]
    values = frames[..., 5:7]
    fields = frames[..., 7:9]
    initial_copy = INITIAL.copy()
    for factor in (0.0, 0.35, 1.0):
        for field_values, suffix in ((None, "no_field"), (fields, "field")):
            actual = compute_weighted_lesion_toads_centroids(
                values,
                memberships,
                labels,
                INITIAL,
                fields=field_values,
                factor=factor,
            )
            assert _bits(actual) == WEIGHTED_REFERENCE["cases"][f"f{factor}_{suffix}"]
    np.testing.assert_array_equal(INITIAL, initial_copy)
    assert all(WEIGHTED_REFERENCE["cases"][f"f{factor}_{suffix}"][0][3] == 0
               for factor in (0.0, 0.35, 1.0)
               for suffix in ("no_field", "field"))
    initial_bits = _bits(INITIAL)
    for suffix in ("no_field", "field"):
        held = WEIGHTED_REFERENCE["cases"][f"f1.0_{suffix}"]
        assert held[0][:3] == initial_bits[0][:3]
        assert held[1][:3] == initial_bits[1][:3]
    with pytest.raises(ValueError, match="factor must be finite"):
        compute_weighted_lesion_toads_centroids(
            values, memberships, labels, INITIAL, factor=np.inf
        )
