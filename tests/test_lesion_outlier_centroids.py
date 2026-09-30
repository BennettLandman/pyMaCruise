from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_statistics import compute_lesion_toads_basic_centroids


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_outlier_centroids"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "LesionOutlierCentroidsReferenceHarness.java"
HARNESS_SHA256 = "527d6a4169261d4c924e7e788b3aebc8aee6faaa5d6e23cfea16743711183b38"
INITIAL = np.asarray(
    [[0.15, 0.55, 1.05, 1.65], [0.30, 0.70, 1.15, 1.80]], dtype=np.float32
)
PRIOR = np.asarray(
    [[0.20, 0.60, 0.95, 1.40], [0.35, 0.75, 1.10, 1.55]], dtype=np.float32
)
LINKS = np.asarray([[0, 1, 1, 2], [0, 0, 1, 1]], dtype=np.int16)
TYPES = ("object", "outlier", "outlier", "outmask")


def _bits(values: np.ndarray) -> list[list[int]]:
    return values.view(np.uint32).astype(np.uint64).tolist()


def _inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    image = nib.load(str(ROOT / "outlier_centroid_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (5, 4, 3, 12)
    return (
        frames[..., 0].astype(np.int16),
        frames[..., 1:5],
        frames[..., 5:7],
        frames[..., 7:9],
        frames[..., 11],
    )


def test_outlier_centroid_modes_match_selected_java_bits() -> None:
    assert REFERENCE["source_sha256"] == "77016c000b8bbbfe55bfbe8d665ddd8be914c8ecbae8563ccce59aaa834bcfe0"
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    labels, memberships, images, fields, outlier = _inputs()
    assert set(np.unique(outlier)) == {0.0, 0.25, 0.5, 0.75, 1.0}
    assert set(np.unique(labels)) == {0, 1, 2}  # class 3 remains absent
    modes = ("none", "linked", "damped", "prior", "fallback")
    for mode in modes:
        for field_key, field_values in (("no_field", None), ("with_field", fields)):
            kwargs: dict[str, object] = {
                "fields": field_values,
                "outlier_membership": outlier,
                "mode": mode,
                "smoothness": 0.35,
            }
            if mode == "linked":
                kwargs["centroid_links"] = LINKS
            if mode == "prior":
                kwargs["prior_centroids"] = PRIOR
            result = compute_lesion_toads_basic_centroids(
                images, memberships, labels, INITIAL, TYPES, **kwargs
            )
            assert _bits(result) == REFERENCE[f"{mode}_{field_key}"]


def test_outlier_suppression_is_type_specific_and_validated() -> None:
    labels, memberships, images, fields, outlier = _inputs()
    without_map = compute_lesion_toads_basic_centroids(
        images, memberships, labels, INITIAL, ("object", "object", "object", "object")
    )
    outmask_control = compute_lesion_toads_basic_centroids(
        images, memberships, labels, INITIAL, ("object", "outmask", "mask", "outmask")
    )
    np.testing.assert_array_equal(outmask_control, without_map)
    with pytest.raises(ValueError, match="required"):
        compute_lesion_toads_basic_centroids(
            images, memberships, labels, INITIAL, TYPES
        )
    with pytest.raises(ValueError, match="match the image grid"):
        compute_lesion_toads_basic_centroids(
            images, memberships, labels, INITIAL, TYPES,
            outlier_membership=outlier[..., :2],
        )
    invalid = outlier.copy()
    invalid[0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        compute_lesion_toads_basic_centroids(
            images, memberships, labels, INITIAL, TYPES, outlier_membership=invalid
        )
    # Zero outlier membership preserves ordinary class accumulation exactly;
    # the all-one map removes every sample from OUTLIER classes.
    zeros = np.zeros_like(outlier)
    ones = np.ones_like(outlier)
    expected = compute_lesion_toads_basic_centroids(
        images, memberships, labels, INITIAL, TYPES, outlier_membership=zeros
    )
    ordinary = compute_lesion_toads_basic_centroids(
        images, memberships, labels, INITIAL, ("object", "object", "object", "outmask")
    )
    np.testing.assert_array_equal(expected, ordinary)
    removed = compute_lesion_toads_basic_centroids(
        images, memberships, labels, INITIAL, TYPES, outlier_membership=ones
    )
    np.testing.assert_array_equal(removed[:, 1:3], np.zeros((2, 2), dtype=np.float32))
