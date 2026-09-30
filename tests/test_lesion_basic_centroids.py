from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_statistics import compute_lesion_toads_basic_centroids


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_variances"
REFERENCE = json.loads((ROOT / "basic_centroids_java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "LesionBasicCentroidsReferenceHarness.java"
HARNESS_SHA256 = "dfe263b4d5600253c89b7fb3d741a4b2fc15ae1761ecf969238c47152303e05f"
INITIAL = np.asarray(
    [[0.15, 0.55, 1.05, 1.65], [0.30, 0.70, 1.15, 1.80]], dtype=np.float32
)
PRIOR = np.asarray(
    [[0.20, 0.60, 0.95, 1.40], [0.35, 0.75, 1.10, 1.55]], dtype=np.float32
)
LINKS = np.asarray([[0, 1, 1, 2], [0, 0, 1, 1]], dtype=np.int16)
TYPES = ("object", "rayleigh", "mask", "object")


def _bits(values: np.ndarray) -> list[list[int]]:
    return values.view(np.uint32).astype(np.uint64).tolist()


def _inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    image = nib.load(str(ROOT / "variance_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    frames = np.asarray(image.dataobj, dtype=np.float32)
    return (
        frames[..., 0].astype(np.int16),
        frames[..., 1:5],
        frames[..., 5:7],
        frames[..., 7:9],
    )


def test_generic_centroid_modes_match_selected_java_bits() -> None:
    assert REFERENCE["source_sha256"] == "77016c000b8bbbfe55bfbe8d665ddd8be914c8ecbae8563ccce59aaa834bcfe0"
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    labels, memberships, images, fields = _inputs()
    assert labels.shape == (5, 4, 3)
    assert set(np.unique(labels)) == {0, 1, 2}  # class 3 is deliberately absent
    modes = ("none", "linked", "damped", "prior", "fallback")
    for mode in modes:
        for field_key, field_values in (("no_field", None), ("with_field", fields)):
            kwargs: dict[str, object] = {"fields": field_values, "mode": mode, "smoothness": 0.35}
            if mode == "linked":
                kwargs["centroid_links"] = LINKS
            if mode == "prior":
                kwargs["prior_centroids"] = PRIOR
            actual = compute_lesion_toads_basic_centroids(
                images, memberships, labels, INITIAL, TYPES, **kwargs
            )
            assert _bits(actual) == REFERENCE[f"{mode}_{field_key}"]
    # In damped mode the absent class preserves its prior centroid; other
    # generic modes follow their selected source zero/linked/prior behavior.
    damped = compute_lesion_toads_basic_centroids(
        images, memberships, labels, INITIAL, TYPES, mode="damped", smoothness=0.35
    )
    np.testing.assert_array_equal(damped[:, 3], INITIAL[:, 3])


def test_generic_centroids_reject_bad_grids_modes_and_special_types() -> None:
    labels, memberships, images, fields = _inputs()
    with pytest.raises(ValueError, match="fields"):
        compute_lesion_toads_basic_centroids(
            images, memberships, labels, INITIAL, TYPES, fields=fields[..., :1]
        )
    invalid_labels = labels.copy()
    invalid_labels[0, 0, 0] = 4
    with pytest.raises(ValueError, match="labels"):
        compute_lesion_toads_basic_centroids(images, memberships, invalid_labels, INITIAL, TYPES)
    with pytest.raises(ValueError, match="centroid_links"):
        compute_lesion_toads_basic_centroids(
            images, memberships, labels, INITIAL, TYPES, mode="linked", centroid_links=np.zeros((2, 3), int)
        )
    with pytest.raises(ValueError, match="prior_centroids"):
        compute_lesion_toads_basic_centroids(
            images, memberships, labels, INITIAL, TYPES, mode="prior"
        )
    for kind in (
        "lesion", "blackhole", "optimized_lesion", "optimized_blackhole",
    ):
        unsupported = ("object", kind, "mask", "object")
        with pytest.raises(ValueError, match="specialized"):
            compute_lesion_toads_basic_centroids(images, memberships, labels, INITIAL, unsupported)
