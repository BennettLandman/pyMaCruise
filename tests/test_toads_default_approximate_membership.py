"""Selected Java comparison for default-atlas TOADS top-k initialization."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_default_approximate_membership import (
    initialize_toads_default_approximate_memberships,
)
from pymacruise.toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_default_approximate_membership"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
DESCRIPTOR = Path(__file__).parents[1] / "data" / "atlases" / "cruise-atlas-10obj-toads2008.txt"


def _bits(value: np.ndarray) -> list[int]:
    return np.asarray(value, dtype=np.float32).ravel().view(np.int32).tolist()


def _fixture():
    image = nib.load(str(ROOT / "images.nii.gz"))
    shapes = nib.load(str(ROOT / "transformed_shapes.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(shapes.affine) == ("R", "A", "S")
    assert image.shape == (7, 6, 5, 2)
    assert shapes.shape == (7, 6, 5, 11)
    np.testing.assert_array_equal(image.affine, shapes.affine)
    np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), (.9, 1.2, 1.45), rtol=0, atol=1e-6)
    return np.asarray(image.dataobj, dtype=np.float32), np.asarray(shapes.dataobj, dtype=np.float32)


def _centroids() -> np.ndarray:
    return np.asarray(
        [[np.float32(.19) + np.float32(.021) * k + np.float32(.013) * c for k in range(11)] for c in range(2)],
        dtype=np.float32,
    )


def test_default_approximate_memberships_match_selected_java_reference() -> None:
    assert REFERENCE["source_sha256"] == "31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7"
    assert hashlib.sha256(DESCRIPTOR.read_bytes()).hexdigest() == REFERENCE["descriptor_sha256"]
    assert REFERENCE["descriptor_class_types"] == list(DEFAULT_ATLAS_CLASS_TYPES)
    assert REFERENCE["approx"] == 3
    assert DESCRIPTOR.is_file()
    evidence = REFERENCE["p18czo_evidence_sha256"]
    for key, path in {
        "harness": Path(__file__).parent / "java_reference/ToadsDefaultApproximateMembershipReferenceHarness.java",
        "numerics": Path(__file__).parent / "java_reference/upstream/edu/jhmi/rad/medic/utilities/Numerics.java",
        "images": ROOT / "images.nii.gz",
        "transformed_shapes": ROOT / "transformed_shapes.nii.gz",
        "descriptor": DESCRIPTOR,
    }.items():
        assert evidence[key] == hashlib.sha256(path.read_bytes()).hexdigest(), key
    images, shapes = _fixture()
    originals = (images.copy(), shapes.copy())
    result = initialize_toads_default_approximate_memberships(
        images, _centroids(), shapes, approx=REFERENCE["approx"],
    )
    assert result.class_types == DEFAULT_ATLAS_CLASS_TYPES
    assert result.memberships.shape == (7, 6, 5, 3)
    assert result.best_class_indices.shape == result.memberships.shape
    assert _bits(result.memberships) == REFERENCE["initial_memberships_bits"]
    np.testing.assert_array_equal(
        result.best_class_indices.ravel().astype(int), REFERENCE["best_class_indices"]
    )
    assert np.all(np.sum(result.memberships, axis=-1) > 0.0)
    assert np.all(result.best_class_indices < len(DEFAULT_ATLAS_CLASS_TYPES))
    assert np.any(result.best_class_indices == 10)
    assert np.all(np.diff(np.sort(result.best_class_indices, axis=-1), axis=-1) > 0)
    assert np.any(result.memberships[0, :, :, :] != 0.0)
    np.testing.assert_array_equal(images, originals[0])
    np.testing.assert_array_equal(shapes, originals[1])


def test_default_approximate_memberships_preserve_unprocessed_voxels() -> None:
    images = np.zeros((3, 3, 3, 1), dtype=np.float32)
    shapes = np.zeros((3, 3, 3, 11), dtype=np.float32)
    centers = np.zeros((1, 11), dtype=np.float32)
    centers[0, 5] = 1.0
    shapes[1, 1, 1, 5] = 1.0
    mask = np.zeros(images.shape[:3], dtype=bool)
    mask[1, 1, 1] = True
    result = initialize_toads_default_approximate_memberships(
        images, centers, shapes, approx=3, process_mask=mask,
    )
    assert _bits(result.memberships) == REFERENCE["process_mask_case"]["memberships_bits"]
    assert result.best_class_indices.ravel().astype(int).tolist() == REFERENCE["process_mask_case"]["best_class_indices"]
    assert np.count_nonzero(result.memberships[~mask]) == 0
    assert np.count_nonzero(result.best_class_indices[~mask]) == 0
    assert np.count_nonzero(result.memberships[mask]) > 0


def test_default_approximate_memberships_keep_shape_seed_and_low_index_ties() -> None:
    image = np.zeros((1, 1, 1, 1), dtype=np.float32)
    centers = np.zeros((1, 11), dtype=np.float32)
    centers[0, 5] = 1.0
    shapes = np.zeros((1, 1, 1, 11), dtype=np.float32)
    shapes[..., 5] = 1.0
    result = initialize_toads_default_approximate_memberships(
        image, centers, shapes, approx=3,
    )
    np.testing.assert_array_equal(
        result.best_class_indices[0, 0, 0], REFERENCE["shape_seeded_membership_tie"]["indices"]
    )
    assert _bits(result.memberships[0, 0, 0]) == REFERENCE["shape_seeded_membership_tie"]["membership_bits"]


def test_default_approximate_memberships_zero_cost_fallback_matches_java() -> None:
    result = initialize_toads_default_approximate_memberships(
        np.zeros((1, 1, 1, 1), dtype=np.float32),
        np.zeros((1, 11), dtype=np.float32),
        np.zeros((1, 1, 1, 11), dtype=np.float32),
        approx=3,
    )
    np.testing.assert_array_equal(
        result.best_class_indices[0, 0, 0], REFERENCE["zero_cost_tie"]["indices"]
    )
    assert _bits(result.memberships[0, 0, 0]) == REFERENCE["zero_cost_tie"]["membership_bits"]


@pytest.mark.parametrize("approx", [0, 1, 11, 12])
def test_default_approximate_memberships_reject_invalid_topk_size(approx: int) -> None:
    images, shapes = _fixture()
    with pytest.raises(ValueError, match="approx"):
        initialize_toads_default_approximate_memberships(
            images, _centroids(), shapes, approx=approx,
        )
