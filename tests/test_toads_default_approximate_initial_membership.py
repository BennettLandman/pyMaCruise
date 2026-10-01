"""Java-pinned no-atlas approximate initialization used by the default pyramid."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_default_approximate_initial_membership import (
    initialize_toads_default_approximate_initial_memberships,
)
from pymacruise.toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_default_approximate_registration"
REFERENCE = json.loads((ROOT / "initial_membership_reference.json").read_text())
DESCRIPTOR = Path(__file__).parents[1] / "data" / "atlases" / "cruise-atlas-10obj-toads2008.txt"
HARNESS = Path(__file__).parent / "java_reference" / "ToadsDefaultApproximateInitialMembershipReferenceHarness.java"
HARNESS_SHA256 = "96f4bc66b121ff773eefe2e02bd85e2a443f4c89572556e37191adab55f67dcb"


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def _image() -> np.ndarray:
    image = nib.load(str(ROOT / "images.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert image.shape == (9, 8, 7, 2)
    np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), (.9, 1.2, 1.45), rtol=0, atol=1e-6)
    return np.asarray(image.dataobj, dtype=np.float32)


def _centroids() -> np.ndarray:
    result = np.asarray(
        [[np.float32(.19) + np.float32(.021) * k + np.float32(.013) * c for k in range(11)] for c in range(2)],
        dtype=np.float32,
    )
    # The selected source groups initial classes by exact first-channel
    # centroid equality; this models its repeated CSF intensity prior.
    result[0, 2] = result[0, 1]
    return result


def test_no_atlas_initialization_matches_selected_java_method() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert REFERENCE["algorithm_source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert REFERENCE["preprocess_source_sha256"] == "31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7"
    assert hashlib.sha256(DESCRIPTOR.read_bytes()).hexdigest() == REFERENCE["descriptor_sha256"]
    assert REFERENCE["class_types"] == list(DEFAULT_ATLAS_CLASS_TYPES)
    assert REFERENCE["approx"] == 3
    image = _image()
    original = image.copy()
    result = initialize_toads_default_approximate_initial_memberships(
        image, _centroids(), approx=REFERENCE["approx"],
    )
    assert result.class_types == DEFAULT_ATLAS_CLASS_TYPES
    assert result.memberships.shape == (9, 8, 7, 3)
    assert _bits(result.memberships) == REFERENCE["initial_memberships_bits"]
    np.testing.assert_array_equal(
        result.best_class_indices.ravel().astype(int), REFERENCE["best_class_indices"],
    )
    assert np.any(result.best_class_indices == 1)
    assert np.any(result.best_class_indices == 2)
    np.testing.assert_array_equal(image, original)


def test_zero_cost_ties_use_lowest_class_indices_without_shape_seed() -> None:
    image = np.zeros((3, 3, 3, 1), dtype=np.float32)
    centroids = np.zeros((1, 11), dtype=np.float32)
    result = initialize_toads_default_approximate_initial_memberships(image, centroids, approx=3)
    np.testing.assert_array_equal(result.best_class_indices[1, 1, 1], [0, 1, 2])
    # Source float32 summation and normalization can produce the next float
    # above one (three identical terms sum to 2.9999998 before division).
    np.testing.assert_allclose(
        result.memberships[1, 1, 1], np.ones(3, dtype=np.float32), rtol=0, atol=2e-7,
    )


@pytest.mark.parametrize("approx", [0, 1, 11, 12])
def test_no_atlas_initialization_rejects_invalid_topk_size(approx: int) -> None:
    image = _image()
    with pytest.raises(ValueError, match="approx"):
        initialize_toads_default_approximate_initial_memberships(
            image, _centroids(), approx=approx,
        )
