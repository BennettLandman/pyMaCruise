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
HARNESS_SHA256 = "49739272c35e61dfa65909c5e7e48e169ee9f4546f2879da0533bdead9336496"


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
    assert REFERENCE["lesion_algorithm_source_sha256"] == "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
    assert REFERENCE["lesion_preprocess_source_sha256"] == "72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08"
    assert REFERENCE["image_sha256"] == hashlib.sha256((ROOT / "images.nii.gz").read_bytes()).hexdigest()
    assert REFERENCE["selected_wrapper"] == "AlgorithmMultipleToads"
    assert REFERENCE["selected_overload"] == "initialMemberships(memberships,best,approx,centroids,classes)"
    assert REFERENCE["wrapper_call_order"][0:2] == [
        "AlgorithmMultipleToads.runAlgorithm",
        "preprocess.initialMemberships(memberships,best,approx,centroids,classes)",
    ]
    assert "separate full-class method, no approx parameter" in REFERENCE["dedicated_lesion_wrapper_call"]
    assert REFERENCE["has_atlas_shape_penalty"] is False
    assert REFERENCE["has_forced_candidate_slot"] is False
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
    # The fixture contains the duplicated first-channel centroid used by
    # VandyCRUISE to group repeated tissue priors. Confirm boundary voxels are
    # included in the full selected-method bit comparison as well.
    assert _bits(result.memberships[0, 0, 0]) == [
        REFERENCE["initial_memberships_bits"][i]
        for i in range(3)
    ]
    tail = (9 * 8 * 7 - 1) * 3
    assert _bits(result.memberships[-1, -1, -1]) == REFERENCE["initial_memberships_bits"][tail:tail + 3]
    np.testing.assert_array_equal(image, original)


def test_zero_cost_ties_use_lowest_class_indices_without_shape_seed() -> None:
    image = np.zeros((3, 3, 3, 1), dtype=np.float32)
    centroids = np.zeros((1, 11), dtype=np.float32)
    result = initialize_toads_default_approximate_initial_memberships(image, centroids, approx=3)
    reference = REFERENCE["zero_cost_tie"]
    np.testing.assert_array_equal(result.best_class_indices[1, 1, 1], [0, 1, 2])
    assert result.best_class_indices[1, 1, 1].tolist() == reference["indices"]
    assert _bits(result.memberships[1, 1, 1]) == reference["membership_bits"]


@pytest.mark.parametrize(
    "case, image_value, centroid_value",
    [("zero_cost_tie", 0.0, 0.0), ("below_cutoff_tie", 1.0e-11, 0.0)],
)
def test_cutoff_and_exact_ties_match_java_reference(case, image_value, centroid_value) -> None:
    image = np.full((1, 1, 1, 1), np.float32(image_value), dtype=np.float32)
    centroids = np.full((1, 11), np.float32(centroid_value), dtype=np.float32)
    result = initialize_toads_default_approximate_initial_memberships(image, centroids, approx=3)
    reference = REFERENCE[case]
    assert result.best_class_indices[0, 0, 0].tolist() == reference["indices"]
    assert _bits(result.memberships[0, 0, 0]) == reference["membership_bits"]


def test_duplicate_centroid_multiplicity_matches_java_reference() -> None:
    centroids = np.zeros((1, 11), dtype=np.float32)
    centroids[0, 3] = np.float32(1.0e-10)
    centroids[0, 4:] = np.float32(2.0e-10)
    result = initialize_toads_default_approximate_initial_memberships(
        np.zeros((1, 1, 1, 1), dtype=np.float32), centroids, approx=3,
    )
    reference = REFERENCE["duplicate_centroid_multiplicity"]
    assert result.best_class_indices[0, 0, 0].tolist() == reference["indices"]
    assert _bits(result.memberships[0, 0, 0]) == reference["membership_bits"]
    assert len(set(reference["indices"])) == 3


@pytest.mark.parametrize("approx", [0, 1, 11, 12])
def test_no_atlas_initialization_rejects_invalid_topk_size(approx: int) -> None:
    image = _image()
    with pytest.raises(ValueError, match="approx"):
        initialize_toads_default_approximate_initial_memberships(
            image, _centroids(), approx=approx,
        )
