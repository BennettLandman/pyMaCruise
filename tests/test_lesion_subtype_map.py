from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_subtype_map import update_lesion_subtype_map


FIXTURE_DIR = Path(__file__).parent / "data" / "synthetic" / "lesion_subtype_map"
HARNESS = Path(__file__).parent / "java_reference" / "LesionSubtypeMapReferenceHarness.java"
REFERENCE = FIXTURE_DIR / "java_reference.json"
HARNESS_SHA256 = "96dadf93f5e5daf0bb1e4d9c89f4efa8e2a313f1a7095fb702fd26ae4b4078fa"
REFERENCE_SHA256 = "0e335e99eb44f485798650c95f50f06b712aca6cf1a584c70a330ab6ab665f24"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
WRAPPER_SHA256 = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"


def _load(name: str) -> np.ndarray:
    image = load_nifti_ras(FIXTURE_DIR / f"{name}.nii.gz")
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert np.array_equal(nib.affines.voxel_sizes(image.affine), (1.0, 1.0, 1.0))
    return np.asarray(image.dataobj)


def test_subtype_map_matches_pinned_java_excerpt_and_nifti_fixtures() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(REFERENCE.read_bytes()).hexdigest() == REFERENCE_SHA256
    ref = json.loads(REFERENCE.read_text())
    assert ref["source_method"] == "LesionToadSegmentation.updateLesionClassifcation"
    assert ref["source_sha256"] == SOURCE_SHA256
    assert ref["wrapper_sha256"] == WRAPPER_SHA256

    classification = _load("classification")
    memberships = _load("memberships")
    lesion_membership = _load("lesion_membership")
    black_hole_membership = _load("black_hole_membership")
    initial_map = _load("initial_map")
    shape = tuple(ref["shape"])
    assert shape == (7, 1, 1)
    assert classification.shape == lesion_membership.shape == shape
    assert black_hole_membership.shape == initial_map.shape == shape
    assert memberships.shape == (*shape, len(ref["template_labels"]))
    np.testing.assert_array_equal(classification.ravel(), ref["classification"])
    np.testing.assert_array_equal(
        memberships.ravel(), np.asarray(ref["memberships"], dtype=np.float32)
    )
    np.testing.assert_array_equal(
        lesion_membership.ravel(), np.asarray(ref["lesion_membership"], dtype=np.float32)
    )
    np.testing.assert_array_equal(
        black_hole_membership.ravel(), np.asarray(ref["black_hole_membership"], dtype=np.float32)
    )
    np.testing.assert_array_equal(initial_map.ravel(), ref["initial_map"])

    for case in ref["cases"]:
        result = initial_map.copy()
        count = update_lesion_subtype_map(
            classification,
            ref["template_labels"],
            memberships,
            lesion_membership,
            black_hole_membership,
            ref["lesion_class"],
            case["lesions_enabled"],
            result,
        )
        expected = _load(f"expected_{case['name']}")
        assert count == case["qualifying_count"]
        np.testing.assert_array_equal(result, expected)
        np.testing.assert_array_equal(result.ravel(), case["map"])

    enabled = ref["cases"][0]
    assert enabled["map"][0] == 1  # Lesion membership is the strict winner.
    assert enabled["map"][1] == 2  # Black-hole membership is the strict winner.
    assert enabled["map"][2] == 1  # Lesion wins an exact subtype tie.
    assert enabled["map"][3] == 0  # Strict threshold equality does not qualify.
    assert enabled["map"][5] == 0  # Wrong input classification label is cleared.
    assert enabled["map"][4] == 0  # Both candidates below ordinary membership.
    assert enabled["map"][6] == 2  # Black-hole can qualify when lesion ties threshold.
    assert enabled["qualifying_count"] == 4
    assert ref["cases"][1]["map"] == [0] * shape[0]  # Disabled lesions clear stale map.


def test_subtype_map_rejects_invalid_inputs() -> None:
    classification = _load("classification")
    memberships = _load("memberships")
    lesion_membership = _load("lesion_membership")
    black_hole_membership = _load("black_hole_membership")
    out = _load("initial_map")
    args = (
        classification,
        [0, 9],
        memberships,
        lesion_membership,
        black_hole_membership,
        1,
        True,
        out,
    )
    with np.testing.assert_raises_regex(ValueError, "outside the membership class axis"):
        update_lesion_subtype_map(*(*args[:5], 2, *args[6:]))
    with np.testing.assert_raises_regex(TypeError, "has_lesion must be boolean"):
        update_lesion_subtype_map(*(*args[:6], 1, *args[7:]))
