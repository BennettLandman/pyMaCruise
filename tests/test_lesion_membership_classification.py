from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_membership_classification import (
    lesion_toads_membership_classification,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_membership_classification"
REFERENCE_PATH = ROOT / "java_reference.json"
HARNESS_PATH = Path(__file__).parent / "java_reference" / "LesionMembershipClassificationReferenceHarness.java"
REFERENCE_SHA256 = "daefc8ce8888cb2a77656761d4d473c74ce38dfaf731b717ba3f50942ace8cf3"
HARNESS_SHA256 = "96d8f809654ab628b4bec3dd468b688f03c3ffe4d2254a1ad07822e765bbdaf9"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
REFERENCE = json.loads(REFERENCE_PATH.read_text())


def _image(name: str) -> np.ndarray:
    image = load_nifti_ras(ROOT / f"{name}.nii.gz")
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)[:3]) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj)


def test_membership_classification_matches_pinned_java_excerpt() -> None:
    assert hashlib.sha256(HARNESS_PATH.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(REFERENCE_PATH.read_bytes()).hexdigest() == REFERENCE_SHA256
    assert REFERENCE["source_sha256"] == SOURCE_SHA256
    assert REFERENCE["source_method"] == (
        "LesionToadSegmentation.updateLesionClassifcationformMemberships"
    )

    classification = _image("classification")
    memberships = _image("memberships")
    lesion = _image("lesion_membership")
    shape = tuple(REFERENCE["shape"])
    assert shape == (3, 2, 2)
    assert tuple(classification.shape) == shape
    assert tuple(memberships.shape) == (*shape, len(REFERENCE["template_labels"]))
    np.testing.assert_array_equal(classification.ravel(), REFERENCE["classification"])
    np.testing.assert_array_equal(
        memberships.ravel(), np.asarray(REFERENCE["memberships"], dtype=np.float32)
    )
    np.testing.assert_array_equal(
        lesion.ravel(), np.asarray(REFERENCE["lesion_membership"], dtype=np.float32)
    )

    for case in REFERENCE["cases"]:
        result = lesion_toads_membership_classification(
            classification,
            REFERENCE["template_labels"],
            case["atlas_names"],
            memberships,
            lesion,
            REFERENCE["lesion_class"],
            has_lesion=case["has_lesion"],
        )
        assert result.gm_class_index == case["gm_class_index"]
        assert result.classification_map.dtype == np.int16
        np.testing.assert_array_equal(result.classification_map.ravel(), case["map"])
        assert result.lesion_volume == case["lesion_volume"]

    duplicate_case = REFERENCE["cases"][0]
    assert duplicate_case["gm_class_index"] == 3  # Last case-insensitive match.
    assert duplicate_case["map"][0] == 0  # Strict equality at the threshold is excluded.
    assert REFERENCE["cases"][1]["gm_class_index"] == 0  # Source default when absent.
    assert REFERENCE["cases"][2]["lesion_volume"] == 0


@pytest.mark.parametrize(
    "changes, error, message",
    [
        ({"classification": np.zeros((2, 2), dtype=np.int8)}, ValueError, "3-D"),
        ({"atlas_names": ["OnlyOne"]}, ValueError, "one string per membership class"),
        ({"lesion_class": 9}, ValueError, "outside the membership class axis"),
        ({"has_lesion": 1}, ValueError, "has_lesion must be boolean"),
    ],
)
def test_invalid_membership_classification_arguments(changes, error, message) -> None:
    args = {
        "classification": _image("classification"),
        "template_labels": REFERENCE["template_labels"],
        "atlas_names": REFERENCE["cases"][0]["atlas_names"],
        "memberships": _image("memberships"),
        "lesion_membership": _image("lesion_membership"),
        "lesion_class": REFERENCE["lesion_class"],
        "has_lesion": True,
    }
    args.update(changes)
    with pytest.raises(error, match=message):
        lesion_toads_membership_classification(**args)
