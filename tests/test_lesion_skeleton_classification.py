from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_skeleton_classification import update_lesion_classification_from_skeleton


FIXTURE_DIR = Path(__file__).parent / "data" / "synthetic" / "lesion_skeleton_classification"
HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsSkeletonClassificationReferenceHarness.java"
REFERENCE = FIXTURE_DIR / "java_reference.json"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
HARNESS_SHA256 = "8f98f9c1306fc09d83c382b07899fb1e96f6e93b8dbe115c1b6cc2e4c2fa7689"
REFERENCE_SHA256 = "2a6568450b2fe269d2c139940ed7b5038eef139e77520cd1a83f7d1a4bb0988a"


def _load(name: str) -> np.ndarray:
    image = load_nifti_ras(FIXTURE_DIR / f"{name}.nii.gz")
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert np.array_equal(nib.affines.voxel_sizes(image.affine), (1.0, 1.0, 1.0))
    return np.asarray(image.dataobj)


def test_skeleton_classification_matches_selected_java_method() -> None:
    ref = json.loads(REFERENCE.read_text())
    assert ref["source_method"] == "LesionToadSegmentation.updateClassificationFromSkeleton"
    assert ref["source_sha256"] == SOURCE_SHA256
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(REFERENCE.read_bytes()).hexdigest() == REFERENCE_SHA256

    classification = _load("initial_classification").copy()
    skeleton = _load("skeleton")
    np.testing.assert_array_equal(classification.ravel(), ref["initial"])
    np.testing.assert_array_equal(skeleton.ravel(), ref["skeleton"])
    changed = update_lesion_classification_from_skeleton(
        classification, ref["template_labels"], skeleton
    )
    assert changed == ref["changed"] == 4
    assert classification.dtype == np.int8
    np.testing.assert_array_equal(classification.ravel(), ref["classification"])
    assert ref["classification"] == [0, -10, 20, -128, 99, 0, -128, -10, 99, 20]


def test_skeleton_classification_rejects_invalid_class_index_without_mutating() -> None:
    classification = np.asarray([9], dtype=np.int8).reshape((1, 1, 1))
    before = classification.copy()
    skeleton = np.asarray([0], dtype=np.int8).reshape((1, 1, 1))
    with np.testing.assert_raises_regex(ValueError, "index template_labels"):
        update_lesion_classification_from_skeleton(classification, [1], np.asarray([2], dtype=np.int8).reshape((1, 1, 1)))
    np.testing.assert_array_equal(classification, before)
