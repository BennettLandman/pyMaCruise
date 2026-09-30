from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_classification_update import (
    update_lesion_classification_from_memberships,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_classification_membership_update"
HARNESS_PATH = Path(__file__).parent / "java_reference" / "LesionClassificationFromMembershipsReferenceHarness.java"
REFERENCE_PATH = ROOT / "java_reference.json"
HARNESS_SHA256 = "48112895f385a6d19326fe2361c648cf104bb6b8fc2204650a263f04c8bd0335"
REFERENCE_SHA256 = "95ecea654ef1a320014f18207249597b2b04ca55fd85c9276758ee781730deca"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
WRAPPER_SHA256 = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
REFERENCE = json.loads(REFERENCE_PATH.read_text())


def _image(case_name: str, field: str) -> np.ndarray:
    image = load_nifti_ras(ROOT / f"{case_name}_{field}.nii.gz")
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)[:3]) == (1.0, 1.0, 1.0)
    return np.asarray(image.dataobj)


def test_matches_pinned_java_update_for_multiple_and_single_class_arrays() -> None:
    assert hashlib.sha256(HARNESS_PATH.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(REFERENCE_PATH.read_bytes()).hexdigest() == REFERENCE_SHA256
    assert REFERENCE["source_sha256"] == SOURCE_SHA256
    assert REFERENCE["wrapper_sha256"] == WRAPPER_SHA256
    assert REFERENCE["source_method"] == (
        "LesionToadSegmentation.updateClassificationFromMemberships"
    )

    for case in REFERENCE["scenarios"]:
        previous = _image(case["name"], "previous")
        memberships = _image(case["name"], "memberships")
        lesion = _image(case["name"], "lesion_membership")
        shape = tuple(case["shape"])
        assert previous.shape == shape
        assert memberships.shape == (*shape, len(case["template_labels"]))
        assert lesion.shape == shape
        assert previous.dtype == np.int8
        assert memberships.dtype == np.float32
        before = previous.copy()

        changed = update_lesion_classification_from_memberships(
            previous,
            np.asarray(case["template_labels"], dtype=np.int16),
            memberships,
            lesion,
            case["lesion_class"],
        )

        assert changed == case["class_diff"]
        assert int(np.count_nonzero(before != previous)) == changed
        np.testing.assert_array_equal(previous.ravel(), case["classification"])


def test_source_edge_cases_are_explicit_in_the_pinned_inputs() -> None:
    multi = REFERENCE["scenarios"][0]
    assert multi["classification"] == [37, 0, 99, -7, 5, 99, -7, -7]
    assert multi["class_diff"] == 6
    # An ordinary tie picks the first class; a lesion tie does not displace it.
    assert multi["classification"][0] == 37
    assert multi["classification"][3] == -7
    # A positive winning class-zero membership maps through template_labels[0].
    assert multi["classification"][4] == 5

    single = REFERENCE["scenarios"][1]
    assert single["classification"] == [-13, -13, -13]
    assert single["class_diff"] == 2
