from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_hard_classification import generate_lesion_toads_classification


FIXTURE_DIR = Path(__file__).parent / "data" / "synthetic" / "lesion_hard_classification"
HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsHardClassificationReferenceHarness.java"
REFERENCE = FIXTURE_DIR / "java_reference.json"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
HARNESS_SHA256 = "885d221e7414fc76fd1000c9dca36a1725643b4ba5aa18c686567611c0d0c10a"
REFERENCE_SHA256 = "4a9deec7593852552c1a1f6c5baa9b0120c4e07f26008b2aaad0fc75c5c37064"
INPUT_SHA256 = "bd351d39da63b3431eff206d3b8c4dc9fe9ae3ab5aa267a78d435907a0e0c0bc"


def test_hard_classification_matches_selected_java_method() -> None:
    ref = json.loads(REFERENCE.read_text())
    assert ref["source_method"] == "LesionToadSegmentation.generateClassification"
    assert ref["source_sha256"] == SOURCE_SHA256
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(REFERENCE.read_bytes()).hexdigest() == REFERENCE_SHA256
    assert hashlib.sha256((FIXTURE_DIR / "memberships_and_lesion.nii.gz").read_bytes()).hexdigest() == INPUT_SHA256

    image = load_nifti_ras(FIXTURE_DIR / "memberships_and_lesion.nii.gz")
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert np.array_equal(nib.affines.voxel_sizes(image.affine), (1.0, 1.0, 1.0))
    values = np.asarray(image.dataobj, dtype=np.float32)
    assert values.shape == (10, 1, 1, 4)
    memberships, lesion = values[..., :3], values[..., 3]
    result = generate_lesion_toads_classification(
        ref["template_labels"], memberships, lesion, ref["lesion_class"]
    )
    assert result.dtype == np.int8
    np.testing.assert_array_equal(result.ravel(), ref["labels"])
    # Strict lesion comparison: ordinary wins an equality, and all-zero
    # memberships preserve the method's zero output even for nonzero labels.
    assert result[:, 0, 0].tolist() == [99, -106, 0, 0, 99, -106, -106, 99, 20, 0]


def test_hard_classification_rejects_shape_and_label_errors() -> None:
    memberships = np.ones((2, 1, 1, 2), dtype=np.float32)
    lesion = np.ones((2, 1, 1), dtype=np.float32)
    with np.testing.assert_raises_regex(ValueError, "spatial shape"):
        generate_lesion_toads_classification([1, 2, 3], memberships, lesion[:1], 2)
    with np.testing.assert_raises_regex(ValueError, "outside template_labels"):
        generate_lesion_toads_classification([1, 2], memberships, lesion, 2)
    with np.testing.assert_raises_regex(ValueError, "signed Java byte"):
        generate_lesion_toads_classification([1, 256], memberships, lesion, 1)
