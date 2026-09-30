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
HARNESS_SHA256 = "c79e53379c7ab8c4aa5b9360551c86e94c8e21c9ec15606e69543b7ecd4dd912"
REFERENCE_SHA256 = "f071a7e45f7be8522b759e2ef4ac2241b9af8c04536eef06c5d1541d45a3c212"


def test_hard_classification_matches_selected_java_method() -> None:
    ref = json.loads(REFERENCE.read_text())
    assert ref["source_method"] == "LesionToadSegmentation.generateClassification"
    assert ref["source_sha256"] == SOURCE_SHA256
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(REFERENCE.read_bytes()).hexdigest() == REFERENCE_SHA256

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
    assert result[:, 0, 0].tolist() == [99, 10, 0, 0, 99, 10, 10, 99, 20, 0]


def test_hard_classification_rejects_shape_and_label_errors() -> None:
    memberships = np.ones((2, 1, 1, 2), dtype=np.float32)
    lesion = np.ones((2, 1, 1), dtype=np.float32)
    with np.testing.assert_raises_regex(ValueError, "spatial shape"):
        generate_lesion_toads_classification([1, 2, 3], memberships, lesion[:1], 2)
    with np.testing.assert_raises_regex(ValueError, "outside template_labels"):
        generate_lesion_toads_classification([1, 2], memberships, lesion, 2)
    with np.testing.assert_raises_regex(ValueError, "signed Java byte"):
        generate_lesion_toads_classification([1, 256], memberships, lesion, 1)
