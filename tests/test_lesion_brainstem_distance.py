from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_brainstem import lesion_toads_brainstem_distance


FIXTURE_DIR = Path(__file__).parent / "data" / "synthetic" / "lesion_brainstem_distance"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsBrainstemDistanceReferenceHarness.java"
HARNESS_SHA256 = "0de0717502f1601244a0d498122bf3b0d15a0df2368d1c76952835300411d97f"
REFERENCE = FIXTURE_DIR / "java_reference.json"
REFERENCE_SHA256 = "b93533e9d4029fd6f15e477b8cf75375927582003a1741bfe1ac237b21d89ab4"


def test_brainstem_mask_and_capped_distance_match_java_method_excerpt() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(REFERENCE.read_bytes()).hexdigest() == REFERENCE_SHA256
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_method"] == "LesionToadSegmentation.computeBstemDistance"
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["names"] == [
        "Background",
        "Brainstem",
        "InterVentricular-WM",
        "FornixWM",
        "brainstem",
        "InterVentricularWM",
        "Fornix-WM",
        "Other",
    ]

    source_images = {
        "asymmetric_seed_labels_cap4": "seeded_labels.nii.gz",
        "asymmetric_seed_labels_cap1": "seeded_labels.nii.gz",
        "alias_only_empty_mask_cap3": "aliases_only_labels.nii.gz",
    }
    for case in reference["cases"]:
        image = load_nifti_ras(FIXTURE_DIR / source_images[case["name"]])
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert tuple(image.shape) == tuple(case["shape"])
        segmentation = np.asarray(image.dataobj, dtype=np.int64)
        expected_segmentation = np.asarray(case["segmentation_xyz_flat"], dtype=np.int64).reshape(
            case["shape"], order="C"
        )
        np.testing.assert_array_equal(segmentation, expected_segmentation)

        mask, distance = lesion_toads_brainstem_distance(
            segmentation, reference["names"], case["max_distance"]
        )
        expected_mask = np.asarray(case["mask_xyz_flat"], dtype=np.uint8).reshape(
            case["shape"], order="C"
        )
        expected_distance = np.asarray(case["distance_xyz_flat"], dtype=np.int8).reshape(
            case["shape"], order="C"
        )
        np.testing.assert_array_equal(mask.astype(np.uint8), expected_mask)
        np.testing.assert_array_equal(distance, expected_distance)

    # Exact spellings seed the mask; case and punctuation aliases do not.
    seed_case = reference["cases"][0]
    assert sum(seed_case["mask_xyz_flat"]) == 3
    empty_case = reference["cases"][2]
    assert not any(empty_case["mask_xyz_flat"])
    assert set(empty_case["distance_xyz_flat"]) == {-1}


def test_brainstem_distance_rejects_invalid_atlas_indices_and_labels() -> None:
    with np.testing.assert_raises_regex(ValueError, "outside atlas_names"):
        lesion_toads_brainstem_distance(np.asarray([[[1]]]), ["Background"], 3)
    with np.testing.assert_raises_regex(TypeError, "integer atlas indices"):
        lesion_toads_brainstem_distance(np.asarray([[[0.0]]]), ["Background"], 3)
