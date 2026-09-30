from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_ventricle import lesion_toads_ventricle_distance


FIXTURE_DIR = Path(__file__).parent / "data" / "synthetic" / "lesion_ventricle_distance"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsVentricleDistanceReferenceHarness.java"
HARNESS_SHA256 = "aa4f8a7ccddd3644bae6b9a675efc08da30aa20269d36f927d3517b34d5da45c"
REFERENCE = FIXTURE_DIR / "java_reference.json"
REFERENCE_SHA256 = "d886c02988386465684ac014b1c444d3a609c129a0cc950972004f794104dd62"


def test_ventricle_mask_and_byte_distance_match_java_method_excerpt() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(REFERENCE.read_bytes()).hexdigest() == REFERENCE_SHA256
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_method"] == "LesionToadSegmentation.computeVentirclsDistance"
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["names"] == [
        "Background",
        "Ventricles",
        "LeftVentricle",
        "RightVentricle",
        "ThirdVentricle",
        "ventricles",
        "Left-Ventricle",
        "Thirdventricle",
        "Other",
    ]

    source_images = {
        "all_exact_names_and_controls_cap4": "seeded_labels.nii.gz",
        "all_exact_names_cap1": "seeded_labels.nii.gz",
        "nonmatching_names_empty_mask_cap3": "nonmatching_labels.nii.gz",
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

        mask, distance = lesion_toads_ventricle_distance(
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

    ordinary, cap_one, empty = reference["cases"]
    assert sum(ordinary["mask_xyz_flat"]) == 4
    assert set(empty["mask_xyz_flat"]) == {0}
    assert set(empty["distance_xyz_flat"]) == {-1}
    # fastDistance performs one scan before testing its byte cap, so distance
    # two is assigned even when max_distance is 1.
    assert 2 in cap_one["distance_xyz_flat"]


def test_ventricle_distance_rejects_invalid_atlas_indices_and_labels() -> None:
    with np.testing.assert_raises_regex(ValueError, "outside atlas_names"):
        lesion_toads_ventricle_distance(np.asarray([[[1]]]), ["Background"], 3)
    with np.testing.assert_raises_regex(TypeError, "integer atlas indices"):
        lesion_toads_ventricle_distance(np.asarray([[[0.0]]]), ["Background"], 3)
