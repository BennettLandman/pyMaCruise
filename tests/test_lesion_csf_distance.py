from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_csf_distance import lesion_toads_csf_distance


FIXTURE_DIR = Path(__file__).parent / "data" / "synthetic" / "lesion_csf_distance"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
WRAPPER_SHA256 = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsCsfDistanceReferenceHarness.java"
HARNESS_SHA256 = "94384a693850c706c10744e8dcf4816c28fe4f7c37f81420117e9ded874f759b"
REFERENCE = FIXTURE_DIR / "java_reference.json"
REFERENCE_SHA256 = "160836366079cf9bfc57a953fa582439161ecf833bd765cf24d6df6599274269"


def test_csf_mask_and_distance_match_selected_java_excerpt() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(REFERENCE.read_bytes()).hexdigest() == REFERENCE_SHA256
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_method"] == "LesionToadSegmentation.computeCSFDistance"
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["wrapper_sha256"] == WRAPPER_SHA256
    assert reference["default_max_distance"] == 5
    assert reference["names"] == [
        "Background", "Caudate", "Thalamus", "Putamen", "Subcortical-GM",
        "SubcorticalGM", "BasalGanglia", "Caudate", "caudate", "Thalamus ",
        "Subcortical_GM", "Subcortical-Gm", "Basal-Ganglia", "Other",
    ]

    image_for_case = {
        "all_exact_names_duplicate_and_controls_cap1": "seeded_atlas_indices.nii.gz",
        "all_exact_names_duplicate_and_controls_default_cap5": "seeded_atlas_indices.nii.gz",
        "near_match_names_only_empty_mask_default_cap5": "near_miss_atlas_indices.nii.gz",
    }
    for case in reference["cases"]:
        image = load_nifti_ras(FIXTURE_DIR / image_for_case[case["name"]])
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert tuple(image.shape) == tuple(case["shape"])
        segmentation = np.asarray(image.dataobj, dtype=np.int64)
        expected_segmentation = np.asarray(case["segmentation_xyz_flat"], dtype=np.int64).reshape(
            case["shape"], order="C"
        )
        np.testing.assert_array_equal(segmentation, expected_segmentation)
        before = segmentation.copy()
        cap = case["max_distance"]
        kwargs = {} if case["name"].endswith("default_cap5") else {"max_distance": cap}
        mask, distance = lesion_toads_csf_distance(segmentation, reference["names"], **kwargs)
        np.testing.assert_array_equal(segmentation, before)
        np.testing.assert_array_equal(
            mask.astype(np.uint8),
            np.asarray(case["mask_xyz_flat"], dtype=np.uint8).reshape(case["shape"], order="C"),
        )
        np.testing.assert_array_equal(
            distance,
            np.asarray(case["distance_xyz_flat"], dtype=np.int8).reshape(case["shape"], order="C"),
        )

    seeded = reference["cases"][0]
    assert sum(seeded["mask_xyz_flat"]) == 7  # Six unique names plus duplicate Caudate.
    assert 2 in seeded["distance_xyz_flat"]  # selected fastDistance cap-one quirk.
    empty = reference["cases"][2]
    assert not any(empty["mask_xyz_flat"])
    assert set(empty["distance_xyz_flat"]) == {-1}


def test_csf_distance_rejects_invalid_labels_and_cap() -> None:
    with np.testing.assert_raises_regex(ValueError, "outside atlas_names"):
        lesion_toads_csf_distance(np.asarray([[[1]]]), ["Background"])
    with np.testing.assert_raises_regex(TypeError, "integer atlas indices"):
        lesion_toads_csf_distance(np.asarray([[[0.0]]]), ["Background"])
    with np.testing.assert_raises_regex(ValueError, "signed-byte range"):
        lesion_toads_csf_distance(np.asarray([[[0]]]), ["Background"], 128)
