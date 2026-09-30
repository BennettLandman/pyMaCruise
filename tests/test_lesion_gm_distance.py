from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_gm_distance import lesion_toads_gm_distance


FIXTURE_DIR = Path(__file__).parent / "data" / "synthetic" / "lesion_gm_distance"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsGmDistanceReferenceHarness.java"
HARNESS_SHA256 = "be18a8292350505187e34f669d4963b47f106359a0f513dd8202626e7e3768e3"
REFERENCE = FIXTURE_DIR / "java_reference.json"
REFERENCE_SHA256 = "08e3aeabaec1588821f3e7b63e0c50117cdc39924f581804c1258bc57e6c559d"


def test_gm_mask_and_byte_distance_match_selected_java_excerpt() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(REFERENCE.read_bytes()).hexdigest() == REFERENCE_SHA256
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_method"] == "LesionToadSegmentation.computeGMDistance"
    assert reference["source_sha256"] == SOURCE_SHA256
    exact_names = [
        "Cerebrum-GM", "CerebralGM", "CorticalGM", "Caudate", "GM",
        "Thalamus", "Putamen", "Sulcal-CSF", "SulcalCSF", "Subcortical-GM",
        "SubcorticalGM", "BasalGanglia", "GlobusPallidus",
    ]
    assert reference["names"][1:14] == exact_names
    assert reference["names"][14:] == [
        "cerebrum-gm", "CerebrumGM", "Cerebrum_GM", "Cerebrum-Gm", "Other"
    ]

    source_images = {
        "all_exact_names_multiple_and_boundary_seeds_cap4": "seeded_atlas_indices.nii.gz",
        "all_exact_names_multiple_and_boundary_seeds_cap1": "seeded_atlas_indices.nii.gz",
        "near_miss_names_only_empty_mask_cap3": "near_miss_atlas_indices.nii.gz",
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

        mask, distance = lesion_toads_gm_distance(
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

    seeded_case = reference["cases"][0]
    assert sum(seeded_case["mask_xyz_flat"]) == 13
    empty_case = reference["cases"][2]
    assert not any(empty_case["mask_xyz_flat"])
    assert set(empty_case["distance_xyz_flat"]) == {-1}
    cap_one = reference["cases"][1]
    # Java's first propagation scan occurs even with cap=1.
    assert 2 in cap_one["distance_xyz_flat"]


def test_gm_distance_rejects_invalid_atlas_indices_and_labels() -> None:
    with np.testing.assert_raises_regex(ValueError, "outside atlas_names"):
        lesion_toads_gm_distance(np.asarray([[[1]]]), ["Background"], 3)
    with np.testing.assert_raises_regex(TypeError, "integer atlas indices"):
        lesion_toads_gm_distance(np.asarray([[[0.0]]]), ["Background"], 3)
