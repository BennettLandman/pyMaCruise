from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_cruise_inputs import lesion_toads_generate_cruise_inputs


FIXTURE_DIR = Path(__file__).parent / "data" / "synthetic" / "lesion_cruise_inputs"
HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsCruiseInputsReferenceHarness.java"
HARNESS_SHA256 = "b8e88b741d8bac82f1cc463abd2c31556a96044506d1f3fc3df04bf3f54a77db"
REFERENCE = FIXTURE_DIR / "java_reference.json"
REFERENCE_SHA256 = "464e51cc508bf77dd48e3f1ada930caf6bf566cc38c54db2eab648b28d1a67d1"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
FIXTURE_HASHES = {
    "memberships.nii.gz": "419206862abb90ac7891614de3eb7abfeff016723b7af36924c193e70877ceeb",
    "lesion_membership.nii.gz": "3a48b16796783cb794baa0c1eeec5932680d6404ac08fafda1c12e8ea09b139b",
    "black_hole_membership.nii.gz": "509d732585976e12591b9d5f655ef8d41c2639f49c5f01450e3817166a58ce13",
    "segmentation.nii.gz": "9da1b6ffd20357cd4819ee347821ff4efb1cfe8886b0abb140781c5452c0e3b8",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _float_bits(values: np.ndarray) -> list[int]:
    return values.astype(np.float32, copy=False).ravel(order="C").view(np.int32).tolist()


def test_generate_cruise_inputs_matches_selected_java_excerpt_on_ras_fixture() -> None:
    assert _sha256(HARNESS) == HARNESS_SHA256
    assert _sha256(REFERENCE) == REFERENCE_SHA256
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_method"] == "LesionToadSegmentation.generateCruiseInputs"
    assert reference["source_sha256"] == SOURCE_SHA256
    names = tuple(reference["atlas_names"])
    for filename, digest in FIXTURE_HASHES.items():
        assert _sha256(FIXTURE_DIR / filename) == digest

    loaded = {
        key: load_nifti_ras(FIXTURE_DIR / filename)
        for key, filename in {
            "memberships": "memberships.nii.gz",
            "lesion": "lesion_membership.nii.gz",
            "black_hole": "black_hole_membership.nii.gz",
            "segmentation": "segmentation.nii.gz",
        }.items()
    }
    for image in loaded.values():
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), (0.8, 1.2, 2.5))
        assert tuple(image.shape[:3]) == tuple(reference["shape_xyz"])
    memberships = np.asarray(loaded["memberships"].dataobj, dtype=np.float32)
    lesion = np.asarray(loaded["lesion"].dataobj, dtype=np.float32)
    black_hole = np.asarray(loaded["black_hole"].dataobj, dtype=np.float32)
    segmentation = np.asarray(loaded["segmentation"].dataobj, dtype=np.int16)

    actual = lesion_toads_generate_cruise_inputs(
        memberships, lesion, black_hole, segmentation, names
    )
    assert actual.dtype == np.float32
    assert actual.shape == (4, *reference["shape_xyz"])
    assert [_float_bits(channel) for channel in actual] == reference[
        "channel_float32_bits_xyz_flat"
    ]

    # Both alias families occur more than once. Java's loops overwrite with
    # the last exact spelling in atlas order rather than summing aliases.
    np.testing.assert_array_equal(actual[0], memberships[..., 2])
    np.testing.assert_array_equal(actual[1], memberships[..., 5])
    # Cerebellar, brainstem, background, and the excluded CSF/GM classes do
    # not set hard-WM; included classes do.
    included = {9, 10, 11, 12}
    np.testing.assert_array_equal(actual[3], np.isin(segmentation, list(included)).astype(np.float32))


def test_cruise_input_shapes_and_atlas_labels_are_validated() -> None:
    memberships = np.zeros((2, 3, 4, 1), dtype=np.float32)
    maps = np.zeros((2, 3, 4), dtype=np.float32)
    labels = np.zeros((2, 3, 4), dtype=np.int16)

    with pytest.raises(ValueError, match="class axis"):
        lesion_toads_generate_cruise_inputs(memberships, maps, maps, labels, ("Background", "Other"))
    with pytest.raises(ValueError, match="outside atlas_names"):
        lesion_toads_generate_cruise_inputs(
            memberships, maps, maps, labels + 1, ("Background",)
        )
    with pytest.raises(TypeError, match="integer atlas indices"):
        lesion_toads_generate_cruise_inputs(
            memberships, maps, maps, labels.astype(np.float32), ("Background",)
        )
