from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_dura_outputs import generate_lesion_dura_removal_outputs


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_dura_outputs"
REFERENCE = ROOT / "java_reference.json"
HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsDuraRemovalOutputsReferenceHarness.java"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
HARNESS_SHA256 = "9b9a395cedc27f0f63126432a90c7c2a022bf8c05fc53c18e8493fd8e8d4fbf7"
REFERENCE_SHA256 = "f9ced8b7269e1fb71d506b16ff1abe8f739b31625d46fc33892dc45d2ee7cb37"


def test_lesion_dura_products_match_pinned_selected_java_excerpt_on_ras_nifti() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(REFERENCE.read_bytes()).hexdigest() == REFERENCE_SHA256
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_method"] == "LesionToadSegmentation.generateDuraRemovalOutputs"
    assert reference["source_sha256"] == SOURCE_SHA256

    names = tuple(reference["names"])
    loaded = {
        name: load_nifti_ras(ROOT / f"{name}.nii.gz")
        for name in ("memberships", "segmentation", "lesion_membership", "black_hole_membership")
    }
    for image in loaded.values():
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert image.shape[:3] == tuple(reference["shape"])

    actual = generate_lesion_dura_removal_outputs(
        np.asarray(loaded["memberships"].dataobj, dtype=np.float32),
        np.asarray(loaded["segmentation"].dataobj, dtype=np.int16),
        names,
        np.asarray(loaded["lesion_membership"].dataobj, dtype=np.float32),
        np.asarray(loaded["black_hole_membership"].dataobj, dtype=np.float32),
    ).as_channels()
    assert actual.shape == (4, *reference["shape"])
    assert actual.dtype == np.float32
    assert all(channel.dtype == np.float32 for channel in actual)
    expected = np.asarray(reference["channels"], dtype=np.float32).reshape(actual.shape)
    np.testing.assert_array_equal(actual, expected)

    # The selected LesionTOADS adapter includes lesion + black-hole maps in WM
    # fill, does not clip values, and retains its exact hard-label behavior.
    assert float(actual[2].max()) > 1.0
    assert actual[1, 0, 0, 0] > 0.0  # GM probabilities are summed.
    assert actual[0, 0, 0, 0] == np.float32(0.014)  # last CSF alias wins.
    assert actual[3, 3, 0, 1] == 0.0  # correctly spelled CerebellarGM is not hard-set.


def test_lesion_dura_output_zero_defaults_and_unclipped_lesion_initialization() -> None:
    names = ("Other", "Background")
    memberships = np.zeros((2, 1, 1, 2), dtype=np.float32)
    labels = np.ones((2, 1, 1), dtype=np.int16)
    lesion = np.asarray([0.8, 0.7], dtype=np.float32).reshape(2, 1, 1)
    black_hole = np.asarray([0.5, 0.6], dtype=np.float32).reshape(2, 1, 1)
    result = generate_lesion_dura_removal_outputs(
        memberships, labels, names, lesion, black_hole
    )
    np.testing.assert_array_equal(result.csf, 0.0)
    np.testing.assert_array_equal(result.gray_matter, 0.0)
    expected_fill = np.add(lesion, black_hole, dtype=np.float32)[:, 0, 0]
    np.testing.assert_array_equal(result.white_matter_fill[:, 0, 0], expected_fill)
    np.testing.assert_array_equal(result.white_matter_mask, 0.0)


def test_lesion_dura_outputs_reject_invalid_shapes_and_class_indices() -> None:
    with np.testing.assert_raises_regex(ValueError, "match"):
        generate_lesion_dura_removal_outputs(
            np.zeros((2, 1, 1, 1), dtype=np.float32),
            np.zeros((3, 1, 1), dtype=np.int16),
            ("Background",),
            np.zeros((2, 1, 1), dtype=np.float32),
            np.zeros((2, 1, 1), dtype=np.float32),
        )
    with np.testing.assert_raises_regex(ValueError, "outside atlas_names"):
        generate_lesion_dura_removal_outputs(
            np.zeros((2, 1, 1, 1), dtype=np.float32),
            np.ones((2, 1, 1), dtype=np.int16),
            ("Background",),
            np.zeros((2, 1, 1), dtype=np.float32),
            np.zeros((2, 1, 1), dtype=np.float32),
        )
