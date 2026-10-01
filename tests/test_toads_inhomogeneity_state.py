"""Selected Java comparison for TOADS per-channel field attachment state."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_inhomogeneity_state import ToadsInhomogeneityFieldState


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_inhomogeneity_field_state"
HARNESS = Path(__file__).parent / "java_reference" / "ToadsInhomogeneityFieldAttachmentReferenceHarness.java"
EXPECTED_HARNESS_SHA256 = "419d0c2c94e0e1fb9b8356ac067defea5e78460a527e917244437b70dc74cc95"
JAVA = json.loads((ROOT / "java_reference.json").read_text())
SOURCE_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"


def test_field_attachment_matches_selected_java_reference() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == EXPECTED_HARNESS_SHA256
    assert JAVA["source_sha256"] == SOURCE_SHA256
    assert JAVA["source_method"] == "MultipleToadSegmentation.addInhomogeneityCorrection"

    field_image = nib.load(str(ROOT / "field.nii.gz"))
    assert nib.aff2axcodes(field_image.affine) == ("R", "A", "S")
    field = np.asarray(field_image.dataobj, dtype=np.float32)
    np.testing.assert_array_equal(field, np.asarray(JAVA["fixture_field"], dtype=np.float32).reshape(field.shape))

    state = ToadsInhomogeneityFieldState(channel_count=2)
    assert state.use_field is False
    assert state.field_by_channel == [None, None]

    state.add_inhomogeneity_correction(field, 1)
    assert state.use_field is JAVA["assigned"]["use_field"]
    assert state.field_by_channel[1] is field
    assert state.field_by_channel[0] is None

    field[0, 0, 0] = np.float32(3.0)
    assert float(state.field_by_channel[1][0, 0, 0]) == JAVA["mutation_visible_value"]

    replacement = np.full(field.shape, 2.0, dtype=np.float32)
    state.add_inhomogeneity_correction(replacement, 1)
    assert state.field_by_channel[1] is replacement
    state.add_inhomogeneity_correction(None, 0)
    assert state.field_by_channel[0] is None
    assert state.use_field is True

    for key in (
        "initial_false",
        "assigned_by_identity",
        "other_channel_untouched",
        "enabled_after_assignment",
        "mutation_visible",
        "replacement_by_identity",
        "null_accepted_and_enabled",
    ):
        assert JAVA["reference_state"][key] is True


def test_field_attachment_channel_index_follows_array_bounds() -> None:
    state = ToadsInhomogeneityFieldState(channel_count=2)
    try:
        state.add_inhomogeneity_correction(np.ones((2, 2, 2), dtype=np.float32), 2)
    except IndexError:
        pass
    else:
        raise AssertionError("out-of-range channel must not be accepted")
    assert state.use_field is False
