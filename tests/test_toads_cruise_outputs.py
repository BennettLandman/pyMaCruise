from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_cruise_outputs import (
    generate_toads_cruise_outputs,
    pack_toads_cruise_outputs,
)


DATA = Path(__file__).parent / "data/synthetic/toads_cruise_outputs"
HARNESS = Path(__file__).parent / "java_reference/ToadsCruiseOutputsReferenceHarness.java"
HARNESS_SHA256 = "30a5f3f793a5d1f9f399409cadefff40dbbc0ab228667e1ed843d6ca6907a2fc"
FIXTURE_SHA256 = "e9805204b7c0a25bc11e3ead5cbe0f4cb5f65a0d5619e42eba06f284441912a9"
SOURCE_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
WRAPPER_SHA256 = "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"


def _load_fixture() -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    reference = json.loads((DATA / "java_reference.json").read_text())
    memberships_img = nib.load(DATA / "memberships.nii.gz")
    segmentation_img = nib.load(DATA / "segmentation.nii.gz")
    best_img = nib.load(DATA / "best_classes.nii.gz")
    assert nib.aff2axcodes(memberships_img.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(segmentation_img.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(best_img.affine) == ("R", "A", "S")
    return (
        np.asarray(memberships_img.dataobj, dtype=np.float32),
        np.asarray(segmentation_img.dataobj, dtype=np.int16),
        np.asarray(best_img.dataobj, dtype=np.int16),
        reference,
    )


def test_harness_and_fixture_provenance_is_pinned() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    fixture = (DATA / "java_reference.json").read_bytes()
    assert hashlib.sha256(fixture).hexdigest() == FIXTURE_SHA256
    metadata = json.loads(fixture)
    assert metadata["source_sha256"] == SOURCE_SHA256
    assert metadata["wrapper_sha256"] == WRAPPER_SHA256


def test_exact_and_approximate_outputs_match_java_excerpt_on_asymmetric_ras_niftis() -> None:
    memberships, labels, best, reference = _load_fixture()
    names = tuple(reference["atlas_names"])
    assert memberships.shape == (*reference["shape_xyz"], len(names))
    for approximate, key in ((False, "exact"), (True, "approx")):
        output = generate_toads_cruise_outputs(
            memberships, labels, names, use_approx=approximate, best=best if approximate else None
        )
        assert output.as_channels().dtype == np.float32
        assert output.as_channels().shape == (4, *labels.shape)
        for channel_index, channel in enumerate(reference["channels"]):
            actual = getattr(output, channel)
            expected = np.asarray(reference["java_reference"][key][channel], dtype=np.float32).reshape(labels.shape)
            np.testing.assert_array_equal(actual, expected)
            assert output.as_channels()[channel_index].tobytes() == expected.tobytes()


def test_approximate_duplicate_matches_use_last_candidate_slot_and_alias_class() -> None:
    # Duplicate aliases are visited in atlas order; candidate slots are then
    # visited in ascending order inside each matching class.
    names = ("SulcalCSF", "Sulcal-CSF", "CorticalGM", "Cerebrum-GM", "WhiteMatter")
    memberships = np.asarray([[[[0.1, 0.2, 0.3, 0.4, 0.5]]]], dtype=np.float32)
    best = np.asarray([[[[0, 1, 1]]]], dtype=np.int16)
    labels = np.asarray([[[4]]], dtype=np.int16)
    result = generate_toads_cruise_outputs(memberships, labels, names, use_approx=True, best=best)
    # CSF alias 1 overwrites alias 0; two slots matching class 1 leave slot 2.
    np.testing.assert_array_equal(result.csf, np.asarray([[[0.3]]], dtype=np.float32))
    # No GM class occurs in a candidate slot.
    np.testing.assert_array_equal(result.gray_matter, np.zeros((1, 1, 1), dtype=np.float32))
    # Assigned WM label contributes a hard mask even though it is absent from best.
    np.testing.assert_array_equal(result.white_matter_mask, np.ones((1, 1, 1), dtype=np.float32))
    np.testing.assert_array_equal(result.white_matter_fill, np.zeros((1, 1, 1), dtype=np.float32))


def test_near_miss_names_are_case_sensitive_and_cerebell_prefix_is_excluded() -> None:
    names = ("Sulcalcsf", "corticalGM", "CerebellumXYZ", "Brainstem-like", "Background ")
    memberships = np.asarray([0.1, 0.2, 0.3, 0.4, 0.5], dtype=np.float32).reshape(1, 1, 1, 5)
    labels = np.asarray([3], dtype=np.int16).reshape(1, 1, 1)
    result = generate_toads_cruise_outputs(memberships, labels, names, use_approx=False)
    assert result.csf.item() == 0.0
    assert result.gray_matter.item() == 0.0
    # Only CerebellumXYZ is excluded by the prefix; Brainstem-like and a
    # whitespace-suffixed Background are eligible WM atlas classes.
    np.testing.assert_array_equal(result.white_matter_fill, np.asarray([[[1.2]]], dtype=np.float32))
    np.testing.assert_array_equal(result.white_matter_mask, np.ones((1, 1, 1), dtype=np.float32))


def test_approximate_mode_requires_valid_best_candidate_indices() -> None:
    memberships = np.zeros((1, 1, 1, 1), dtype=np.float32)
    labels = np.zeros((1, 1, 1), dtype=np.int16)
    names = ("Background",)
    with np.testing.assert_raises_regex(ValueError, "required"):
        generate_toads_cruise_outputs(memberships, labels, names, use_approx=True)
    with np.testing.assert_raises_regex(ValueError, "outside"):
        generate_toads_cruise_outputs(
            memberships,
            labels,
            names,
            use_approx=True,
            best=np.ones((1, 1, 1, 1), dtype=np.int16),
        )


def test_wrapper_uncrop_packs_four_images_separately_with_zero_background() -> None:
    memberships, labels, _best, reference = _load_fixture()
    products = generate_toads_cruise_outputs(
        memberships, labels, tuple(reference["atlas_names"]), use_approx=False
    )
    bounds = (1, 3, 1, 4, 1, 5)
    original_shape = (5, 6, 7)
    buffers = pack_toads_cruise_outputs(products, original_shape, bounds, padding=0)
    assert len(buffers) == 4
    for channel, buffer in zip(products.as_channels(), buffers):
        assert buffer.dtype == np.float32
        assert buffer.size == int(np.prod(original_shape))
        restored = np.zeros(original_shape, dtype=np.float32)
        restored[1:4, 1:5, 1:6] = channel
        expected = restored.transpose(2, 1, 0).reshape(-1)
        np.testing.assert_array_equal(buffer, expected)
