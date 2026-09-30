from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_crop import uncrop_toads_buffer
from pymacruise.toads_thalamus_output import generate_thalamus_output


FIXTURE_DIR = Path(__file__).parent / "data/synthetic/toads_thalamus_output"
HARNESS = Path(__file__).parent / "java_reference/ToadsThalamusOutputReferenceHarness.java"
SOURCE_SHA256 = "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
WRAPPER_SHA256 = "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"


def _load(name: str) -> np.ndarray:
    image = nib.load(FIXTURE_DIR / f"{name}.nii.gz")
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(nib.affines.voxel_sizes(image.affine), (1.0, 1.0, 1.0))
    return np.asarray(image.dataobj)


def _same_float32_bits(actual: np.ndarray, expected: np.ndarray) -> None:
    actual32 = np.ascontiguousarray(actual, dtype=np.float32)
    expected32 = np.ascontiguousarray(expected, dtype=np.float32)
    assert actual32.shape == expected32.shape
    assert actual32.tobytes() == expected32.tobytes()


def test_fixture_and_excerpt_provenance_are_pinned() -> None:
    reference_path = FIXTURE_DIR / "java_reference.json"
    reference = json.loads(reference_path.read_text())
    assert reference["source_method"] == "MultipleToadSegmentation.generateThalamusOutput"
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["wrapper_sha256"] == WRAPPER_SHA256
    assert reference["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    for filename, expected_sha256 in reference["fixture_sha256"].items():
        assert hashlib.sha256((FIXTURE_DIR / filename).read_bytes()).hexdigest() == expected_sha256
    assert tuple(reference["shape_xyz"]) == (4, 3, 5)
    assert reference["affine_ras"] == np.eye(4).tolist()


def test_exact_and_approximate_modes_match_selected_java_excerpt_bitwise() -> None:
    reference = json.loads((FIXTURE_DIR / "java_reference.json").read_text())
    memberships = _load("memberships")
    best = _load("best_classes")
    names = tuple(reference["atlas_names"])
    assert memberships.shape == (*reference["shape_xyz"], len(names))
    assert best.shape == (*reference["shape_xyz"], reference["approx_slots"])
    assert np.count_nonzero(memberships == 0) > 0
    assert np.count_nonzero(memberships) > 0

    for mode in ("exact", "approx"):
        output = generate_thalamus_output(
            memberships,
            names,
            use_approx=(mode == "approx"),
            best=best if mode == "approx" else None,
        )
        expected_reference = np.asarray(
            reference["java_reference"][mode], dtype=np.float32
        ).reshape(tuple(reference["shape_xyz"]))
        expected_image = _load(f"expected_{mode}")
        assert output.shape == (1, *reference["shape_xyz"])
        assert output.dtype == np.float32
        _same_float32_bits(output[0], expected_reference)
        _same_float32_bits(output[0], expected_image)


def test_exact_name_duplicates_near_matches_and_zero_initialization() -> None:
    memberships = np.asarray([[[[0.25, 0.5, 0.75, 1.0]]]], dtype=np.float32)
    # Only an exact, case-sensitive and whitespace-sensitive name is selected;
    # the later duplicate exact name overwrites the first one.
    names = ("Thalamus", "thalamus", "Thalamus ", "Thalamus")
    exact = generate_thalamus_output(memberships, names, use_approx=False)
    np.testing.assert_array_equal(exact, np.asarray([[[[1.0]]]], dtype=np.float32))

    near_matches = generate_thalamus_output(
        memberships[..., :3], ("thalamus", "Thalamus ", "THALAMUS"), use_approx=False
    )
    np.testing.assert_array_equal(near_matches, np.zeros((1, 1, 1, 1), dtype=np.float32))
    approximate_near_matches = generate_thalamus_output(
        memberships[..., :3],
        ("thalamus", "Thalamus ", "THALAMUS"),
        use_approx=True,
        best=np.asarray([[[[0, 1, 2]]]], dtype=np.int16),
    )
    np.testing.assert_array_equal(
        approximate_near_matches, np.zeros((1, 1, 1, 1), dtype=np.float32)
    )


def test_approximate_duplicate_and_candidate_slot_overwrite_order() -> None:
    reference = json.loads((FIXTURE_DIR / "java_reference.json").read_text())
    result = generate_thalamus_output(
        _load("memberships"),
        tuple(reference["atlas_names"]),
        use_approx=True,
        best=_load("best_classes"),
    )[0]
    # At [0,0,0], class 0 is followed by duplicate class 4, whose matching
    # slots 1 and 3 overwrite in order. Slot 3 membership is 0.34.
    assert result[0, 0, 0] == np.float32(0.34)
    # At [0,0,1], duplicate class 4 matches slots 0 and 1; slot 1 wins (0.19).
    assert result[0, 0, 1] == np.float32(0.19)
    # At [0,0,2], no exact Thalamus class is selected, so the initialized zero remains.
    assert result[0, 0, 2] == np.float32(0.0)
    # At [0,0,3], only the first exact class is selected and slot 0 is copied.
    assert result[0, 0, 3] == np.float32(0.22)


def test_wrapper_uncrop_helper_places_channel_and_zero_fills_outside_crop() -> None:
    cropped = _load("expected_approx")
    original_shape = (6, 5, 7)
    bounds = (1, 4, 1, 3, 1, 5)
    restored = uncrop_toads_buffer(cropped, original_shape, bounds, 0, 0.0)
    expected = np.zeros(original_shape, dtype=np.float32)
    expected[1:5, 1:4, 1:6] = cropped
    np.testing.assert_array_equal(restored, expected.transpose(2, 1, 0).reshape(-1))


def test_invalid_approximate_candidate_arrays_are_rejected() -> None:
    memberships = np.zeros((1, 1, 1, 1), dtype=np.float32)
    with np.testing.assert_raises_regex(ValueError, "required"):
        generate_thalamus_output(memberships, ("Thalamus",), use_approx=True)
    with np.testing.assert_raises_regex(ValueError, "outside"):
        generate_thalamus_output(
            memberships,
            ("Thalamus",),
            use_approx=True,
            best=np.ones((1, 1, 1, 1), dtype=np.int16),
        )
    with np.testing.assert_raises_regex(ValueError, "exceed"):
        generate_thalamus_output(
            np.zeros((1, 1, 1, 1), dtype=np.float32),
            ("Thalamus",),
            use_approx=True,
            best=np.zeros((1, 1, 1, 2), dtype=np.int16),
        )


def test_complex_membership_arrays_are_rejected() -> None:
    with np.testing.assert_raises_regex(TypeError, "real numeric"):
        generate_thalamus_output(
            np.zeros((1, 1, 1, 1), dtype=np.complex64),
            ("Thalamus",),
            use_approx=False,
        )
