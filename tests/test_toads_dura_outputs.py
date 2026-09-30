from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.toads_dura_outputs import (
    adapt_toads_wrapper_outputs,
    generate_toads_dura_removal_outputs,
    pack_toads_dura_removal_outputs,
)


DATA = Path(__file__).parent / "data/synthetic/toads_dura_removal_outputs"
FIXTURE = DATA / "java_reference.json"
HARNESS = Path(__file__).parent / "java_reference/ToadsDuraRemovalOutputsReferenceHarness.java"
HARNESS_SHA256 = "d2af265c37a632fafa99c9061851c2eb4ffa9f930054a6b5d43624d434039b7d"
FIXTURE_SHA256 = "85dc342d530d82c07edbc2bfabf6c5243ceedd7472717a99d2dca62bb4a73437"
NAMES = (
    "SulcalCSF",
    "Sulcal-CSF",
    "CorticalGM",
    "Cerebellum-GM",
    "CerebellarGM",
    "CerebellarGM)",
    "WhiteMatter",
    "Background",
)


def test_selected_java_dura_output_harness_and_fixture_are_pinned() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == FIXTURE_SHA256


def _synthetic_inputs(case: dict) -> tuple[np.ndarray, np.ndarray]:
    cx, cy, cz = case["crop_shape"]
    mems = np.empty((cx, cy, cz, len(NAMES)), dtype=np.float32)
    segmentation = np.empty((cx, cy, cz), dtype=np.int32)
    for x in range(cx):
        for y in range(cy):
            for z in range(cz):
                segmentation[x, y, z] = (x + 2 * y + z) % len(NAMES)
                for k in range(len(NAMES)):
                    mems[x, y, z, k] = np.float32(0.01 * (1 + x + 10*y + 100*z + 1000*k))
    return mems, segmentation


def test_dura_products_and_four_image_uncrop_match_selected_java_harness() -> None:
    case = json.loads(FIXTURE.read_text())
    mems, segmentation = _synthetic_inputs(case)
    products = generate_toads_dura_removal_outputs(mems, segmentation, NAMES)
    buffers = pack_toads_dura_removal_outputs(
        products,
        tuple(case["original_shape"]),
        (-1, 1, 0, 1, -1, 0),
        1,
    )
    for actual, key in zip(
        buffers,
        ("raw_csf", "raw_gray_matter", "raw_white_matter_fill", "raw_white_matter_mask"),
    ):
        assert actual.dtype == np.float32
        np.testing.assert_array_equal(actual, np.asarray(case[key], dtype=np.float32))


def test_wrapper_restores_byte_classification_and_default_255_scaled_float_outputs() -> None:
    case = json.loads(FIXTURE.read_text())
    mems, segmentation = _synthetic_inputs(case)
    products = generate_toads_dura_removal_outputs(mems, segmentation, NAMES)
    result = adapt_toads_wrapper_outputs(
        products,
        segmentation,
        tuple(case["original_shape"]),
        (-1, 1, 0, 1, -1, 0),
        1,
    )

    assert result.classification.dtype == np.int8
    np.testing.assert_array_equal(result.classification, np.asarray(case["classification"], dtype=np.int8))
    raw_names = ("raw_csf", "raw_gray_matter", "raw_white_matter_fill", "raw_white_matter_mask")
    scaled_names = ("scaled_csf", "scaled_gray_matter", "scaled_white_matter_fill", "scaled_white_matter_mask")
    for raw, scaled, raw_key, scaled_key in zip(
        result.raw_float_buffers,
        result.wrapper_float_buffers,
        raw_names,
        scaled_names,
    ):
        np.testing.assert_array_equal(raw, np.asarray(case[raw_key], dtype=np.float32))
        np.testing.assert_array_equal(scaled, np.asarray(case[scaled_key], dtype=np.float32))
        np.testing.assert_array_equal(scaled, np.multiply(raw, np.float32(255.0), dtype=np.float32))


def test_cerebellar_gm_spelling_bug_and_csf_alias_order_are_preserved() -> None:
    names = ("SulcalCSF", "Sulcal-CSF", "Cerebellum-GM", "CerebellarGM", "CerebellarGM)", "Background")
    mems = np.zeros((5, 1, 1, len(names)), dtype=np.float32)
    mems[..., 0] = 0.125
    mems[..., 1] = 0.25
    mems[..., 2] = 0.5
    mems[..., 3] = 0.75
    mems[..., 4] = 0.875
    labels = np.asarray([2, 3, 4, 5, 5], dtype=np.int32).reshape(5, 1, 1)

    result = generate_toads_dura_removal_outputs(mems, labels, names)

    # The second CSF alias overwrites the first, as in the Java loop.
    np.testing.assert_array_equal(result.csf[:, 0, 0], np.full(5, 0.25, dtype=np.float32))
    # GM probabilities include both correct CerebellarGM and typo-labelled class.
    np.testing.assert_array_equal(result.gray_matter[:, 0, 0], np.full(5, 1.25, dtype=np.float32))
    # Correct CerebellarGM is intentionally omitted from this hard-mask branch;
    # Cerebellum-GM and CerebellarGM) are included.
    np.testing.assert_array_equal(result.white_matter_mask[:, 0, 0], [1, 0, 1, 0, 0])


def test_dura_output_rejects_invalid_class_indices_and_shapes() -> None:
    with np.testing.assert_raises_regex(ValueError, "match"):
        generate_toads_dura_removal_outputs(
            np.zeros((2, 2, 2, 1), dtype=np.float32),
            np.zeros((2, 2, 3), dtype=np.int32),
            ("Background",),
        )
    with np.testing.assert_raises_regex(ValueError, "outside"):
        generate_toads_dura_removal_outputs(
            np.zeros((2, 2, 2, 1), dtype=np.float32),
            np.full((2, 2, 2), 1, dtype=np.int32),
            ("Background",),
        )
