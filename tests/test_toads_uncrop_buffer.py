from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.toads_crop import uncrop_toads_buffer


DATA = Path(__file__).parent / "data/synthetic/toads_uncrop_buffer"
FIXTURE = DATA / "java_reference.json"
HARNESS = Path(__file__).parent / "java_reference/ToadsUncropBufferReferenceHarness.java"
HARNESS_SHA256 = "e3957825edb2b2a8233ea0744edf1f457299afd2dfaa567f06d673abc1ec2e8b"
FIXTURE_SHA256 = "ddff3cda6186906b2881b45ac37a9325950cc98816fbaaec2feea9b230c8367a"


def test_selected_java_uncrop_buffer_harness_and_fixture_are_pinned() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == FIXTURE_SHA256


def _synthetic_crops(case: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cx, cy, cz = case["crop_dims"]
    scalar = np.empty((cx, cy, cz), dtype=np.float32)
    labels = np.empty((cx, cy, cz), dtype=np.int8)
    channels = np.empty((2, cx, cy, cz), dtype=np.float32)
    for x in range(cx):
        for y in range(cy):
            for z in range(cz):
                scalar[x, y, z] = np.float32(x * 100 + y * 10 + z + 0.25)
                raw_label = 100 + 20 * x + 5 * y + z
                labels[x, y, z] = np.int8((raw_label + 128) % 256 - 128)
                channels[0, x, y, z] = np.float32(x * 1000 + y * 100 + z + 0.5)
                channels[1, x, y, z] = np.float32(-x * 1000 - y * 100 - z - 0.75)
    return scalar, labels, channels


def test_scalar_byte_and_multichannel_uncrop_buffers_match_selected_java_order() -> None:
    case = json.loads(FIXTURE.read_text())
    scalar, labels, channels = _synthetic_crops(case)
    shape = tuple(case["dims"])
    bounds = tuple(case["bounds"])
    padding = case["padding"]

    scalar_buffer = uncrop_toads_buffer(scalar, shape, bounds, padding, np.float32(-7.5))
    label_buffer = uncrop_toads_buffer(labels, shape, bounds, padding, np.int8(9))
    channels_buffer = uncrop_toads_buffer(channels, shape, bounds, padding, np.float32(2.25))

    assert scalar_buffer.dtype == np.float32
    assert label_buffer.dtype == np.int8
    assert channels_buffer.dtype == np.float32
    np.testing.assert_array_equal(scalar_buffer, np.asarray(case["scalar"], dtype=np.float32))
    np.testing.assert_array_equal(label_buffer, np.asarray(case["labels"], dtype=np.int8))
    np.testing.assert_array_equal(channels_buffer, np.asarray(case["channels"], dtype=np.float32))
    assert np.count_nonzero(scalar_buffer == np.float32(-7.5)) > 0
    assert np.count_nonzero(label_buffer == np.int8(9)) > 0
    assert np.count_nonzero(channels_buffer == np.float32(2.25)) > 0
    assert np.any(label_buffer < 0)  # Java byte narrowing keeps the signed byte representation.


def test_uncrop_buffer_rejects_crop_shape_that_disagrees_with_bounds() -> None:
    with np.testing.assert_raises_regex(ValueError, "does not match"):
        uncrop_toads_buffer(np.zeros((2, 2, 2), dtype=np.float32), (4, 5, 3), (-1, 1, 1, 2, -1, 0), 1, 0)
