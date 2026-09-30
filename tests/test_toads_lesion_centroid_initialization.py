from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_lesion_centroid_initialization import initial_toads_lesion_centroid

ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_lesion_centroid_initialization"
REFERENCE = ROOT / "java_reference.json"
HARNESS = Path(__file__).parent / "java_reference" / "ToadsLesionCentroidInitializationReferenceHarness.java"
HARNESS_SHA256 = "abe57ca7e1adc6ec4e4f2685e5534edad87a86f111ee5557eb0d892aba5a4cc1"


def _bits(values: np.ndarray) -> list[int]:
    return [int(value) for value in np.asarray(values, dtype=np.float32).view(np.uint32)]


def test_initial_lesion_centroid_matches_selected_java_expression() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert reference["source_sha256"] == {
        "AlgorithmMultipleToads.java": "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d",
        "MultipleToadPreprocess.java": "31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7",
    }
    for case in reference["cases"]:
        low = np.asarray(case["low"], dtype=np.float32)
        high = np.asarray(case["high"], dtype=np.float32)
        prior = np.asarray(case["prior"], dtype=np.float32)
        saved = tuple(array.copy() for array in (low, high, prior))
        result = initial_toads_lesion_centroid(low, high, prior)
        assert result.shape == low.shape
        assert _bits(result) == [int(value) & 0xFFFFFFFF for value in case["centroid_bits"]]
        for actual, before in zip((low, high, prior), saved, strict=True):
            np.testing.assert_array_equal(actual, before)

    # Source adds the scalar prior after low/span and does not clamp the result.
    last = reference["cases"][-1]
    assert initial_toads_lesion_centroid(last["low"], last["high"], last["prior"])[-1] > 1.0


@pytest.mark.parametrize(
    ("low", "high", "prior"),
    [
        ([1.0, 2.0], [3.0, 4.0], [0.1]),
        ([1.0, 2.0], [3.0, 4.0], [[0.1, 0.2]]),
        ([1.0], [1.0], [0.1]),
        ([1.0], [2.0], [np.nan]),
    ],
)
def test_initial_lesion_centroid_rejects_invalid_inputs(low, high, prior) -> None:
    with pytest.raises(ValueError):
        initial_toads_lesion_centroid(low, high, prior)
