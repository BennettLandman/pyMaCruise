from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_special_centroid_handoff import (
    LesionToadsSpecialCentroidState,
    initialize_lesion_toads_special_centroids,
)


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "lesion_special_centroid_handoff"
HARNESS = ROOT / "java_reference" / "LesionToadSegmentationSpecialCentroidSettersReferenceHarness.java"
REFERENCE = FIXTURE / "java_reference.json"
SOURCE_HASHES = {
    "segmentation": "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99",
    "wrapper": "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb",
}


def _bits(values: np.ndarray) -> list[int]:
    return values.view(np.uint32).astype(np.uint64).tolist()


def test_special_centroid_setters_match_java_reference_identity_and_order() -> None:
    ref = json.loads(REFERENCE.read_text())
    fixture = json.loads((FIXTURE / "input.json").read_text())
    assert ref["source_sha256"] == SOURCE_HASHES
    assert ref["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    lesion = np.asarray(fixture["lesion_bits"], dtype=np.uint32).view(np.float32)
    black_hole = np.asarray(fixture["black_hole_bits"], dtype=np.uint32).view(np.float32)
    lesion_before, black_hole_before = lesion.copy(), black_hole.copy()
    state = LesionToadsSpecialCentroidState(
        lesion_centroid=np.zeros_like(lesion),
        black_hole_centroid=np.ones_like(black_hole),
    )

    initialize_lesion_toads_special_centroids(state, lesion, black_hole)

    assert ref["source_methods"] == [
        "LesionToadSegmentation.setLesionCentroid",
        "LesionToadSegmentation.setBlackHoleCentroid",
    ]
    assert ref["call_order"] == ["lesion", "black_hole"]
    assert ref["lesion_input_bits"] == fixture["lesion_bits"]
    assert ref["black_hole_input_bits"] == fixture["black_hole_bits"]
    assert state.lesion_centroid is lesion is not black_hole
    assert state.black_hole_centroid is black_hole
    assert _bits(state.lesion_centroid) == ref["lesion_bits"] == fixture["lesion_bits"]
    assert _bits(state.black_hole_centroid) == ref["black_hole_bits"] == fixture["black_hole_bits"]
    np.testing.assert_array_equal(lesion.view(np.uint32), lesion_before.view(np.uint32))
    np.testing.assert_array_equal(black_hole.view(np.uint32), black_hole_before.view(np.uint32))
