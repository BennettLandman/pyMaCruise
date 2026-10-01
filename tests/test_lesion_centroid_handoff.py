from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_centroid_handoff import (
    LesionToadsCentroidHandoffState,
    initialize_lesion_toads_segmentation_centroids,
)


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "lesion_segmentation_init_centroids"
HARNESS = ROOT / "java_reference" / "LesionToadSegmentationInitCentroidsReferenceHarness.java"
REFERENCE = FIXTURE / "java_reference.json"
SOURCE_HASHES = {
    "segmentation": "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99",
    "wrapper": "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb",
}


def _bits(values: np.ndarray) -> list[list[int]]:
    return values.view(np.uint32).astype(np.uint64).tolist()


def test_init_centroids_matches_selected_java_copy_handoff() -> None:
    ref = json.loads(REFERENCE.read_text())
    fixture = json.loads((FIXTURE / "input.json").read_text())
    assert ref["source_method"] == "LesionToadSegmentation.initCentroids"
    assert ref["source_sha256"] == SOURCE_HASHES
    assert ref["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    assert fixture["shape"] == ref["shape"] == [2, 4]
    values = np.asarray(fixture["input_bits"], dtype=np.uint32).view(np.float32)
    before = values.copy()
    state = LesionToadsCentroidHandoffState.empty(2, 4)

    initialize_lesion_toads_segmentation_centroids(state, values)

    assert _bits(state.centroid) == ref["centroid_bits"] == fixture["input_bits"]
    assert _bits(state.prior_centroid) == ref["prior_centroid_bits"] == fixture["input_bits"]
    assert not np.shares_memory(state.centroid, values)
    assert not np.shares_memory(state.prior_centroid, values)
    np.testing.assert_array_equal(values.view(np.uint32), before.view(np.uint32))
