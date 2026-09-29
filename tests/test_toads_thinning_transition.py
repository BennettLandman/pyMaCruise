from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.toads_thinning import toads_thinning_voxel_transition


REFERENCE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_thinning_transition"
    / "java_reference.json"
)


def _labels(value: str) -> np.ndarray:
    return np.fromiter((int(label) for label in value), dtype=np.int16, count=27).reshape(3, 3, 3)


def test_thinning_transition_matches_java_source_and_lut_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source"] == "VandyCRUISE MultipleToadSegmentation.propagateThinning local queue item"
    assert reference["source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert reference["lut_source_sha256"] == "50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f"
    assert {case["kind"] for case in reference["cases"]} == {
        "same", "low_volume", "not_skeleton", "no_connection", "direct_allowed",
        "direct_rejected", "related_allowed", "related_rejected", "related_tie",
    }

    for case in reference["cases"]:
        result = toads_thinning_voxel_transition(
            _labels(case["segmentation"]),
            _labels(case["skeleton"]),
            np.asarray(case["memberships"], dtype=np.float32),
            np.asarray(case["relations"], dtype=np.bool_),
            np.asarray(case["shape"], dtype=np.float32).reshape(3, 3, 3, 4),
            np.asarray(case["volumes"], dtype=np.float32),
            np.asarray(case["minimums"], dtype=np.float32),
            available_at_center=np.asarray(case["available"], dtype=np.bool_),
        )
        expected = case["expected"]
        assert result.segmentation_label == expected["segmentation_label"], case["kind"]
        assert result.skeleton_label == expected["skeleton_label"], case["kind"]
        np.testing.assert_array_equal(result.class_volumes, np.asarray(expected["volumes"], dtype=np.float32))
        np.testing.assert_array_equal(result.available_at_center, np.asarray(expected["available"], dtype=np.bool_))
        assert result.candidate_label == expected["candidate_label"], case["kind"]
        assert result.topology_accepted is expected["topology_accepted"], case["kind"]
        assert result.eligible is expected["eligible"], case["kind"]
