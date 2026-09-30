from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_relations import toads_update_boundary_distances


REFERENCE = Path(__file__).parent / "data/synthetic/toads_boundary_distances/java_reference.json"


def test_boundary_distance_update_matches_selected_java_method() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    for case in reference["cases"]:
        shape = tuple(case["shape"])
        labels = np.asarray(case["segmentation"], dtype=np.int8).reshape(shape)
        initial = np.asarray(case["initial"], dtype=np.int8).reshape((*shape, case["classes"]))
        actual = toads_update_boundary_distances(labels, initial, case["limit"])
        expected = np.asarray(case["boundary_distances"], dtype=np.int8).reshape((*shape, case["classes"]))
        np.testing.assert_array_equal(actual, expected)
        if case["limit"] == 0:
            assert np.count_nonzero(actual) > 0  # Seeding precedes the optional diffusion loop.


@pytest.mark.parametrize(
    ("labels", "distances", "limit"),
    [
        (np.zeros((2, 3, 3), dtype=np.int8), np.zeros((2, 3, 3, 2), dtype=np.int8), 1),
        (np.zeros((3, 3, 3), dtype=np.float32), np.zeros((3, 3, 3, 2), dtype=np.int8), 1),
        (np.zeros((3, 3, 3), dtype=np.int8), np.zeros((2, 3, 3, 2), dtype=np.int8), 1),
        (np.zeros((3, 3, 3), dtype=np.int8), np.full((3, 3, 3, 2), 128), 1),
        (np.zeros((3, 3, 3), dtype=np.int8), np.zeros((3, 3, 3, 2), dtype=np.int8), -1),
    ],
)
def test_boundary_distance_update_rejects_invalid_inputs(labels, distances, limit) -> None:
    with pytest.raises(ValueError):
        toads_update_boundary_distances(labels, distances, limit)
