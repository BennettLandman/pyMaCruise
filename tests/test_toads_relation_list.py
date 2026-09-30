from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_relations import toads_build_relation_list


REFERENCE = Path(__file__).parent / "data/synthetic/toads_relation_list/java_reference.json"


def test_relation_list_construction_matches_selected_java_method() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    for case in reference["cases"]:
        shape = tuple(case["shape"])
        labels = np.asarray(case["segmentation"], dtype=np.int8).reshape(shape)
        actual = toads_build_relation_list(labels, case["classes"])
        expected_map = np.asarray(case["relation_map"], dtype=np.int16).reshape(shape)
        expected_relations = np.asarray(case["relations"], dtype=bool).reshape(
            case["relation_count"], case["classes"]
        )
        np.testing.assert_array_equal(actual.relation_map, expected_map)
        np.testing.assert_array_equal(actual.relation_membership, expected_relations)
        assert actual.relation_membership.shape[0] == case["relation_count"]
        assert actual.class_count == case["classes"]
        assert np.all(actual.class_adjacency.diagonal())
        for axis in range(3):
            low = [slice(None)] * 3
            high = [slice(None)] * 3
            low[axis] = 0
            high[axis] = -1
            np.testing.assert_array_equal(actual.relation_map[tuple(low)], labels[tuple(low)])
            np.testing.assert_array_equal(actual.relation_map[tuple(high)], labels[tuple(high)])


@pytest.mark.parametrize(
    ("labels", "classes"),
    [
        (np.zeros((2, 3, 3), dtype=np.int8), 2),
        (np.zeros((3, 3, 3), dtype=np.float32), 2),
        (np.zeros((3, 3, 3), dtype=np.int8), 0),
        (np.zeros((3, 3, 3), dtype=np.int8), 64),
        (np.full((3, 3, 3), -1, dtype=np.int8), 2),
        (np.full((3, 3, 3), 2, dtype=np.int8), 2),
    ],
)
def test_relation_list_rejects_invalid_inputs(labels, classes) -> None:
    with pytest.raises(ValueError):
        toads_build_relation_list(labels, classes)
