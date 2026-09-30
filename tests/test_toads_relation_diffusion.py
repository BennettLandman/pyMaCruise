from __future__ import annotations

import json
import struct
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_relations import toads_compute_relations


REFERENCE = Path(__file__).parent / "data/synthetic/toads_relation_diffusion/java_reference.json"


def _float32_from_bits(values: list[int]) -> np.ndarray:
    packed = b"".join(struct.pack("!I", value) for value in values)
    return np.frombuffer(packed, dtype=np.dtype(">f4")).astype(np.float32)


def test_relation_diffusion_matches_selected_java_exact_membership_method() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    for case in reference["cases"]:
        shape = tuple(case["shape"])
        labels = np.asarray(case["segmentation"], dtype=np.int8).reshape(shape)
        memberships = _float32_from_bits(case["membership_values"]).reshape((*shape, case["classes"]))
        result = toads_compute_relations(labels, memberships, case["limit"])
        expected = np.asarray(case["relation_map"], dtype=np.int16).reshape(shape)
        np.testing.assert_array_equal(result.relation_map, expected)
        assert result.diffusion_iterations == case["diffusion_iterations"]
        assert result.relation_membership.shape[0] == case["initial_relation_count"]
        for axis in range(3):
            low = [slice(None)] * 3
            high = [slice(None)] * 3
            low[axis] = 0
            high[axis] = -1
            np.testing.assert_array_equal(result.relation_map[tuple(low)], labels[tuple(low)])
            np.testing.assert_array_equal(result.relation_map[tuple(high)], labels[tuple(high)])


@pytest.mark.parametrize(
    ("labels", "memberships", "limit"),
    [
        (np.zeros((3, 3, 3), dtype=np.int8), np.zeros((3, 3, 3), dtype=np.float32), 1),
        (np.zeros((3, 3, 3), dtype=np.int8), np.zeros((2, 3, 3, 2), dtype=np.float32), 1),
        (np.zeros((3, 3, 3), dtype=np.int8), np.full((3, 3, 3, 2), np.nan), 1),
        (np.zeros((3, 3, 3), dtype=np.int8), np.zeros((3, 3, 3, 1)), -1),
    ],
)
def test_relation_diffusion_rejects_invalid_inputs(labels, memberships, limit) -> None:
    with pytest.raises(ValueError):
        toads_compute_relations(labels, memberships, limit)
