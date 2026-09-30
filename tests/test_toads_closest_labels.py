from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_relations import toads_closest_labels


REFERENCE = Path(__file__).parent / "data/synthetic/toads_closest_labels/java_reference.json"


def test_closest_label_seed_and_diffusion_match_selected_java_method() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    for case in reference["cases"]:
        shape = tuple(case["shape"])
        segmentation = np.asarray(case["segmentation"], dtype=np.int8).reshape(shape)
        actual = toads_closest_labels(segmentation, case["limit"])
        expected = np.asarray(case["closest_labels"], dtype=np.int8).reshape((*shape, 4))
        np.testing.assert_array_equal(actual, expected)
        assert actual.dtype == np.int8
        # The source initializes all boundary slots but never diffuses them.
        for axis in range(3):
            low = [slice(None)] * 4
            high = [slice(None)] * 4
            low[axis] = 0
            high[axis] = -1
            assert np.all(actual[tuple(low)][..., 0] == segmentation[tuple(low[:3])])
            assert np.all(actual[tuple(high)][..., 0] == segmentation[tuple(high[:3])])
            assert np.all(actual[tuple(low)][..., 1:] == -1)
            assert np.all(actual[tuple(high)][..., 1:] == -1)


@pytest.mark.parametrize(
    ("segmentation", "limit"),
    [
        (np.zeros((2, 3, 3), dtype=np.int8), 1),
        (np.zeros((3, 3, 3), dtype=np.float32), 1),
        (np.full((3, 3, 3), -1, dtype=np.int8), 1),
        (np.full((3, 3, 3), 128, dtype=np.int16), 1),
        (np.zeros((3, 3, 3), dtype=np.int8), -1),
        (np.zeros((3, 3, 3), dtype=np.int8), True),
    ],
)
def test_closest_label_rejects_invalid_inputs(segmentation, limit) -> None:
    with pytest.raises(ValueError):
        toads_closest_labels(segmentation, limit)
