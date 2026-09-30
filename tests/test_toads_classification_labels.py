from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_control import toads_classification_from_labels


REFERENCE = Path(__file__).parent / "data/synthetic/toads_classification_labels/java_reference.json"


def test_classification_transfer_matches_selected_java_method() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    for case in reference["cases"]:
        shape = tuple(case["shape"])
        segmentation = np.asarray(case["segmentation"], dtype=np.int32).reshape(shape)
        previous = np.asarray(case["previous"], dtype=np.int8).reshape(shape)
        output, changed = toads_classification_from_labels(
            segmentation, case["template_labels"], previous
        )
        np.testing.assert_array_equal(output.ravel().tolist(), case["classification"])
        assert output.dtype == np.int8
        assert changed == case["class_diff"]


@pytest.mark.parametrize(
    ("segmentation", "template", "previous"),
    [
        (np.zeros((2, 2), dtype=np.int8), [0], np.zeros((2, 2), dtype=np.int8)),
        (np.array([[[1]]]), [0], np.zeros((1, 1, 1), dtype=np.int8)),
        (np.zeros((1, 1, 1), dtype=np.float32), [0], np.zeros((1, 1, 1), dtype=np.int8)),
        (np.zeros((1, 1, 1), dtype=np.int8), [0], np.zeros((2, 1, 1), dtype=np.int8)),
        (np.zeros((1, 1, 1), dtype=np.int8), [128], np.zeros((1, 1, 1), dtype=np.int8)),
    ],
)
def test_classification_transfer_rejects_invalid_arrays(segmentation, template, previous) -> None:
    with pytest.raises(ValueError):
        toads_classification_from_labels(segmentation, template, previous)
