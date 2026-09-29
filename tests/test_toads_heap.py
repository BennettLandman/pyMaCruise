from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_heap import ToadsPriorityHeap4D


REFERENCE = Path(__file__).parent / "data" / "synthetic" / "toads_binary_heap" / "java_reference.json"


def _signed_bits(value: float) -> int:
    bits = int(np.asarray(np.float32(value)).view(np.int32))
    return bits


def test_binary_heap_matches_java_min_max_tie_growth_and_reuse() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source"] == "VandyCRUISE BinaryHeap4D"
    assert reference["source_sha256"] == "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
    assert {case["name"] for case in reference["cases"]} == {
        "min_growth_ties", "max_growth_ties", "change_mode_with_items", "reset_then_reuse"
    }

    for case in reference["cases"]:
        heap = ToadsPriorityHeap4D(case["capacity"], case["initial_mode"])
        for operation in case["operations"]:
            if operation[0] == "add":
                heap.add(*operation[1:])
            elif operation[0] == "reset":
                heap.reset()
            elif operation[0] == "set_min":
                heap.set_min_tree()
            elif operation[0] == "set_max":
                heap.set_max_tree()
            else:  # pragma: no cover - fixture integrity guard
                raise AssertionError(f"unknown heap operation {operation[0]}")

        actual = []
        while heap:
            first = heap.peek_first()
            item = heap.pop_first()
            assert item == first
            actual.append([_signed_bits(item.value), item.x, item.y, item.z, item.label])
        assert actual == case["popped"], case["name"]
        assert len(heap) == 0


def test_heap_rejects_invalid_configuration_and_empty_access() -> None:
    with pytest.raises(ValueError, match="capacity"):
        ToadsPriorityHeap4D(0)
    with pytest.raises(ValueError, match="mode"):
        ToadsPriorityHeap4D(2, "middle")
    heap = ToadsPriorityHeap4D(2)
    with pytest.raises(IndexError, match="empty"):
        heap.peek_first()
    with pytest.raises(IndexError, match="empty"):
        heap.pop_first()
