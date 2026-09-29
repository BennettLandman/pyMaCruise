from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.toads_heap import ToadsPriorityHeap4D
from pymacruise.toads_thinning import toads_schedule_thinning_neighbors


REFERENCE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_thinning_neighbor_schedule"
    / "java_reference.json"
)


def _float_from_bits(value: int) -> np.float32:
    return np.asarray(np.int32(value)).view(np.float32).item()


def _int_bits(value: np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.int32))


def _memberships() -> np.ndarray:
    out = np.empty((5, 5, 5, 3), dtype=np.float32)
    for x in range(5):
        for y in range(5):
            for z in range(5):
                for k in range(3):
                    numerator = (x * 5 + y * 7 + z * 11 + k * 3) % 29 + 1
                    out[x, y, z, k] = np.float32(numerator) / np.float32(30.0)
    return out


def test_thinning_neighbor_scheduling_matches_java_source_and_heap() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source"] == "VandyCRUISE MultipleToadSegmentation.propagateThinning neighbor scheduling"
    assert reference["source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert reference["heap_sha256"] == "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
    assert {c["name"] for c in reference["cases"]} == {
        "all_six", "single_count_tie", "volume_equal", "limit_equal", "limit_above", "availability_gate"
    }
    mems = _memberships()
    for case in reference["cases"]:
        segmentation = np.asarray(case["segmentation"], dtype=np.int16).reshape(5, 5, 5)
        available = np.asarray(case["available_before"], dtype=np.bool_).reshape(5, 5, 5, 3)
        queue = ToadsPriorityHeap4D(2, "max")
        for bits, x, y, z, label in case["seeds"]:
            queue.add(_float_from_bits(bits), x, y, z, label)
        updated, inserted = toads_schedule_thinning_neighbors(
            segmentation,
            mems,
            available,
            queue,
            (2, 2, 2),
            removed_label=0,
            popped_priority=_float_from_bits(case["value_bits"]),
            class_volume=_float_from_bits(case["volume_bits"]),
            minimum_class_volume=_float_from_bits(case["minimum_bits"]),
            limit=_float_from_bits(case["limit_bits"]),
        )
        assert updated.ravel().tolist() == case["available_after"], case["name"]
        inserted_ref = [[_int_bits(item.value), item.x, item.y, item.z, item.label] for item in inserted]
        assert inserted_ref == case["inserted"], case["name"]
        popped = []
        while queue:
            item = queue.pop_first()
            popped.append([_int_bits(item.value), item.x, item.y, item.z, item.label])
        assert popped == case["popped"], case["name"]
