from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.toads_heap import ToadsPriorityHeap4D
from pymacruise.toads_thinning import toads_process_next_thinning_item


REFERENCE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_thinning_queue_step"
    / "java_reference.json"
)


def _from_bits(bits: int) -> np.float32:
    return np.asarray(np.int32(bits)).view(np.float32).item()


def _bits(value: np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.int32))


def _skeleton_text(skeleton: np.ndarray) -> str:
    return "".join("a" if value < 0 else str(int(value)) for value in skeleton.ravel())


def test_one_queue_event_matches_java_lut_heap_and_source_order() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source"] == "VandyCRUISE MultipleToadSegmentation.propagateThinning one queue event"
    assert reference["source_sha256"] == "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert reference["lut_source_sha256"] == "50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f"
    assert reference["heap_sha256"] == "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
    assert {case["name"] for case in reference["cases"]} == {
        "direct_allowed", "direct_rejected", "related_allowed", "no_connection",
        "below_floor", "stale_skeleton", "queued_label_differs_from_segmentation",
    }

    for case in reference["cases"]:
        initial = case["input"]
        segmentation = np.fromiter(
            (int(v) for v in initial["segmentation"]), dtype=np.int16, count=125
        ).reshape(5, 5, 5)
        skeleton = np.fromiter(
            (-1 if v == "a" else int(v) for v in initial["skeleton"]),
            dtype=np.int16,
            count=125,
        ).reshape(5, 5, 5)
        available = np.fromiter(
            (v == "1" for v in initial["available"]), dtype=np.bool_, count=500
        ).reshape(5, 5, 5, 4)
        memberships = np.empty((5, 5, 5, 4), dtype=np.float32)
        memberships[...] = np.asarray([0.1, 0.8, 0.05, 0.05], dtype=np.float32)
        memberships[2, 2, 2] = np.asarray(initial["center_membership"], dtype=np.float32)
        relations = np.zeros((4, 4), dtype=np.bool_)
        shapes = np.zeros((5, 5, 5, 4), dtype=np.float32)
        if case["name"] == "related_allowed":
            relations[2, 1] = True
            shapes[1, 2, 2, 1] = np.float32(0.7)
        volumes = np.asarray(initial["volumes"], dtype=np.float32)
        minimums = np.asarray(initial["minimums"], dtype=np.float32)
        queue = ToadsPriorityHeap4D(2, "max")
        for value_bits, x, y, z, label in initial["queue"]:
            queue.add(_from_bits(value_bits), x, y, z, label)

        result = toads_process_next_thinning_item(
            segmentation,
            skeleton,
            memberships,
            available,
            relations,
            shapes,
            volumes,
            minimums,
            queue,
            initial["limit"],
        )

        assert "".join(str(int(v)) for v in result.segmentation.ravel()) == case["segmentation"], case["name"]
        assert _skeleton_text(result.skeleton) == case["skeleton"], case["name"]
        assert "".join("1" if v else "0" for v in result.available.ravel()) == case["available"], case["name"]
        np.testing.assert_array_equal(result.class_volumes, np.asarray(case["volumes"], dtype=np.float32))
        popped = result.popped_item
        assert [_bits(popped.value), popped.x, popped.y, popped.z, popped.label] == initial["queue"][0]
        remaining = []
        while queue:
            item = queue.pop_first()
            remaining.append([_bits(item.value), item.x, item.y, item.z, item.label])
        assert remaining == case["queue_after"], case["name"]

    by_name = {case["name"]: case for case in reference["cases"]}
    assert by_name["direct_allowed"]["segmentation"][62] == "1"
    assert by_name["direct_rejected"]["segmentation"][62] == "1"
    assert by_name["related_allowed"]["segmentation"][62] == "2"
    assert by_name["no_connection"]["segmentation"][62] == "0"
