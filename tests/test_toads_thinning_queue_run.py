from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.toads_heap import ToadsPriorityHeap4D
from pymacruise.toads_thinning import toads_process_next_thinning_item, toads_run_thinning_queue


ROOT = Path(__file__).parent
STEP_REFERENCE = ROOT / "data/synthetic/toads_thinning_queue_step/java_reference.json"
RUN_REFERENCE = ROOT / "data/synthetic/toads_thinning_queue_run/java_reference.json"


def _from_bits(bits: int) -> np.float32:
    return np.asarray(np.int32(bits)).view(np.float32).item()


def _bits(value: np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.int32))


def _labels_text(array: np.ndarray) -> str:
    return "".join(str(int(v)) for v in array.ravel())


def _skeleton_text(array: np.ndarray) -> str:
    return "".join("a" if v < 0 else str(int(v)) for v in array.ravel())


def test_repeated_queue_run_matches_java_source_loop() -> None:
    step_reference = json.loads(STEP_REFERENCE.read_text())
    direct = next(case for case in step_reference["cases"] if case["name"] == "direct_allowed")["input"]
    reference = json.loads(RUN_REFERENCE.read_text())
    assert reference["source"] == "VandyCRUISE MultipleToadSegmentation.propagateThinning outer queue loop"
    assert reference["source_sha256"] == "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert {case["name"] for case in reference["cases"]} == {
        "drain_with_rescheduling", "cross_threshold_at_second_item",
        "top_priority_equals_threshold", "zero_limit_processes_nothing",
    }

    for case in reference["cases"]:
        segmentation = np.fromiter((int(v) for v in direct["segmentation"]), dtype=np.int16, count=125).reshape(5, 5, 5)
        skeleton = np.fromiter((-1 if v == "a" else int(v) for v in direct["skeleton"]), dtype=np.int16, count=125).reshape(5, 5, 5)
        available = np.fromiter((v == "1" for v in direct["available"]), dtype=np.bool_, count=500).reshape(5, 5, 5, 4)
        if case["suppress_scheduling"]:
            available.fill(False)
        memberships = np.full((5, 5, 5, 4), [0.1, 0.8, 0.05, 0.05], dtype=np.float32)
        memberships[2, 2, 2] = np.asarray(direct["center_membership"], dtype=np.float32)
        relations = np.zeros((4, 4), dtype=np.bool_)
        shapes = np.zeros((5, 5, 5, 4), dtype=np.float32)
        volumes = np.asarray(direct["volumes"], dtype=np.float32)
        minimums = np.asarray(direct["minimums"], dtype=np.float32)
        queue = ToadsPriorityHeap4D(2, "max")
        for value_bits, x, y, z, label in direct["queue"]:
            queue.add(_from_bits(value_bits), x, y, z, label)

        snapshot_labels, snapshot_skeleton = segmentation.copy(), skeleton.copy()
        snapshot_available, snapshot_volumes = available.copy(), volumes.copy()
        snapshot_queue = ToadsPriorityHeap4D(2, "max")
        for value_bits, x, y, z, label in direct["queue"]:
            snapshot_queue.add(_from_bits(value_bits), x, y, z, label)
        snapshots = case["snapshots"]
        for expected_item, expected_state in zip(case["processed"], snapshots, strict=True):
            step = toads_process_next_thinning_item(
                snapshot_labels, snapshot_skeleton, memberships,
                snapshot_available, relations, shapes, snapshot_volumes,
                minimums, snapshot_queue, case["limit"],
            )
            popped_step = step.popped_item
            assert [_bits(popped_step.value), popped_step.x, popped_step.y,
                    popped_step.z, popped_step.label] == expected_item
            snapshot_labels, snapshot_skeleton = step.segmentation, step.skeleton
            snapshot_available, snapshot_volumes = step.available, step.class_volumes
            assert _labels_text(snapshot_labels) == expected_state["segmentation"], case["name"]
            assert _skeleton_text(snapshot_skeleton) == expected_state["skeleton"], case["name"]

        result = toads_run_thinning_queue(
            segmentation, skeleton, memberships, available, relations, shapes,
            volumes, minimums, queue, case["limit"],
        )

        popped = [[_bits(item.value), item.x, item.y, item.z, item.label] for item in result.processed_items]
        assert popped == case["processed"], case["name"]
        assert _labels_text(result.segmentation) == case["segmentation"], case["name"]
        assert _skeleton_text(result.skeleton) == case["skeleton"], case["name"]
        assert "".join("1" if value else "0" for value in result.available.ravel()) == case["available"], case["name"]
        np.testing.assert_array_equal(result.class_volumes, np.asarray(case["volumes"], dtype=np.float32))
        remaining = []
        while queue:
            item = queue.pop_first()
            remaining.append([_bits(item.value), item.x, item.y, item.z, item.label])
        assert remaining == case["queue_after"], case["name"]

    by_name = {case["name"]: case for case in reference["cases"]}
    assert len(by_name["drain_with_rescheduling"]["processed"]) == 8
    assert len(by_name["cross_threshold_at_second_item"]["processed"]) == 2
    assert len(by_name["top_priority_equals_threshold"]["processed"]) == 1
    assert by_name["zero_limit_processes_nothing"]["processed"] == []
