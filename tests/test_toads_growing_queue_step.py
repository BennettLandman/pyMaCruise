from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.toads_heap import ToadsPriorityHeap4D
from pymacruise.toads_thinning import toads_process_next_growing_item


REFERENCE = Path(__file__).parent / "data/synthetic/toads_growing_queue_step/java_reference.json"


def _bits(value: np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.int32))


def _from_bits(bits: int) -> np.float32:
    return np.asarray(np.int32(bits)).view(np.float32).item()


def _arrays(name: str):
    size, classes = 5, 4
    segmentation = np.zeros((size, size, size), dtype=np.int16)
    skeleton = np.zeros_like(segmentation)
    available = np.zeros((size, size, size, classes), dtype=np.bool_)
    if name == "direct_allowed":
        segmentation[1, 2, 2] = skeleton[1, 2, 2] = 1
        skeleton[2, 2, 2] = -1
        for delta in ((-1, 0, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (0, 1, 0), (1, 0, 0)):
            x, y, z = (2 + value for value in delta)
            if skeleton[x, y, z] != 1:
                skeleton[x, y, z] = -1
                available[x, y, z, 1] = True
        available[1, 2, 2, 1] = False
    elif name == "direct_rejected":
        encoded = "020033121011213000023100323"
        patch = np.fromiter((int(value) for value in encoded), dtype=np.int16).reshape(3, 3, 3)
        segmentation[1:4, 1:4, 1:4] = patch
        skeleton[1:4, 1:4, 1:4] = patch
        skeleton[2, 2, 2] = -1
    elif name == "stale_skeleton":
        segmentation[2, 2, 2] = 1
        skeleton[2, 2, 2] = 0
    else:
        raise AssertionError(name)

    memberships = np.empty((size, size, size, classes), dtype=np.float32)
    for x in range(size):
        for y in range(size):
            for z in range(size):
                for k in range(classes):
                    memberships[x, y, z, k] = np.float32(
                        ((x * 5 + y * 7 + z * 11 + k * 3) % 29 + 1) / 30.0
                    )
    ordering = np.zeros((size, size, size), dtype=np.float32)
    volumes = np.asarray([10, 20, 30, 40], dtype=np.float32)
    queued_label = {"direct_allowed": 1, "direct_rejected": 3, "stale_skeleton": 3}[name]
    queue = ToadsPriorityHeap4D(8, "min")
    queue.add(0.25, 2, 2, 2, queued_label)
    queue.add(0.5, 3, 3, 3, 2)
    return segmentation, skeleton, memberships, available, ordering, volumes, queue


def _labels(array: np.ndarray, *, skeleton: bool = False) -> str:
    if skeleton:
        return "".join("a" if value < 0 else str(int(value)) for value in array.ravel())
    return "".join(str(int(value)) for value in array.ravel())


def test_one_growing_queue_event_matches_java_lut_and_heap_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source"] == "VandyCRUISE MultipleToadSegmentation.propagateGrowing one queue item"
    assert reference["source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert reference["lut_sha256"] == "50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f"
    assert reference["heap_sha256"] == "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
    assert {case["name"] for case in reference["cases"]} == {
        "direct_allowed", "direct_rejected", "stale_skeleton"
    }

    for case in reference["cases"]:
        arrays = _arrays(case["name"])
        result = toads_process_next_growing_item(*arrays)
        segmentation, skeleton, memberships, available, ordering, volumes, queue = arrays
        assert [_bits(result.popped_item.value), result.popped_item.x, result.popped_item.y,
                result.popped_item.z, result.popped_item.label] == case["popped"]
        assert result.topology_accepted is case["topology_accepted"], case["name"]
        assert _labels(result.segmentation) == case["segmentation"], case["name"]
        assert _labels(result.skeleton, skeleton=True) == case["skeleton"], case["name"]
        assert "".join("1" if value else "0" for value in result.available.ravel()) == case["available"]
        assert [_bits(value) for value in result.ordering[2, 2, 2 : 3]] == case["ordering_bits"]
        assert [_bits(value) for value in result.class_volumes] == case["volumes_bits"]
        inserted = [[_bits(item.value), item.x, item.y, item.z, item.label]
                    for item in result.inserted_neighbors]
        assert inserted == case["inserted"], case["name"]
        remaining = []
        while queue:
            item = queue.pop_first()
            remaining.append([_bits(item.value), item.x, item.y, item.z, item.label])
        assert remaining == case["queue_after"], case["name"]

    by_name = {case["name"]: case for case in reference["cases"]}
    assert by_name["direct_allowed"]["topology_accepted"] is True
    assert len(by_name["direct_allowed"]["inserted"]) == 5
    assert by_name["direct_rejected"]["topology_accepted"] is False
    assert by_name["direct_rejected"]["inserted"] == []
    assert by_name["stale_skeleton"]["topology_accepted"] is None
