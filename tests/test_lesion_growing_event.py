from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_growth import process_lesion_toads_growing_event
from pymacruise.toads_heap import ToadsPriorityHeap4D


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_growing_event/java_reference.json"
FACE = ((-1, 0, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (0, 1, 0), (1, 0, 0))


def _bits(value: np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.int32))


def _arrays(case: dict[str, object]):
    name = str(case["name"])
    n, classes = 5, 4
    segmentation = np.zeros((n, n, n), dtype=np.int8)
    skeleton = np.zeros_like(segmentation)
    available = np.zeros((n, n, n, classes), dtype=np.bool_)
    if name.startswith("accepted_"):
        segmentation[1, 2, 2] = 1
        skeleton[1, 2, 2] = 1
        skeleton[2, 2, 2] = -1
        for dx, dy, dz in FACE:
            x, y, z = 2 + dx, 2 + dy, 2 + dz
            if skeleton[x, y, z] != 1:
                skeleton[x, y, z] = -1
                available[x, y, z, 1] = True
    elif name == "rejected":
        encoded = "020033121011213000023100323"
        patch = np.fromiter((int(value) for value in encoded), dtype=np.int8).reshape(3, 3, 3)
        segmentation[1:4, 1:4, 1:4] = patch
        skeleton[1:4, 1:4, 1:4] = patch
        skeleton[2, 2, 2] = -1
    elif name == "stale":
        segmentation[2, 2, 2] = 1
        skeleton[2, 2, 2] = 0
    else:
        raise AssertionError(name)

    memberships = np.empty((n, n, n, classes), dtype=np.float32)
    lesion = np.empty((n, n, n), dtype=np.float32)
    black_hole = np.empty((n, n, n), dtype=np.float32)
    for x in range(n):
        for y in range(n):
            for z in range(n):
                for k in range(classes):
                    memberships[x, y, z, k] = np.float32(
                        ((x * 5 + y * 7 + z * 11 + k * 3) % 29 + 1) / 30.0
                    )
                lesion[x, y, z] = np.float32(((x + 3 * y + 5 * z) % 7 + 1) / 10.0)
                black_hole[x, y, z] = np.float32(((2 * x + y + 3 * z) % 9 + 1) / 12.0)
    ordering = np.zeros((n, n, n), dtype=np.float32)
    volume = np.asarray([10, 20, 30, 40], dtype=np.float32)
    queue = ToadsPriorityHeap4D(32, "min")
    label = 1 if name.startswith("accepted_") else 3
    queue.add(0.25, 2, 2, 2, label)
    queue.add(0.5, 3, 3, 3, 2)
    object_types = (
        ("obj", "obj", "obj", "obj")
        if name == "accepted_nonlesion"
        else ("obj", "lesion", "blackhole", "optimized_lesion")
    )
    return (
        segmentation,
        skeleton,
        memberships,
        available,
        ordering,
        volume,
        object_types,
        lesion,
        black_hole,
        queue,
    )


def _encode_labels(array: np.ndarray, *, empty: bool = False) -> str:
    if empty:
        return "".join("a" if value < 0 else str(int(value)) for value in array.ravel())
    return "".join(str(int(value)) for value in array.ravel())


def test_one_lesion_growing_event_matches_java_lut_heap_and_source_formulas() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert reference["lut_sha256"] == "50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f"
    assert reference["heap_sha256"] == "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
    assert reference["numerics_sha256"] == "13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb"
    assert {case["name"] for case in reference["cases"]} == {
        "accepted_lesion_sum",
        "accepted_lesion_max",
        "accepted_nonlesion",
        "rejected",
        "stale",
    }

    for case in reference["cases"]:
        (
            segmentation,
            skeleton,
            memberships,
            available,
            ordering,
            volume,
            object_types,
            lesion,
            black_hole,
            queue,
        ) = _arrays(case)
        result = process_lesion_toads_growing_event(
            segmentation,
            skeleton,
            memberships,
            available,
            ordering,
            volume,
            object_types,
            lesion,
            black_hole,
            queue,
            lesion_class=int(case["lesion_class"]),
            use_lesion_max=bool(case["use_lesion_max"]),
        )
        popped = result.popped_item
        assert [_bits(popped.value), popped.x, popped.y, popped.z, popped.label] == [
            1048576000,
            2,
            2,
            2,
            int(case["label"]),
        ]
        assert result.topology_accepted is case["topology_accepted"], case["name"]
        assert _encode_labels(result.segmentation) == case["segmentation"], case["name"]
        assert _encode_labels(result.skeleton, empty=True) == case["skeleton"], case["name"]
        assert "".join("1" if value else "0" for value in result.available.ravel()) == case["available"]
        assert [_bits(value) for value in result.ordering.ravel()] == case["ordering_bits"]
        assert [_bits(value) for value in result.class_volume] == case["volume_bits"]
        inserted = [
            [_bits(item.value), item.x, item.y, item.z, item.label]
            for item in result.inserted_neighbors
        ]
        assert inserted == case["inserted"], case["name"]
        queue_after = []
        while queue:
            item = queue.pop_first()
            queue_after.append([_bits(item.value), item.x, item.y, item.z, item.label])
        assert queue_after == case["queue_after"], case["name"]

    cases = {case["name"]: case for case in reference["cases"]}
    assert cases["accepted_lesion_sum"]["topology_accepted"] is True
    assert cases["accepted_lesion_max"]["topology_accepted"] is True
    assert cases["accepted_nonlesion"]["topology_accepted"] is True
    assert cases["rejected"]["topology_accepted"] is False
    assert cases["stale"]["topology_accepted"] is None


def test_growing_event_rejects_empty_queue() -> None:
    volume = np.zeros(2, dtype=np.float32)
    image = np.zeros((3, 3, 3), dtype=np.float32)
    memberships = np.zeros((3, 3, 3, 2), dtype=np.float32)
    queue = ToadsPriorityHeap4D(1, "min")
    with pytest.raises(ValueError, match="empty"):
        process_lesion_toads_growing_event(
            np.zeros((3, 3, 3), dtype=np.int8),
            np.full((3, 3, 3), -1, dtype=np.int8),
            memberships,
            np.zeros((3, 3, 3, 2), dtype=np.bool_),
            image,
            volume,
            ("obj", "lesion"),
            image,
            image,
            queue,
            lesion_class=1,
        )
