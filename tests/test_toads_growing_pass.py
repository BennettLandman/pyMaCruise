from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.toads_thinning import toads_propagate_growing


REFERENCE = Path(__file__).parent / "data/synthetic/toads_growing_pass/java_reference.json"


def _bits(values) -> list[int]:
    return [int(np.asarray(value, dtype=np.float32).view(np.int32)) for value in np.asarray(values).ravel()]


def _inputs():
    size, classes = 9, 4
    segmentation = np.zeros((size, size, size), dtype=np.int16)
    skeleton = np.full_like(segmentation, -1)
    ordering = np.empty((size, size, size), dtype=np.float32)
    memberships = np.empty((size, size, size, classes), dtype=np.float32)
    for x in range(size):
        for y in range(size):
            for z in range(size):
                ordering[x, y, z] = np.float32(-((x * 13 + y * 7 + z * 3) % 19)) / np.float32(10.0)
                for k in range(classes):
                    memberships[x, y, z, k] = np.float32(
                        ((x * 5 + y * 7 + z * 11 + k * 3) % 29 + 1)
                    ) / np.float32(30.0)
                center = 2 + ((y + z) % 3)
                if 2 <= y <= 6 and 2 <= z <= 6 and abs(x - center) <= 1:
                    segmentation[x, y, z] = 1
                if 5 <= x <= 6 and 2 <= y <= 3 and 5 <= z <= 6:
                    segmentation[x, y, z] = 2
                if x == 1 and 2 <= y <= 6 and z == 7:
                    segmentation[x, y, z] = 3
    center = 2 + ((4 + 3) % 3)
    skeleton[center, 4, 3] = 1
    skeleton[5, 2, 5] = 2
    volumes = np.asarray([1, 2, 3, 4], dtype=np.float32)
    return segmentation, skeleton, memberships, ordering, volumes


def _labels(values, *, skeleton: bool = False) -> str:
    if skeleton:
        return "".join("a" if value < 0 else str(int(value)) for value in values.ravel())
    return "".join(str(int(value)) for value in values.ravel())


def test_initialized_growing_pass_matches_java_lut_and_heap_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source"] == "VandyCRUISE MultipleToadSegmentation.initGrowingLabels + propagateGrowing queue exhaustion"
    assert reference["source_sha256"] == "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert reference["lut_sha256"] == "50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f"
    assert reference["heap_sha256"] == "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
    assert reference["shape"] == [9, 9, 9, 4]

    segmentation, skeleton, memberships, ordering, volumes = _inputs()
    result = toads_propagate_growing(segmentation, skeleton, memberships, ordering, volumes)
    python_events = []
    for event in result.events:
        item = event.popped_item
        inserted = [
            [_bits([neighbor.value])[0], neighbor.x, neighbor.y, neighbor.z, neighbor.label]
            for neighbor in event.inserted_neighbors
        ]
        code = -1 if event.topology_accepted is None else int(event.topology_accepted)
        python_events.append(
            [[_bits([item.value])[0], item.x, item.y, item.z, item.label], code, inserted]
        )

    assert python_events == reference["events"]
    assert _labels(result.segmentation) == reference["segmentation"]
    assert _labels(result.skeleton, skeleton=True) == reference["skeleton"]
    assert "".join("1" if value else "0" for value in result.available.ravel()) == reference["available"]
    assert _bits(result.ordering) == reference["ordering_bits"]
    assert _bits(result.class_volumes) == reference["volumes_bits"]
    assert [
        [_bits([item.value])[0], item.x, item.y, item.z, item.label]
        for item in result.remaining_queue
    ] == reference["queue_after"]
    assert result.remaining_queue == ()

    topology_codes = [event[1] for event in reference["events"]]
    assert any(code == 1 for code in topology_codes)
    assert any(code == 0 for code in topology_codes)
    assert any(code == -1 for code in topology_codes)
    assert {code: topology_codes.count(code) for code in (-1, 0, 1)} == {-1: 51, 0: 60, 1: 337}
    assert len(reference["events"]) == 448
