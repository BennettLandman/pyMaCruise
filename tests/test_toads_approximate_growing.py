from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.toads_approximate_growing import toads_approximate_propagate_growing


ROOT = Path(__file__).parent
SOURCE = ROOT.parent.parent / "toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java"
HARNESS = ROOT / "java_reference/ToadsApproximateGrowingPassReferenceHarness.java"
REFERENCE = ROOT / "data/synthetic/toads_approximate_growing/java_reference.json"
SOURCE_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
HARNESS_SHA256 = "fe33597307cd64570f61d95bb0fad58f5cf9357055f728d13441d6e9f4b866a5"


def _inputs() -> tuple[np.ndarray, ...]:
    size, classes, slots = 9, 3, 2
    segmentation = np.empty((size, size, size), dtype=np.int8)
    skeleton = np.full_like(segmentation, -1)
    memberships = np.empty((size, size, size, slots), dtype=np.float32)
    best = np.empty((size, size, size, slots), dtype=np.int8)
    ordering = np.empty((size, size, size), dtype=np.float32)
    for x in range(size):
        for y in range(size):
            for z in range(size):
                folded = (x * 3 + y * 2) % 3 - 1
                dx, dy, dz = x - 4, y - 4, z - 4
                radius = 2 + folded
                if dx * dx + dy * dy + dz * dz < radius * radius:
                    current = 1
                else:
                    fold_z = 4 + (x + 2 * y) % 3 - 1
                    current = int(abs(z - fold_z) <= 1 and abs(y - 4) <= 3 and abs(x - 4) <= 3) * 2
                segmentation[x, y, z] = current
                ordering[x, y, z] = np.float32(-((x * 13 + y * 7 + z * 3) % 19)) / np.float32(10.0)
                for slot in range(slots):
                    best[x, y, z, slot] = current if slot == 0 else (current + 1) % classes
                    numerator = (x * 17 + y * 11 + z * 7 + slot * 13) % 39 + 1
                    memberships[x, y, z, slot] = np.float32(numerator) / np.float32(100.0)
                if 0 < x < size - 1 and 0 < y < size - 1 and 0 < z < size - 1 and (x * 2 + y * 3 + z * 5) % 7 == 0:
                    skeleton[x, y, z] = current
    volumes = np.asarray([1.0, 2.0, 3.0], dtype=np.float32)
    return segmentation, skeleton, memberships, best, ordering, volumes


def _bits(values) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.int32).astype(int).ravel().tolist()


def _label_string(values: np.ndarray, *, skeleton: bool = False) -> str:
    return "".join("a" if skeleton and value < 0 else str(int(value)) for value in values.ravel())


def test_approximate_growing_pass_matches_java_events_and_terminal_state() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    segmentation, skeleton, memberships, best, ordering, volumes = _inputs()
    result = toads_approximate_propagate_growing(
        segmentation, skeleton, memberships, best, ordering, volumes,
        class_count=3,
    )
    actual_events = []
    for event in result.events:
        item = event.popped_item
        inserted = [[int(v.value.view(np.int32)), v.x, v.y, v.z, v.label] for v in event.inserted_neighbors]
        topology = -1 if event.topology_accepted is None else int(event.topology_accepted)
        actual_events.append([[int(item.value.view(np.int32)), item.x, item.y, item.z, item.label], topology, inserted])
    assert actual_events == reference["events"]
    assert _label_string(result.segmentation) == reference["segmentation"]
    assert _label_string(result.skeleton, skeleton=True) == reference["skeleton"]
    assert "".join("1" if v else "0" for v in result.available.ravel()) == reference["available"]
    assert _bits(result.ordering) == reference["ordering_bits"]
    assert _bits(result.class_volumes) == reference["volumes_bits"]
    assert result.remaining_queue == ()
    assert any(event.topology_accepted is True for event in result.events)
    assert any(event.topology_accepted is False for event in result.events)


def test_approximate_growing_production_mode_does_not_retain_event_snapshots() -> None:
    segmentation, skeleton, memberships, best, ordering, volumes = _inputs()
    result = toads_approximate_propagate_growing(
        segmentation, skeleton, memberships, best, ordering, volumes,
        class_count=3, record_events=False,
    )
    assert result.events == ()
    assert result.remaining_queue == ()
    assert np.count_nonzero(result.skeleton[1:-1, 1:-1, 1:-1] >= 0) > np.count_nonzero(
        skeleton[1:-1, 1:-1, 1:-1] >= 0
    )
