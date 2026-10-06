from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.toads_approximate_thinning_startup import initialize_toads_approximate_thinning
from pymacruise.toads_approximate_thinning_step import process_next_approximate_thinning_item
from pymacruise.toads_approximate_neighbor_schedule import schedule_approximate_thinning_neighbors


ROOT = Path(__file__).parent
REFERENCE = ROOT / "data/synthetic/toads_approximate_thinning_one_step/java_reference.json"
HARNESS = ROOT / "java_reference/ToadsApproximateThinningOneStepReferenceHarness.java"
SOURCE = ROOT.parent.parent / "toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java"
SOURCE_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
HARNESS_SHA256 = "44a21b9b7fd5c35e4e2bde88eb45110363d4295f54e4d983325243d414a3aef3"
NEIGHBOR_REFERENCE = ROOT / "data/synthetic/toads_approximate_thinning_one_step/neighbor_schedule_java_reference.json"
NEIGHBOR_HARNESS = ROOT / "java_reference/ToadsApproximateNeighborScheduleReferenceHarness.java"
NEIGHBOR_HARNESS_SHA256 = "483a3de5f045d093cd5a84978fd0c425415d2847b15daed0d6e8462a13cace9e"


def _arrays() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n, classes, slots = 9, 3, 2
    segmentation = np.empty((n, n, n), dtype=np.int8)
    best = np.empty((n, n, n, slots), dtype=np.int8)
    mems = np.empty((n, n, n, slots), dtype=np.float32)
    for x in range(n):
        for y in range(n):
            for z in range(n):
                folded = (x * 3 + y * 2) % 3 - 1
                dx, dy, dz = x - 4, y - 4, z - 4
                radius = 2 + folded
                radius2 = dx * dx + dy * dy + dz * dz
                if radius2 < radius * radius:
                    current = 1
                else:
                    fold_z = 4 + (x + 2 * y) % 3 - 1
                    current = 2 if abs(z - fold_z) <= 1 and abs(y - 4) <= 3 and abs(x - 4) <= 3 else 0
                segmentation[x, y, z] = current
                best[x, y, z, 0] = current
                best[x, y, z, 1] = (current + 1) % classes
                for slot in range(slots):
                    numerator = (x * 17 + y * 11 + z * 7 + slot * 13) % 39 + 1
                    mems[x, y, z, slot] = np.float32(numerator) / np.float32(100.0)
    return segmentation, mems, best


def test_sparse_approximate_thinning_initial_state_processes_one_java_pinned_event() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["shape"] == [9, 9, 9]
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256

    segmentation, mems, best = _arrays()
    startup = initialize_toads_approximate_thinning(
        segmentation,
        mems,
        best,
        class_count=3,
    )
    before_volumes = startup.class_volumes.copy()
    result = process_next_approximate_thinning_item(
        segmentation,
        mems,
        best,
        startup,
        np.zeros((3, 3), dtype=bool),
        np.zeros(segmentation.shape + (3,), dtype=np.float32),
        np.zeros(3, dtype=np.float32),
    )

    item = result.popped_item
    assert [int(item.value.view(np.int32)), item.x, item.y, item.z, item.label] == reference["popped"]
    assert result.transition is not None and result.transition.eligible
    assert result.segmentation[item.x, item.y, item.z] == reference["center_label"]
    assert result.skeleton[item.x, item.y, item.z] == reference["center_skeleton"]
    assert result.available[item.x, item.y, item.z].tolist() == reference["center_available"]
    assert result.class_volumes.view(np.int32).astype(int).tolist() == reference["volume_bits"]
    assert result.remaining_queue.size == reference["queue_remaining"]
    assert int(result.class_volumes[item.label].view(np.int32)) != int(before_volumes[item.label].view(np.int32))

    # This fixture selects slot 0 equal to the current tissue: it exercises
    # sparse volume removal and skeleton update while preserving segmentation.
    assert np.count_nonzero(result.segmentation != segmentation) == 0
    assert result.skeleton[item.x, item.y, item.z] == -1


def test_approximate_thinning_reschedules_face_neighbors_like_java() -> None:
    reference = json.loads(NEIGHBOR_REFERENCE.read_text())
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256
    assert hashlib.sha256(NEIGHBOR_HARNESS.read_bytes()).hexdigest() == NEIGHBOR_HARNESS_SHA256
    segmentation, mems, best = _arrays()
    startup = initialize_toads_approximate_thinning(
        segmentation, mems, best, class_count=3,
    )
    event = process_next_approximate_thinning_item(
        segmentation, mems, best, startup,
        np.zeros((3, 3), dtype=bool),
        np.zeros(segmentation.shape + (3,), dtype=np.float32),
        np.zeros(3, dtype=np.float32),
    )
    item = event.popped_item
    assert [int(item.value.view(np.int32)), item.x, item.y, item.z, item.label] == reference["popped"]
    queue_size_before_schedule = event.remaining_queue.size
    scheduled = schedule_approximate_thinning_neighbors(
        event.segmentation,
        mems,
        best,
        event.available,
        event.remaining_queue,
        (item.x, item.y, item.z),
        item.label,
        item.value,
        event.class_volumes[item.label],
        0.0,
        1000.0,
        copy_state=False,
    )
    actual = [
        [int(entry.value.view(np.int32)), entry.x, entry.y, entry.z, entry.label]
        for entry in scheduled.inserted
    ]
    assert actual == reference["inserted"]
    assert event.remaining_queue.size == queue_size_before_schedule + len(actual)
    for _, x, y, z, label in actual:
        assert not scheduled.available[x, y, z, label]
