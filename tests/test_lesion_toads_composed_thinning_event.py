from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_growth import process_next_lesion_toads_thinning_item
from pymacruise.toads_heap import ToadsPriorityHeap4D


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_toads_composed_thinning_event"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())


def _signed_bits(value: np.float32) -> int:
    bits = int(np.asarray(value, dtype=np.float32).view(np.uint32))
    return bits if bits < 0x80000000 else bits - 0x100000000


def test_composed_lesion_toads_queue_event_matches_java_and_heap() -> None:
    assert REFERENCE["popped"] == [-1098907648, 2, 2, 2, 0]
    image = nib.load(str(ROOT / "event_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (5, 5, 5, 14)
    labels = frames[..., 0].astype(np.int16)
    skeleton = frames[..., 1].astype(np.int16)
    memberships = np.moveaxis(frames[..., 2:5], -1, -1)
    lesion = frames[..., 5]
    black_hole = frames[..., 6]
    shapes = frames[..., 7:10]
    available = np.zeros((5, 5, 5, 3), dtype=bool)
    available[..., 0] = frames[..., 10] > 0.5
    ordering = frames[..., 13]
    volumes = np.array([12.0, 17.0, 9.0], dtype=np.float32)
    minimums = np.array([2.0, 2.0, 2.0], dtype=np.float32)
    queue = ToadsPriorityHeap4D(16, mode="max")
    queue.add(-0.25, 2, 2, 2, 0)
    queue.add(-0.5, 4, 4, 4, 2)

    result = process_next_lesion_toads_thinning_item(
        labels,
        skeleton,
        memberships,
        ("object", "lesion", "object"),
        lesion,
        black_hole,
        shapes,
        available,
        ordering,
        volumes,
        minimums,
        queue,
        100.0,
        lesion_class=1,
    )
    assert [_signed_bits(result.popped_item.value), result.popped_item.x,
            result.popped_item.y, result.popped_item.z, result.popped_item.label] == REFERENCE["popped"]
    assert result.transition is not None
    assert result.transition.candidate_label == REFERENCE["candidate"]
    assert result.transition.topology_accepted is REFERENCE["topology"]
    assert int(result.segmentation[2, 2, 2]) == REFERENCE["seg"]
    assert int(result.skeleton[2, 2, 2]) == REFERENCE["skeleton"]
    assert result.available[2, 2, 2].tolist() == REFERENCE["available_center"]
    assert _signed_bits(result.ordering[2, 2, 2]) == REFERENCE["ordering_bits"]
    assert [_signed_bits(v) for v in result.class_volume] == REFERENCE["volume_bits"]
    inserted = [
        [_signed_bits(item.value), item.x, item.y, item.z, item.label]
        for item in result.inserted_neighbors
    ]
    assert inserted == REFERENCE["inserted"]
    remaining = []
    while queue:
        item = queue.pop_first()
        remaining.append([_signed_bits(item.value), item.x, item.y, item.z, item.label])
    assert remaining == REFERENCE["remaining"]


def test_composed_event_discards_stale_and_volume_blocked_items() -> None:
    image = nib.load(str(ROOT / "event_inputs.nii.gz"))
    frames = np.asarray(image.dataobj, dtype=np.float32)
    labels = frames[..., 0].astype(np.int16)
    skeleton = frames[..., 1].astype(np.int16)
    memberships = frames[..., 2:5]
    lesion = frames[..., 5]
    black_hole = frames[..., 6]
    shapes = frames[..., 7:10]
    availability = np.zeros((5, 5, 5, 3), dtype=bool)
    availability[..., 0] = frames[..., 10] > 0.5
    ordering = frames[..., 13]
    minimums = np.array([2.0, 2.0, 2.0], dtype=np.float32)

    for processing_class, volumes in (
        (2, np.array([12.0, 17.0, 9.0], dtype=np.float32)),  # skeleton at center is class 0
        (0, np.array([2.0, 17.0, 9.0], dtype=np.float32)),  # strict volume > minimum fails
    ):
        queue = ToadsPriorityHeap4D(4, mode="max")
        queue.add(0.25, 2, 2, 2, processing_class)
        result = process_next_lesion_toads_thinning_item(
            labels, skeleton, memberships, ("object", "lesion", "object"),
            lesion, black_hole, shapes, availability, ordering, volumes,
            minimums, queue, 100.0, lesion_class=1,
        )
        assert result.transition is None
        assert result.inserted_neighbors == ()
        np.testing.assert_array_equal(result.segmentation, labels)
        np.testing.assert_array_equal(result.skeleton, skeleton)
        np.testing.assert_array_equal(result.available, availability)
        np.testing.assert_array_equal(result.ordering, ordering)
        np.testing.assert_array_equal(result.class_volume, volumes)
        assert not queue
