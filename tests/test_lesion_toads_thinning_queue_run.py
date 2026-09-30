from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_growth import run_lesion_toads_thinning_queue
from pymacruise.toads_heap import ToadsPriorityHeap4D


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_toads_thinning_queue_run"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())


def _signed_bits(value: np.float32 | float) -> int:
    bits = int(np.asarray(value, dtype=np.float32).view(np.uint32))
    return bits if bits < 0x80000000 else bits - 0x100000000


def test_lesion_toads_outer_queue_loop_matches_java_reference() -> None:
    assert REFERENCE["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    image = nib.load(str(ROOT / "queue_run_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (5, 5, 5, 14)
    segmentation = frames[..., 0].astype(np.int16)
    skeleton = frames[..., 1].astype(np.int16)
    memberships = frames[..., 2:5]
    lesion = frames[..., 5]
    black_hole = frames[..., 6]
    shapes = frames[..., 7:10]
    available = np.zeros((5, 5, 5, 3), dtype=bool)
    available[..., 0] = frames[..., 10] > 0.5
    ordering = frames[..., 13]
    volumes = np.array([20.0, 1.0, 0.0], dtype=np.float32)
    minimums = np.zeros(3, dtype=np.float32)
    queue = ToadsPriorityHeap4D(256, mode="max")
    queue.add(0.5, 2, 2, 2, 0)
    queue.add(-1.0, 3, 3, 3, 1)  # stale: skeleton label differs
    queue.add(-1.25, 3, 1, 3, 2)  # volume is at its minimum
    queue.add(-2.0, 1, 1, 1, 0)  # source processes this crossing item once
    queue.add(-2.5, 1, 3, 1, 0)  # retained because the next loop check stops

    result = run_lesion_toads_thinning_queue(
        segmentation, skeleton, memberships, ("object", "object", "object"),
        lesion, black_hole, shapes, available, ordering, volumes, minimums,
        queue, 2.0, lesion_class=1,
    )
    processed = [
        [_signed_bits(event.popped_item.value), event.popped_item.x,
         event.popped_item.y, event.popped_item.z, event.popped_item.label,
         event.processed]
        for event in result.events
    ]
    assert processed == REFERENCE["processed"]
    assert any(event[0] == _signed_bits(-2.0) and event[-1] for event in processed)
    np.testing.assert_array_equal(result.segmentation.ravel().tolist(), REFERENCE["segmentation"])
    np.testing.assert_array_equal(result.skeleton.ravel().tolist(), REFERENCE["skeleton"])
    np.testing.assert_array_equal(result.available.ravel().tolist(), REFERENCE["available"])
    assert [_signed_bits(v) for v in result.ordering.ravel()] == REFERENCE["ordering_bits"]
    assert [_signed_bits(v) for v in result.class_volume] == REFERENCE["volume_bits"]
    remaining = []
    while result.queue:
        item = result.queue.pop_first()
        remaining.append([_signed_bits(item.value), item.x, item.y, item.z, item.label])
    assert remaining == REFERENCE["remaining"]
