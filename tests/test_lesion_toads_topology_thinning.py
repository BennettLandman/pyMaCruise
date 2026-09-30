from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_growth import propagate_lesion_toads_thinning


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_toads_topology_thinning"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())


def _signed_bits(value: np.float32 | float) -> int:
    bits = int(np.asarray(value, dtype=np.float32).view(np.uint32))
    return bits if bits < 0x80000000 else bits - 0x100000000


def test_topology_changing_lesion_thinning_pass_matches_java() -> None:
    assert REFERENCE["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    image = nib.load(str(ROOT / "topology_thinning_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (11, 11, 11, 14)
    segmentation = frames[..., 0].astype(np.int16)
    result = propagate_lesion_toads_thinning(
        segmentation,
        frames[..., 2:5],
        ("object", "object", "object"),
        frames[..., 5],
        frames[..., 6],
        frames[..., 7:10],
        frames[..., 1].astype(np.int16),
        frames[..., 10:13] > 0.5,
        frames[..., 13],
        lesion_class=1,
        volume_ratio=0.4,
        limit=10.0,
    )

    assert [_signed_bits(v) for v in result.initialization.class_volume] == REFERENCE["initial_volume_bits"]
    assert [_signed_bits(v) for v in result.minimum_class_volume] == REFERENCE["minimum_bits"]
    events = []
    for event in result.run.events:
        events.append([
            _signed_bits(event.popped_item.value), event.popped_item.x,
            event.popped_item.y, event.popped_item.z, event.popped_item.label,
            event.processed, event.candidate_label,
            event.topology_accepted,
            [[_signed_bits(item.value), item.x, item.y, item.z, item.label]
             for item in event.inserted_neighbors],
        ])
    assert events == REFERENCE["events"]
    np.testing.assert_array_equal(result.run.segmentation.ravel(), REFERENCE["segmentation"])
    np.testing.assert_array_equal(result.run.skeleton.ravel(), REFERENCE["skeleton"])
    np.testing.assert_array_equal(result.run.available.ravel(), REFERENCE["available"])
    assert [_signed_bits(v) for v in result.run.ordering.ravel()] == REFERENCE["ordering_bits"]
    assert [_signed_bits(v) for v in result.run.class_volume] == REFERENCE["volume_bits"]
    remaining = []
    while result.run.queue:
        item = result.run.queue.pop_first()
        remaining.append([_signed_bits(item.value), item.x, item.y, item.z, item.label])
    assert remaining == REFERENCE["remaining"]
    assert sum(event.topology_accepted is True for event in result.run.events) == REFERENCE["topology_accepted_count"] > 0
    assert sum(event.topology_accepted is False for event in result.run.events) == REFERENCE["topology_rejected_count"] > 0
    assert np.count_nonzero(result.run.segmentation != segmentation) == REFERENCE["changed_label_count"] > 0
    assert sum(len(event.inserted_neighbors) for event in result.run.events) == REFERENCE["inserted_neighbor_count"] > 0
