from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_growth import propagate_lesion_toads_thinning


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_toads_thinning_pass"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())


def _signed_bits(value: np.float32 | float) -> int:
    bits = int(np.asarray(value, dtype=np.float32).view(np.uint32))
    return bits if bits < 0x80000000 else bits - 0x100000000


def test_selected_lesion_toads_thinning_pass_matches_java() -> None:
    assert REFERENCE["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    image = nib.load(str(ROOT / "thinning_pass_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (11, 11, 11, 14)
    segmentation = frames[..., 0].astype(np.int16)
    skeleton = frames[..., 1].astype(np.int16)
    memberships = frames[..., 2:5]
    lesion = frames[..., 5]
    black_hole = frames[..., 6]
    shapes = frames[..., 7:10]
    available = frames[..., 10:13] > 0.5
    ordering = frames[..., 13]
    result = propagate_lesion_toads_thinning(
        segmentation,
        memberships,
        ("object", "lesion", "object"),
        lesion,
        black_hole,
        shapes,
        skeleton,
        available,
        ordering,
        lesion_class=1,
        volume_ratio=0.4,
        limit=1000.0,
    )
    assert [_signed_bits(v) for v in result.initialization.class_volume] == REFERENCE["initial_volume_bits"]
    assert [_signed_bits(v) for v in result.minimum_class_volume] == REFERENCE["minimum_bits"]
    assert result.minimum_class_volume[0] > np.float32(100.0)
    np.testing.assert_array_equal(result.minimum_class_volume[1:], np.array([100.0, 100.0], dtype=np.float32))
    inserted = [
        [_signed_bits(item.value), item.x, item.y, item.z, item.label]
        for item in result.initialization.inserted_items
    ]
    assert inserted == REFERENCE["initial_inserted"]
    events = [
        [_signed_bits(item.popped_item.value), item.popped_item.x,
         item.popped_item.y, item.popped_item.z, item.popped_item.label, item.processed]
        for item in result.run.events
    ]
    assert events == REFERENCE["events"]
    np.testing.assert_array_equal(result.run.segmentation.ravel().tolist(), REFERENCE["segmentation"])
    np.testing.assert_array_equal(result.run.skeleton.ravel().tolist(), REFERENCE["skeleton"])
    np.testing.assert_array_equal(result.run.available.ravel().tolist(), REFERENCE["available"])
    assert [_signed_bits(v) for v in result.run.ordering.ravel()] == REFERENCE["ordering_bits"]
    assert [_signed_bits(v) for v in result.run.class_volume] == REFERENCE["volume_bits"]
    remaining = []
    while result.run.queue:
        item = result.run.queue.pop_first()
        remaining.append([_signed_bits(item.value), item.x, item.y, item.z, item.label])
    assert remaining == REFERENCE["remaining"]
