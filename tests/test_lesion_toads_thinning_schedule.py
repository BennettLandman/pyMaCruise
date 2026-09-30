from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_growth import schedule_lesion_toads_thinning_neighbors
from pymacruise.toads_heap import ToadsPriorityHeap4D


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_toads_thinning_schedule"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
CASES = (
    ("all_faces_ordinary", 0, False, -0.25, 8.0, 2.0, 100.0),
    ("lesion_max", 1, True, -0.125, 8.5, 2.25, 100.0),
    ("lesion_sum", 1, False, -0.125, 8.5, 2.25, 100.0),
    ("ratio_tie_limit_equal", 0, False, 0.0, 8.0, 2.0, 1.0),
    ("ratio_tie_limit_above", 0, False, 0.0, 8.0, 2.0, 1.01),
    ("volume_equal", 0, False, -0.25, 4.0, 4.0, 100.0),
    ("volume_below", 0, False, -0.25, 3.99, 4.0, 100.0),
    ("availability_gate", 0, False, -0.25, 8.0, 2.0, 100.0),
)


def _float_bits(values) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def test_lesion_toads_thinning_neighbor_loop_matches_java_and_heap() -> None:
    assert REFERENCE["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert REFERENCE["heap_sha256"] == "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
    image = nib.load(str(ROOT / "schedule_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert tuple(nib.affines.voxel_sizes(image.affine)) == (1.0, 1.0, 1.0)
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (5, 5, 5, len(CASES) * 6)
    assert [item["name"] for item in REFERENCE["cases"]] == [case[0] for case in CASES]

    for case_index, (name, removed_class, use_max, value, volume, minimum, limit) in enumerate(CASES):
        base = case_index * 6
        memberships = np.moveaxis(frames[..., base : base + 3], -1, -1)
        lesion = frames[..., base + 3]
        black_hole = frames[..., base + 4]
        available = np.zeros((*frames.shape[:3], 3), dtype=bool)
        available[..., removed_class] = frames[..., base + 5] > 0.5
        queue = ToadsPriorityHeap4D(16, mode="max")
        result = schedule_lesion_toads_thinning_neighbors(
            memberships,
            ("object", "lesion", "blackhole"),
            lesion,
            black_hole,
            available,
            queue,
            (2, 2, 2),
            removed_class,
            value,
            volume,
            minimum,
            limit,
            use_lesion_max=use_max,
        )
        reference = REFERENCE["cases"][case_index]
        np.testing.assert_array_equal(result.available.ravel().tolist(), reference["available_after"])
        inserted = [
            [int(item.value.view(np.uint32)) if int(item.value.view(np.uint32)) < 0x80000000
             else int(item.value.view(np.uint32)) - 0x100000000,
             item.x, item.y, item.z, item.label]
            for item in result.inserted_items
        ]
        assert inserted == reference["inserted"]
        popped = []
        while queue:
            item = queue.pop_first()
            bits = int(item.value.view(np.uint32))
            popped.append([bits if bits < 0x80000000 else bits - 0x100000000,
                           item.x, item.y, item.z, item.label])
        assert popped == reference["popped"]
