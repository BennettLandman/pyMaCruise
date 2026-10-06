from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.toads_thinning import toads_initialize_thinning_queue, toads_propagate_thinning


REFERENCE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_initialized_thinning_pass"
    / "java_reference.json"
)


def _bits(value: np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.int32))


def _from_bits(value: int) -> np.float32:
    return np.asarray(np.int32(value)).view(np.float32).item()


def _phantom() -> tuple[np.ndarray, np.ndarray]:
    size, classes = 12, 4
    segmentation = np.empty((size, size, size), dtype=np.int16)
    memberships = np.empty((size, size, size, classes), dtype=np.float32)
    for x in range(size):
        for y in range(size):
            for z in range(size):
                segmentation[x, y, z] = (x + 2 * y + 3 * z) % classes
                for k in range(classes):
                    numerator = (x * 5 + y * 7 + z * 11 + k * 3) % 29 + 1
                    memberships[x, y, z, k] = np.float32(numerator) / np.float32(30.0)
    return segmentation, memberships


def test_initialized_full_thinning_pass_matches_java_source_harness() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source"] == "VandyCRUISE MultipleToadSegmentation.initThinningLabels + propagateThinning"
    assert reference["source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert reference["lut_sha256"] == "50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f"
    assert reference["heap_sha256"] == "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
    assert reference["shape"] == [12, 12, 12, 4]
    assert len(reference["processed"]) > 1

    segmentation, memberships = _phantom()
    initialized = toads_initialize_thinning_queue(segmentation, memberships)
    np.testing.assert_array_equal(
        initialized.class_volumes.view(np.int32),
        np.asarray(reference["initial_volumes_bits"], dtype=np.int32),
    )
    assert np.all(initialized.class_volumes > np.float32(100.0))
    result = toads_propagate_thinning(
        segmentation,
        memberships,
        relation_matrix=np.zeros((4, 4), dtype=np.bool_),
        transformed_shape_values=np.zeros_like(memberships),
        volume_ratio=_from_bits(reference["volume_ratio_bits"]),
        limit=_from_bits(reference["limit_bits"]),
    )

    np.testing.assert_array_equal(
        result.class_volumes.view(np.int32), np.asarray(reference["final_volumes_bits"], dtype=np.int32)
    )
    np.testing.assert_array_equal(
        result.minimum_class_volumes.view(np.int32), np.asarray(reference["minimum_volumes_bits"], dtype=np.int32)
    )
    assert result.boundary_counts.tolist() == reference["boundary_counts"]
    processed = [[_bits(item.value), item.x, item.y, item.z, item.label] for item in result.processed_items]
    assert processed == reference["processed"]
    assert "".join(str(int(v)) for v in result.segmentation.ravel()) == reference["segmentation"]
    assert "".join("a" if v < 0 else str(int(v)) for v in result.skeleton.ravel()) == reference["skeleton"]
    assert "".join("1" if v else "0" for v in result.available.ravel()) == reference["available"]
    queue_after = []
    while result.queue:
        item = result.queue.pop_first()
        queue_after.append([_bits(item.value), item.x, item.y, item.z, item.label])
    assert queue_after == reference["queue_after"]

    assert np.all(result.class_volumes <= result.minimum_class_volumes)
    assert np.any(result.class_volumes < result.minimum_class_volumes)
    assert np.any(result.minimum_class_volumes == np.float32(100.0))
    assert np.any(result.minimum_class_volumes > np.float32(100.0))

    compiled = toads_propagate_thinning(
        segmentation,
        memberships,
        relation_matrix=np.zeros((4, 4), dtype=np.bool_),
        transformed_shape_values=np.zeros_like(memberships),
        volume_ratio=_from_bits(reference["volume_ratio_bits"]),
        limit=_from_bits(reference["limit_bits"]),
        record_processed_items=False,
    )
    assert compiled.processed_items == ()
    assert compiled.processed_count == len(reference["processed"])
    np.testing.assert_array_equal(compiled.segmentation, result.segmentation)
    np.testing.assert_array_equal(compiled.skeleton, result.skeleton)
    np.testing.assert_array_equal(compiled.available, result.available)
    np.testing.assert_array_equal(compiled.class_volumes, result.class_volumes)
    np.testing.assert_array_equal(compiled.ordering, result.ordering)
    compiled_queue_after = []
    while compiled.queue:
        item = compiled.queue.pop_first()
        compiled_queue_after.append([_bits(item.value), item.x, item.y, item.z, item.label])
    assert compiled_queue_after == reference["queue_after"]
