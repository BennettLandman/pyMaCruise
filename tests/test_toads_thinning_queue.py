from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_thinning import toads_initialize_thinning_queue


REFERENCE = Path(__file__).parent / "data" / "synthetic" / "toads_thinning_queue" / "java_reference.json"


def _signed_bits(value: np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.int32))


def _phantom() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    nx = ny = nz = 7
    classes = 4
    labels = np.empty((nx, ny, nz), dtype=np.int16)
    memberships = np.empty((nx, ny, nz, classes), dtype=np.float32)
    available = np.zeros((nx, ny, nz, classes), dtype=np.bool_)
    skeleton = np.full((nx, ny, nz), -1, dtype=np.int16)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                labels[x, y, z] = (x + 2 * y + 3 * z) % classes
                for k in range(classes):
                    numerator = (x * 5 + y * 7 + z * 11 + k * 3) % 29 + 1
                    memberships[x, y, z, k] = np.float32(numerator) / np.float32(30.0)
                    border = x in (0, nx - 1) or y in (0, ny - 1) or z in (0, nz - 1)
                    available[x, y, z, k] = border and (x + y + z + k) % 2 == 0
    return labels, memberships, skeleton, available


def test_thinning_queue_initialization_matches_java_method_and_heap() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source"] == "VandyCRUISE MultipleToadSegmentation.initThinningLabels"
    assert reference["source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert reference["heap_source_sha256"] == "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
    labels, memberships, skeleton, available = _phantom()
    result = toads_initialize_thinning_queue(
        labels,
        memberships,
        initial_skeleton=skeleton,
        initial_available=available,
    )

    np.testing.assert_array_equal(
        result.class_volumes.view(np.int32), np.asarray(reference["volumes_bits"], dtype=np.int32)
    )
    assert result.boundary_counts.tolist() == reference["boundary_counts"]
    assert result.skeleton.ravel().tolist() == reference["skeleton"]
    assert result.available.ravel().tolist() == reference["available"]

    popped = []
    while result.queue:
        item = result.queue.pop_first()
        popped.append([_signed_bits(item.value), item.x, item.y, item.z, item.label])
    assert popped == reference["popped"]
    # The selected source's endpoint counter counts all 20 diagonals, so its
    # `Nngb == 1` branch cannot be reached for an interior neighborhood.
    assert not any(item[0] == _signed_bits(np.float32(-0.0)) for item in popped)

    # No source initialization writes the border state.
    for axis in range(3):
        np.testing.assert_array_equal(np.take(result.skeleton, [0, -1], axis=axis), -1)


def test_thinning_queue_rejects_bad_shapes_and_labels() -> None:
    labels, memberships, _, _ = _phantom()
    with pytest.raises(ValueError, match="at least 3"):
        toads_initialize_thinning_queue(labels[:, :, :2], memberships[:, :, :2])

    bad_labels = labels.astype(np.int16)
    bad_labels[2, 2, 2] = 4
    with pytest.raises(ValueError, match="valid class indices"):
        toads_initialize_thinning_queue(bad_labels, memberships)
