from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.toads_approximate_thinning_startup import initialize_toads_approximate_thinning


ROOT = Path(__file__).parent
REFERENCE = ROOT / "data/synthetic/toads_approximate_thinning_startup/java_reference.json"
HARNESS = ROOT / "java_reference/ToadsApproximateThinningStartupReferenceHarness.java"
SOURCE = ROOT.parent.parent / "toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java"
SOURCE_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
HARNESS_SHA256 = "6369e2a3289e0295804f512c7728327345f37660a0078517b26b60a7897f796a"


def _fixture_arrays() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n, classes, slots = 5, 3, 2
    segmentation = np.empty((n, n, n), dtype=np.int8)
    best = np.empty((n, n, n, slots), dtype=np.int8)
    mems = np.empty((n, n, n, slots), dtype=np.float32)
    skeleton = np.full((n, n, n), -9, dtype=np.int16)
    available = np.zeros((n, n, n, classes), dtype=bool)
    for x in range(n):
        for y in range(n):
            for z in range(n):
                segmentation[x, y, z] = (x + 2 * y + z) % classes
                for slot in range(slots):
                    best[x, y, z, slot] = (x + y + z + slot) % classes
                    numerator = (x * 17 + y * 11 + z * 7 + slot * 13) % 79 + 1
                    mems[x, y, z, slot] = np.float32(numerator) / np.float32(100.0)
                if x in (0, n - 1) or y in (0, n - 1) or z in (0, n - 1):
                    available[x, y, z, 1] = True
    return segmentation, mems, best, skeleton, available


def _float_bits(values: np.ndarray) -> list[int]:
    return values.astype(np.float32, copy=False).view(np.int32).astype(int).tolist()


def test_sparse_approximate_thinning_startup_matches_java_method_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["shape"] == [5, 5, 5]
    assert reference["classes"] == 3 and reference["slots"] == 2
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256

    segmentation, mems, best, skeleton, available = _fixture_arrays()
    result = initialize_toads_approximate_thinning(
        segmentation,
        mems,
        best,
        class_count=3,
        initial_skeleton=skeleton,
        initial_available=available,
    )

    assert result.queue.mode == "max"
    assert _float_bits(result.class_volumes) == reference["volumes_bits"]
    assert result.boundary_counts.tolist() == reference["boundary_counts"]
    assert "".join(str(int(value)) for value in result.skeleton.ravel()) == reference["skeleton"]
    assert "".join("1" if value else "0" for value in result.available.ravel()) == reference["available"]

    observed = []
    while result.queue:
        item = result.queue.pop_first()
        observed.append([int(np.asarray(item.value, dtype=np.float32).view(np.int32)), item.x, item.y, item.z, item.label])
    # Compare the full source insertion set and exact float32 priorities; heap
    # tie order is checked by the shared heap tests, not duplicated here.
    assert sorted(observed) == sorted(reference["priority_events"])


def test_sparse_approximate_thinning_rejects_repeated_class_slots() -> None:
    segmentation, mems, best, _, _ = _fixture_arrays()
    best[2, 2, 2, 1] = best[2, 2, 2, 0]
    try:
        initialize_toads_approximate_thinning(segmentation, mems, best, class_count=3)
    except ValueError as error:
        assert "unique per voxel" in str(error)
    else:
        raise AssertionError("duplicate sparse labels must be rejected")
