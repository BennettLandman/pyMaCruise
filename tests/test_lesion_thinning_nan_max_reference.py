from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np

from pymacruise.lesion_growth import propagate_lesion_toads_thinning


ROOT = Path(__file__).parent
JAVA = ROOT / "java_reference" / "LesionThinningNanMaxReferenceHarness.java"
HEAP_JAVA = ROOT / "java_reference/upstream/edu/jhmi/rad/medic/structures/BinaryHeap4D.java"
NUMERICS_JAVA = ROOT / "java_reference/upstream/edu/jhmi/rad/medic/utilities/Numerics.java"
LESION_SOURCE_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
HEAP_SOURCE_SHA = "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
NUMERICS_SOURCE_SHA = "13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb"
HARNESS_SHA = "4cf1774116345c406be30393c9e668cfbcc9361c4e65d74f03e804a819df4f79"
TYPES = tuple("lesion" if index == 10 else "object" for index in range(13))


def _inputs():
    shape = (5, 5, 5)
    segmentation = np.empty(shape, dtype=np.int8)
    memberships = np.empty((*shape, 13), dtype=np.float32)
    lesion = np.zeros(shape, dtype=np.float32)
    black_hole = np.zeros(shape, dtype=np.float32)
    skeleton = np.empty(shape, dtype=np.int8)
    available = np.empty((*shape, 13), dtype=np.bool_)
    for x in range(shape[0]):
        for y in range(shape[1]):
            for z in range(shape[2]):
                segmentation[x, y, z] = 0 if x < 2 else 10
                border = x in (0, 4) or y in (0, 4) or z in (0, 4)
                skeleton[x, y, z] = 7 if border else -1
                for k in range(13):
                    available[x, y, z, k] = (x + y + z + k) % 2 == 0
                    if k not in (0, 10) or (x, y, z, k) == (2, 2, 2, 10):
                        memberships[x, y, z, k] = np.float32(np.nan)
                    elif k == 0:
                        memberships[x, y, z, k] = np.float32((8 + (x + 2 * y + 3 * z) % 5) / 32.0)
                    else:
                        memberships[x, y, z, k] = np.float32((18 + (2 * x + y + z) % 7) / 32.0)
    return segmentation, memberships, lesion, black_hole, skeleton, available


def _bits(values) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def _expected_bits(values) -> np.ndarray:
    return np.asarray(values, dtype=np.int64) & 0xFFFFFFFF


def _matches_bits(actual, expected) -> bool:
    return _bits(actual) == _expected_bits(expected).astype(np.uint64).tolist()


def _compile_and_run(tmp_path: Path) -> dict:
    assert LESION_SOURCE_SHA in JAVA.read_text(encoding="utf-8")
    assert hashlib.sha256(JAVA.read_bytes()).hexdigest() == HARNESS_SHA
    assert hashlib.sha256(HEAP_JAVA.read_bytes()).hexdigest() == HEAP_SOURCE_SHA
    assert hashlib.sha256(NUMERICS_JAVA.read_bytes()).hexdigest() == NUMERICS_SOURCE_SHA
    javac, java = shutil.which("javac"), shutil.which("java")
    assert javac is not None and java is not None, "a JDK is required for pinned Java references"
    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run(
        [javac, "-d", str(classes), str(HEAP_JAVA), str(NUMERICS_JAVA), str(JAVA)],
        check=True,
        capture_output=True,
        text=True,
    )
    completed = subprocess.run(
        [java, "-cp", str(classes), "LesionThinningNanMaxReferenceHarness"],
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(completed.stdout)


def test_use_lesion_max_masks_nan_before_thinning_minimum_and_queue_gates(tmp_path: Path) -> None:
    segmentation, memberships, lesion, black_hole, skeleton, available = _inputs()
    shapes = np.zeros((*segmentation.shape, memberships.shape[-1]), dtype=np.float32)
    ordering = np.zeros(segmentation.shape, dtype=np.float32)
    original = tuple(
        array.copy()
        for array in (segmentation, memberships, lesion, black_hole, skeleton, available, shapes, ordering)
    )
    result = propagate_lesion_toads_thinning(
        segmentation,
        memberships,
        TYPES,
        lesion,
        black_hole,
        shapes,
        skeleton,
        available,
        ordering,
        lesion_class=10,
        volume_ratio=0.4,
        limit=1000.0,
        use_lesion_max=True,
        source_faithful_nan=True,
    )
    expected = _compile_and_run(tmp_path)

    # Java Numerics.max(NaN, 0, 0) returns the finite inner max. The single
    # NaN therefore contributes zero to class 10's volume instead of poisoning
    # the minimum and keeping a NaN queue priority.
    assert np.isnan(memberships[2, 2, 2, 10])
    assert np.isfinite(result.initialization.class_volume).all()
    assert result.initialization.class_volume[10] < 100.0
    assert _matches_bits(result.initialization.class_volume, expected["volume_bits"])
    assert _matches_bits(result.minimum_class_volume, expected["minimum_bits"])
    assert result.minimum_class_volume.tolist() == [100.0] * 13
    assert result.initialization.boundary_counts.tolist() == expected["boundary_counts"]
    assert result.initialization.skeleton.ravel(order="C").astype(int).tolist() == expected["skeleton"]
    assert result.initialization.available.ravel(order="C").astype(int).tolist() == expected["available"]

    inserted = result.initialization.inserted_items
    assert _matches_bits([item.value for item in inserted], expected["inserted_bits"])
    assert [part for item in inserted for part in (item.x, item.y, item.z, item.label)] == expected[
        "inserted_payloads"
    ]
    events = result.run.events
    assert _matches_bits([event.popped_item.value for event in events], expected["pop_bits"])
    assert [
        part
        for event in events
        for part in (event.popped_item.x, event.popped_item.y, event.popped_item.z, event.popped_item.label)
    ] == expected["pop_payloads"]
    assert [int(event.processed) for event in events] == expected["processed"]
    assert all(event.candidate_label is None for event in events)
    assert all(not event.inserted_neighbors for event in events)

    remaining = []
    while result.run.queue:
        item = result.run.queue.pop_first()
        remaining.append(item)
    assert _matches_bits([item.value for item in remaining], expected["remaining_bits"])
    assert [part for item in remaining for part in (item.x, item.y, item.z, item.label)] == expected[
        "remaining_payloads"
    ]
    np.testing.assert_array_equal(result.run.segmentation, segmentation)
    np.testing.assert_array_equal(result.run.skeleton, result.initialization.skeleton)
    np.testing.assert_array_equal(result.run.available, result.initialization.available)
    np.testing.assert_array_equal(result.run.ordering, ordering)
    for actual, before in zip(
        (segmentation, memberships, lesion, black_hole, skeleton, available, shapes, ordering),
        original,
        strict=True,
    ):
        np.testing.assert_array_equal(actual, before)
