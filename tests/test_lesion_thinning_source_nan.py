from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from pymacruise.lesion_growth import initialize_lesion_toads_thinning_frontier
from pymacruise.lesion_growth import propagate_lesion_toads_thinning


ROOT = Path(__file__).parent
JAVA = ROOT / "java_reference" / "LesionThinningSourceNanHarness.java"
HEAP_JAVA = ROOT / "java_reference/upstream/edu/jhmi/rad/medic/structures/BinaryHeap4D.java"
NUMERICS_JAVA = ROOT / "java_reference/upstream/edu/jhmi/rad/medic/utilities/Numerics.java"
LESION_SOURCE_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
HEAP_SOURCE_SHA = "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
NUMERICS_SOURCE_SHA = "13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb"
HARNESS_SHA = "a27368df3269e251219dccdb571110412248275c1d787b8ffb2591ed4e495cd4"
TYPES = tuple("lesion" if index == 10 else "object" for index in range(13))


def _inputs():
    """Build a standard-13 class grid with Java-style absent-class NaNs."""
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
                border = x in (0, shape[0] - 1) or y in (0, shape[1] - 1) or z in (0, shape[2] - 1)
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


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def test_selected_thinning_initializer_consumes_java_nan_memberships(tmp_path: Path) -> None:
    assert LESION_SOURCE_SHA in JAVA.read_text(encoding="utf-8")
    assert hashlib.sha256(HEAP_JAVA.read_bytes()).hexdigest() == HEAP_SOURCE_SHA
    assert hashlib.sha256(NUMERICS_JAVA.read_bytes()).hexdigest() == NUMERICS_SOURCE_SHA
    assert hashlib.sha256(JAVA.read_bytes()).hexdigest() == HARNESS_SHA
    javac, java = shutil.which("javac"), shutil.which("java")
    assert javac is not None and java is not None, "a JDK is required for pinned Java references"
    segmentation, memberships, lesion, black_hole, skeleton, available = _inputs()
    kwargs = dict(
        lesion_class=10,
        object_types=TYPES,
    )
    with pytest.raises(ValueError, match="only finite"):
        initialize_lesion_toads_thinning_frontier(
            segmentation, memberships, lesion_membership=lesion,
            black_hole_membership=black_hole, initial_skeleton=skeleton,
            initial_available=available, **kwargs,
        )

    result = initialize_lesion_toads_thinning_frontier(
        segmentation, memberships, lesion_membership=lesion,
        black_hole_membership=black_hole, initial_skeleton=skeleton,
        initial_available=available, source_faithful_nan=True, **kwargs,
    )
    assert np.isnan(result.class_volume[10])
    assert result.class_volume[2] == 0.0
    assert any(np.isnan(item.value) for item in result.inserted_items)

    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run(
        [javac, "-d", str(classes), str(HEAP_JAVA), str(NUMERICS_JAVA), str(JAVA)],
        check=True, capture_output=True, text=True,
    )
    completed = subprocess.run(
        [java, "-cp", str(classes), "LesionThinningSourceNanHarness"],
        check=True, capture_output=True, text=True,
    )
    expected = json.loads(completed.stdout)

    assert _bits(result.class_volume) == expected["volume_bits"]
    assert result.boundary_counts.tolist() == expected["boundary_counts"]
    assert result.skeleton.ravel(order="C").astype(int).tolist() == expected["skeleton"]
    assert result.available.ravel(order="C").astype(int).tolist() == expected["available"]
    inserted_bits = [int(value) for value in _bits(np.asarray([item.value for item in result.inserted_items]))]
    inserted_payloads = [part for item in result.inserted_items for part in (item.x, item.y, item.z, item.label)]
    assert inserted_bits == [value & 0xFFFFFFFF for value in expected["inserted_priority_bits"]]
    assert inserted_payloads == expected["inserted_payloads"]

    removed = [result.queue.pop_first() for _ in range(len(result.queue))]
    assert [int(value) for value in _bits(np.asarray([item.value for item in removed]))] == [
        value & 0xFFFFFFFF for value in expected["removed_priority_bits"]
    ]
    assert [part for item in removed for part in (item.x, item.y, item.z, item.label)] == expected["removed_payloads"]


def test_selected_nan_bearing_thinning_pass_stops_after_popping_nan_priority(tmp_path: Path) -> None:
    """Trace only the source branch whose volume gate skips every topology edit."""
    javac, java = shutil.which("javac"), shutil.which("java")
    assert javac is not None and java is not None, "a JDK is required for pinned Java references"
    segmentation, memberships, lesion, black_hole, skeleton, available = _inputs()
    shapes = np.zeros((*segmentation.shape, memberships.shape[-1]), dtype=np.float32)
    ordering = np.zeros(segmentation.shape, dtype=np.float32)
    pass_kwargs = dict(
        segmentation=segmentation,
        memberships=memberships,
        object_types=TYPES,
        lesion_membership=lesion,
        black_hole_membership=black_hole,
        transformed_shape=shapes,
        initial_skeleton=skeleton,
        initial_available=available,
        initial_ordering=ordering,
        lesion_class=10,
        volume_ratio=0.4,
        limit=1000.0,
    )
    with pytest.raises(ValueError, match="only finite"):
        propagate_lesion_toads_thinning(**pass_kwargs)
    result = propagate_lesion_toads_thinning(**pass_kwargs, source_faithful_nan=True)

    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run(
        [javac, "-d", str(classes), str(HEAP_JAVA), str(NUMERICS_JAVA), str(JAVA)],
        check=True, capture_output=True, text=True,
    )
    completed = subprocess.run(
        [java, "-cp", str(classes), "LesionThinningSourceNanHarness", "pass"],
        check=True, capture_output=True, text=True,
    )
    expected = json.loads(completed.stdout)
    assert np.isnan(result.minimum_class_volume[10])
    assert _bits(result.minimum_class_volume) == [
        value & 0xFFFFFFFF for value in expected["minimum_bits"]
    ]
    assert [int(value) for value in _bits(np.asarray([event.popped_item.value for event in result.run.events]))] == [
        value & 0xFFFFFFFF for value in expected["pop_bits"]
    ]
    assert [
        part
        for event in result.run.events
        for part in (event.popped_item.x, event.popped_item.y, event.popped_item.z, event.popped_item.label)
    ] == expected["pop_payloads"]
    assert [event.processed for event in result.run.events] == expected["processed"]
    assert all(event.candidate_label is None for event in result.run.events)
    assert all(not event.inserted_neighbors for event in result.run.events)
    remaining = []
    while result.run.queue:
        item = result.run.queue.pop_first()
        remaining.append(item)
    assert [int(value) for value in _bits(np.asarray([item.value for item in remaining]))] == [
        value & 0xFFFFFFFF for value in expected["remaining_bits"]
    ]
    assert [
        part for item in remaining for part in (item.x, item.y, item.z, item.label)
    ] == expected["remaining_payloads"]
    np.testing.assert_array_equal(result.run.segmentation, segmentation)
    np.testing.assert_array_equal(result.run.skeleton, result.initialization.skeleton)
    np.testing.assert_array_equal(result.run.available, result.initialization.available)
    np.testing.assert_array_equal(result.run.ordering, ordering)
    assert np.isnan(result.run.class_volume[10])


def test_source_faithful_nan_mode_rejects_infinity_and_nonboolean_flag() -> None:
    segmentation, memberships, lesion, black_hole, skeleton, available = _inputs()
    memberships[1, 1, 1, 2] = np.float32(np.inf)
    kwargs = dict(
        segmentation=segmentation,
        memberships=memberships,
        lesion_membership=lesion,
        black_hole_membership=black_hole,
        initial_skeleton=skeleton,
        initial_available=available,
        lesion_class=10,
        object_types=TYPES,
        source_faithful_nan=True,
    )
    with pytest.raises(ValueError, match="not infinities"):
        initialize_lesion_toads_thinning_frontier(**kwargs)
    memberships[1, 1, 1, 2] = np.float32(np.nan)
    kwargs["source_faithful_nan"] = 1
    with pytest.raises(TypeError, match="must be boolean"):
        initialize_lesion_toads_thinning_frontier(**kwargs)
