from __future__ import annotations

import hashlib
import json
import struct
from pathlib import Path

import numpy as np

from pymacruise.toads_approximate_thinning_pass import run_approximate_thinning_pass


ROOT = Path(__file__).parent
SOURCE = ROOT.parent.parent / "toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java"
HARNESS = ROOT / "java_reference/ToadsApproximateThinningPassReferenceHarness.java"
REFERENCE = ROOT / "data/synthetic/toads_approximate_thinning_pass/java_reference.json"
SOURCE_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
HARNESS_SHA256 = "e15f0898cffb762632ce162b6a2bb1230deff8a5dea9018afac0556aa5536ec6"


def _folded_case() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n, classes, slots = 9, 3, 2
    labels = np.empty((n, n, n), dtype=np.int8)
    candidates = np.empty((n, n, n, slots), dtype=np.int8)
    memberships = np.empty((n, n, n, slots), dtype=np.float32)
    for x in range(n):
        for y in range(n):
            for z in range(n):
                folded = (x * 3 + y * 2) % 3 - 1
                dx, dy, dz = x - 4, y - 4, z - 4
                radius = 2 + folded
                radius2 = dx * dx + dy * dy + dz * dz
                if radius2 < radius * radius:
                    current = 1
                else:
                    fold_z = 4 + (x + 2 * y) % 3 - 1
                    current = 2 if abs(z - fold_z) <= 1 and abs(y - 4) <= 3 and abs(x - 4) <= 3 else 0
                labels[x, y, z] = current
                candidates[x, y, z, 0] = current
                candidates[x, y, z, 1] = (current + 1) % classes
                for slot in range(slots):
                    numerator = (x * 17 + y * 11 + z * 7 + slot * 13) % 39 + 1
                    memberships[x, y, z, slot] = np.float32(numerator) / np.float32(100.0)
    return labels, memberships, candidates


def test_approximate_thinning_queue_pass_matches_java_event_trace_and_terminal_state() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    labels, memberships, candidates = _folded_case()
    result = run_approximate_thinning_pass(
        labels,
        memberships,
        candidates,
        np.zeros((3, 3), dtype=bool),
        np.zeros(labels.shape + (3,), dtype=np.float32),
        class_count=3,
        limit=1000.0,
        explicit_minimum_class_volumes=np.zeros(3, dtype=np.float32),
        record_processed_items=True,
    )
    trace = hashlib.sha256()
    for item in result.processed_items:
        trace.update(struct.pack(">iiiii", int(item.value.view(np.int32)), item.x, item.y, item.z, item.label))
    assert result.processed_count == reference["processed_count"] == len(result.processed_items)
    assert trace.hexdigest() == reference["event_trace_sha256"]
    assert len(result.queue) == reference["remaining_queue"] == 0
    assert result.stop_reason == reference["stop_reason"] == "queue_empty"
    assert result.class_volumes.view(np.int32).astype(int).tolist() == reference["class_volume_bits"]
    assert np.count_nonzero(result.skeleton[1:-1, 1:-1, 1:-1] >= 0) == reference["skeleton_voxels_remaining"] == 0
    np.testing.assert_array_equal(result.segmentation, labels)


def test_approximate_thinning_pass_applies_source_default_volume_floor() -> None:
    labels, memberships, candidates = _folded_case()
    result = run_approximate_thinning_pass(
        labels,
        memberships,
        candidates,
        np.zeros((3, 3), dtype=bool),
        np.zeros(labels.shape + (3,), dtype=np.float32),
        class_count=3,
        limit=1000.0,
    )
    assert np.all(result.minimum_class_volumes == np.float32(100.0))
    assert result.processed_count == int(result.boundary_counts.sum())
    assert np.count_nonzero(result.skeleton[1:-1, 1:-1, 1:-1] >= 0) == 343
    assert np.all(result.class_volumes > 0.0)


def test_approximate_queue_validates_full_volume_inputs_once(monkeypatch) -> None:
    labels, memberships, candidates = _folded_case()
    shapes = np.zeros(labels.shape + (3,), dtype=np.float32)
    original = np.isfinite
    whole_volume_checks = {memberships.shape: 0, shapes.shape: 0}

    def counted_isfinite(values, *args, **kwargs):
        shape = np.shape(values)
        if shape in whole_volume_checks:
            whole_volume_checks[shape] += 1
        return original(values, *args, **kwargs)

    monkeypatch.setattr(np, "isfinite", counted_isfinite)
    run_approximate_thinning_pass(
        labels,
        memberships,
        candidates,
        np.zeros((3, 3), dtype=bool),
        shapes,
        class_count=3,
        limit=1000.0,
        explicit_minimum_class_volumes=np.zeros(3, dtype=np.float32),
    )
    assert whole_volume_checks == {memberships.shape: 1, shapes.shape: 1}
