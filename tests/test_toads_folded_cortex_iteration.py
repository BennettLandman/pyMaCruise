"""Larger folded-cortex integration regression through one TOADS loop seam."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_iteration_slice import toads_continuing_iteration_slice
from pymacruise.toads_relations import toads_closest_labels
import pymacruise.toads_thinning as thinning_module


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_folded_cortex_iteration"
REFERENCE = json.loads((ROOT / "python_integration_reference.json").read_text())


def _load(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / name))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    return np.asarray(image.dataobj)


def _digest(values: np.ndarray) -> str:
    array = np.asarray(values)
    return hashlib.sha256(array.tobytes(order="C")).hexdigest()


def _run():
    segmentation = _load("segmentation.nii.gz").astype(np.int16)
    images = _load("image_values.nii.gz").astype(np.float32)
    memberships = _load("initial_memberships.nii.gz").astype(np.float32)
    shapes = _load("transformed_shapes.nii.gz").astype(np.float32)
    previous_closest = toads_closest_labels(segmentation, 2)
    nx, ny, nz = segmentation.shape
    ordering = np.empty(segmentation.shape, dtype=np.float32)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                ordering[x, y, z] = np.float32(-0.006 * (x * 7 + y * 5 + z * 3))

    result = toads_continuing_iteration_slice(
        segmentation,
        memberships,
        np.ones((3, 3), dtype=np.bool_),
        shapes,
        0.16,
        0.8,
        images,
        np.asarray([[0.28, 0.47, 0.69]], dtype=np.float32),
        np.ones(1, dtype=np.float32),
        np.ones((1, 3), dtype=np.float32),
        np.asarray([0.08], dtype=np.float32),
        shapes,
        previous_closest,
        np.asarray([0.25, 0.3, 0.35], dtype=np.float32),
        np.asarray([[0.0, 0.2, 0.1], [0.2, 0.0, 0.15], [0.1, 0.15, 0.0]], dtype=np.float32),
        np.asarray([0.2, 0.3, 0.25], dtype=np.float32),
        0.2,
        0.6,
        0.1,
        ordering,
        np.zeros((*segmentation.shape, 3), dtype=np.int8),
        1.0,
        1,
        4,
        1e-4,
        relation_limit=2,
    )
    return segmentation, result


def test_larger_folded_cortex_iteration_matches_committed_python_snapshot(monkeypatch) -> None:
    topology_counts = {"accepted": 0, "rejected": 0}
    source_gate = thinning_module.toads_homeomorphic_label_change

    def count_gate(neighborhood, candidate_label):
        accepted = source_gate(neighborhood, candidate_label)
        topology_counts["accepted" if accepted else "rejected"] += 1
        return accepted

    monkeypatch.setattr(thinning_module, "toads_homeomorphic_label_change", count_gate)
    segmentation, result = _run()
    snapshot = {
        "shape": list(segmentation.shape),
        "initial_label_counts": np.bincount(segmentation.ravel(), minlength=3).tolist(),
        "thinning_processed": len(result.thinning.processed_items),
        "topology_gate_calls": topology_counts,
        "thinning_segmentation_sha256": _digest(result.thinning.segmentation),
        "thinning_skeleton_sha256": _digest(result.thinning.skeleton),
        "membership_sha256": _digest(result.memberships),
        "membership_energy_float32_bits": int(np.asarray(result.membership_energy, dtype=np.float32).view(np.uint32)),
        "control_energy_ratio_float32_bits": int(np.asarray(result.control.energy_ratio, dtype=np.float32).view(np.uint32)),
        "control_flags": [
            result.control.stop,
            result.control.run_growing,
            result.control.update_classification,
            result.control.refresh_relations,
        ],
        "growing_processed": len(result.growing.events),
        "growing_segmentation_sha256": _digest(result.growing.segmentation),
        "growing_skeleton_sha256": _digest(result.growing.skeleton),
        "relation_map_sha256": _digest(result.relations.relation.relation_map),
        "boundary_distance_sha256": _digest(result.relations.boundary_distances),
        "closest_labels_sha256": _digest(result.relations.closest_labels),
    }
    assert snapshot == REFERENCE
    assert result.control.stop is False
    assert result.control.run_growing is True
    assert snapshot["thinning_processed"] > 500
    assert snapshot["growing_processed"] > 100
    assert topology_counts["accepted"] > 0
    assert topology_counts["rejected"] > 0
