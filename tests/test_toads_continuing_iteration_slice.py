"""Exact field snapshots for a bounded folded-phantom TOADS iteration slice.

Individual source kernels have their own Java comparisons. This integration
snapshot checks their ordering and state handoff together; it is not a claim
that this synthetic composition is a complete TOADS run.
"""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_iteration_slice import toads_continuing_iteration_slice
from pymacruise.toads_relations import toads_closest_labels


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_continuing_folded_slice"
REFERENCE = ROOT / "composition_snapshot.json"


def _float_bits(values: np.ndarray) -> list[int]:
    return values.astype(np.float32, copy=False).ravel().view(np.uint32).astype(np.uint64).tolist()


def _load(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / name))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    return np.asarray(image.dataobj)


def _run():
    segmentation = _load("segmentation.nii.gz").astype(np.int16)
    images = _load("image_values.nii.gz").astype(np.float32)
    memberships = _load("initial_memberships.nii.gz").astype(np.float32)
    shapes = _load("transformed_shapes.nii.gz").astype(np.float32)
    previous_closest = toads_closest_labels(segmentation, 1)
    nx, ny, nz = segmentation.shape
    classes = memberships.shape[-1]
    ordering = np.empty(segmentation.shape, dtype=np.float32)
    boundary = np.zeros((*segmentation.shape, classes), dtype=np.int8)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                ordering[x, y, z] = np.float32(-0.01 * (x * 7 + y * 5 + z * 3))

    return toads_continuing_iteration_slice(
        segmentation,
        memberships,
        np.ones((classes, classes), dtype=np.bool_),
        shapes,
        0.1,
        0.75,
        images,
        np.asarray([[0.16, 0.48, 0.78]], dtype=np.float32),
        np.ones(1, dtype=np.float32),
        np.ones((1, classes), dtype=np.float32),
        np.asarray([0.08], dtype=np.float32),
        shapes,
        previous_closest,
        np.asarray([0.25, 0.3, 0.35], dtype=np.float32),
        np.asarray(
            [[0.0, 0.2, 0.1], [0.2, 0.0, 0.15], [0.1, 0.15, 0.0]],
            dtype=np.float32,
        ),
        np.asarray([0.2, 0.3, 0.25], dtype=np.float32),
        0.2,
        0.6,
        0.1,
        ordering,
        boundary,
        1.0,
        1,
        4,
        1e-4,
        relation_limit=1,
    )


def _actual_snapshot(result) -> dict[str, object]:
    return {
        "stage_order": ["thinning", "object_membership", "decision", "growing", "relations"],
        "thinning_segmentation": result.thinning.segmentation.ravel().astype(int).tolist(),
        "thinning_skeleton": result.thinning.skeleton.ravel().astype(int).tolist(),
        "thinning_available": result.thinning.available.ravel().astype(int).tolist(),
        "thinning_volumes_bits": _float_bits(result.thinning.class_volumes),
        "thinning_processed": len(result.thinning.processed_items),
        "memberships_bits": _float_bits(result.memberships),
        "membership_energy_bits": _float_bits(np.asarray([result.membership_energy])),
        "control_ratio_bits": _float_bits(np.asarray([result.control.energy_ratio])),
        "control_flags": [
            result.control.stop,
            result.control.run_growing,
            result.control.update_classification,
            result.control.refresh_relations,
        ],
        "growing_segmentation": result.growing.segmentation.ravel().astype(int).tolist(),
        "growing_skeleton": result.growing.skeleton.ravel().astype(int).tolist(),
        "growing_available": result.growing.available.ravel().astype(int).tolist(),
        "growing_ordering_bits": _float_bits(result.growing.ordering),
        "growing_volumes_bits": _float_bits(result.growing.class_volumes),
        "growing_processed": len(result.growing.events),
        "relation_map": result.relations.relation.relation_map.ravel().astype(int).tolist(),
        "relation_membership": result.relations.relation.relation_membership.ravel().astype(int).tolist(),
        "boundary_distances": result.relations.boundary_distances.ravel().astype(int).tolist(),
        "closest_labels": result.relations.closest_labels.ravel().astype(int).tolist(),
        "relation_passes": [
            result.relations.relation.diffusion_iterations,
            result.relations.boundary_diffusion_passes,
            result.relations.closest_label_passes,
        ],
    }


def test_folded_asymmetric_continuing_slice_matches_fieldwise_snapshot() -> None:
    reference = json.loads(REFERENCE.read_text())
    result = _run()
    actual = _actual_snapshot(result)
    assert actual == reference
    assert result.control.stop is False
    assert result.control.run_growing is True
    assert len(result.thinning.processed_items) > 0
    assert len(result.growing.events) > 0
