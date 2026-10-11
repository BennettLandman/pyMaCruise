from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.toads_iteration import toads_continuing_iteration
from pymacruise.toads_default_continue_thinning import toads_default_continue_thinning
from pymacruise.toads_thinning import toads_propagate_thinning


ROOT = Path(__file__).parent / "data" / "synthetic"
GROWTH = ROOT / "toads_growing_pass" / "java_reference.json"
REFERENCE = ROOT / "toads_default_continue_thinning" / "java_reference.json"


def _inputs():
    size, classes = 9, 4
    segmentation = np.zeros((size, size, size), dtype=np.int16)
    skeleton = np.full_like(segmentation, -1)
    ordering = np.empty((size, size, size), dtype=np.float32)
    memberships = np.empty((size, size, size, classes), dtype=np.float32)
    shapes = np.zeros_like(memberships)
    for x in range(size):
        for y in range(size):
            for z in range(size):
                ordering[x, y, z] = np.float32(-((x * 13 + y * 7 + z * 3) % 19)) / np.float32(10.0)
                for k in range(classes):
                    memberships[x, y, z, k] = np.float32(
                        ((x * 5 + y * 7 + z * 11 + k * 3) % 29 + 1)
                    ) / np.float32(30.0)
                center = 2 + ((y + z) % 3)
                if 2 <= y <= 6 and 2 <= z <= 6 and abs(x - center) <= 1:
                    segmentation[x, y, z] = 1
                if 5 <= x <= 6 and 2 <= y <= 3 and 5 <= z <= 6:
                    segmentation[x, y, z] = 2
                if x == 1 and 2 <= y <= 6 and z == 7:
                    segmentation[x, y, z] = 3
    center = 2 + ((4 + 3) % 3)
    skeleton[center, 4, 3] = 1
    skeleton[5, 2, 5] = 2
    return segmentation, skeleton, memberships, ordering, shapes


def _bits(values: np.ndarray) -> list[int]:
    return values.astype(np.float32, copy=False).ravel().view(np.int32).astype(np.int64).tolist()


def test_next_thinning_pass_consumes_java_referenced_continuation_state() -> None:
    growth_ref = json.loads(GROWTH.read_text())
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_order"] == [
        "continue.growing.segmentation",
        "continue.relation_state.class_adjacency",
        "exact_memberships",
        "fixed_transformed_shapes",
        "initThinningLabels",
        "propagateThinning",
    ]
    seg, skeleton, memberships, ordering, shapes = _inputs()
    growth_state = growth_ref["relation_state_after_growth"]
    boundary = np.asarray(growth_state["initial_boundary_distances"], dtype=np.int8).reshape(
        (*seg.shape, memberships.shape[-1])
    )
    transition = toads_continuing_iteration(
        1.0,
        0.5,
        1,
        4,
        1e-4,
        seg,
        skeleton,
        memberships,
        ordering,
        np.asarray([1, 2, 3, 4], dtype=np.float32),
        boundary,
        relation_limit=2,
        do_growing=False,
    )
    assert "".join(str(int(v)) for v in transition.growing.segmentation.ravel()) == growth_ref["segmentation"]
    assert transition.relation_state.relation.relation_membership.astype(int).ravel().tolist() == growth_state["relation_membership"]
    assert "".join("1" if value else "0" for value in transition.relation_state.relation.class_adjacency.ravel()) == "".join(
        str(value) for value in reference["relation_adjacency_bits"] if value in "01"
    )
    assert not transition.control.stop  # selected continuation grows regardless of doGrowing

    labels_before = transition.growing.segmentation.copy()
    memberships_before = memberships.copy()
    shapes_before = shapes.copy()
    actual = toads_default_continue_thinning(transition, memberships, shapes)
    # The selected wrapper applies a 100-voxel floor; this chained fixture
    # reaches it for one class and exercises a topology-gated label change.
    assert _bits(actual.minimum_class_volumes) == reference["minimum_volumes_bits"]
    assert [int(value) for value in actual.boundary_counts] == reference["boundary_counts"]
    np.testing.assert_array_equal(transition.growing.segmentation, labels_before)
    np.testing.assert_array_equal(memberships, memberships_before)
    np.testing.assert_array_equal(shapes, shapes_before)

    direct = toads_propagate_thinning(
        transition.growing.segmentation,
        memberships,
        transition.relation_state.relation.class_adjacency,
        shapes,
        volume_ratio=0.5,
        limit=3.0,
        initial_skeleton=transition.growing.skeleton,
    )
    np.testing.assert_array_equal(actual.segmentation, direct.segmentation)
    np.testing.assert_array_equal(actual.skeleton, direct.skeleton)
    np.testing.assert_array_equal(actual.available, direct.available)
    np.testing.assert_array_equal(actual.class_volumes.view(np.uint32), direct.class_volumes.view(np.uint32))
    assert len(actual.processed_items) == reference["processed_count"]
    processed = [
        [int(np.asarray(item.value, dtype=np.float32).view(np.int32)), item.x, item.y, item.z, item.label]
        for item in actual.processed_items
    ]
    assert processed == reference["processed"]
    assert _bits(actual.class_volumes) == reference["final_volumes_bits"]
    assert "".join(str(int(v)) for v in actual.segmentation.ravel()) == reference["segmentation"]
    assert reference["segmentation"] != growth_ref["segmentation"]
    assert "".join("a" if v < 0 else str(int(v)) for v in actual.skeleton.ravel()) == reference["skeleton"]
    assert "".join("1" if v else "0" for v in actual.available.ravel()) == reference["available"]


def test_default_continue_thinning_rejects_wrong_shapes() -> None:
    seg, skeleton, memberships, ordering, _ = _inputs()
    boundary = np.zeros((*seg.shape, memberships.shape[-1]), dtype=np.int8)
    transition = toads_continuing_iteration(
        1.0, 0.5, 1, 4, 1e-4, seg, skeleton, memberships, ordering,
        np.asarray([1, 2, 3, 4], dtype=np.float32), boundary,
        relation_limit=1, do_growing=True,
    )
    with np.testing.assert_raises_regex(ValueError, "transformed_shape_values"):
        toads_default_continue_thinning(transition, memberships, np.zeros(seg.shape, dtype=np.float32))
