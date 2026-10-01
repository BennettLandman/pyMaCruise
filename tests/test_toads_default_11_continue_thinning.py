from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_default_continue_thinning import toads_default_continue_thinning
from pymacruise.toads_default_atlas_membership import DEFAULT_ATLAS_CLASS_TYPES
from pymacruise.toads_iteration import toads_continuing_iteration


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_default_11_continue_thinning"
REFERENCE = ROOT / "java_reference.json"


def _load(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / name))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    return np.asarray(image.dataobj)


def _signed_bits(value: np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.int32))


def _labels(values: np.ndarray, *, empty: bool = False) -> str:
    return "".join(
        "." if empty and int(value) == -1 else (str(int(value)) if int(value) < 10 else chr(ord("a") + int(value) - 10))
        for value in values.ravel()
    )


def test_selected_11_class_continue_growth_relation_and_next_thinning_match_java() -> None:
    ref = json.loads(REFERENCE.read_text())
    assert ref["source_order"] == [
        "continue.growing",
        "computeRelations.class_adjacency",
        "exact_memberships",
        "fixed_zero_transformed_shapes",
        "initThinningLabels",
        "propagateThinning",
    ]
    seg = _load("initial_segmentation.nii.gz").astype(np.int16)
    skeleton = _load("initial_skeleton.nii.gz").astype(np.int16)
    memberships = _load("exact_memberships.nii.gz").astype(np.float32)
    ordering = _load("growing_ordering.nii.gz").astype(np.float32)
    shapes = _load("transformed_shapes.nii.gz").astype(np.float32)
    assert seg.shape == tuple(ref["shape"][:3])
    assert memberships.shape == tuple(ref["shape"])
    assert len(np.unique(seg)) == 11
    assert DEFAULT_ATLAS_CLASS_TYPES == (
        "MASK",
        "OPTIMIZED_CSF",
        "OPTIMIZED_CSF",
        "OPTIMIZED_GM",
        "OPTIMIZED_GM",
        "OPTIMIZED_GM",
        "OPTIMIZED_GM",
        "OPTIMIZED_GM",
        "OPTIMIZED_WM",
        "OPTIMIZED_WM",
        "OPTIMIZED_WM",
    )
    assert np.allclose(memberships.sum(axis=-1), 1.0, rtol=0.0, atol=2e-7)
    loaded = [nib.load(str(ROOT / name)) for name in (
        "initial_segmentation.nii.gz",
        "initial_skeleton.nii.gz",
        "exact_memberships.nii.gz",
        "growing_ordering.nii.gz",
        "transformed_shapes.nii.gz",
    )]
    for image in loaded[1:]:
        np.testing.assert_array_equal(image.affine, loaded[0].affine)
    fixed_inputs = [seg.copy(), skeleton.copy(), memberships.copy(), ordering.copy(), shapes.copy()]

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
        np.zeros(11, dtype=np.float32),
        np.zeros((*seg.shape, 11), dtype=np.int8),
        relation_limit=2,
        do_growing=False,
    )
    assert not transition.control.stop
    assert len(transition.growing.events) == ref["growth_event_count"]
    assert _labels(transition.growing.segmentation) == ref["grown_segmentation"]
    assert _labels(transition.growing.skeleton, empty=True) == ref["grown_skeleton"]
    growth_events = [
        [
            _signed_bits(event.popped_item.value),
            event.popped_item.x,
            event.popped_item.y,
            event.popped_item.z,
            event.popped_item.label,
            -1 if event.topology_accepted is None else int(event.topology_accepted),
            [
                [
                    _signed_bits(item.value),
                    item.x,
                    item.y,
                    item.z,
                    item.label,
                ]
                for item in event.inserted_neighbors
            ],
        ]
        for event in transition.growing.events
    ]
    assert growth_events == ref["growth_events"]

    adjacency = transition.relation_state.relation.class_adjacency
    assert "".join("1" if value else "0" for value in adjacency.ravel()) == ref["relation_adjacency_bits"]
    thinning = toads_default_continue_thinning(transition, memberships, shapes)
    for actual_input, original_input in zip(
        (seg, skeleton, memberships, ordering, shapes), fixed_inputs, strict=True
    ):
        np.testing.assert_array_equal(actual_input, original_input)
    assert _signed_bits(np.float32(0.5)) == ref["volume_ratio_bits"]
    assert _signed_bits(np.float32(3.0)) == ref["limit_bits"]
    assert thinning.class_volumes.view(np.int32).astype(int).tolist() == ref["final_volumes_bits"]
    assert thinning.minimum_class_volumes.view(np.int32).astype(int).tolist() == ref["minimum_volumes_bits"]
    assert thinning.boundary_counts.astype(int).tolist() == ref["boundary_counts"]
    processed = [
        [_signed_bits(item.value), item.x, item.y, item.z, item.label]
        for item in thinning.processed_items
    ]
    assert processed == ref["processed"]
    assert _labels(thinning.segmentation) == ref["segmentation"]
    assert _labels(thinning.skeleton, empty=True) == ref["skeleton"]
    assert "".join("1" if v else "0" for v in thinning.available.ravel()) == ref["available"]
