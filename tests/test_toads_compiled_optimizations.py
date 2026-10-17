"""Exact regression of accelerated TOADS against the retained reference paths."""
from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_inhomogeneity_exact_fit import fit_toads_exact_chebyshev_field
from pymacruise.toads_relations import (
    toads_closest_labels, toads_compute_relations, toads_update_boundary_distances,
)
from pymacruise.toads_thinning import toads_propagate_growing


def _exact(a, b):
    assert a.shape == b.shape and a.dtype == b.dtype
    assert a.tobytes() == b.tobytes()


@pytest.mark.parametrize("degree", [1, 2, 3])
@pytest.mark.parametrize("strided", [False, True])
def test_field_kernels_preserve_every_intermediate_and_strided_field_mutation(degree, strided):
    fixture = Path(__file__).parent / "data/synthetic/toads_exact_polynomial_field_degrees23"
    image, mem, seg, field = (
        np.asarray(nib.load(fixture / f"{name}.nii.gz").dataobj)
        for name in ("image", "memberships", "segmentation", "initial_field")
    )
    centers = np.asarray(json.loads((fixture / "centroids.json").read_text()), dtype=np.float32)
    if strided:
        image, mem, seg, field = (a[::-1] for a in (image, mem, seg, field))
    reference = fit_toads_exact_chebyshev_field(image, mem, seg, centers, field.copy(), degree=degree, backend="python")
    # Channel-last fields are strided in the multi-channel production path.
    storage = np.stack((field, field), axis=-1)
    actual = fit_toads_exact_chebyshev_field(image, mem, seg, centers, storage[..., 0], degree=degree)
    for name in ("field", "coefficients", "normal_matrix", "right_hand_side", "correction_mask"):
        _exact(getattr(actual, name), getattr(reference, name))
    _exact(storage[..., 0], reference.field)
    _exact(storage[..., 1], field)


@pytest.mark.parametrize("limit", [0, 1, 3])
@pytest.mark.parametrize("kind", ["blocks", "random", "uniform"])
def test_relation_kernels_preserve_diffusion_order_borders_and_carried_bytes(limit, kind):
    rng = np.random.default_rng(6934)
    shape = (11, 10, 9)
    xyz = np.indices(shape)
    if kind == "random":
        labels = rng.integers(0, 4, shape, dtype=np.int16)
    elif kind == "blocks":
        labels = ((xyz[0] // 4 + xyz[1] // 4 + xyz[2] // 4) % 4).astype(np.int16)
    else:
        labels = np.zeros(shape, dtype=np.int16)
    mem = rng.random((*shape, 4), dtype=np.float32)
    mem[::2, ::2, ::2] = 0  # Ties and zero-support candidates.
    boundary = rng.integers(-128, 128, (*shape, 4), dtype=np.int8)
    original = boundary.copy()
    a = toads_compute_relations(labels, mem, limit)
    b = toads_compute_relations(labels, mem, limit, backend="python")
    _exact(a.relation_map, b.relation_map)
    _exact(a.relation_membership, b.relation_membership)
    assert a.diffusion_iterations == b.diffusion_iterations
    _exact(toads_closest_labels(labels, limit), toads_closest_labels(labels, limit, backend="python"))
    _exact(toads_update_boundary_distances(labels, boundary, limit),
           toads_update_boundary_distances(labels, boundary, limit, backend="python"))
    _exact(boundary, original)


@pytest.mark.parametrize("tied", [False, True])
def test_compiled_growing_matches_full_event_reference_and_progress(tied):
    from test_toads_growing_pass import _inputs
    labels, skeleton, mem, ordering, volumes = _inputs()
    if tied:
        mem[:] = np.float32(.25)
        ordering[:] = 0
    available = np.zeros_like(mem, dtype=np.bool_)
    original = tuple(a.copy() for a in (labels, skeleton, mem, ordering, volumes, available))
    expected = toads_propagate_growing(labels, skeleton, mem, ordering, volumes,
                                      initial_available=available, backend="python")
    messages = []
    actual = toads_propagate_growing(labels, skeleton, mem, ordering, volumes,
                                    initial_available=available, record_events=False,
                                    progress_callback=messages.append, progress_interval=17)
    for name in ("segmentation", "skeleton", "available", "ordering", "class_volumes"):
        _exact(getattr(actual, name), getattr(expected, name))
    assert actual.events == actual.remaining_queue == ()
    assert messages[-1] == f"growing queue complete after {len(expected.events)} events"
    for before, after in zip(original, (labels, skeleton, mem, ordering, volumes, available), strict=True):
        _exact(before, after)


def test_growing_preserves_zero_threshold_nan_priority_behavior():
    from test_toads_growing_pass import _inputs
    labels, skeleton, mem, ordering, volumes = _inputs()
    mem[:] = 0
    with np.errstate(invalid="ignore", divide="ignore"):
        expected = toads_propagate_growing(labels, skeleton, mem, ordering, volumes,
                                          zero_threshold=0, record_events=False, backend="python")
        actual = toads_propagate_growing(labels, skeleton, mem, ordering, volumes,
                                        zero_threshold=0, record_events=False)
    for name in ("segmentation", "skeleton", "available", "ordering", "class_volumes"):
        _exact(getattr(actual, name), getattr(expected, name))
