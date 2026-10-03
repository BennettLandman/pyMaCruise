from __future__ import annotations

import numpy as np

from pymacruise.toads_heap import ToadsPriorityHeap4D
from pymacruise.toads_thinning import toads_initialize_thinning_queue


def _source_order_reference(
    labels: np.ndarray,
    memberships: np.ndarray,
    initial_skeleton: np.ndarray,
    initial_available: np.ndarray,
    zero_threshold: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[tuple[int, int, int, int, int]]]:
    """Literal, small-grid form of the Java initialization loops."""
    nx, ny, nz = labels.shape
    classes = memberships.shape[-1]
    skeleton = initial_skeleton.copy()
    available = initial_available.copy()
    volumes = np.zeros(classes, dtype=np.float32)
    boundaries = np.zeros(classes, dtype=np.int32)
    queue = ToadsPriorityHeap4D(nx * ny + ny * nz + nz * nx, "max")

    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                current = int(labels[x, y, z])
                skeleton[x, y, z] = np.int16(current)
                volumes[current] = np.float32(volumes[current] + memberships[x, y, z, current])
                for candidate in range(classes):
                    available[x, y, z, candidate] = candidate == current

    for x in range(1, nx - 1):
        for y in range(1, ny - 1):
            for z in range(1, nz - 1):
                candidate = int(labels[x, y, z])
                boundary = False
                neighbor_count = 0
                ni = nj = nk = 0
                for i in (-1, 0, 1):
                    for j in (-1, 0, 1):
                        for k in (-1, 0, 1):
                            distance2 = i * i + j * j + k * k
                            neighbor_label = int(labels[x + i, y + j, z + k])
                            if distance2 == 1 and neighbor_label != candidate:
                                boundary = True
                            elif distance2 > 0:
                                neighbor_count += 1
                                ni, nj, nk = x + i, y + j, z + k
                if not boundary:
                    continue
                available[x, y, z, candidate] = False
                endpoint = False
                if neighbor_count == 1:
                    neighbor_count = 0
                    for i in (-1, 0, 1):
                        for j in (-1, 0, 1):
                            for k in (-1, 0, 1):
                                if i * i + j * j + k * k > 0 and labels[ni + i, nj + j, nk + k] == candidate:
                                    neighbor_count += 1
                    endpoint = neighbor_count <= 2
                if endpoint:
                    order = np.float32(-0.0)
                else:
                    best = 0
                    maximum = np.float32(0.0)
                    for other in range(classes):
                        if memberships[x, y, z, other] > maximum:
                            best = other
                            maximum = np.float32(memberships[x, y, z, other])
                    if candidate == best:
                        ratio = np.float32(1.0)
                    else:
                        ratio = np.float32(
                            np.float32(memberships[x, y, z, best] + np.float32(zero_threshold))
                            / np.float32(memberships[x, y, z, candidate] + np.float32(zero_threshold))
                        )
                    order = np.float32(np.float32(-1.0) / ratio)
                queue.add(order, x, y, z, candidate)
                boundaries[candidate] += 1

    popped = []
    while queue:
        item = queue.pop_first()
        popped.append((int(item.value.view(np.int32)), item.x, item.y, item.z, item.label))
    return skeleton, available, volumes, boundaries, popped


def test_vectorized_initialization_matches_literal_java_loop_order() -> None:
    rng = np.random.default_rng(79031)
    labels = rng.integers(0, 5, size=(8, 7, 9), dtype=np.int16)
    memberships = rng.normal(0.1, 0.4, size=(*labels.shape, 5)).astype(np.float32)
    memberships[2, 3, 4, :] = np.float32(0.0)
    memberships[4, 2, 5, :] = np.float32(-0.25)
    memberships[5, 1, 3, :] = np.float32(0.5)
    initial_skeleton = rng.integers(-2, 5, size=labels.shape, dtype=np.int16)
    initial_available = rng.integers(0, 2, size=(*labels.shape, 5), dtype=np.uint8).astype(bool)
    zero_threshold = 1e-10

    expected = _source_order_reference(
        labels, memberships, initial_skeleton, initial_available, zero_threshold
    )
    actual = toads_initialize_thinning_queue(
        labels,
        memberships,
        initial_skeleton=initial_skeleton,
        initial_available=initial_available,
        zero_threshold=zero_threshold,
    )
    popped = []
    while actual.queue:
        item = actual.queue.pop_first()
        popped.append((int(item.value.view(np.int32)), item.x, item.y, item.z, item.label))

    np.testing.assert_array_equal(actual.skeleton, expected[0])
    np.testing.assert_array_equal(actual.available, expected[1])
    np.testing.assert_array_equal(actual.class_volumes.view(np.uint32), expected[2].view(np.uint32))
    np.testing.assert_array_equal(actual.boundary_counts, expected[3])
    assert popped == expected[4]
