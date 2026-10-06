"""Numba execution path for the selected source-ordered TOADS thinning loop.

The kernel preserves the Python reference implementation's heap tie rules,
float32 priority arithmetic, x/y/z neighborhood order, and 18/6 LUT checks.
It is an execution optimization only; ``toads_thinning.py`` remains the
readable reference path and supplies validation and event-trace behavior.
"""

from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def _heap_push(values, payloads, size, priority, x, y, z, label):
    if size + 1 >= values.shape[0]:
        capacity = max(16, values.shape[0] * 2)
        grown_values = np.empty(capacity, dtype=np.float32)
        grown_payloads = np.empty((capacity, 4), dtype=np.int16)
        grown_values[: size + 1] = values[: size + 1]
        grown_payloads[: size + 1, :] = payloads[: size + 1, :]
        values = grown_values
        payloads = grown_payloads

    size += 1
    hole = size
    while hole > 1:
        parent = hole // 2
        if not priority > values[parent]:
            break
        values[hole] = values[parent]
        payloads[hole, :] = payloads[parent, :]
        hole = parent
    values[hole] = priority
    payloads[hole, 0] = x
    payloads[hole, 1] = y
    payloads[hole, 2] = z
    payloads[hole, 3] = label
    return values, payloads, size


@njit(cache=True)
def _heap_pop(values, payloads, size):
    priority = values[1]
    x = payloads[1, 0]
    y = payloads[1, 1]
    z = payloads[1, 2]
    label = payloads[1, 3]
    last_priority = values[size]
    last_payload = payloads[size, :].copy()
    size -= 1
    if size == 0:
        return values, payloads, size, priority, x, y, z, label

    hole = 1
    while hole * 2 <= size:
        child = hole * 2
        if child != size and values[child + 1] > values[child]:
            child += 1
        if not values[child] > last_priority:
            break
        values[hole] = values[child]
        payloads[hole, :] = payloads[child, :]
        hole = child
    values[hole] = last_priority
    payloads[hole, :] = last_payload
    return values, payloads, size, priority, x, y, z, label


@njit(cache=True, inline="always")
def _regular(labels, x, y, z, first, second, third, lut):
    key = 0
    bit = 0
    for dx in range(-1, 2):
        for dy in range(-1, 2):
            for dz in range(-1, 2):
                if dx == 0 and dy == 0 and dz == 0:
                    continue
                value = labels[x + dx, y + dy, z + dz]
                if value == first or value == second or value == third:
                    key |= 1 << bit
                bit += 1
    return lut[key] == 1


@njit(cache=True)
def _homeomorphic(labels, x, y, z, new_label, lut):
    old_label = labels[x, y, z]
    if old_label == new_label:
        return True
    if not _regular(labels, x, y, z, new_label, -999999, -999999, lut):
        return False
    if not _regular(labels, x, y, z, old_label, -999999, -999999, lut):
        return False

    related = np.empty(26, dtype=np.int32)
    related_count = 0
    for dx in range(-1, 2):
        for dy in range(-1, 2):
            for dz in range(-1, 2):
                if dx == 0 and dy == 0 and dz == 0:
                    continue
                label = labels[x + dx, y + dy, z + dz]
                if label == new_label or label == old_label:
                    continue
                seen = False
                for index in range(related_count):
                    if related[index] == label:
                        seen = True
                        break
                if not seen:
                    related[related_count] = label
                    related_count += 1

    for first_index in range(related_count):
        relation = related[first_index]
        if not _regular(labels, x, y, z, old_label, relation, -999999, lut):
            return False
        if not _regular(labels, x, y, z, new_label, relation, -999999, lut):
            return False

    for first_index in range(related_count):
        for second_index in range(first_index + 1, related_count):
            first_relation = related[first_index]
            second_relation = related[second_index]
            if not _regular(
                labels, x, y, z, old_label, first_relation, second_relation, lut
            ):
                return False
            if not _regular(
                labels, x, y, z, new_label, first_relation, second_relation, lut
            ):
                return False
    return True


@njit(cache=True)
def _run_thinning(
    labels,
    skeleton,
    memberships,
    available,
    relations,
    shapes,
    volumes,
    minimums,
    ordering,
    heap_values,
    heap_payloads,
    heap_size,
    limit,
    empty_label,
    zero,
    lut,
):
    previous_priority = np.float32(0.0)
    processed = 0
    nx, ny, nz = labels.shape
    classes = memberships.shape[3]
    neighbor_x = np.array([-1, 0, 0, 0, 0, 1], dtype=np.int64)
    neighbor_y = np.array([0, -1, 0, 0, 1, 0], dtype=np.int64)
    neighbor_z = np.array([0, 0, -1, 1, 0, 0], dtype=np.int64)

    while heap_size > 0 and previous_priority > np.float32(-limit):
        (
            heap_values,
            heap_payloads,
            heap_size,
            popped_priority,
            x_short,
            y_short,
            z_short,
            label_short,
        ) = _heap_pop(heap_values, heap_payloads, heap_size)
        x = int(x_short)
        y = int(y_short)
        z = int(z_short)
        label = int(label_short)
        previous_priority = popped_priority
        processed += 1

        if not (
            1 <= x < nx - 1
            and 1 <= y < ny - 1
            and 1 <= z < nz - 1
            and 0 <= label < classes
            and skeleton[x, y, z] == label
            and volumes[label] > minimums[label]
        ):
            continue

        target = 0
        maximum = np.float32(0.0)
        for candidate in range(classes):
            score = memberships[x, y, z, candidate]
            if score > maximum:
                target = candidate
                maximum = score

        candidate_label = target
        change_boundary = False
        if target != label:
            direct = False
            next_label = -1
            best_prior = np.float32(-0.1)
            for neighbor in range(6):
                nx1 = x + neighbor_x[neighbor]
                ny1 = y + neighbor_y[neighbor]
                nz1 = z + neighbor_z[neighbor]
                neighbor_label = int(labels[nx1, ny1, nz1])
                if neighbor_label == target:
                    direct = True
                elif neighbor_label != label and relations[neighbor_label, target]:
                    shape_prior = shapes[nx1, ny1, nz1, target]
                    if shape_prior > best_prior:
                        next_label = neighbor_label
                        best_prior = shape_prior
            if direct:
                change_boundary = True
            elif next_label >= 0:
                candidate_label = next_label
                change_boundary = True
            if change_boundary and not _homeomorphic(
                labels, x, y, z, candidate_label, lut
            ):
                change_boundary = False

        # The selected Python transition computes the local topology predicate
        # before this volume gate, but does not mutate state when it is false.
        if volumes[label] <= minimums[label]:
            continue

        if target != label and change_boundary:
            labels[x, y, z] = candidate_label
        skeleton[x, y, z] = empty_label
        volumes[label] = np.float32(volumes[label] - memberships[x, y, z, label])
        available[x, y, z, label] = False
        ordering[x, y, z] = popped_priority

        for neighbor in range(6):
            nx1 = x + neighbor_x[neighbor]
            ny1 = y + neighbor_y[neighbor]
            nz1 = z + neighbor_z[neighbor]
            if not available[nx1, ny1, nz1, label]:
                continue

            neighbor_count = 0
            for dx in range(-1, 2):
                for dy in range(-1, 2):
                    for dz in range(-1, 2):
                        if labels[nx1 + dx, ny1 + dy, nz1 + dz] == label:
                            neighbor_count += 1

            if neighbor_count == 1:
                order = popped_priority
            else:
                best = 0
                best_membership = np.float32(0.0)
                for class_index in range(classes):
                    score = memberships[nx1, ny1, nz1, class_index]
                    if score > best_membership:
                        best = class_index
                        best_membership = score
                if label == best:
                    ratio = np.float32(1.0)
                else:
                    numerator = np.float32(best_membership + zero)
                    denominator = np.float32(memberships[nx1, ny1, nz1, label] + zero)
                    ratio = np.float32(numerator / denominator)
                order = np.float32(
                    popped_priority - np.float32(np.float32(1.0) / ratio)
                )

            if volumes[label] > minimums[label] and order > np.float32(-limit):
                available[nx1, ny1, nz1, label] = False
                (
                    heap_values,
                    heap_payloads,
                    heap_size,
                ) = _heap_push(
                    heap_values,
                    heap_payloads,
                    heap_size,
                    order,
                    nx1,
                    ny1,
                    nz1,
                    label,
                )

    return labels, skeleton, available, volumes, ordering, heap_values, heap_payloads, heap_size, processed


def run_compiled_thinning_queue(
    labels: np.ndarray,
    skeleton: np.ndarray,
    memberships: np.ndarray,
    available: np.ndarray,
    relations: np.ndarray,
    shapes: np.ndarray,
    volumes: np.ndarray,
    minimums: np.ndarray,
    ordering: np.ndarray,
    queue,
    limit: np.float32,
    empty_label: int,
    zero: np.float32,
    lut: np.ndarray,
):
    """Execute a no-event-trace thinning pass in the compiled source loop."""
    values = np.asarray(queue._values, dtype=np.float32)
    payloads = np.zeros((len(queue._payloads), 4), dtype=np.int16)
    for index, payload in enumerate(queue._payloads):
        if payload is not None:
            payloads[index, :] = payload
    (
        labels,
        skeleton,
        available,
        volumes,
        ordering,
        values,
        payloads,
        size,
        processed,
    ) = _run_thinning(
        labels,
        skeleton,
        memberships,
        available,
        relations,
        shapes,
        volumes,
        minimums,
        ordering,
        values,
        payloads,
        queue.size,
        limit,
        empty_label,
        zero,
        lut,
    )
    queue._values = values[: size + 1].astype(np.float32).tolist()
    queue._payloads = [None] + [
        tuple(int(value) for value in payloads[index, :])
        for index in range(1, size + 1)
    ]
    return labels, skeleton, available, volumes, ordering, processed
