"""Compiled unit-speed FMM used to reinitialize CRUISE TGDM level sets.

This is a Numba translation of the ordered JIST ``DistanceField`` heap path.
The general weighted ACE FMM remains in :mod:`pymacruise.ace`; this kernel is
specific to unit-speed TGDM reinitialization.
"""

from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True, fastmath=False)
def _heap_set(
    nodes: np.ndarray,
    priorities: np.ndarray,
    positions: np.ndarray,
    size: int,
    node: int,
    priority_value: float,
) -> int:
    """Indexed binary heap with JIST's strict-less-than tie behavior."""
    priority = np.float32(priority_value)
    index = int(positions[node])
    if index < 0:
        index = size
        size += 1
        nodes[index] = node
        priorities[index] = priority
        positions[node] = index
        item_node = node
        item_priority = priority
        while index > 0:
            parent = (index - 1) // 2
            if not item_priority < priorities[parent]:
                break
            nodes[index] = nodes[parent]
            priorities[index] = priorities[parent]
            positions[nodes[index]] = index
            index = parent
        nodes[index] = item_node
        priorities[index] = item_priority
        positions[item_node] = index
        return size

    previous = priorities[index]
    priorities[index] = priority
    if priority < previous:
        item_node = node
        item_priority = priority
        while index > 0:
            parent = (index - 1) // 2
            if not item_priority < priorities[parent]:
                break
            nodes[index] = nodes[parent]
            priorities[index] = priorities[parent]
            positions[nodes[index]] = index
            index = parent
        nodes[index] = item_node
        priorities[index] = item_priority
        positions[item_node] = index
        return size

    item_node = node
    item_priority = priority
    while 2 * index + 1 < size:
        child = 2 * index + 1
        right = child + 1
        # JIST keeps the left child when child priorities are equal.
        if right < size and priorities[right] < priorities[child]:
            child = right
        if not priorities[child] < item_priority:
            break
        nodes[index] = nodes[child]
        priorities[index] = priorities[child]
        positions[nodes[index]] = index
        index = child
    nodes[index] = item_node
    priorities[index] = item_priority
    positions[item_node] = index
    return size


@njit(cache=True, fastmath=False)
def _heap_pop(
    nodes: np.ndarray,
    priorities: np.ndarray,
    positions: np.ndarray,
    size: int,
) -> tuple[np.float32, int, int]:
    priority = priorities[0]
    node = int(nodes[0])
    size -= 1
    positions[node] = -1
    if size > 0:
        item_node = int(nodes[size])
        item_priority = priorities[size]
        index = 0
        while 2 * index + 1 < size:
            child = 2 * index + 1
            right = child + 1
            if right < size and priorities[right] < priorities[child]:
                child = right
            if not priorities[child] < item_priority:
                break
            nodes[index] = nodes[child]
            priorities[index] = priorities[child]
            positions[nodes[index]] = index
            index = child
        nodes[index] = item_node
        priorities[index] = item_priority
        positions[item_node] = index
    return priority, node, size


@njit(cache=True, fastmath=False)
def _crossing_seed(
    values: np.ndarray,
    node: int,
    i: int,
    j: int,
    k: int,
    nx: int,
    ny: int,
    nz: int,
) -> tuple[bool, np.float64]:
    """Return the source-compatible subvoxel distance seed for one voxel."""
    center = np.float64(values[node])
    inv_square_sum = np.float64(0.0)
    found = False
    # Keep the source axis order and its asymmetric equal-magnitude decisions:
    # north/south chooses south on a tie, east/west chooses west, front/back
    # chooses back.
    for axis_index in range(3):
        if axis_index == 0:
            axis = 1
            first_delta = -1
            second_delta = 1
            tie_to_second = True
        elif axis_index == 1:
            axis = 0
            first_delta = -1
            second_delta = 1
            tie_to_second = False
        else:
            axis = 2
            first_delta = 1
            second_delta = -1
            tie_to_second = True

        first_node = -1
        second_node = -1
        first_value = np.float64(0.0)
        second_value = np.float64(0.0)
        if axis == 0:
            if 0 <= i + first_delta < nx:
                first_node = node + first_delta * ny * nz
                first_value = np.float64(values[first_node])
            if 0 <= i + second_delta < nx:
                second_node = node + second_delta * ny * nz
                second_value = np.float64(values[second_node])
        elif axis == 1:
            if 0 <= j + first_delta < ny:
                first_node = node + first_delta * nz
                first_value = np.float64(values[first_node])
            if 0 <= j + second_delta < ny:
                second_node = node + second_delta * nz
                second_value = np.float64(values[second_node])
        else:
            if 0 <= k + first_delta < nz:
                first_node = node + first_delta
                first_value = np.float64(values[first_node])
            if 0 <= k + second_delta < nz:
                second_node = node + second_delta
                second_value = np.float64(values[second_node])

        first_opposite = first_node >= 0 and center * first_value < 0.0
        second_opposite = second_node >= 0 and center * second_value < 0.0
        if first_opposite or second_opposite:
            found = True
            adjacent = first_value if first_opposite else second_value
            if first_opposite and second_opposite:
                first_abs = abs(first_value)
                second_abs = abs(second_value)
                if second_abs > first_abs or (tie_to_second and second_abs == first_abs):
                    adjacent = second_value
            offset = abs(center / (center - adjacent))
            inv_square_sum += 1.0 / (offset * offset)
    if not found:
        return False, np.float64(0.0)
    return True, 1.0 / np.sqrt(inv_square_sum)


@njit(cache=True, fastmath=False)
def _march_candidate(
    node: int,
    i: int,
    j: int,
    k: int,
    nx: int,
    ny: int,
    nz: int,
    state: np.ndarray,
    time: np.ndarray,
    characteristic: np.ndarray,
) -> tuple[np.float64, np.float64]:
    """Compute one FMM update with the JIST neighbor and arithmetic order."""
    min_time = np.float64(np.inf)
    min_char = np.float64(0.0)
    # Source order used to select the first minimum characteristic.
    if j > 0:
        n = node - nz
        if state[n] == 2 and np.float64(time[n]) < min_time:
            min_time = np.float64(time[n]); min_char = np.float64(characteristic[n])
    if j + 1 < ny:
        n = node + nz
        if state[n] == 2 and np.float64(time[n]) < min_time:
            min_time = np.float64(time[n]); min_char = np.float64(characteristic[n])
    if k + 1 < nz:
        n = node + 1
        if state[n] == 2 and np.float64(time[n]) < min_time:
            min_time = np.float64(time[n]); min_char = np.float64(characteristic[n])
    if k > 0:
        n = node - 1
        if state[n] == 2 and np.float64(time[n]) < min_time:
            min_time = np.float64(time[n]); min_char = np.float64(characteristic[n])
    if i + 1 < nx:
        n = node + ny * nz
        if state[n] == 2 and np.float64(time[n]) < min_time:
            min_time = np.float64(time[n]); min_char = np.float64(characteristic[n])
    if i > 0:
        n = node - ny * nz
        if state[n] == 2 and np.float64(time[n]) < min_time:
            min_time = np.float64(time[n]); min_char = np.float64(characteristic[n])

    # Unit speed: Python's ``min_time += 1.0 / local_speed`` is double math.
    min_time += 1.0
    summation = np.float64(0.0)
    squares = np.float64(0.0)
    char_sum = np.float64(0.0)
    char_time_sum = np.float64(0.0)
    count = 0

    # Group order: N/S, E/W, F/B. A tie within a pair selects the first.
    for group in range(3):
        first = -1
        second = -1
        if group == 0:
            if j > 0: first = node - nz
            if j + 1 < ny: second = node + nz
        elif group == 1:
            if i + 1 < nx: first = node + ny * nz
            if i > 0: second = node - ny * nz
        else:
            if k + 1 < nz: first = node + 1
            if k > 0: second = node - 1

        chosen = -1
        if first >= 0 and state[first] == 2:
            chosen = first
        if second >= 0 and state[second] == 2:
            if chosen < 0 or np.float64(time[second]) < np.float64(time[chosen]):
                chosen = second
        if chosen >= 0:
            chosen_time = np.float64(time[chosen])
            if chosen_time <= min_time:
                chosen_char = np.float64(characteristic[chosen])
                summation += chosen_time
                squares += chosen_time * chosen_time
                char_sum += chosen_char
                char_time_sum += chosen_char * chosen_time
                count += 1

    if count == 0:
        return min_time, min_char + 1.0
    discriminant = summation * summation - count * (squares - 1.0)
    if discriminant < 0.0:
        return min_time, min_char + 1.0
    result = (summation + np.sqrt(discriminant)) / count
    denominator = count * result - summation
    char = (1.0 + char_sum * result - char_time_sum) / denominator
    return result, char


@njit(cache=True, fastmath=False)
def _unit_speed_fmm_kernel(
    values: np.ndarray, nx: int, ny: int, nz: int, maximum: np.float64
) -> np.ndarray:
    count = nx * ny * nz
    time = np.full(count, np.inf, dtype=np.float32)
    characteristic = np.zeros(count, dtype=np.float32)
    state = np.zeros(count, dtype=np.uint8)  # far, narrow band, accepted

    # Java's source scans z/y/x with x innermost.
    accepted_count = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                node = (i * ny + j) * nz + k
                value = np.float64(values[node])
                if value == 0.0:
                    time[node] = np.float32(0.0)
                    state[node] = 2
                    accepted_count += 1
                    continue
                found, char = _crossing_seed(values, node, i, j, k, nx, ny, nz)
                if found:
                    characteristic[node] = np.float32(char)
                    time[node] = np.float32(char)
                    state[node] = 2
                    accepted_count += 1

    if accepted_count == 0:
        return np.zeros((nx, ny, nz), dtype=np.float32)

    heap_nodes = np.empty(count, dtype=np.int64)
    heap_priorities = np.empty(count, dtype=np.float32)
    positions = np.full(count, -1, dtype=np.int64)
    heap_size = 0
    offsets = (ny * nz, nz, -ny * nz, -nz, 1, -1)

    # Initial narrow band: accepted scan order then MaskVolume6 neighbor order.
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                accepted = (i * ny + j) * nz + k
                if state[accepted] != 2:
                    continue
                for index in range(6):
                    candidate = accepted + offsets[index]
                    if index == 0 and i + 1 >= nx: continue
                    if index == 1 and j + 1 >= ny: continue
                    if index == 2 and i == 0: continue
                    if index == 3 and j == 0: continue
                    if index == 4 and k + 1 >= nz: continue
                    if index == 5 and k == 0: continue
                    if state[candidate] != 0:
                        continue
                    state[candidate] = 1
                    value, char = _march_candidate(
                        candidate,
                        candidate // (ny * nz),
                        (candidate // nz) % ny,
                        candidate % nz,
                        nx, ny, nz, state, time, characteristic,
                    )
                    time[candidate] = np.float32(value)
                    characteristic[candidate] = np.float32(char)
                    heap_size = _heap_set(
                        heap_nodes, heap_priorities, positions,
                        heap_size, candidate, value,
                    )

    while heap_size > 0:
        value32, node, heap_size = _heap_pop(
            heap_nodes, heap_priorities, positions, heap_size
        )
        value = np.float64(value32)
        if value > maximum:
            break
        time[node] = value32
        state[node] = 2
        i = node // (ny * nz)
        j = (node // nz) % ny
        k = node % nz
        for index in range(6):
            candidate = node + offsets[index]
            if index == 0 and i + 1 >= nx: continue
            if index == 1 and j + 1 >= ny: continue
            if index == 2 and i == 0: continue
            if index == 3 and j == 0: continue
            if index == 4 and k + 1 >= nz: continue
            if index == 5 and k == 0: continue
            if state[candidate] == 2:
                continue
            ci = candidate // (ny * nz)
            cj = (candidate // nz) % ny
            ck = candidate % nz
            tentative, char = _march_candidate(
                candidate, ci, cj, ck, nx, ny, nz, state, time, characteristic
            )
            time[candidate] = np.float32(tentative)
            characteristic[candidate] = np.float32(char)
            if state[candidate] == 1:
                heap_size = _heap_set(
                    heap_nodes, heap_priorities, positions,
                    heap_size, candidate, tentative,
                )
            else:
                state[candidate] = 1
                heap_size = _heap_set(
                    heap_nodes, heap_priorities, positions,
                    heap_size, candidate, tentative,
                )

    for node in range(count):
        if state[node] != 2:
            time[node] = np.float32(maximum)
        if values[node] < 0.0:
            time[node] *= np.float32(-1.0)
    return time.reshape((nx, ny, nz))


def tgdm_signed_distance_fmm(level_set: np.ndarray, max_distance: float) -> np.ndarray:
    """Reinitialize a TGDM level set with the ordered unit-speed JIST FMM."""
    values = np.asarray(level_set, dtype=np.float32)
    if not np.isfinite(max_distance) or max_distance <= 0:
        raise ValueError("TGDM signed-distance maximum must be finite and positive")
    if values.ndim != 3 or any(size < 2 for size in values.shape):
        raise ValueError("TGDM signed distance requires a 3D grid with dimensions >= 2")
    if not np.isfinite(values).all():
        raise ValueError("TGDM signed-distance input must be finite")
    return _unit_speed_fmm_kernel(
        values.ravel(order="C"),
        values.shape[0],
        values.shape[1],
        values.shape[2],
        np.float64(max_distance),
    )
