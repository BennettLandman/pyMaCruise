"""Source-ordered TOADS growing without per-event Python objects.

The diagnostic path in toads_thinning retains event snapshots. This kernel
uses the same min-heap ties, float32 priorities, and 18/6 topology LUT.
"""
from __future__ import annotations
import numpy as np
from numba import njit
from .toads_thinning_compiled import _homeomorphic

_FACES = ((-1,0,0), (0,-1,0), (0,0,-1), (0,0,1), (0,1,0), (1,0,0))

@njit(cache=True)
def _min_push(values, payloads, size, priority, x, y, z, label):
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
        if not priority < values[parent]:
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
def _min_pop(values, payloads, size):
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
        if child != size and values[child + 1] < values[child]:
            child += 1
        if not values[child] < last_priority:
            break
        values[hole] = values[child]
        payloads[hole, :] = payloads[child, :]
        hole = child
    values[hole] = last_priority
    payloads[hole, :] = last_payload
    return values, payloads, size, priority, x, y, z, label



@njit(cache=True, fastmath=False)
def _ratio(mems, x, y, z, label, zero):
    best = 0
    maximum = np.float32(0.0)
    for candidate in range(mems.shape[3]):
        if mems[x,y,z,candidate] > maximum:
            best = candidate
            maximum = mems[x,y,z,candidate]
    if label == best:
        return np.float32(1.0)
    return np.float32(np.float32(mems[x,y,z,best] + zero) / np.float32(mems[x,y,z,label] + zero))


@njit(cache=True, fastmath=False)
def initialize_growing(skel, mems, orderings, available, empty_label, zero):
    nx, ny, nz = skel.shape
    classes = mems.shape[3]
    capacity = max(16, nx*ny + ny*nz + nz*nx)
    values = np.empty(capacity, dtype=np.float32)
    payloads = np.empty((capacity,4), dtype=np.int16)
    size = 0
    for x in range(1,nx-1):
        for y in range(1,ny-1):
            for z in range(1,nz-1):
                for label in range(classes):
                    available[x,y,z,label] = skel[x,y,z] == empty_label
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for label in range(classes):
                    if not available[x,y,z,label]:
                        continue
                    boundary = False
                    for dx,dy,dz in _FACES:
                        if skel[x+dx,y+dy,z+dz] == label:
                            boundary = True
                            break
                    if not boundary:
                        continue
                    available[x,y,z,label] = False
                    priority = np.float32(orderings[x,y,z] + _ratio(mems,x,y,z,label,zero))
                    values,payloads,size = _min_push(values,payloads,size,priority,x,y,z,label)
    return values,payloads,size


@njit(cache=True, fastmath=False)
def grow_queue(labels, skel, mems, available, orderings, volumes, values, payloads, size, empty_label, zero, lut, max_events):
    processed = 0
    while size and processed < max_events:
        values,payloads,size,priority,x,y,z,label = _min_pop(values,payloads,size)
        processed += 1
        if skel[x,y,z] != empty_label:
            continue
        if _homeomorphic(labels,x,y,z,label,lut):
            skel[x,y,z] = label
            available[x,y,z,:] = False
            orderings[x,y,z] = priority
            volumes[label] = np.float32(volumes[label] + mems[x,y,z,label])
            labels[x,y,z] = label
            for dx,dy,dz in _FACES:
                xi,yj,zl = x+dx,y+dy,z+dz
                if not available[xi,yj,zl,label]:
                    continue
                ratio = _ratio(mems,xi,yj,zl,label,zero)
                next_priority = np.float32(priority + ratio)
                available[xi,yj,zl,label] = False
                values,payloads,size = _min_push(values,payloads,size,next_priority,xi,yj,zl,label)
        else:
            available[x,y,z,label] = True
    return values,payloads,size,processed
