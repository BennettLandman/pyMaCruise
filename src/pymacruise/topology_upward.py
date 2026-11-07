"""V3.2 intensity-initialized background-to-object topology propagation.

Port of AlgorithmTopologyCorrection and TopologyPropagation's upward exact
smoothing. Source notices are retained under licenses/upstream/. The paint
input is deliberately unused in this branch, as in the selected JIST wrapper.
"""
from __future__ import annotations

import numpy as np
from numba import njit

from .toads_growing_compiled import _min_pop, _min_push


@njit(cache=True)
def _components(cube, connectivity, face_only):
    seen = np.zeros((3, 3, 3), dtype=np.bool_)
    queue = np.empty((27, 3), dtype=np.int32)
    count = 0
    for x in range(3):
        for y in range(3):
            for z in range(3):
                if not cube[x, y, z] or seen[x, y, z]:
                    continue
                seen[x, y, z] = True
                queue[0, 0], queue[0, 1], queue[0, 2] = x, y, z
                head, size = 0, 1
                touches = False
                while head < size:
                    a, b, c = queue[head, 0], queue[head, 1], queue[head, 2]
                    head += 1
                    if abs(a-1) + abs(b-1) + abs(c-1) == 1:
                        touches = True
                    for dx in range(-1, 2):
                        for dy in range(-1, 2):
                            for dz in range(-1, 2):
                                norm = abs(dx) + abs(dy) + abs(dz)
                                if norm == 0 or norm > connectivity:
                                    continue
                                i, j, k = a+dx, b+dy, c+dz
                                if 0 <= i < 3 and 0 <= j < 3 and 0 <= k < 3:
                                    if cube[i,j,k] and not seen[i,j,k]:
                                        seen[i,j,k] = True
                                        queue[size,0],queue[size,1],queue[size,2] = i,j,k
                                        size += 1
                if not face_only or touches:
                    count += 1
                    if count > 1:
                        return count
    return count


@njit(cache=True)
def _regular_background(labels, x, y, z):
    positive = np.zeros((3,3,3), dtype=np.bool_)
    negative = np.zeros((3,3,3), dtype=np.bool_)
    for dx in range(-1,2):
        for dy in range(-1,2):
            for dz in range(-1,2):
                norm = abs(dx)+abs(dy)+abs(dz)
                if norm == 0:
                    continue
                value = labels[x+dx,y+dy,z+dz]
                positive[dx+1,dy+1,dz+1] = value > 70
                if norm <= 2:
                    negative[dx+1,dy+1,dz+1] = value <= 70
    return _components(positive, 2, False) == 1 and _components(negative, 1, True) == 1


@njit(cache=True)
def upward_exact_smoothing(image, low, high, minimum_distance):
    """Execute the source float32 min-tree sequence on an indexing lattice."""
    nx,ny,nz = image.shape
    labels = np.full((nx+2,ny+2,nz+2), 30, dtype=np.uint8)
    padded = np.zeros(labels.shape, dtype=np.float32)
    padded[1:-1,1:-1,1:-1] = image
    distance = np.zeros(labels.shape, dtype=np.float32)
    # The source expands by one zero plane, then seeds two bounding planes.
    # Only the interior plane of each pair is in its processing mask.
    for x in range(nx+2):
        for y in range(ny+2):
            for z in range(nz+2):
                if x==0 or y==0 or z==0 or x==nx+1 or y==ny+1 or z==nz+1:
                    labels[x,y,z] = 110  # BACKFRONTIER participates in topology.
                elif x==1 or y==1 or z==1 or x==nx or y==ny or z==nz:
                    labels[x,y,z] = 80
                    distance[x,y,z] = padded[x,y,z]
    capacity = max(16, nx*ny + ny*nz + nz*nx)
    values = np.empty(capacity, dtype=np.float32)
    payloads = np.empty((capacity,4), dtype=np.int16)
    size = 0
    for x in range(1,nx+1):
        for y in range(1,ny+1):
            for z in range(1,nz+1):
                if labels[x,y,z] == 80:
                    continue
                adjacent = False
                for dx in range(-1,2):
                    for dy in range(-1,2):
                        for dz in range(-1,2):
                            if dx*dx+dy*dy+dz*dz < 3 and labels[x+dx,y+dy,z+dz] == 80:
                                adjacent = True
                if adjacent:
                    labels[x,y,z] = 40
                    distance[x,y,z] = padded[x,y,z]
                    values,payloads,size = _min_push(values,payloads,size,padded[x,y,z],x,y,z,0)
    value = np.float32(0)
    while size:
        values,payloads,size,value,x,y,z,_ = _min_pop(values,payloads,size)
        higher = value
        for dx in range(-1,2):
            for dy in range(-1,2):
                for dz in range(-1,2):
                    if labels[x+dx,y+dy,z+dz] == 80:
                        higher = max(higher, distance[x+dx,y+dy,z+dz])
        if higher > value:
            values,payloads,size = _min_push(values,payloads,size,higher,x,y,z,0)
        else:  # The v3.2 wrapper selects -inStop false.
            labels[x,y,z] = 70
            if _regular_background(labels,x,y,z):
                labels[x,y,z] = 80
                distance[x,y,z] = value
                for dx in range(-1,2):
                    for dy in range(-1,2):
                        for dz in range(-1,2):
                            i,j,k = x+dx,y+dy,z+dz
                            if dx*dx+dy*dy+dz*dz < 3 and labels[i,j,k] == 30:
                                labels[i,j,k] = 40
                                proposed = max(padded[i,j,k], np.float32(value+minimum_distance))
                                distance[i,j,k] = proposed
                                values,payloads,size = _min_push(values,payloads,size,proposed,i,j,k,0)
                            if labels[i,j,k] == 20:
                                labels[i,j,k] = 40
                                values,payloads,size = _min_push(values,payloads,size,np.float32(value+minimum_distance),i,j,k,0)
            else:
                labels[x,y,z] = 20
    for x in range(1,nx+1):
        for y in range(1,ny+1):
            for z in range(1,nz+1):
                if labels[x,y,z] == 20 or labels[x,y,z] == 30:
                    distance[x,y,z] = value
    return distance[1:-1,1:-1,1:-1].copy()
