"""TOADS relation stencils with source-order synchronous double buffering."""
from __future__ import annotations

import numpy as np
from numba import njit

_FACES = ((-1,0,0), (0,-1,0), (0,0,-1), (0,0,1), (0,1,0), (1,0,0))


@njit(cache=True, fastmath=False)
def diffuse_relations(relation_map, table, memberships, limit):
    current = relation_map.copy()
    updated = current.copy()
    nx, ny, nz = current.shape
    classes = memberships.shape[3]
    iterations = 0
    change = True
    for _ in range(1, 10 * limit):
        if not change:
            break
        change = False
        for x in range(1,nx-1):
            for y in range(1,ny-1):
                for z in range(1,nz-1):
                    best = current[x,y,z]
                    if best >= classes:
                        updated[x,y,z] = best
                        continue
                    change = True
                    best_support = np.float32(0.0)
                    for dx,dy,dz in _FACES:
                        xi,yj,zl = x+dx,y+dy,z+dz
                        candidate = current[xi,yj,zl]
                        if candidate < classes:
                            continue
                        support = np.float32(0.0)
                        for k in range(classes):
                            term = memberships[xi,yj,zl,k] if table[candidate,k] else np.float32(0.0)
                            support = np.float32(support + term)
                        if support > best_support:
                            best_support = support
                            best = candidate
                    updated[x,y,z] = best
        current,updated = updated,current
        iterations += 1
    return current,iterations


@njit(cache=True, fastmath=False)
def closest_labels(labels, limit):
    nx, ny, nz = labels.shape
    current = np.full((nx,ny,nz,4), np.int8(-1), dtype=np.int8)
    current[:,:,:,0] = labels
    for x in range(1,nx-1):
        for y in range(1,ny-1):
            for z in range(1,nz-1):
                own = labels[x,y,z]
                for dx,dy,dz in _FACES:
                    candidate = labels[x+dx,y+dy,z+dz]
                    if candidate == own:
                        continue
                    if current[x,y,z,1] == -1:
                        current[x,y,z,1] = candidate
                    elif current[x,y,z,1] != candidate and current[x,y,z,2] == -1:
                        current[x,y,z,2] = candidate
                    elif current[x,y,z,1] != candidate and current[x,y,z,2] != candidate and current[x,y,z,3] == -1:
                        current[x,y,z,3] = candidate
    updated = current.copy()
    for _ in range(1,5*limit):
        for x in range(1,nx-1):
            for y in range(1,ny-1):
                for z in range(1,nz-1):
                    for k in range(4):
                        updated[x,y,z,k] = current[x,y,z,k]
                    missing = 0
                    for k in range(1,4):
                        if current[x,y,z,k] == -1:
                            missing = k
                            break
                    if missing == 0:
                        continue
                    # The source chooses the missing slot once per pass. Later
                    # faces may overwrite it; they do not fill another slot.
                    for dx,dy,dz in _FACES:
                        for k in range(1,2 if missing == 1 else 4):
                            candidate = current[x+dx,y+dy,z+dz,k]
                            if candidate == -1:
                                continue
                            duplicate = False
                            for previous in range(missing):
                                if candidate == current[x,y,z,previous]:
                                    duplicate = True
                            if not duplicate:
                                updated[x,y,z,missing] = candidate
                                break
        current,updated = updated,current
    return current


@njit(cache=True, fastmath=False)
def boundary_distances(labels, distances):
    output = distances.copy()
    nx,ny,nz = labels.shape
    for x in range(1,nx-1):
        for y in range(1,ny-1):
            for z in range(1,nz-1):
                own = labels[x,y,z]
                for dx,dy,dz in _FACES:
                    xi,yj,zl = x+dx,y+dy,z+dz
                    neighbor = labels[xi,yj,zl]
                    if own != neighbor:
                        output[x,y,z,neighbor] = np.int8(1)
                        output[xi,yj,zl,own] = np.int8(1)
    # The selected subsequent diffusion is a no-op: only zero entries can
    # change and min(0, positive_neighbor + 1) remains zero. Keep other values.
    return output
