"""Python reproduction of the selected JIST ``KdTree`` traversal ordering.

This implements the selected JIST ``KdTree`` surface-area split, Java
priority-queue ordering, float32 ``KdTriangle.distance`` arithmetic, and
float32 projected points used by the closest-point thickness operation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math

import numpy as np


def _jist_triangle_distance_point(point, p1, p2, p3):
    """Float32 transcription of the selected Java ``KdTriangle.distance``."""
    zero = np.float32(0.0)
    one = np.float32(1.0)
    two = np.float32(2.0)
    diff = np.empty(3, dtype=np.float32)
    edge0 = np.empty(3, dtype=np.float32)
    edge1 = np.empty(3, dtype=np.float32)
    for axis in range(3):
        diff[axis] = np.float32(p1[axis] - point[axis])
        edge0[axis] = np.float32(p2[axis] - p1[axis])
        edge1[axis] = np.float32(p3[axis] - p1[axis])
    a00 = np.float32(np.float32(edge0[0] * edge0[0]) + np.float32(edge0[1] * edge0[1]))
    a00 = np.float32(a00 + np.float32(edge0[2] * edge0[2]))
    a01 = np.float32(np.float32(edge0[0] * edge1[0]) + np.float32(edge0[1] * edge1[1]))
    a01 = np.float32(a01 + np.float32(edge0[2] * edge1[2]))
    a11 = np.float32(np.float32(edge1[0] * edge1[0]) + np.float32(edge1[1] * edge1[1]))
    a11 = np.float32(a11 + np.float32(edge1[2] * edge1[2]))
    b0 = np.float32(np.float32(diff[0] * edge0[0]) + np.float32(diff[1] * edge0[1]))
    b0 = np.float32(b0 + np.float32(diff[2] * edge0[2]))
    b1 = np.float32(np.float32(diff[0] * edge1[0]) + np.float32(diff[1] * edge1[1]))
    b1 = np.float32(b1 + np.float32(diff[2] * edge1[2]))
    c = np.float32(np.float32(diff[0] * diff[0]) + np.float32(diff[1] * diff[1]))
    c = np.float32(c + np.float32(diff[2] * diff[2]))
    det = np.abs(np.float32(np.float32(a00 * a11) - np.float32(a01 * a01)))
    s = np.float32(np.float32(a01 * b1) - np.float32(a11 * b0))
    t = np.float32(np.float32(a01 * b0) - np.float32(a00 * b1))
    if np.float32(s + t) <= det:
        if s < zero:
            if t < zero:  # region 4
                if b0 < zero:
                    t = zero
                    if -b0 >= a00:
                        s = one
                        sq = np.float32(a00 + np.float32(two * b0) + c)
                    else:
                        s = np.float32(-b0 / a00)
                        sq = np.float32(np.float32(b0 * s) + c)
                else:
                    s = zero
                    if b1 >= zero:
                        t = zero
                        sq = c
                    elif -b1 >= a11:
                        t = one
                        sq = np.float32(a11 + np.float32(two * b1) + c)
                    else:
                        t = np.float32(-b1 / a11)
                        sq = np.float32(np.float32(b1 * t) + c)
            else:  # region 3
                s = zero
                if b1 >= zero:
                    t = zero
                    sq = c
                elif -b1 >= a11:
                    t = one
                    sq = np.float32(a11 + np.float32(two * b1) + c)
                else:
                    t = np.float32(-b1 / a11)
                    sq = np.float32(np.float32(b1 * t) + c)
        elif t < zero:  # region 5
            t = zero
            if b0 >= zero:
                s = zero
                sq = c
            elif -b0 >= a00:
                s = one
                sq = np.float32(a00 + np.float32(two * b0) + c)
            else:
                s = np.float32(-b0 / a00)
                sq = np.float32(np.float32(b0 * s) + c)
        else:  # region 0
            inv_det = np.float32(one / det)
            s = np.float32(s * inv_det)
            t = np.float32(t * inv_det)
            term0 = np.float32(np.float32(a00 * s) + np.float32(a01 * t) + np.float32(two * b0))
            term1 = np.float32(np.float32(a01 * s) + np.float32(a11 * t) + np.float32(two * b1))
            sq = np.float32(np.float32(s * term0) + np.float32(t * term1) + c)
    else:
        if s < zero:  # region 2
            tmp0 = np.float32(a01 + b0)
            tmp1 = np.float32(a11 + b1)
            if tmp1 > tmp0:
                numer = np.float32(tmp1 - tmp0)
                denom = np.float32(a00 - np.float32(two * a01) + a11)
                if numer >= denom:
                    s, t = one, zero
                    sq = np.float32(a00 + np.float32(two * b0) + c)
                else:
                    s = np.float32(numer / denom)
                    t = np.float32(one - s)
                    term0 = np.float32(np.float32(a00 * s) + np.float32(a01 * t) + np.float32(two * b0))
                    term1 = np.float32(np.float32(a01 * s) + np.float32(a11 * t) + np.float32(two * b1))
                    sq = np.float32(np.float32(s * term0) + np.float32(t * term1) + c)
            else:
                s = zero
                if tmp1 <= zero:
                    t = one
                    sq = np.float32(a11 + np.float32(two * b1) + c)
                elif b1 >= zero:
                    t = zero
                    sq = c
                else:
                    t = np.float32(-b1 / a11)
                    sq = np.float32(np.float32(b1 * t) + c)
        elif t < zero:  # region 6
            tmp0 = np.float32(a01 + b1)
            tmp1 = np.float32(a00 + b0)
            if tmp1 > tmp0:
                numer = np.float32(tmp1 - tmp0)
                denom = np.float32(a00 - np.float32(two * a01) + a11)
                if numer >= denom:
                    t, s = one, zero
                    sq = np.float32(a11 + np.float32(two * b1) + c)
                else:
                    t = np.float32(numer / denom)
                    s = np.float32(one - t)
                    term0 = np.float32(np.float32(a00 * s) + np.float32(a01 * t) + np.float32(two * b0))
                    term1 = np.float32(np.float32(a01 * s) + np.float32(a11 * t) + np.float32(two * b1))
                    sq = np.float32(np.float32(s * term0) + np.float32(t * term1) + c)
            else:
                t = zero
                if tmp1 <= zero:
                    s = one
                    sq = np.float32(a00 + np.float32(two * b0) + c)
                elif b0 >= zero:
                    s = zero
                    sq = c
                else:
                    s = np.float32(-b0 / a00)
                    sq = np.float32(np.float32(b0 * s) + c)
        else:  # region 1
            numer = np.float32(a11 + b1 - a01 - b0)
            if numer <= zero:
                s, t = zero, one
                sq = np.float32(a11 + np.float32(two * b1) + c)
            else:
                denom = np.float32(a00 - np.float32(two * a01) + a11)
                if numer >= denom:
                    s, t = one, zero
                    sq = np.float32(a00 + np.float32(two * b0) + c)
                else:
                    s = np.float32(numer / denom)
                    t = np.float32(one - s)
                    term0 = np.float32(np.float32(a00 * s) + np.float32(a01 * t) + np.float32(two * b0))
                    term1 = np.float32(np.float32(a01 * s) + np.float32(a11 * t) + np.float32(two * b1))
                    sq = np.float32(np.float32(s * term0) + np.float32(t * term1) + c)
    if sq < zero:
        sq = zero
    closest = np.empty(3, dtype=np.float32)
    for axis in range(3):
        e0 = np.float32(edge0[axis] * s)
        e1 = np.float32(edge1[axis] * t)
        closest[axis] = np.float32(np.float32(p1[axis] + e0) + e1)
    return math.sqrt(float(sq)), closest


@dataclass(slots=True)
class _Face:
    face_id: int
    minimum: np.ndarray
    maximum: np.ndarray
    depth: int = 1


@dataclass(slots=True)
class _Box:
    depth: int
    children: list[_Box | _Face] = field(default_factory=list)
    minimum: np.ndarray | None = None
    maximum: np.ndarray | None = None

    def update(self) -> None:
        if not self.children:
            self.minimum = np.full(3, 1e10, dtype=np.float32)
            self.maximum = np.full(3, -1e10, dtype=np.float32)
            return
        minima = []
        maxima = []
        for child in self.children:
            if isinstance(child, _Box):
                child.update()
            minima.append(child.minimum)
            maxima.append(child.maximum)
        self.minimum = np.min(np.stack(minima), axis=0).astype(np.float32)
        self.maximum = np.max(np.stack(maxima), axis=0).astype(np.float32)


class JistSurfaceTree:
    """Selected JIST SAH triangle tree and Java-compatible closest-point order."""

    def __init__(self, vertices: np.ndarray, faces: np.ndarray, *, max_depth: int = 13):
        points = np.asarray(vertices, dtype=np.float32)
        indices = np.asarray(faces, dtype=np.int64)
        triangles = points[indices]
        self.triangles = triangles
        eps = np.float32(5e-4)
        self.root = _Box(0, [
            _Face(
                i,
                np.subtract(triangle.min(axis=0).astype(np.float32), eps, dtype=np.float32),
                np.add(triangle.max(axis=0).astype(np.float32), eps, dtype=np.float32),
            )
            for i, triangle in enumerate(triangles)
        ])
        self.root.update()
        self._build(max_depth)

    @staticmethod
    def _split_position(
        children: list[_Box | _Face], box: _Box
    ) -> tuple[int, list[tuple[float, bool, _Box | _Face]]]:
        extent = np.subtract(box.maximum, box.minimum, dtype=np.float32)
        if extent[0] > extent[1] and extent[0] > extent[2]:
            dim = 0
        elif extent[1] > extent[2]:
            dim = 1
        else:
            dim = 2
        voxel = 1.0 / float(np.float32(np.float32(extent[0] + extent[1]) + extent[2]))
        retry = 3
        while retry:
            if dim == 0:
                summed = float(np.float32(extent[1] + extent[2]))
                product = float(np.float32(extent[1] * extent[2]))
            elif dim == 1:
                summed = float(np.float32(extent[0] + extent[2]))
                product = float(np.float32(extent[0] * extent[2]))
            else:
                summed = float(np.float32(extent[0] + extent[1]))
                product = float(np.float32(extent[0] * extent[1]))
            edge_list: list[tuple[float, bool, _Box | _Face]] = []
            for child in children:
                edge_list.append((float(child.maximum[dim]), False, child))
                edge_list.append((float(child.minimum[dim]), True, child))
            edge_list.sort(key=lambda edge: (edge[0], not edge[1]))
            lower, upper = float(box.minimum[dim]), float(box.maximum[dim])
            n_below, n_above = 0, len(children)
            min_cost, best = 1e30, -1
            for i, (value, is_min, _) in enumerate(edge_list):
                if not is_min:
                    n_above -= 1
                if value > lower and upper > value:
                    p_below = (summed * (value - lower) + product) * voxel
                    p_above = (summed * (upper - value) + product) * voxel
                    empty = n_above == 0 or n_below == 0
                    cost = 1.0 + 80.0 * (0.8 if empty else 1.0) * (
                        p_below * n_below + p_above * n_above
                    )
                    if cost < min_cost:
                        min_cost, best = cost, i
                if is_min:
                    n_below += 1
            if best >= 0:
                return best, edge_list
            dim = (dim + 1) % 3
            retry -= 1
        return -1, []

    def _build(self, max_depth: int) -> None:
        pending = [self.root]
        while pending:
            box = pending.pop(0)
            if box.depth > max_depth or len(box.children) < 2:
                continue
            split, edges = self._split_position(box.children, box)
            if split <= 0 or split == len(edges) - 1:
                continue
            left = _Box(box.depth + 1)
            right = _Box(box.depth + 1)
            for _, is_min, child in edges[:split]:
                if is_min:
                    child.depth = left.depth + 1
                    left.children.append(child)
            for _, is_min, child in edges[split + 1 :]:
                if not is_min:
                    child.depth = right.depth + 1
                    right.children.append(child)
            if len(left.children) >= len(box.children) or len(right.children) >= len(box.children):
                continue
            left.update()
            right.update()
            box.children = [left, right]
            box.update()
            pending.extend((left, right))

    @staticmethod
    def _box_distance(point: np.ndarray, box: _Box | _Face) -> float:
        point = np.asarray(point, dtype=np.float32)
        minimum, maximum = box.minimum, box.maximum
        if np.all(point >= minimum) and np.all(point <= maximum):
            return -1.0
        reference = np.maximum(np.minimum(point, maximum), minimum).astype(np.float32)
        delta = np.subtract(point, reference, dtype=np.float32)
        squared = np.float32(np.float32(delta[0] * delta[0]) + np.float32(delta[1] * delta[1]))
        squared = np.float32(squared + np.float32(delta[2] * delta[2]))
        return float(np.sqrt(squared, dtype=np.float32))

    @staticmethod
    def _heap_push(heap: list[tuple[float, int, _Box | _Face]], item: tuple[float, int, _Box | _Face]) -> None:
        # java.util.PriorityQueue leaves equal keys in their insertion position.
        heap.append(item)
        index = len(heap) - 1
        while index:
            parent = (index - 1) >> 1
            if heap[index][0] >= heap[parent][0]:
                break
            heap[index], heap[parent] = heap[parent], heap[index]
            index = parent

    @staticmethod
    def _heap_pop(heap: list[tuple[float, int, _Box | _Face]]) -> tuple[float, int, _Box | _Face]:
        result = heap[0]
        last = heap.pop()
        if heap:
            heap[0] = last
            index = 0
            half = len(heap) >> 1
            while index < half:
                child = (index << 1) + 1
                right = child + 1
                if right < len(heap) and heap[right][0] < heap[child][0]:
                    child = right
                if heap[index][0] <= heap[child][0]:
                    break
                heap[index], heap[child] = heap[child], heap[index]
                index = child
        return result

    def query(self, point: np.ndarray) -> tuple[np.ndarray, float, int]:
        query = np.asarray(point, dtype=np.float32)
        heap: list[tuple[float, int, _Box | _Face]] = []
        sequence = 0
        self._heap_push(heap, (self._box_distance(query, self.root), sequence, self.root))
        best_distance = 1e30
        best_point: np.ndarray | None = None
        best_face = -1
        while heap:
            queued_distance, _, box = self._heap_pop(heap)
            if not queued_distance < best_distance:
                continue
            if isinstance(box, _Face):
                triangle = self.triangles[box.face_id]
                distance, closest = _jist_triangle_distance_point(
                    query, triangle[0], triangle[1], triangle[2]
                )
                if distance < best_distance:
                    best_distance = distance
                    best_point = closest.astype(np.float64)
                    best_face = box.face_id
            else:
                for child in box.children:
                    sequence += 1
                    self._heap_push(heap, (self._box_distance(query, child), sequence, child))
        if best_point is None:
            raise RuntimeError("JIST surface tree found no closest triangle")
        return best_point, best_distance, best_face

    def query_many(self, points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        values = np.asarray(points, dtype=np.float32)
        closest = np.empty((len(values), 3), dtype=np.float64)
        distances = np.empty(len(values), dtype=np.float64)
        face_ids = np.empty(len(values), dtype=np.int64)
        for index, point in enumerate(values):
            closest[index], distances[index], face_ids[index] = self.query(point)
        return closest, distances, face_ids

    def structure_signature(self) -> list[tuple[object, ...]]:
        """Return path-ordered nodes for direct comparison with the Java tree."""
        signature: list[tuple[object, ...]] = []

        def visit(node: _Box | _Face, path: str) -> None:
            if isinstance(node, _Face):
                signature.append((path, "T", node.face_id))
                return
            signature.append((
                path,
                "B",
                node.depth,
                len(node.children),
                tuple(float(value) for value in node.minimum),
                tuple(float(value) for value in node.maximum),
            ))
            for index, child in enumerate(node.children):
                visit(child, path + str(index))

        visit(self.root, "r")
        return signature
