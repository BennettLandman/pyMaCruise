"""Tie-faithful port of VandyCRUISE's four-coordinate binary priority heap.

Behavior reference: ``BinaryHeap4D.java`` (July 2004, Pierre-Louis Bazin);
the original Java source is retained as a test-only reference.

Paper context: the source-compatible priority queue schedules the fast-marching
style multi-object front used by TOADS (Bazin & Pham, 2007, 2008). The heap is
infrastructure only; it does not encode tissue likelihoods or topology rules.
See ``docs/papers/toads.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ToadsHeapItem:
    value: np.float32
    x: int
    y: int
    z: int
    label: int


def _short(value: int) -> int:
    return ((int(value) + 32768) % 65536) - 32768


class ToadsPriorityHeap4D:
    """Source-order-compatible min/max queue for ``(x, y, z, label)`` items.

    Equal priorities retain the selected Java heap's insertion/removal path;
    this is not a stable FIFO queue. Values are rounded to float32 and integer
    coordinates/labels are narrowed to signed 16-bit values like Java shorts.
    """

    def __init__(self, capacity: int, mode: str = "min") -> None:
        if not isinstance(capacity, int) or capacity < 1:
            raise ValueError("capacity must be a positive integer")
        if mode not in {"min", "max"}:
            raise ValueError("mode must be 'min' or 'max'")
        self.capacity = capacity
        self.mode = mode
        self._values: list[np.float32] = [np.float32(-1e12 if mode == "min" else 1e12)]
        self._payloads: list[tuple[int, int, int, int] | None] = [None]

    @property
    def size(self) -> int:
        return len(self._values) - 1

    def __len__(self) -> int:
        return self.size

    def __bool__(self) -> bool:
        return self.size > 0

    def reset(self) -> None:
        self._values = self._values[:1]
        self._payloads = self._payloads[:1]

    def set_min_tree(self) -> None:
        self.mode = "min"
        self._values[0] = np.float32(-1e12)

    def set_max_tree(self) -> None:
        self.mode = "max"
        self._values[0] = np.float32(1e12)

    def add(self, value: float, x: int, y: int, z: int, label: int) -> None:
        priority = np.float32(value)
        payload = (_short(x), _short(y), _short(z), _short(label))
        hole = len(self._values)
        self._values.append(priority)
        self._payloads.append(payload)
        while hole > 1:
            parent = hole // 2
            parent_value = self._values[parent]
            moves = priority < parent_value if self.mode == "min" else priority > parent_value
            if not moves:
                break
            self._values[hole] = parent_value
            self._payloads[hole] = self._payloads[parent]
            hole = parent
        self._values[hole] = priority
        self._payloads[hole] = payload

    def peek_first(self) -> ToadsHeapItem:
        if not self:
            raise IndexError("cannot peek an empty TOADS heap")
        return self._item_at(1)

    def pop_first(self) -> ToadsHeapItem:
        if not self:
            raise IndexError("cannot pop an empty TOADS heap")
        first = self._item_at(1)
        last_value = self._values.pop()
        last_payload = self._payloads.pop()
        if not self:
            return first

        hole = 1
        current_size = self.size
        while hole * 2 <= current_size:
            child = hole * 2
            if child != current_size:
                left, right = self._values[child], self._values[child + 1]
                prefer_right = right < left if self.mode == "min" else right > left
                if prefer_right:
                    child += 1
            child_value = self._values[child]
            moves = child_value < last_value if self.mode == "min" else child_value > last_value
            if not moves:
                break
            self._values[hole] = child_value
            self._payloads[hole] = self._payloads[child]
            hole = child
        self._values[hole] = last_value
        self._payloads[hole] = last_payload
        return first

    def _item_at(self, index: int) -> ToadsHeapItem:
        payload = self._payloads[index]
        assert payload is not None
        return ToadsHeapItem(self._values[index], *payload)
