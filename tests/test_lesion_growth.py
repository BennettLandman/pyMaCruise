from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_growth import initialize_lesion_toads_growing_frontier


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_growing_frontier/java_reference.json"


def _case_arrays(name: str) -> tuple[np.ndarray, ...]:
    size, classes = 7, 4
    skeleton = np.full((size, size, size), -1, dtype=np.int8)
    memberships = np.empty((size, size, size, classes), dtype=np.float32)
    ordering = np.empty((size, size, size), dtype=np.float32)
    lesion = np.empty((size, size, size), dtype=np.float32)
    black_hole = np.empty((size, size, size), dtype=np.float32)
    for x in range(size):
        for y in range(size):
            for z in range(size):
                order_numerator = (x * 13 + y * 7 + z * 3) % 19
                ordering[x, y, z] = np.float32(-order_numerator) / np.float32(10.0)
                for k in range(classes):
                    numerator = (x * 5 + y * 7 + z * 11 + k * 3) % 29 + 1
                    memberships[x, y, z, k] = (
                        np.float32(0.0)
                        if name == "zero_ties"
                        else np.float32(numerator) / np.float32(30.0)
                    )
                lesion[x, y, z] = (
                    np.float32(0.0)
                    if name == "zero_ties"
                    else np.float32((x + 3 * y + 5 * z) % 7 + 1) / np.float32(10.0)
                )
                black_hole[x, y, z] = (
                    np.float32(0.0)
                    if name == "zero_ties"
                    else np.float32((2 * x + y + 3 * z) % 9 + 1) / np.float32(12.0)
                )
    skeleton[3, 3, 2] = 0
    skeleton[3, 3, 4] = 1
    skeleton[2, 3, 3] = 2
    skeleton[4, 3, 3] = 3
    return skeleton, memberships, ordering, lesion, black_hole


def _bits(value: np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.int32))


def test_growing_frontier_matches_selected_java_formula_and_heap() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert reference["numerics_sha256"] == "13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb"
    assert reference["heap_sha256"] == "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
    assert {case["name"] for case in reference["cases"]} == {
        "additive",
        "max",
        "zero_ties",
        "optimized_types",
    }

    for case in reference["cases"]:
        skeleton, memberships, ordering, lesion, black_hole = _case_arrays(case["name"])
        result = initialize_lesion_toads_growing_frontier(
            skeleton,
            memberships,
            ordering,
            tuple(case["object_types"]),
            lesion,
            black_hole,
            use_lesion_max=case["use_lesion_max"],
        )
        available = "".join("1" if value else "0" for value in result.available.ravel())
        assert available == case["available"], case["name"]
        popped = []
        while result.queue:
            item = result.queue.pop_first()
            popped.append([_bits(item.value), item.x, item.y, item.z, item.label])
        assert popped == case["queue"], case["name"]


def test_growing_frontier_preserves_interior_reset_and_border_state() -> None:
    skeleton = np.full((5, 5, 5), -1, dtype=np.int8)
    skeleton[2, 2, 1] = 0
    available = np.zeros((5, 5, 5, 2), dtype=np.bool_)
    available[0, 2, 2, 0] = False
    memberships = np.ones((5, 5, 5, 2), dtype=np.float32)
    zeros = np.zeros((5, 5, 5), dtype=np.float32)
    result = initialize_lesion_toads_growing_frontier(
        skeleton, memberships, zeros, ("obj", "lesion"), zeros, zeros,
        initial_available=available,
    )
    assert not result.available[0, 2, 2, 0]
    assert not result.available[2, 2, 2, 0]
    assert result.available[2, 2, 2, 1]
    assert len(result.queue) == 5


def test_growing_frontier_rejects_true_border_availability() -> None:
    skeleton = np.full((3, 3, 3), -1, dtype=np.int8)
    memberships = np.ones((3, 3, 3, 2), dtype=np.float32)
    zeros = np.zeros((3, 3, 3), dtype=np.float32)
    available = np.zeros((3, 3, 3, 2), dtype=np.bool_)
    available[0, 1, 1, 0] = True
    with pytest.raises(ValueError, match="border availability"):
        initialize_lesion_toads_growing_frontier(
            skeleton, memberships, zeros, ("obj", "lesion"), zeros, zeros,
            initial_available=available,
        )


@pytest.mark.parametrize(
    ("skeleton", "memberships", "types", "ordering", "lesion", "black_hole", "message"),
    [
        (np.zeros((2, 3, 3), dtype=np.int16), np.zeros((2, 3, 3, 2)), ("obj", "lesion"), np.zeros((2, 3, 3)), np.zeros((2, 3, 3)), np.zeros((2, 3, 3)), "every dimension"),
        (np.zeros((3, 3, 3), dtype=np.int16), np.zeros((3, 3, 3, 2)), ("obj",), np.zeros((3, 3, 3)), np.zeros((3, 3, 3)), np.zeros((3, 3, 3)), "one nonempty object type"),
        (np.zeros((3, 3, 3), dtype=np.int16), np.zeros((3, 3, 3, 2)), ("obj", "lesion"), np.zeros((3, 3, 2)), np.zeros((3, 3, 3)), np.zeros((3, 3, 3)), "ordering must match"),
        (np.zeros((3, 3, 3), dtype=np.int16), np.zeros((3, 3, 3, 2)), ("obj", "lesion"), np.zeros((3, 3, 3)), np.zeros((2, 3, 3)), np.zeros((3, 3, 3)), "memberships must match skeleton"),
    ],
)
def test_growing_frontier_rejects_invalid_inputs(
    skeleton: np.ndarray,
    memberships: np.ndarray,
    types: tuple[str, ...],
    ordering: np.ndarray,
    lesion: np.ndarray,
    black_hole: np.ndarray,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        initialize_lesion_toads_growing_frontier(
            skeleton, memberships, ordering, types, lesion, black_hole
        )
