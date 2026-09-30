from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.toads_thinning import toads_initialize_growing_queue


REFERENCE = Path(__file__).parent / "data/synthetic/toads_growing_queue/java_reference.json"


def _bits(value: np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.int32))


def _case_arrays(name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    size, classes = 7, 4
    skeleton = np.full((size, size, size), -1, dtype=np.int16)
    if name == "mixed":
        skeleton[3, 3, 2] = 0
        skeleton[3, 3, 4] = 1
        skeleton[2, 3, 3] = 2
        skeleton[4, 3, 3] = 3
    elif name == "single_seed":
        skeleton[3, 3, 3] = 0
    elif name == "no_boundary":
        skeleton[1:-1, 1:-1, 1:-1] = 0

    memberships = np.empty((size, size, size, classes), dtype=np.float32)
    ordering = np.empty((size, size, size), dtype=np.float32)
    for x in range(size):
        for y in range(size):
            for z in range(size):
                numerator = (x * 13 + y * 7 + z * 3) % 19
                ordering[x, y, z] = np.float32(-numerator) / np.float32(10.0)
                for k in range(classes):
                    numerator = (x * 5 + y * 7 + z * 11 + k * 3) % 29 + 1
                    memberships[x, y, z, k] = np.float32(numerator) / np.float32(30.0)
    return skeleton, memberships, ordering


def test_growing_queue_initialization_matches_java_min_heap_reference() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source"] == "VandyCRUISE MultipleToadSegmentation.initGrowingLabels"
    assert reference["source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert reference["heap_sha256"] == "9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a"
    assert reference["shape"] == [7, 7, 7, 4]
    assert {case["name"] for case in reference["cases"]} == {"mixed", "single_seed", "no_boundary"}

    for case in reference["cases"]:
        skeleton, memberships, ordering = _case_arrays(case["name"])
        available = np.zeros(skeleton.shape + (4,), dtype=np.bool_)
        result = toads_initialize_growing_queue(
            skeleton,
            memberships,
            ordering,
            initial_available=available,
        )
        assert result.queue.mode == "min"
        assert "".join("a" if v < 0 else str(int(v)) for v in result.skeleton.ravel()) == case["skeleton"]
        assert "".join("1" if v else "0" for v in result.available.ravel()) == case["available"]
        popped = []
        while result.queue:
            item = result.queue.pop_first()
            popped.append([_bits(item.value), item.x, item.y, item.z, item.label])
        assert popped == case["queue"], case["name"]

    cases = {case["name"]: case for case in reference["cases"]}
    assert len(cases["mixed"]["queue"]) >= 4
    assert len(cases["single_seed"]["queue"]) == 6
    assert cases["no_boundary"]["queue"] == []
