from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_distance import lesion_toads_distance_to_mask


FIXTURES = Path(__file__).parent / "data" / "synthetic" / "lesion_toads_distance"


def _mask_from_spec(case: dict[str, object]) -> np.ndarray:
    mask = np.zeros(tuple(case["shape"]), dtype=bool)
    points = case["object_voxels"]
    if points == "all":
        mask[...] = True
    else:
        for x, y, z in points:
            mask[x, y, z] = True
    return mask


def test_lesion_toads_capped_distance_matches_java_reference() -> None:
    spec = json.loads((FIXTURES / "spec.json").read_text())
    reference = json.loads((FIXTURES / "java_reference.json").read_text())
    assert reference["source"].endswith(
        "LesionToadSegmentation.fastDistance(boolean[][][], byte, int, int, int)"
    )
    assert [case["name"] for case in spec["cases"]] == [
        case["name"] for case in reference["cases"]
    ]

    for source_case, java_case in zip(spec["cases"], reference["cases"], strict=True):
        assert java_case["shape"] == source_case["shape"]
        assert java_case["max_distance"] == source_case["max_distance"]
        assert java_case["object_voxels"] == (
            [
                [x, y, z]
                for x in range(source_case["shape"][0])
                for y in range(source_case["shape"][1])
                for z in range(source_case["shape"][2])
            ]
            if source_case["object_voxels"] == "all"
            else source_case["object_voxels"]
        )
        actual = lesion_toads_distance_to_mask(
            _mask_from_spec(source_case), source_case["max_distance"]
        )
        expected = np.asarray(java_case["distance_xyz_flat"], dtype=np.int8).reshape(
            source_case["shape"], order="C"
        )
        np.testing.assert_array_equal(actual, expected)


def test_lesion_toads_distance_uses_six_connected_voxel_steps_and_inclusive_cap() -> None:
    diagonal = np.zeros((3, 3, 3), dtype=bool)
    diagonal[0, 0, 0] = True
    distance = lesion_toads_distance_to_mask(diagonal, 3)
    assert distance[1, 0, 0] == 1
    assert distance[1, 1, 0] == 2
    assert distance[1, 1, 1] == 3

    capped = np.zeros((7, 5, 5), dtype=bool)
    capped[3, 2, 2] = True
    distance = lesion_toads_distance_to_mask(capped, 1)
    assert distance[3, 2, 2] == 0
    assert distance[2, 2, 2] == 1
    assert distance[1, 2, 2] == 2
    assert distance[0, 2, 2] == -1


@pytest.mark.parametrize(
    ("mask", "max_distance", "exception", "message"),
    [
        (np.zeros((2, 3), dtype=bool), 1, ValueError, "3D"),
        (np.zeros((2, 3, 0), dtype=bool), 1, ValueError, "nonempty"),
        (np.zeros((2, 3, 4), dtype=np.uint8), 1, TypeError, "boolean"),
        (np.zeros((2, 3, 4), dtype=bool), 1.5, TypeError, "integer"),
        (np.zeros((2, 3, 4), dtype=bool), 0, ValueError, "1..127"),
        (np.zeros((2, 3, 4), dtype=bool), 128, ValueError, "1..127"),
    ],
)
def test_lesion_toads_distance_validates_inputs(
    mask: np.ndarray, max_distance: object, exception: type[Exception], message: str
) -> None:
    with pytest.raises(exception, match=message):
        lesion_toads_distance_to_mask(mask, max_distance)
