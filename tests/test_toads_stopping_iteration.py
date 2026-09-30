from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_iteration import toads_stopping_iteration


GROWING_REFERENCE = (
    Path(__file__).parent / "data/synthetic/toads_growing_pass/java_reference.json"
)
CONTROL_REFERENCE = (
    Path(__file__).parent / "data/synthetic/toads_iteration_control/java_reference.json"
)


def _inputs():
    size, classes = 9, 4
    segmentation = np.zeros((size, size, size), dtype=np.int16)
    skeleton = np.full_like(segmentation, -1)
    ordering = np.empty((size, size, size), dtype=np.float32)
    memberships = np.empty((size, size, size, classes), dtype=np.float32)
    for x in range(size):
        for y in range(size):
            for z in range(size):
                ordering[x, y, z] = np.float32(-((x * 13 + y * 7 + z * 3) % 19)) / np.float32(10.0)
                for k in range(classes):
                    memberships[x, y, z, k] = np.float32(
                        ((x * 5 + y * 7 + z * 11 + k * 3) % 29 + 1)
                    ) / np.float32(30.0)
                center = 2 + ((y + z) % 3)
                if 2 <= y <= 6 and 2 <= z <= 6 and abs(x - center) <= 1:
                    segmentation[x, y, z] = 1
                if 5 <= x <= 6 and 2 <= y <= 3 and 5 <= z <= 6:
                    segmentation[x, y, z] = 2
                if x == 1 and 2 <= y <= 6 and z == 7:
                    segmentation[x, y, z] = 3
    center = 2 + ((4 + 3) % 3)
    skeleton[center, 4, 3] = 1
    skeleton[5, 2, 5] = 2
    volumes = np.asarray([1, 2, 3, 4], dtype=np.float32)
    return segmentation, skeleton, memberships, ordering, volumes


def test_stopping_iteration_composes_java_referenced_final_growth_and_label_transfer() -> None:
    reference = json.loads(GROWING_REFERENCE.read_text())
    control_reference = json.loads(CONTROL_REFERENCE.read_text())
    assert reference["source_sha256"] == "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert control_reference["source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    segmentation, skeleton, memberships, ordering, volumes = _inputs()
    result = toads_stopping_iteration(
        0.9,
        0.8,
        3,
        3,
        1e-5,
        segmentation,
        skeleton,
        memberships,
        ordering,
        volumes,
        (-10, 20, 30, 40),
        np.full(segmentation.shape, 99, dtype=np.int8),
        do_growing=True,
    )

    assert result.control.stop and result.control.run_growing
    assert result.control.update_classification and not result.control.refresh_relations
    assert result.growing is not None
    assert "".join(str(int(value)) for value in result.growing.segmentation.ravel()) == reference["segmentation"]
    assert "".join("a" if value < 0 else str(int(value)) for value in result.growing.skeleton.ravel()) == reference["skeleton"]
    expected_classification = np.asarray(
        [-10, 20, 30, 40], dtype=np.int8
    )[result.growing.segmentation]
    np.testing.assert_array_equal(result.classification, expected_classification)
    assert result.changed_classification_voxels == segmentation.size


def test_stopping_iteration_skips_final_growth_when_disabled() -> None:
    segmentation, skeleton, memberships, ordering, volumes = _inputs()
    previous = np.zeros(segmentation.shape, dtype=np.int8)
    result = toads_stopping_iteration(
        0.9,
        0.8,
        3,
        3,
        1e-5,
        segmentation,
        skeleton,
        memberships,
        ordering,
        volumes,
        (0, 11, 22, 33),
        previous,
        do_growing=False,
    )

    assert result.control.stop and not result.control.run_growing
    assert result.growing is None
    expected = np.asarray((0, 11, 22, 33), dtype=np.int8)[segmentation]
    np.testing.assert_array_equal(result.classification, expected)
    assert result.changed_classification_voxels == int(np.count_nonzero(expected))


def test_stopping_iteration_rejects_continuation_branch() -> None:
    segmentation, skeleton, memberships, ordering, volumes = _inputs()
    with pytest.raises(ValueError, match="requires a selected stop condition"):
        toads_stopping_iteration(
            1.0,
            0.5,
            1,
            5,
            1e-5,
            segmentation,
            skeleton,
            memberships,
            ordering,
            volumes,
            (0, 1, 2, 3),
            np.zeros(segmentation.shape, dtype=np.int8),
            do_growing=True,
        )
