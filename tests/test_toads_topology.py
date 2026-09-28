from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_topology import toads_homeomorphic_label_change


REFERENCE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_homeomorphic_change"
    / "java_reference.json"
)


def test_homeomorphic_change_matches_selected_java_source_excerpt() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source"] == (
        "VandyCRUISE MultipleToadSegmentation.isHomeomorphicSegmentation excerpt"
    )
    assert reference["source_sha256"] == (
        "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    )
    assert reference["source_lines"] == "1554-1673"
    assert reference["lut_source_sha256"] == (
        "50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f"
    )
    assert reference["check_composed"] is False
    assert len(reference["cases"]) == 130

    for case in reference["cases"]:
        neighborhood = np.fromiter(
            (int(value) for value in case["neighborhood"]),
            dtype=np.int16,
            count=27,
        ).reshape((3, 3, 3))
        assert toads_homeomorphic_label_change(
            neighborhood, case["candidate"]
        ) is case["allowed"]


def test_homeomorphic_change_accepts_unchanged_center_before_topology_checks() -> None:
    # The selected source returns immediately when the center already has k.
    neighborhood = np.zeros((3, 3, 3), dtype=np.int16)
    neighborhood[1, 1, 1] = 1
    assert toads_homeomorphic_label_change(neighborhood, 1)


@pytest.mark.parametrize(
    ("neighborhood", "candidate", "message"),
    [
        (np.zeros((3, 3), dtype=np.int16), 1, "shape"),
        (np.zeros((3, 3, 3), dtype=np.float32), 1, "integers"),
        (np.full((3, 3, 3), 32768, dtype=np.int32), 1, "signed short"),
        (np.zeros((3, 3, 3), dtype=np.int16), 32768, "signed short"),
    ],
)
def test_homeomorphic_change_rejects_invalid_input(
    neighborhood: np.ndarray, candidate: int, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        toads_homeomorphic_label_change(neighborhood, candidate)
