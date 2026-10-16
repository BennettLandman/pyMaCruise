from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.atlas_em_registration import build_atlas_em_rigid_registration_pyramid


ROOT = Path(__file__).parent
HARNESS = ROOT / "java_reference" / "AtlasEMRigidPyramidReferenceHarness.java"
UPSTREAM = ROOT / "java_reference" / "upstream"


def _field(salt: int) -> np.ndarray:
    values = np.empty((8, 8, 8, 3), dtype=np.float32)
    for x in range(8):
        for y in range(8):
            for z in range(8):
                for klass in range(3):
                    values[x, y, z, klass] = np.float32(
                        np.float32((x * 17 + y * 11 + z * 7 + klass * 13 + salt) % 97)
                        / np.float32(32.0)
                    )
    return values


def _bits(values: np.ndarray) -> list[int]:
    return [int(value) for value in np.asarray(values, dtype=np.float32).ravel(order="C").view(np.uint32)]


def test_atlas_em_rigid_registration_pyramid_matches_java_subsampling_excerpt(
    tmp_path: Path,
) -> None:
    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        pytest.skip("JDK is required for pinned Java method references")
    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run(
        [javac, "-sourcepath", str(UPSTREAM), "-d", str(classes), str(HARNESS)],
        check=True,
        capture_output=True,
        text=True,
    )
    result = subprocess.run(
        [java, "-cp", str(classes), "AtlasEMRigidPyramidReferenceHarness"],
        check=True,
        capture_output=True,
        text=True,
    )
    reference = json.loads(result.stdout)
    assert reference["demon_atlas_source_sha256"] == (
        "96cb5321f9cc95a84d9115b6ea638e9d84aded422171c8792906fc2130bd415a"
    )
    assert reference["image_functions_sha256"] == (
        "be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892"
    )

    actual = build_atlas_em_rigid_registration_pyramid(
        _field(3), _field(19), levels=3
    )
    assert [level.scale for level in actual] == [level["scale"] for level in reference["levels"]]
    for level, expected in zip(actual, reference["levels"], strict=True):
        assert list(level.memberships.shape[:3]) == expected["dims"]
        assert _bits(level.memberships) == expected["membership_bits"]
        assert _bits(level.atlas_shapes) == expected["shape_bits"]
