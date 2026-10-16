from __future__ import annotations

import hashlib
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.atlas_em_preprocessing import transform_atlas_em_rigid_priors


ROOT = Path(__file__).parent
HARNESS = ROOT / "java_reference" / "AtlasEMRigidPriorReferenceHarness.java"
UPSTREAM = ROOT / "java_reference" / "upstream"
ATLAS_SOURCE_SHA256 = "96cb5321f9cc95a84d9115b6ea638e9d84aded422171c8792906fc2130bd415a"
PARAMETRIC_TRANSFORM_SHA256 = "22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac"
IMAGE_FUNCTIONS_SHA256 = "be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892"


def _shape_maps() -> np.ndarray:
    values = np.zeros((13, 11, 9, 2), dtype=np.float32)
    for k in range(2):
        for x in range(13):
            for y in range(11):
                for z in range(9):
                    if 2 <= x <= 10 and 1 <= y <= 9 and 1 <= z <= 7:
                        code = (x * 37 + y * 23 + z * 17 + k * 29) % 101
                        values[x, y, z, k] = np.float32(
                            np.float32(0.08) + np.float32(code) / np.float32(125.0)
                        )
    return values


def test_rigid_atlas_priors_match_selected_java_source_excerpt(tmp_path: Path) -> None:
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
        [java, "-cp", str(classes), "AtlasEMRigidPriorReferenceHarness"],
        check=True,
        capture_output=True,
        text=True,
    )
    import json

    reference = json.loads(result.stdout.splitlines()[-1])
    assert reference["source"] == (
        "DemonToadDeformableAtlas rigid transformed shape excerpt + "
        "AtlasEM computeAtlas normalization"
    )
    assert reference["atlas_source_sha256"] == ATLAS_SOURCE_SHA256
    assert reference["parametric_transform_sha256"] == PARAMETRIC_TRANSFORM_SHA256
    assert reference["image_functions_sha256"] == IMAGE_FUNCTIONS_SHA256
    assert reference["active_voxels"] == sum(
        (x * 3 + y * 5 + z * 7) % 9 != 0
        for x in range(15)
        for y in range(13)
        for z in range(11)
    )

    mask = np.fromfunction(
        lambda x, y, z: (x * 3 + y * 5 + z * 7) % 9 != 0,
        (15, 13, 11),
        dtype=int,
    ).astype(bool)
    actual = transform_atlas_em_rigid_priors(
        _shape_maps(),
        mask,
        image_dimensions_xyz=(15, 13, 11),
        image_resolutions_mm=(0.9, 1.1, 1.3),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        transform=(0.07, -0.04, 0.03, 0.25, -0.30, 0.15),
    )
    actual_sha = hashlib.sha256(actual.astype(">f4", copy=False).tobytes(order="C")).hexdigest()
    assert actual_sha == reference["output_sha256"]

    assert hashlib.sha256((UPSTREAM / "edu/jhmi/rad/medic/methods/ParametricTransform.java").read_bytes()).hexdigest() == PARAMETRIC_TRANSFORM_SHA256
    assert hashlib.sha256((UPSTREAM / "edu/jhmi/rad/medic/libraries/ImageFunctions.java").read_bytes()).hexdigest() == IMAGE_FUNCTIONS_SHA256
