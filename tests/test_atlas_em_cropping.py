from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.atlas_em_preprocessing import (
    crop_atlas_em_images,
    find_atlas_em_crop_bounds,
    uncrop_atlas_em_images,
)


ROOT = Path(__file__).parent
HARNESS = ROOT / "java_reference" / "AtlasEMCropReferenceHarness.java"
UPSTREAM = ROOT / "java_reference" / "upstream"


def _synthetic_images() -> np.ndarray:
    values = np.zeros((13, 11, 9, 2), dtype=np.float32)
    for x in range(13):
        for y in range(11):
            for z in range(9):
                values[x, y, z, 0] = np.float32(
                    np.float32(0.02) + np.float32((x * 7 + y * 5 + z * 3) % 13) * np.float32(0.001)
                )
                values[x, y, z, 1] = np.float32(
                    np.float32(0.01) + np.float32((x * 3 + y * 11 + z * 2) % 9) * np.float32(0.001)
                )
                if 3 <= x <= 9 and 2 <= y <= 8 and 1 <= z <= 7:
                    values[x, y, z, 0] = np.float32(
                        np.float32(0.4) + np.float32(x + y + z) * np.float32(0.01)
                    )
                if 1 <= x <= 4 and 5 <= y <= 9 and 3 <= z <= 7:
                    values[x, y, z, 1] = np.float32(
                        np.float32(0.35) + np.float32(x + y + z) * np.float32(0.012)
                    )
    return values


def _sha(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype=">f4").tobytes(order="C")).hexdigest()


def _bits(values: np.ndarray) -> list[int]:
    return [int(value) for value in np.asarray(values, dtype=np.float32).view(np.uint32)]


def test_atlas_em_crop_and_uncrop_match_java_excerpt(tmp_path: Path) -> None:
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
        [java, "-cp", str(classes), "AtlasEMCropReferenceHarness"],
        check=True,
        capture_output=True,
        text=True,
    )
    reference = json.loads(result.stdout)
    assert reference["atlas_em_preprocess_sha256"] == "acec23f916b3134e53381e83ff3de7ca569f785e6b866bfe135fb409c2926d2b"
    rotation_path = UPSTREAM / "edu/jhmi/rad/medic/structures/RotationMatrix.java"
    assert hashlib.sha256(rotation_path.read_bytes()).hexdigest() == reference["rotation_matrix_sha256"]

    images = _synthetic_images()
    bounds = find_atlas_em_crop_bounds(
        images,
        [0.1, 0.1],
        [1.1, 1.1],
        background_ratio=0.15,
    )
    assert list(bounds[0]) == reference["lower"]
    assert list(bounds[1]) == reference["upper"]
    cropped = crop_atlas_em_images(
        images,
        *bounds,
        border=2,
        signed_resolutions_mm=(0.8, -1.1, 1.4),
        transform=(0.12, -0.08, 0.05, 0.4, -0.3, 0.2),
    )
    assert list(cropped.images.shape[:3]) == reference["cropped_dims"]
    assert _sha(cropped.images) == reference["crop_sha256"]
    assert _bits(cropped.transform) == reference["crop_transform_bits"]

    uncropped, restored_transform = uncrop_atlas_em_images(cropped)
    assert _sha(uncropped) == reference["uncrop_sha256"]
    assert _bits(restored_transform) == reference["restored_transform_bits"]
