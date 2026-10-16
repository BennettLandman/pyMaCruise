from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.mipav_morphology import mipav_fill_holes_25d

ROOT = Path(__file__).parent
HARNESS = ROOT / "java_reference" / "MipavMorphology25DFillHolesReferenceHarness.java"
SOURCE_SHA = "d3e271151bce7e635777764648ef21c0253485699e577897c4ff4821f75c70f6"


def _phantom() -> np.ndarray:
    x, y, z = np.indices((8, 7, 3))
    image = (((x * 7 + y * 11 + z * 13) % 9) < 5).astype(np.uint8)
    image[1:-1, 1:-1, 0] = 1
    image[2:6, 2:5, 1] = 1
    image[3, 3, 1] = 0
    image[4, 3, 1] = 0
    image[1, 1, 2] = 1
    return image


def test_mipav_fill_holes_25d_matches_selected_java_method_excerpt(tmp_path: Path) -> None:
    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        pytest.skip("JDK is required for pinned Java method references")
    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run([javac, "-d", str(classes), str(HARNESS)], check=True, capture_output=True, text=True)
    result = subprocess.run([java, "-cp", str(classes), "MipavMorphology25DFillHolesReferenceHarness"], check=True, capture_output=True, text=True)
    reference = json.loads(result.stdout)
    assert reference["source_sha256"] == SOURCE_SHA
    expected = np.asarray(reference["data"], dtype=bool).reshape((8, 7, 3))
    actual = mipav_fill_holes_25d(_phantom())
    np.testing.assert_array_equal(actual, expected)
    # MIPAV clears all 2D slice borders before flooding, including border foreground.
    assert not np.any(actual[[0, -1], :, :])
    assert not np.any(actual[:, [0, -1], :])


def test_mipav_fill_holes_25d_rejects_bad_grids_and_nan() -> None:
    with pytest.raises(ValueError, match="3D binary grid"):
        mipav_fill_holes_25d(np.ones((4, 4)))
    bad = np.ones((3, 3, 3), np.float32)
    bad[1, 1, 1] = np.nan
    with pytest.raises(ValueError, match="finite"):
        mipav_fill_holes_25d(bad)
