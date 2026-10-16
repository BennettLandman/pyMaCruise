from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.atlas_em_statistics import normalize_atlas_em_priors


ROOT = Path(__file__).parent
HARNESS = ROOT / "java_reference" / "AtlasEMPriorNormalizationReferenceHarness.java"


def test_atlas_prior_normalization_matches_java_excerpt(tmp_path: Path) -> None:
    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        pytest.skip("JDK is required for pinned Java method references")
    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run([javac, "-d", str(classes), str(HARNESS)], check=True, capture_output=True, text=True)
    output = subprocess.run([java, "-cp", str(classes), "AtlasEMPriorNormalizationReferenceHarness"], check=True, capture_output=True, text=True)
    reference = json.loads(output.stdout)

    priors = np.empty((3, 2, 2, 3), np.float32)
    mask = np.empty((3, 2, 2), bool)
    for x in range(3):
        for y in range(2):
            for z in range(2):
                mask[x, y, z] = (x + y + z) % 2 == 0
                for k in range(3):
                    priors[x, y, z, k] = np.float32(x * 3 + y * 2 + z + k) / np.float32(7)
    priors[0, 0, 0] = 0.0
    result = normalize_atlas_em_priors(priors, mask)

    assert reference["source_sha256"] == "b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a"
    actual = result.ravel().view(np.uint32).astype(np.uint64).tolist()
    assert actual == [int(value) & 0xFFFFFFFF for value in reference["priors"]]
    np.testing.assert_array_equal(result[~mask], priors[~mask])
    np.testing.assert_array_equal(result[0, 0, 0], [0.0, 0.0, 0.0])
