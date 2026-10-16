from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.atlas_em_outputs import generate_atlas_em_classification


ROOT = Path(__file__).parent
HARNESS = ROOT / "java_reference" / "AtlasEMClassificationReferenceHarness.java"


def test_atlas_em_classification_matches_java_excerpt(tmp_path: Path) -> None:
    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        pytest.skip("JDK is required for pinned Java method references")
    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run([javac, "-d", str(classes), str(HARNESS)], check=True, capture_output=True, text=True)
    output = subprocess.run([java, "-cp", str(classes), "AtlasEMClassificationReferenceHarness"], check=True, capture_output=True, text=True)
    reference = json.loads(output.stdout)

    memberships = np.empty((3, 2, 2, 3), np.float32)
    for x in range(3):
        for y in range(2):
            for z in range(2):
                for k in range(3):
                    value = 0.0 if (x + y + z + k) % 5 == 0 else np.float32(1 + x * 2 + y + z + k) / np.float32(9)
                    memberships[x, y, z, k] = value
    memberships[0, 0, 0] = (0.5, 0.5, 0.1)  # strict greater-than keeps class 0 on a tie
    memberships[1, 1, 1] = 0.0  # source emits label zero when all memberships are zero

    labels = generate_atlas_em_classification(memberships, np.asarray([4, 17, 91], np.int8))
    assert reference["source_sha256"] == "b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a"
    assert labels.ravel().tolist() == reference["labels"]


def test_atlas_em_classification_rejects_non_byte_template_labels() -> None:
    with pytest.raises(ValueError, match="signed-byte range"):
        generate_atlas_em_classification(np.ones((2, 2, 2, 1), np.float32), [256])
