from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.atlas_em_outputs import generate_atlas_em_dura_inputs

ROOT = Path(__file__).parent
HARNESS = ROOT / "java_reference" / "AtlasEMDuraOutputsReferenceHarness.java"
SOURCE_SHA = "b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a"
NAMES = ["Background", "SulcalCSF", "CorticalGM", "OtherA", "Cerebellum-GM", "OtherB"]


def _inputs():
    x, y, z = np.indices((7, 7, 2))
    mask = (x + y + z) % 5 != 0
    labels = ((x * 3 + y + z * 2) % len(NAMES)).astype(np.int8)
    memberships = np.empty((7, 7, 2, len(NAMES)), np.float32)
    for k in range(len(NAMES)):
        memberships[..., k] = np.float32(1 + x * 3 + y * 5 + z * 7 + k * 11) / np.float32(101.0)
    return memberships, labels, mask


def test_atlas_em_dura_inputs_match_selected_java_output_method(tmp_path: Path) -> None:
    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        pytest.skip("JDK is required for pinned Java method references")
    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run([javac, "-d", str(classes), str(HARNESS)], check=True, capture_output=True, text=True)
    result = subprocess.run([java, "-cp", str(classes), "AtlasEMDuraOutputsReferenceHarness"], check=True, capture_output=True, text=True)
    reference = json.loads(result.stdout)
    assert reference["source_sha256"] == SOURCE_SHA
    memberships, labels, mask = _inputs()
    before = memberships.copy()
    products = generate_atlas_em_dura_inputs(memberships, labels, mask, NAMES)
    actual_bits = products.ravel().view(np.uint32).astype(np.uint64).tolist()
    assert actual_bits == reference["data"]
    assert products.shape == (7, 7, 2, 3)
    assert not np.shares_memory(products, memberships)
    np.testing.assert_array_equal(memberships, before)


def test_atlas_em_dura_inputs_validate_class_count_and_grids() -> None:
    memberships, labels, mask = _inputs()
    with pytest.raises(ValueError, match="one string per membership class"):
        generate_atlas_em_dura_inputs(memberships, labels, mask, NAMES[:-1])
    with pytest.raises(ValueError, match="integer grid"):
        generate_atlas_em_dura_inputs(memberships, labels.astype(np.float32), mask, NAMES)
