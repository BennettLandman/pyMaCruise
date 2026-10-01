"""Pin the synthetic identity-grid contract and Java oracle for P39zzzzzzag."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "lesion_standard13_voxel"
HARNESS = ROOT / "java_reference" / "LesionToadsStandard13VoxelReferenceHarness.java"
HARNESS_SHA256 = "a48def0d650e94ec53c327f4663534096f5243ea48bc9d7aa39d8619d2914d3c"


def test_standard13_voxel_fixture_is_identity_ras_and_java_harness_is_pinned():
    image = nib.load(str(FIXTURE / "inputs_ras.nii.gz"))
    assert image.shape == (3, 3, 3, 2)
    np.testing.assert_array_equal(image.affine, np.eye(4))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")

    values = np.asarray(image.dataobj, dtype=np.float32)
    reference = json.loads((FIXTURE / "java_reference.json").read_text())
    image_sha = hashlib.sha256(values.astype("<f4", copy=False).tobytes()).hexdigest()
    assert image_sha == reference["image_sha256"]
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert reference["source_sha256"] == (
        "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    )
