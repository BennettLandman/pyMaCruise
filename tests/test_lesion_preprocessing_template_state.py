from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_preprocessing_template_state import (
    LesionToadsTemplateState,
    update_lesion_transformed_template,
)


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "lesion_preprocess_template_state"
HARNESS = ROOT / "java_reference" / "LesionToadPreprocessUpdateTransformedTemplateReferenceHarness.java"
REFERENCE = FIXTURE / "java_reference.json"
SOURCE_HASHES = {
    "preprocess": "72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08",
    "wrapper": "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb",
}


def test_update_transformed_template_matches_java_reference() -> None:
    ref = json.loads(REFERENCE.read_text())
    assert ref["source_method"] == "LesionToadPreprocess.updateTransformedTemplate"
    assert ref["source_sha256"] == SOURCE_HASHES
    assert ref["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    image = load_nifti_ras(FIXTURE / "transformed_template.nii.gz")
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), (0.8, 1.0, 1.2), atol=1e-6)
    template = np.asarray(image.dataobj, dtype=np.int8)
    assert template.shape == tuple(ref["shape"])
    before = template.copy()

    state = LesionToadsTemplateState(topology=np.zeros_like(template), transformed=False)
    update_lesion_transformed_template(state, template)

    assert ref["transformed"] is True
    assert ref["same_reference"] is True
    assert ref["same_contents"] is True
    assert state.topology is template
    assert state.transformed is True
    np.testing.assert_array_equal(template.ravel(order="C"), np.asarray(ref["labels"], dtype=np.int8))
    np.testing.assert_array_equal(template, before)
