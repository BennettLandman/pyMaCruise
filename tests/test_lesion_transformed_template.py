from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_transformed_template import transform_lesion_toads_template


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "lesion_transformed_template"
HARNESS = ROOT / "java_reference" / "LesionToadAtlasTransformedTemplateReferenceHarness.java"
REFERENCE = FIXTURE / "java_reference.json"
SOURCE_HASHES = {
    "atlas": "718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e",
    "transform": "22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac",
    "image_functions": "be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892",
    "numerics": "13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb",
    "rotation_matrix": "144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0",
}


def test_transformed_template_matches_selected_java_rigid_mapping() -> None:
    ref = json.loads(REFERENCE.read_text())
    assert ref["source_method"] == "LesionToadAtlas.computeTransformedTemplate"
    assert ref["source_sha256"] == SOURCE_HASHES
    assert ref["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    image = load_nifti_ras(FIXTURE / "template_labels.nii.gz")
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), (0.9, 1.0, 1.1), atol=1e-6)
    atlas = np.asarray(image.dataobj)
    assert atlas.shape == tuple(ref["template_shape"])
    assert atlas[0, 0, 0] == 42 and ref["background_label"] == 1
    before = atlas.copy()

    output = transform_lesion_toads_template(
        atlas,
        (0.02, -0.03, 0.01, 0.3, -0.15, 0.4),
        tuple(ref["image_shape"]),
        image_resolution_mm=(1.1, 1.3, 0.9),
        template_resolution_mm=(0.9, 1.0, 1.1),
        background_label=1,
    )
    assert output.shape == tuple(ref["image_shape"])
    np.testing.assert_array_equal(output.ravel(order="C"), np.asarray(ref["output_labels"], dtype=np.int8))
    np.testing.assert_array_equal(atlas, before)
    assert np.count_nonzero(output == 1) > 0
    assert not np.any(output == 42)  # the corner sentinel is not the Java fill value
    assert np.count_nonzero(output != 1) > 0


def test_transformed_template_validates_rigid_inputs() -> None:
    with np.testing.assert_raises_regex(ValueError, "six finite float32"):
        transform_lesion_toads_template(
            np.ones((3, 3, 3), dtype=np.uint8),
            np.zeros(5, dtype=np.float32),
            (3, 3, 3),
            image_resolution_mm=(1, 1, 1),
            template_resolution_mm=(1, 1, 1),
            background_label=1,
        )
