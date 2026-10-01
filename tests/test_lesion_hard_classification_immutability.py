"""Input-preservation assertion for the selected hard-label map adapter."""

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_hard_classification import generate_lesion_toads_classification


FIXTURE_DIR = Path(__file__).parent / "data" / "synthetic" / "lesion_hard_classification"


def test_hard_classification_leaves_synthetic_ras_inputs_unchanged():
    reference = json.loads((FIXTURE_DIR / "java_reference.json").read_text())
    image = load_nifti_ras(FIXTURE_DIR / "memberships_and_lesion.nii.gz")
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    values = np.asarray(image.dataobj, dtype=np.float32)
    memberships = values[..., :3]
    lesion = values[..., 3]
    memberships_before = memberships.copy()
    lesion_before = lesion.copy()

    result = generate_lesion_toads_classification(
        reference["template_labels"], memberships, lesion, reference["lesion_class"]
    )

    np.testing.assert_array_equal(result.ravel(), reference["labels"])
    np.testing.assert_array_equal(memberships, memberships_before)
    np.testing.assert_array_equal(lesion, lesion_before)
