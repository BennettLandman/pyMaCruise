from pathlib import Path
import json

import nibabel as nib
import numpy as np
import pytest

from pymacruise.base_preprocessing import axial_windowed_sinc_ras, prepare_base_stripped_t1
from pymacruise.cli import main


FIXTURE = Path(__file__).parent / "data/synthetic/base_axial_sinc"


def test_base_axial_sinc_matches_independent_java_transform():
    image = nib.load(FIXTURE / "masked_t1.nii.gz")
    expected = np.load(FIXTURE / "java_expected_ras.npz")["values"]
    observed = axial_windowed_sinc_ras(np.asarray(image.dataobj), image.header.get_zooms())
    np.testing.assert_array_equal(observed, expected)
    assert np.any(observed < 0)  # The source explicitly disables clipping.
    assert observed.max() > np.asarray(image.dataobj).max()


def test_base_masking_precedes_sinc_and_keeps_its_external_support(tmp_path):
    image = nib.load(FIXTURE / "masked_t1.nii.gz")
    mask = nib.load(FIXTURE / "brain_mask.nii.gz")
    values = np.asarray(image.dataobj).copy()
    values[np.asarray(mask.dataobj) == 0] = 9999
    unstripped = nib.Nifti1Image(values, image.affine)
    result, toads_mask = prepare_base_stripped_t1(unstripped, mask)
    expected = np.load(FIXTURE / "java_expected_ras.npz")["values"]
    np.testing.assert_array_equal(result.dataobj, expected)
    assert np.count_nonzero(expected[np.asarray(mask.dataobj) == 0]) > 0
    np.testing.assert_array_equal(toads_mask.dataobj, np.ones(image.shape, dtype=np.uint8))
    np.testing.assert_array_equal(result.affine, image.affine)
    raw_path = tmp_path / "unstripped.nii.gz"
    nib.save(unstripped, raw_path)
    output = tmp_path / "output"
    provenance = output / "provenance.json"
    assert main(["base-preprocessing", "--t1-n4", str(raw_path), "--brain-mask",
                 str(FIXTURE / "brain_mask.nii.gz"), "--output-dir", str(output),
                 "--provenance", str(provenance)]) == 0
    np.testing.assert_array_equal(nib.load(output / "t1-n4.nii.gz").dataobj, expected)
    assert json.loads(provenance.read_text())["parameters"]["source_clip"] is False


def test_base_preprocessing_preserves_physical_orientation():
    image = nib.load(FIXTURE / "masked_t1.nii.gz")
    mask = nib.load(FIXTURE / "brain_mask.nii.gz")
    flip = np.array([[0, -1], [1, -1], [2, 1]])
    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        result, _ = prepare_base_stripped_t1(image.as_reoriented(flip), mask.as_reoriented(flip))
    np.testing.assert_array_equal(result.dataobj, np.load(FIXTURE / "java_expected_ras.npz")["values"])
    np.testing.assert_allclose(result.affine, image.affine, rtol=0, atol=1e-6)


def test_base_preprocessing_rejects_invalid_mask():
    image = nib.load(FIXTURE / "masked_t1.nii.gz")
    invalid = nib.Nifti1Image(np.full(image.shape, 0.5, dtype=np.float32), image.affine)
    with pytest.raises(ValueError, match="binary"):
        prepare_base_stripped_t1(image, invalid)
