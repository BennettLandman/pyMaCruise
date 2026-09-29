from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.gvf import generalized_gvf_image
from pymacruise.plus_pipeline import PlusPipelineOptions, run_plus_gvf_stage


PHANTOM = Path(__file__).parent / "data/synthetic/plus_gvf/scalar.nii.gz"


def test_plus_gvf_stage_forwards_source_lambda_and_returns_ras_vector_field() -> None:
    scalar = nib.load(PHANTOM)
    options = PlusPipelineOptions()
    actual = run_plus_gvf_stage(scalar, options).vector_field
    expected = generalized_gvf_image(
        scalar, lambda_value=0.2, max_iterations=200
    ).vector_field
    base = generalized_gvf_image(scalar).vector_field

    assert options.gvf_lambda == 0.2
    assert actual.shape == (*scalar.shape, 3)
    assert nib.aff2axcodes(actual.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(np.asarray(actual.dataobj), np.asarray(expected.dataobj))
    assert not np.allclose(actual.get_fdata(), base.get_fdata(), rtol=1e-5, atol=1e-7)
