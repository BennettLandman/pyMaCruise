import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise import n4 as n4_module
from pymacruise.n4 import N4Config, n4_bias_correct
from pymacruise.n4_legacy_otsu import legacy_n4_otsu_mask


FIXTURE = Path(__file__).parent / "data/synthetic/n4_legacy_otsu"


@pytest.mark.parametrize("case", ["bimodal", "discrete", "zero_background"])
def test_historical_otsu_mask_matches_independent_itk(case):
    reference = np.load(FIXTURE / "reference.npz")
    metadata = json.loads((FIXTURE / "provenance.json").read_text())
    mask, threshold = legacy_n4_otsu_mask(reference[case])
    np.testing.assert_array_equal(mask, reference[case + "_mask"])
    assert threshold == metadata["cases"][case]["threshold"]


def test_legacy_n4_estimation_mask_does_not_limit_output_correction(monkeypatch):
    values = np.arange(1, 121, dtype=np.float32).reshape(4, 5, 6)
    # NIfTI-1 stores its affine as float32; use that representable grid here.
    affine = np.diag([0.7, 0.9, 1.2, 1]).astype(np.float32).astype(np.float64)
    image = nib.Nifti1Image(values, affine)
    captured = {}

    class Converted:
        def new_image_like(self, array):
            return array

    class BiasField:
        def numpy(self):
            return np.full(values.shape, 2, dtype=np.float32)

    class Backend:
        def n4_bias_field_correction(self, converted, **parameters):
            captured.update(parameters)
            return BiasField()

    monkeypatch.setattr(n4_module, "_nib_to_ants", lambda _: (Backend(), Converted()))
    result = n4_bias_correct(image, config=N4Config(estimation_mask="legacy-otsu"))
    assert captured["return_bias_field"] is True
    assert 0 < np.count_nonzero(captured["mask"]) < values.size
    np.testing.assert_array_equal(result.dataobj, values / 2)
    np.testing.assert_array_equal(result.affine, image.affine)
    assert result.extra["n4_estimation_mask"]["output_correction"] == "entire-image"


def test_legacy_n4_rejects_unsupported_rescaling():
    with pytest.raises(ValueError, match="does not rescale"):
        N4Config(estimation_mask="legacy-otsu", rescale_intensities=True)
