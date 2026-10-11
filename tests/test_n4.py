from __future__ import annotations

import nibabel as nib
import numpy as np
import pytest

from pymacruise.n4 import N4Config, n4_bias_correct


def test_n4_reduces_a_smooth_synthetic_bias_and_preserves_ras_grid() -> None:
    pytest.importorskip("ants")
    size = 32
    axis = np.linspace(-1.0, 1.0, size, dtype=np.float32)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    known_bias = np.exp(0.25 * x + 0.16 * y - 0.12 * z)
    data = (100.0 * known_bias).astype(np.float32)
    affine = np.diag([1.1, 1.2, 1.3, 1.0])
    source = nib.Nifti1Image(data, affine)

    # Keep this algorithm-quality check independent of the source wrapper's
    # CLI output-range rescaling, which is a separate intensity post-step.
    corrected = n4_bias_correct(source, config=N4Config(rescale_intensities=False))
    result = np.asarray(corrected.dataobj)
    before_cv = float(np.std(data) / np.mean(data))
    after_cv = float(np.std(result) / np.mean(result))

    assert before_cv > 0.1
    assert after_cv < before_cv * 0.2
    assert corrected.shape == source.shape
    assert nib.aff2axcodes(corrected.affine) == ("R", "A", "S")
    np.testing.assert_allclose(corrected.affine, affine, atol=1e-6)
    assert np.isfinite(result).all()
    repeated = np.asarray(
        n4_bias_correct(source, config=N4Config(rescale_intensities=False)).dataobj
    )
    np.testing.assert_allclose(repeated, result, rtol=0, atol=1e-5)


def test_n4_config_and_input_contract_errors() -> None:
    assert N4Config().rescale_intensities is True
    with pytest.raises(ValueError, match="shrink factor"):
        N4Config(shrink_factor=0)
    with pytest.raises(ValueError, match="iteration counts"):
        N4Config(iterations=(0,))

    negative = nib.Nifti1Image(-np.ones((4, 4, 4), dtype=np.float32), np.eye(4))
    with pytest.raises(ValueError, match="negative"):
        n4_bias_correct(negative)

    four_d = nib.Nifti1Image(np.ones((4, 4, 4, 2), dtype=np.float32), np.eye(4))
    with pytest.raises(ValueError, match="3D"):
        n4_bias_correct(four_d)
