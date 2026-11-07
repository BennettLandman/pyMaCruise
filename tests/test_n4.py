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
    assert N4Config().rescale_intensities is False
    with pytest.raises(ValueError, match="shrink factor"):
        N4Config(shrink_factor=0)
    with pytest.raises(ValueError, match="iteration counts"):
        N4Config(iterations=(0,))

    negative = nib.Nifti1Image(-np.ones((4, 4, 4), dtype=np.float32), np.eye(4))
    with pytest.raises(ValueError, match="at least one positive"):
        n4_bias_correct(negative)

    four_d = nib.Nifti1Image(np.ones((4, 4, 4, 2), dtype=np.float32), np.eye(4))
    with pytest.raises(ValueError, match="3D"):
        n4_bias_correct(four_d)


def test_n4_oblique_nifti_matches_native_ants_file_reader(tmp_path) -> None:
    ants = pytest.importorskip('ants')
    shape = (20, 21, 22)
    x, y, z = np.indices(shape, dtype=np.float32)
    data = (100 * np.exp(.005*x + .007*y - .003*z)).astype(np.float32)
    angle = np.deg2rad(11)
    affine = np.eye(4)
    affine[:3, :3] = np.array([[np.cos(angle), -np.sin(angle), 0],
                              [np.sin(angle), np.cos(angle), 0], [0, 0, 1]]) @ np.diag([.89007622, .91, 1.2])
    affine[:3, 3] = [-113.5, -101, -67.5]
    source = nib.Nifti1Image(data, affine)
    path = tmp_path / 'oblique.nii.gz'
    nib.save(source, path)
    source = nib.load(path)
    config = N4Config(shrink_factor=1, iterations=(3,))
    expected = ants.n4_bias_field_correction(
        ants.image_read(str(path)), mask=None, rescale_intensities=False,
        shrink_factor=1, convergence={'iters': [3], 'tol': 0.0}, spline_param=[1,1,1],
    )
    corrected = n4_bias_correct(source, config=config)
    np.testing.assert_array_equal(np.asarray(corrected.dataobj), expected.numpy())
    np.testing.assert_array_equal(corrected.affine, source.affine)


def test_n4_accepts_interpolation_induced_negative_values() -> None:
    pytest.importorskip("ants")
    size = 24
    axis = np.linspace(-1.0, 1.0, size, dtype=np.float32)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    values = (100.0 * np.exp(0.12 * x - 0.08 * y + 0.05 * z)).astype(np.float32)
    # Cubic interpolation near the image boundary can overshoot below zero;
    # ANTsPy's N4 accepts this mixed-valued image and the source wrapper clips
    # interpolation negatives after N4/masking.
    values[:2] = -2.0
    image = nib.Nifti1Image(values, np.eye(4))
    corrected = n4_bias_correct(
        image,
        config=N4Config(shrink_factor=1, iterations=(2,), spline_mesh=(1, 1, 1)),
    )
    result = np.asarray(corrected.dataobj)
    assert np.isfinite(result).all()
    assert np.any(result < 0)
