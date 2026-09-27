from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.gvf import (
    _source_gradient,
    _source_laplacian,
    generalized_gvf,
    generalized_gvf_image,
)

GVF_PHANTOM = Path(__file__).parent / "data" / "synthetic" / "gvf"


def test_source_gradient_matches_analytic_axis_ramp() -> None:
    ramp = np.broadcast_to(np.arange(8, dtype=np.float32)[:, None, None], (8, 7, 6))
    gradient = _source_gradient(ramp / 7.0)

    np.testing.assert_allclose(gradient[..., 0], 1.0 / 7.0, atol=1e-7)
    np.testing.assert_array_equal(gradient[..., 1:], 0.0)


def test_generalized_gvf_preserves_constant_ramp_gradient() -> None:
    ramp = np.broadcast_to(np.arange(8, dtype=np.float32)[:, None, None], (8, 7, 6))

    vectors = generalized_gvf(ramp, lambda_value=1.0)

    assert vectors.shape == (*ramp.shape, 3)
    assert vectors.dtype == np.float32
    assert np.isfinite(vectors).all()
    np.testing.assert_allclose(vectors[..., 0], 1.0 / 7.0, atol=1e-6)
    np.testing.assert_allclose(vectors[..., 1:], 0.0, atol=1e-7)


def test_generalized_gvf_constant_and_zero_iteration_cases() -> None:
    constant = np.full((5, 6, 7), 12.0, dtype=np.float32)
    np.testing.assert_array_equal(generalized_gvf(constant), 0.0)

    ramp = np.broadcast_to(np.arange(5, dtype=np.float32)[:, None, None], (5, 4, 3))
    no_solve = generalized_gvf(ramp, max_iterations=0)
    np.testing.assert_allclose(no_solve[..., 0], 0.25, atol=1e-7)


def test_generalized_gvf_solves_each_source_poisson_system_and_respects_lambda() -> None:
    axes = np.indices((9, 9, 9), dtype=np.float32)
    scalar = np.exp(-sum((axis - 4) ** 2 for axis in axes) / 8.0).astype(np.float32)
    normalized = (scalar - scalar.min()) / (scalar.max() - scalar.min())
    gradient = _source_gradient(normalized)
    edge_weight = np.exp(-np.sum(gradient * gradient, axis=-1) / (0.5**2))
    g = edge_weight / 6.0
    h = 1.0 - edge_weight

    vectors = generalized_gvf(scalar, lambda_value=0.5)
    relative_residuals = []
    for component in range(3):
        residual = (
            h * vectors[..., component]
            - g * _source_laplacian(vectors[..., component])
            - h * gradient[..., component]
        )
        relative_residuals.append(
            np.linalg.norm(residual) / np.linalg.norm(h * gradient[..., component])
        )
    assert max(relative_residuals) < 2e-6

    low_lambda = generalized_gvf(scalar, lambda_value=0.1)
    high_lambda = generalized_gvf(scalar, lambda_value=2.0)
    assert np.linalg.norm(low_lambda - high_lambda) > 0.1


def test_generalized_gvf_rejects_invalid_inputs() -> None:
    with pytest.raises(ValueError, match="3D"):
        generalized_gvf(np.ones((5, 6)))
    with pytest.raises(ValueError, match="positive"):
        generalized_gvf(np.ones((4, 4, 4)), lambda_value=0)
    with pytest.raises(ValueError, match="finite"):
        generalized_gvf(np.full((4, 4, 4), np.nan))


def test_gvf_nifti_output_is_component_last_and_ras_with_warning() -> None:
    image = nib.load(GVF_PHANTOM / "gradient_ramp.nii.gz")
    affine = image.affine.copy()
    affine[0, 0] *= -1
    affine[0, 3] = image.shape[0] - 1
    non_ras = nib.Nifti1Image(np.asarray(image.dataobj), affine)

    with pytest.warns(UserWarning, match="reorienting"):
        result = generalized_gvf_image(non_ras)

    assert result.vector_field.shape == (*image.shape, 3)
    assert nib.aff2axcodes(result.vector_field.affine) == ("R", "A", "S")
    assert result.vector_field.header.get_intent()[0] == "vector"
    assert np.isfinite(result.vector_field.get_fdata()).all()
