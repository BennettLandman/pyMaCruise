from __future__ import annotations

import numpy as np
import nibabel as nib

import pymacruise.plus_pipeline as plus_pipeline
from pymacruise.plus_pipeline import PlusPipelineOptions, run_plus_pre_tgdm_stages


def _folded_ras_phantom() -> tuple[nib.Nifti1Image, ...]:
    """Build asymmetric folded-cortex-like stage inputs on one RAS grid."""
    size = 25
    grid = np.indices((size, size, size), dtype=np.float32)
    x, y, z = (grid[axis] - 12.0 for axis in range(3))
    radius = np.sqrt(x * x + y * y + z * z)
    theta = np.arctan2(y, x)
    phi = np.arctan2(z, np.sqrt(x * x + y * y))
    folds = 0.85 * np.sin(3.0 * theta + 0.4) * np.cos(2.0 * phi)
    folds += 0.35 * np.sin(5.0 * theta - 0.7) * np.cos(phi)
    folded_radius = radius + folds
    wm = folded_radius < 5.2
    gm = (folded_radius >= 5.2) & (folded_radius < 8.6)
    brain = folded_radius < 9.7

    # Vary labels across the cortical band so the selected 26-neighbor
    # boundary operation sees folded, asymmetric class interfaces.
    sectors = (np.floor((theta + np.pi) * 2.5).astype(np.int32) % 5)
    segmentation = np.zeros((size, size, size), dtype=np.int16)
    segmentation[brain] = (200 + sectors[brain]).astype(np.int16)
    segmentation[gm] = (210 + sectors[gm]).astype(np.int16)
    # Include one selected overlook label in an offset patch.
    segmentation[gm & (x > 1.0) & (y < -1.0)] = 202

    atlas_gm = np.where(brain, 180.0, 0.0).astype(np.float32)
    atlas_wm = np.where(wm, 255.0, 0.0).astype(np.float32)
    toads_gm = np.where(gm, 200.0, 0.0).astype(np.float32)
    toads_wm = np.where(wm, 255.0, 0.0).astype(np.float32)
    mask = brain.astype(np.uint8) * 255
    boundary_wm = toads_wm.copy()
    scalar = (80.0 + 5.0 * x + 3.0 * y + 2.0 * z + 2.0 * folds).astype(np.float32)

    affine = np.diag((0.8, 0.8, 1.2, 1.0))
    arrays = (
        segmentation,
        boundary_wm,
        atlas_gm,
        atlas_wm,
        mask,
        toads_gm,
        toads_wm,
        mask,
        scalar,
    )
    return tuple(nib.Nifti1Image(data, affine) for data in arrays)


def test_plus_pre_tgdm_composition_routes_boundary_tissue_ace_gvf_in_order(
    monkeypatch,
) -> None:
    events: list[str] = []
    originals = {
        name: getattr(plus_pipeline, name)
        for name in (
            "run_plus_boundary_stage",
            "run_plus_tissue_update_stage",
            "run_plus_ace_stage",
            "run_plus_gvf_stage",
        )
    }
    for name in originals:
        original = originals[name]

        def traced(*args, _name=name, _original=original, **kwargs):
            events.append(_name)
            return _original(*args, **kwargs)

        monkeypatch.setattr(plus_pipeline, name, traced)

    inputs = _folded_ras_phantom()
    before = [np.asarray(image.dataobj).copy() for image in inputs]
    result = run_plus_pre_tgdm_stages(*inputs, options=PlusPipelineOptions())

    assert events == [
        "run_plus_boundary_stage",
        "run_plus_tissue_update_stage",
        "run_plus_ace_stage",
        "run_plus_gvf_stage",
    ]
    for image in (
        result.boundary.raw_boundary,
        result.boundary.cleaned_boundary,
        result.tissue.gm,
        result.tissue.wm,
        result.tissue.mask,
        result.ace.enhanced_gray_matter,
        result.ace.skeleton,
        result.ace.thinned_skeleton,
        result.gvf.vector_field,
    ):
        assert image.shape[:3] == inputs[0].shape
        np.testing.assert_allclose(image.affine, inputs[0].affine, rtol=0, atol=1e-6)
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")

    assert result.gvf.vector_field.shape == (*inputs[0].shape, 3)
    assert np.any(np.asarray(result.boundary.raw_boundary.dataobj))
    assert np.any(np.asarray(result.boundary.cleaned_boundary.dataobj))
    assert np.any(np.asarray(result.ace.skeleton.dataobj))
    assert np.isfinite(np.asarray(result.gvf.vector_field.dataobj)).all()
    for image, original in zip(inputs, before, strict=True):
        np.testing.assert_array_equal(image.dataobj, original)
