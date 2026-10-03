from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

import pymacruise.plus_topology_composition as composition
from pymacruise.cropping import CropParameters
from pymacruise.plus_pipeline import PlusPipelineOptions, run_plus_post_tgdm_stages


FIXTURE = Path(__file__).parent / "data/synthetic/plus_pre_tgdm_composition"
INPUT_NAMES = (
    "segmentation.nii.gz",
    "toads_wm_for_boundary.nii.gz",
    "multi_atlas_gm.nii.gz",
    "multi_atlas_wm.nii.gz",
    "multi_atlas_brain_mask.nii.gz",
    "toads_gm.nii.gz",
    "toads_wm.nii.gz",
    "toads_mask.nii.gz",
)


def _crop_ready_products() -> tuple[nib.Nifti1Image, ...]:
    """Load explicit precomputed upstream TOADS/atlas products.

    These are synthetic products supplied at the crop-ready boundary. This
    test does not run or claim TOADS, LesionTOADS, registration, or cropping.
    """
    provenance = json.loads((FIXTURE / "provenance.json").read_text())
    images = []
    for name in INPUT_NAMES:
        path = FIXTURE / name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == provenance["sha256"][name]
        image = nib.load(str(path))
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        images.append(image)
    return tuple(images)


def test_plus_crop_ready_products_run_through_tgdm_surfaces_and_thickness() -> None:
    inputs = _crop_ready_products()
    before = [np.asarray(image.dataobj).copy() for image in inputs]
    result = composition.run_plus_topology_and_tgdm(
        *inputs,
        options=PlusPipelineOptions(),
    )

    reference = inputs[0]
    # A full-grid crop descriptor keeps this integration focused on available
    # post-TGDM APIs while still validating their grid contract and JIST uncrop.
    crop = CropParameters(
        original_shape=tuple(int(n) for n in reference.shape),
        bounds=tuple((0, int(n) - 1) for n in reference.shape),
        background=0.0,
        original_affine=reference.affine.copy(),
    )
    post = run_plus_post_tgdm_stages(result.tgdm, crop)

    assert result.topology.atlas_refined is not None
    assert result.tgdm_initial_membership is result.topology.atlas_refined
    for image in (
        result.boundary.raw_boundary,
        result.boundary.cleaned_boundary,
        result.tissue.gm,
        result.tissue.wm,
        result.tissue.mask,
        result.tgdm_initial_membership,
        result.ace.enhanced_gray_matter,
        result.gvf.vector_field,
        result.tgdm.inner,
        result.tgdm.central,
        result.tgdm.outer,
    ):
        assert image.shape[:3] == reference.shape
        np.testing.assert_allclose(image.affine, reference.affine, rtol=0, atol=1e-6)
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert np.isfinite(np.asarray(image.dataobj)).all()

    for level_set in (post.level_sets.inner, post.level_sets.central, post.level_sets.outer):
        assert level_set.shape == reference.shape
        assert nib.aff2axcodes(level_set.affine) == ("R", "A", "S")
        assert np.isfinite(np.asarray(level_set.dataobj)).all()
    for mesh in (post.surfaces.inner, post.surfaces.central, post.surfaces.outer):
        assert mesh.coordinate_frame == "RAS"
        assert mesh.units == "mm"
        assert mesh.vertices.shape[1] == 3
        assert mesh.faces.shape[1] == 3
        assert np.isfinite(mesh.vertices).all()
    assert post.thickness.volume.image.shape == reference.shape
    assert nib.aff2axcodes(post.thickness.volume.image.affine) == ("R", "A", "S")
    assert np.isfinite(np.asarray(post.thickness.volume.image.dataobj)).all()
    assert np.isfinite(post.thickness.volume.mean_mm)
    assert np.isfinite(post.thickness.surface.mean_mm)

    for image, original in zip(inputs, before, strict=True):
        np.testing.assert_array_equal(image.dataobj, original)
