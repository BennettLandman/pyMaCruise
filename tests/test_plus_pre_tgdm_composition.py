from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import nibabel as nib

import pymacruise.plus_pipeline as plus_pipeline
from pymacruise.plus_pipeline import PlusPipelineOptions, run_plus_pre_tgdm_stages


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


def _folded_ras_phantom() -> tuple[nib.Nifti1Image, ...]:
    """Load the committed asymmetric folded-cortex RAS stage inputs."""
    provenance = json.loads((FIXTURE / "provenance.json").read_text())
    images = []
    for name in INPUT_NAMES:
        path = FIXTURE / name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == provenance["sha256"][name]
        image = nib.load(str(path))
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        images.append(image)
    return tuple(images)


def test_plus_pre_tgdm_composition_routes_boundary_tissue_ace_gvf_in_order(
    monkeypatch,
) -> None:
    events: list[str] = []
    gvf_inputs: list[nib.spatialimages.SpatialImage] = []
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
            if _name == "run_plus_gvf_stage":
                gvf_inputs.append(args[0])
            return _original(*args, **kwargs)

        monkeypatch.setattr(plus_pipeline, name, traced)

    inputs = _folded_ras_phantom()
    before = [np.asarray(image.dataobj).copy() for image in inputs]
    result = run_plus_pre_tgdm_stages(*inputs, options=PlusPipelineOptions())
    assert gvf_inputs == [result.ace.enhanced_gray_matter]
    assert gvf_inputs[0] is result.ace.enhanced_gray_matter

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
