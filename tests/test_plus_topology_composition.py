from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

import pymacruise.plus_topology_composition as composition
from pymacruise.plus_pipeline import PlusPipelineOptions
from pymacruise.topology import topology_correct


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


def _inputs() -> tuple[nib.Nifti1Image, ...]:
    provenance = json.loads((FIXTURE / "provenance.json").read_text())
    images = []
    for name in INPUT_NAMES:
        path = FIXTURE / name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == provenance["sha256"][name]
        image = nib.load(str(path))
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        images.append(image)
    return tuple(images)


def test_topology_refinement_precedes_ace_gvf_and_tgdm(monkeypatch) -> None:
    events: list[str] = []
    tgdm_inputs: list[tuple[object, ...]] = []
    traced_names = (
        "run_plus_boundary_stage",
        "run_plus_tissue_update_stage",
        "run_plus_topology_stage",
        "run_plus_ace_stage",
        "run_plus_gvf_stage",
        "run_plus_tgdm_stage",
    )
    originals = {name: getattr(composition, name) for name in traced_names}
    for name, original in originals.items():

        def traced(*args, _name=name, _original=original, **kwargs):
            events.append(_name)
            if _name == "run_plus_tgdm_stage":
                tgdm_inputs.append(args)
            return _original(*args, **kwargs)

        monkeypatch.setattr(composition, name, traced)

    inputs = _inputs()
    before = [np.asarray(image.dataobj).copy() for image in inputs]
    result = composition.run_plus_topology_and_tgdm(
        *inputs,
        options=PlusPipelineOptions(),
    )

    assert events == list(traced_names)
    assert result.topology.atlas_refined is not None
    assert result.tgdm_initial_membership is result.topology.atlas_refined
    expected_topology = topology_correct(
        result.tissue.wm,
        result.tissue.mask,
        multi_atlas_wm=inputs[3],
    )
    np.testing.assert_array_equal(
        result.topology.corrected.dataobj,
        expected_topology.corrected.dataobj,
    )
    np.testing.assert_array_equal(
        result.tgdm_initial_membership.dataobj,
        expected_topology.atlas_refined.dataobj,
    )

    assert len(tgdm_inputs) == 1
    routed = tgdm_inputs[0]
    assert routed[0] is result.tgdm_initial_membership
    assert routed[1] is result.tissue.wm
    assert routed[2] is result.tissue.gm
    assert routed[3] is result.ace.enhanced_gray_matter
    assert routed[4] is result.ace.thinned_skeleton
    assert routed[5] is result.gvf.vector_field

    for image in (
        result.boundary.raw_boundary,
        result.boundary.cleaned_boundary,
        result.tissue.gm,
        result.tissue.wm,
        result.tissue.mask,
        result.tgdm_initial_membership,
        result.ace.enhanced_gray_matter,
        result.gvf.vector_field,
        result.tgdm.outer,
    ):
        assert image.shape[:3] == inputs[0].shape
        np.testing.assert_allclose(image.affine, inputs[0].affine, rtol=0, atol=1e-6)
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert np.isfinite(np.asarray(image.dataobj)).all()
    for image, original in zip(inputs, before, strict=True):
        np.testing.assert_array_equal(image.dataobj, original)
