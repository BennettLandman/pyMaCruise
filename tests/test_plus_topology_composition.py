from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

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
    tgdm_inputs: list[tuple[tuple[object, ...], dict[str, object]]] = []
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
                tgdm_inputs.append((args, kwargs))
            return _original(*args, **kwargs)

        monkeypatch.setattr(composition, name, traced)

    inputs = _inputs()
    before = [np.asarray(image.dataobj).copy() for image in inputs]
    result = composition.run_plus_topology_and_tgdm(
        *inputs,
        options=PlusPipelineOptions(),
        brain_mask=inputs[4],
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
    routed, routed_kwargs = tgdm_inputs[0]
    assert routed[0] is result.tgdm_initial_membership
    assert routed[1] is result.tissue.wm
    assert routed[2] is result.tissue.gm
    assert routed[3] is result.ace.enhanced_gray_matter
    assert routed[4] is result.ace.thinned_skeleton
    assert routed[5] is result.gvf.vector_field
    assert routed_kwargs["brain_mask"] is inputs[4]

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


def test_las_segmentation_warns_and_composition_outputs_remain_ras(monkeypatch) -> None:
    ras_inputs = _inputs()
    ras_segmentation = ras_inputs[0]
    las_data = np.flip(np.asarray(ras_segmentation.dataobj), axis=0).copy()
    las_affine = ras_segmentation.affine.copy()
    las_affine[0, 0] *= -1
    las_affine[0, 3] = ras_segmentation.affine[0, 0] * (
        ras_segmentation.shape[0] - 1
    )
    las_segmentation = nib.Nifti1Image(las_data, las_affine)
    assert nib.aff2axcodes(las_segmentation.affine) == ("L", "A", "S")

    received_segmentation: list[nib.spatialimages.SpatialImage] = []
    original_boundary = composition.run_plus_boundary_stage

    def traced_boundary(segmentation, *args, **kwargs):
        received_segmentation.append(segmentation)
        return original_boundary(segmentation, *args, **kwargs)

    monkeypatch.setattr(composition, "run_plus_boundary_stage", traced_boundary)
    with pytest.warns(UserWarning, match="reorienting the image to RAS") as caught:
        result = composition.run_plus_topology_and_tgdm(
            las_segmentation,
            *ras_inputs[1:],
            options=PlusPipelineOptions(),
        )

    assert len(caught) == 1
    assert len(received_segmentation) == 1
    assert nib.aff2axcodes(received_segmentation[0].affine) == ("R", "A", "S")
    np.testing.assert_array_equal(
        received_segmentation[0].dataobj,
        ras_segmentation.dataobj,
    )

    output_images = (
        result.boundary.raw_boundary,
        result.boundary.cleaned_boundary,
        result.tissue.gm,
        result.tissue.wm,
        result.tissue.mask,
        result.topology.corrected,
        result.topology.atlas_refined,
        result.tgdm_initial_membership,
        result.ace.enhanced_gray_matter,
        result.ace.skeleton,
        result.ace.thinned_skeleton,
        result.gvf.vector_field,
        result.tgdm.initial_level_set,
        result.tgdm.inner,
        result.tgdm.central,
        result.tgdm.outer,
    )
    assert all(image is not None for image in output_images)
    for image in output_images:
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
