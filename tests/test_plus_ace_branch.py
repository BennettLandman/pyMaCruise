from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.ace import ace_enhance
from pymacruise.boundary import (
    clean_csf_boundary_image,
    combine_multi_atlas_toads_images,
    multi_atlas_csf_boundary_image,
)
from pymacruise.plus_pipeline import (
    PlusPipelineOptions,
    run_plus_ace_stage,
    run_plus_boundary_stage,
    run_plus_tissue_update_stage,
)


FIXTURE = Path(__file__).parent / "data/synthetic/plus_ace_branch"
BASE_ACE = Path(__file__).parent / "data/synthetic/ace"
TISSUE_INPUTS = (
    "multi_atlas_gm",
    "multi_atlas_wm",
    "multi_atlas_brain_mask",
    "toads_gm",
    "toads_wm",
    "toads_mask",
)


def _load_inputs(lesion: bool):
    segmentation_name = "segmentation_lesion" if lesion else "segmentation_no_lesion"
    inputs = {
        "segmentation": nib.load(FIXTURE / f"{segmentation_name}.nii.gz"),
        "toads_wm_boundary": nib.load(FIXTURE / "toads_wm_boundary.nii.gz"),
    }
    inputs.update(
        {name: nib.load(FIXTURE / f"{name}.nii.gz") for name in TISSUE_INPUTS}
    )
    return inputs


def _snapshot(inputs) -> dict[str, np.ndarray]:
    return {
        name: np.asarray(image.dataobj).copy()
        for name, image in inputs.items()
    }


def test_plus_boundary_tissue_ace_branch_with_and_without_lesion_label() -> None:
    options = PlusPipelineOptions()
    branch_outputs = {}
    saved_inputs = []

    for lesion in (True, False):
        inputs = _load_inputs(lesion)
        before = _snapshot(inputs)
        saved_inputs.append((inputs, before))

        boundary = run_plus_boundary_stage(
            inputs["segmentation"], inputs["toads_wm_boundary"], options
        )
        expected_raw = multi_atlas_csf_boundary_image(inputs["segmentation"])
        expected_cleaned = clean_csf_boundary_image(
            expected_raw, inputs["toads_wm_boundary"], mode="toads_wm"
        )
        np.testing.assert_array_equal(boundary.raw_boundary.dataobj, expected_raw.dataobj)
        np.testing.assert_array_equal(
            boundary.cleaned_boundary.dataobj, expected_cleaned.dataobj
        )
        for image in (boundary.raw_boundary, boundary.cleaned_boundary):
            assert nib.aff2axcodes(image.affine) == ("R", "A", "S")

        tissue_args = [inputs[name] for name in TISSUE_INPUTS]
        tissue = run_plus_tissue_update_stage(*tissue_args)
        direct_tissue = combine_multi_atlas_toads_images(
            *tissue_args, gm_confident_level=0.001
        )
        for output, direct in zip(
            (tissue.gm, tissue.wm, tissue.mask), direct_tissue, strict=True
        ):
            np.testing.assert_array_equal(output.dataobj, direct.dataobj)
            np.testing.assert_array_equal(output.affine, direct.affine)

        ace = run_plus_ace_stage(
            tissue.gm, tissue.wm, boundary.cleaned_boundary, options
        )
        direct_ace = ace_enhance(
            tissue.gm,
            tissue.wm,
            multi_atlas_csf=boundary.cleaned_boundary,
            use_updated_csf_in_ace=3,
            merge_method="Merge",
        )
        for name in ("enhanced_gray_matter", "skeleton", "thinned_skeleton"):
            actual_image = getattr(ace, name)
            direct_image = getattr(direct_ace, name)
            np.testing.assert_array_equal(actual_image.dataobj, direct_image.dataobj)
            np.testing.assert_array_equal(actual_image.affine, direct_image.affine)
            assert nib.aff2axcodes(actual_image.affine) == ("R", "A", "S")
        branch_outputs[lesion] = (boundary, tissue, ace)

    # The selected boundary ignores lesion label 202. The paired branch keeps
    # the same source geometry and upstream tissue products but changes the
    # generated boundary where that label occupied the cortical region.
    lesion_boundary = np.asarray(branch_outputs[True][0].raw_boundary.dataobj)
    ordinary_boundary = np.asarray(branch_outputs[False][0].raw_boundary.dataobj)
    assert np.any(lesion_boundary != ordinary_boundary)
    assert np.any(lesion_boundary)
    assert np.any(ordinary_boundary)

    for inputs, before in saved_inputs:
        for name, image in inputs.items():
            np.testing.assert_array_equal(image.dataobj, before[name])


def test_default_base_ace_call_still_matches_established_java_reference() -> None:
    gray = nib.load(BASE_ACE / "gray_matter.nii.gz")
    white = nib.load(BASE_ACE / "white_matter.nii.gz")
    atlas_csf = nib.load(BASE_ACE / "multi_atlas_csf.nii.gz")
    # No plus options or plus merge placement: preserve the default base path.
    actual = ace_enhance(gray, white, multi_atlas_csf=atlas_csf)
    reference_dir = BASE_ACE / "java_reference"
    for name, filename in (
        ("enhanced_gray_matter", "enhanced_gm_java.nii.gz"),
        ("skeleton", "csf_skeleton_java.nii.gz"),
        ("thinned_skeleton", "thinned_skeleton_java.nii.gz"),
    ):
        np.testing.assert_array_equal(
            getattr(actual, name).dataobj,
            nib.load(reference_dir / filename).dataobj,
        )
