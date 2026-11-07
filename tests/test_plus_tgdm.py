from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.plus_pipeline import (
    PlusPipelineOptions,
    plus_tgdm_parameters,
    run_plus_tgdm_stage,
)
from pymacruise.tgdm import tgdm_nested_images
from pymacruise.tgdm import TGDMBaselineParameters, tgdm_initial_membership_level_set


FIXTURE = Path(__file__).parent / "data/synthetic/tgdm"


def _inputs() -> tuple[nib.spatialimages.SpatialImage, ...]:
    images = {
        path.name.removesuffix(".nii.gz"): nib.load(path)
        for path in FIXTURE.glob("*.nii.gz")
    }
    reference = images["initial_membership"]
    zero_gvf = nib.Nifti1Image(
        np.zeros((*reference.shape, 3), dtype=np.float32), reference.affine
    )
    return (
        reference,
        images["white_matter"],
        images["gray_matter"],
        images["ace_gray_matter"],
        images["ace_skeleton"],
        zero_gvf,
    )


def test_plus_tgdm_stage_routes_cropped_inputs_and_wm_iso_level() -> None:
    inputs = _inputs()
    options = PlusPipelineOptions(wm_iso_level=0.43)
    actual = run_plus_tgdm_stage(*inputs, options=options)
    direct = tgdm_nested_images(
        inputs[0],
        inputs[1],
        inputs[2],
        inputs[3],
        ace_skeleton=inputs[4],
        gvf=inputs[5],
        parameters=plus_tgdm_parameters(options),
    )

    for name in (
        "initial_level_set",
        "tissue_sign",
        "initial_inner",
        "inner",
        "central",
        "outer",
    ):
        actual_image = getattr(actual, name)
        direct_image = getattr(direct, name)
        np.testing.assert_array_equal(actual_image.dataobj, direct_image.dataobj)
        np.testing.assert_array_equal(actual_image.affine, inputs[0].affine)
        assert nib.aff2axcodes(actual_image.affine) == ("R", "A", "S")
    assert actual.execution.outer.completed_steps == direct.execution.outer.completed_steps


def test_plus_tgdm_stage_rejects_a_gvf_on_another_physical_grid() -> None:
    inputs = list(_inputs())
    shifted = inputs[-1].affine.copy()
    shifted[0, 3] += 2.0
    inputs[-1] = nib.Nifti1Image(np.asarray(inputs[-1].dataobj), shifted)
    with pytest.raises(ValueError, match="GVF must match the initial-membership RAS grid"):
        run_plus_tgdm_stage(*inputs)


def test_plus_tgdm_stage_accepts_hash_replayed_precomputed_initial_level_set() -> None:
    inputs = _inputs()
    seed = tgdm_initial_membership_level_set(
        np.asarray(inputs[0].dataobj), iso_level=0.5, max_distance=8.0
    )
    seed_image = nib.Nifti1Image(seed, inputs[0].affine)
    zero_steps = TGDMBaselineParameters(
        initial_iterations=0, inner_iterations=0, central_shrink_iterations=0,
        central_iterations=0, outer_iterations=0,
    )

    actual = run_plus_tgdm_stage(
        *inputs,
        options=PlusPipelineOptions(),
        precomputed_initial_level_set=seed_image,
    )

    np.testing.assert_array_equal(actual.initial_level_set.dataobj, seed)
