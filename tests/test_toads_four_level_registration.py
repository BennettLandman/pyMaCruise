from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_alignment_startup import (
    prepare_toads_alignment_startup,
    register_toads_alignment_startup,
)


FIXTURE = Path(__file__).parent / "data" / "synthetic" / "toads_four_level_registration"


def test_registered_startup_runs_all_three_source_scales_at_four_levels() -> None:
    t1 = nib.load(FIXTURE / "t1_ras.nii.gz")
    brain_mask = nib.load(FIXTURE / "brain_mask_ras.nii.gz")
    topology = np.asarray(nib.load(FIXTURE / "atlas_topology_ras.nii.gz").dataobj)
    atlas_shapes = np.asarray(nib.load(FIXTURE / "atlas_shape_priors_ras.nii.gz").dataobj)
    priors = np.asarray(((0.10, 0.36, 0.63, 0.88),), dtype=np.float32)

    startup = prepare_toads_alignment_startup(
        t1,
        brain_mask,
        topology,
        priors,
        (1.0, 1.0, 1.0),
        padding=2,
    )
    # Every source scale used here must remain sampleable by the trilinear
    # interpolation kernel on both independent input and atlas grids.
    assert min(startup.initial_memberships.shape[:3]) // 8 >= 2
    assert min(atlas_shapes.shape[:3]) // 8 >= 2
    assert nib.aff2axcodes(startup.intensity.image.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(brain_mask.affine) == ("R", "A", "S")

    result = register_toads_alignment_startup(
        startup,
        atlas_shapes,
        levels=4,
        max_iterations=1,
    )

    assert [level.scale for level in result.registration.levels] == [8, 4, 2]
    coarsest, middle, finest = result.registration.levels
    # The selected-source line search keeps the initial transform when its
    # scale-8 energy and relative change are nonfinite; later scales still run.
    assert coarsest.result.transform.shape == (6,)
    np.testing.assert_array_equal(coarsest.result.transform, coarsest.initial_transform)
    assert np.isnan(coarsest.result.energy)
    assert coarsest.result.stop_reason == "nonfinite_relative_change"
    assert len(coarsest.result.iterations) == 1
    rejected = coarsest.result.iterations[0]
    assert not rejected.accepted
    assert np.isnan(rejected.energy)
    assert np.isnan(rejected.relative_change)
    assert rejected.rejected_trials == 11
    for level in (middle, finest):
        assert level.result.transform.shape == (6,)
        assert np.isfinite(level.initial_transform).all()
        assert np.isfinite(level.result.transform).all()
        assert np.isfinite(level.result.energy)
    np.testing.assert_array_equal(
        middle.initial_transform,
        coarsest.result.transform,
    )
    np.testing.assert_array_equal(
        finest.initial_transform,
        middle.result.transform,
    )
