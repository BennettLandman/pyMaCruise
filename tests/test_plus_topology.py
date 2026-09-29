from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.plus_pipeline import PlusPipelineOptions, run_plus_topology_stage
from pymacruise.topology import topology_correct


FIXTURE = Path(__file__).parent / "data/synthetic/plus_topology"


def _inputs():
    return tuple(
        nib.load(FIXTURE / filename)
        for filename in ("topology_volume.nii.gz", "paint_mask.nii.gz", "ma_wm.nii.gz")
    )


def test_plus_topology_routes_mawm_into_existing_refinement() -> None:
    volume, paint, ma_wm = _inputs()
    actual = run_plus_topology_stage(volume, paint, ma_wm, PlusPipelineOptions())
    expected = topology_correct(volume, paint, multi_atlas_wm=ma_wm)

    np.testing.assert_array_equal(actual.corrected.dataobj, expected.corrected.dataobj)
    assert actual.atlas_refined is not None and expected.atlas_refined is not None
    np.testing.assert_array_equal(actual.atlas_refined.dataobj, expected.atlas_refined.dataobj)
    assert nib.aff2axcodes(actual.atlas_refined.affine) == ("R", "A", "S")
    assert np.any(
        np.asarray(actual.atlas_refined.dataobj) != np.asarray(actual.corrected.dataobj)
    )


def test_plus_topology_disabled_refinement_preserves_base_result() -> None:
    volume, paint, ma_wm = _inputs()
    options = PlusPipelineOptions(refine_topology_by_maseg=False)
    actual = run_plus_topology_stage(volume, paint, ma_wm, options)
    expected = topology_correct(volume, paint)

    np.testing.assert_array_equal(actual.corrected.dataobj, expected.corrected.dataobj)
    assert actual.atlas_refined is None


def test_plus_topology_requires_mawm_when_refinement_is_enabled() -> None:
    volume, paint, _ = _inputs()
    with pytest.raises(ValueError, match="requires the selected MaWM"):
        run_plus_topology_stage(volume, paint, None, PlusPipelineOptions())
