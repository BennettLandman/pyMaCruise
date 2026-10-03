from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_atlas import load_toads_atlas, load_toads_atlas_ras_data


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/data/synthetic/toads_downsampled_vandy_atlas"


def test_downsampled_bundled_atlas_fixture_is_folded_ras_and_complete() -> None:
    atlas = load_toads_atlas(FIXTURE / "atlas.txt")
    atlas_data = load_toads_atlas_ras_data(atlas)
    t1 = nib.load(str(FIXTURE / "t1.nii.gz"))
    brain_mask = nib.load(str(FIXTURE / "brain-mask.nii.gz"))
    maseg = nib.load(str(FIXTURE / "maseg.nii.gz"))
    topology_image = nib.load(str(atlas.topology.path))
    topology = np.asanyarray(topology_image.dataobj)
    mask = np.asanyarray(brain_mask.dataobj) > 0

    assert atlas.shape_dimensions == t1.shape == brain_mask.shape == maseg.shape == topology_image.shape
    assert atlas.shape_samples == 8
    assert nib.aff2axcodes(t1.affine) == nib.aff2axcodes(brain_mask.affine) == ("R", "A", "S")
    assert np.array_equal(np.asanyarray(t1.dataobj) == 0, ~mask)
    assert np.array_equal(t1.affine, brain_mask.affine)
    assert np.array_equal(t1.affine, maseg.affine)
    assert set(np.unique(topology)) == {item.label for item in atlas.structures}
    assert set(atlas_data.shapes) == {item.name for item in atlas.structures}
    assert all(np.isfinite(values).all() and np.any(values > 0) for values in atlas_data.shapes.values())
    # Multiple cortical and sulcal transitions across adjacent voxels establish
    # that this downsized atlas retains folded label boundaries.
    assert np.count_nonzero(topology[:-1, :, :] != topology[1:, :, :]) > 1000
    assert np.count_nonzero(topology[:, :-1, :] != topology[:, 1:, :]) > 1000
