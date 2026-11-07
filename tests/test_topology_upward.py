"""Independent selected-Java fixtures for v3.2 upward topology correction."""
from pathlib import Path
import nibabel as nib
import numpy as np
import pytest
from pymacruise.topology import topology_correct

FIXTURES = Path(__file__).parent/'data/synthetic/topology_upward'

@pytest.mark.parametrize('name', ['spatial_hash','plateau','boundary','narrow_branch','disconnected_seeds'])
def test_upward_matches_saved_java(name):
    with np.load(FIXTURES/(name+'.npz')) as f:
        image = nib.Nifti1Image(f['input'],np.eye(4))
        empty = nib.Nifti1Image(np.zeros(image.shape,dtype=np.uint8),image.affine)
        result = topology_correct(image,empty,propagation_direction='background->object')
        np.testing.assert_array_equal(result.corrected.dataobj,f['reference'])
        np.testing.assert_array_equal(result.corrected.affine,image.affine)


def test_lps_order_and_unused_paint_are_restored_to_public_ras():
    with np.load(FIXTURES/'plateau.npz') as f:
        values = f['input'][::-1,::-1,:].copy()
        affine = np.diag([.8,1.2,1.7,1.]);affine[:3,3] = [-55,21,-90]
        image = nib.Nifti1Image(values,affine)
        arbitrary = nib.Nifti1Image(np.ones(values.shape,dtype=np.uint8),affine)
        result = topology_correct(image,arbitrary,propagation_direction='background->object',source_axis_order='LPS')
        np.testing.assert_array_equal(result.corrected.dataobj,f['reference'][::-1,::-1,:])
        np.testing.assert_array_equal(result.corrected.affine,affine)
