from __future__ import annotations

import nibabel as nib
import numpy as np
import pytest

from pymacruise.cortical_separation import cortical_separation


def test_cortical_separation_matches_source_label_sets_and_values(tmp_path) -> None:
    labels = [
        0, 4, 11, 23, 30, 31, 32, 35, 36, 37, 38, 39, 40, 41, 44, 45, 47,
        48, 49, 50, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62, 75, 76, 99,
        100, 109, 110, 112, 125, 126, 127, 128, 129, 130, 131, 132, 157,
        158, 159, 160, 187, 188, 189, 190, 207, 208,
    ]
    data = np.asarray(labels, dtype=np.int16).reshape((len(labels), 1, 1))
    image = nib.Nifti1Image(data, np.eye(4))
    outputs = cortical_separation(image, tmp_path)

    gm_selected = {
        100, 109, 112, 125, 128, 129, 132, 157, 160, 187, 190, 207,
    }
    wm_selected = {4, 23, 30, 36, 37, 44, 45, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62}
    brain_selected = gm_selected | wm_selected | {31, 32, 47, 48, 49, 50, 75, 76}
    expected = {
        "Atlas_GM": gm_selected,
        "Atlas_WM": wm_selected,
        "Atlas_Mask": wm_selected,
        "Atlas_BrainTissueMask": brain_selected,
    }
    for name, selected in expected.items():
        result = outputs[name]
        values = np.asanyarray(result.dataobj).reshape(-1)
        np.testing.assert_array_equal(values, [255 if label in selected else 0 for label in labels])
        assert result.get_data_dtype() == np.dtype(np.uint8)
        assert nib.aff2axcodes(result.affine) == ("R", "A", "S")
        persisted = nib.load(tmp_path / f"{name}.nii.gz")
        np.testing.assert_array_equal(persisted.dataobj, result.dataobj)


def test_cortical_separation_warns_and_reorients_lps_input() -> None:
    data = np.zeros((5, 6, 7), dtype=np.int16)
    data[1, 2, 3] = 100
    flip = np.diag([-1.0, -1.0, 1.0, 1.0])
    flip[0, 3] = data.shape[0] - 1
    flip[1, 3] = data.shape[1] - 1
    image = nib.Nifti1Image(data[::-1, ::-1, :], flip)
    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        outputs = cortical_separation(image)
    assert nib.aff2axcodes(outputs["Atlas_GM"].affine) == ("R", "A", "S")
    assert np.argwhere(np.asanyarray(outputs["Atlas_GM"].dataobj) == 255).tolist() == [[1, 2, 3]]
    assert not np.any(np.asanyarray(outputs["Atlas_WM"].dataobj))
