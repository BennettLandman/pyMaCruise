from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.cli import main
from pymacruise.toads_tissue_crop import crop_toads_tissue_maps


FIXTURE = Path(__file__).parent / "data" / "synthetic" / "toads_tissue_crop"
EXPECTED_BOUNDS = ((1, 11), (0, 10), (1, 9))
EXPECTED_SHAPE = (11, 11, 9)


def _inputs() -> tuple[nib.spatialimages.SpatialImage, ...]:
    return tuple(
        nib.load(str(FIXTURE / filename))
        for filename in ("gm.nii.gz", "wmfill.nii.gz", "wmmask.nii.gz")
    )


def test_crop_toads_maps_uses_one_source_order_gm_box() -> None:
    gm, wmfill, wmmask = _inputs()
    result = crop_toads_tissue_maps(gm, wmfill, wmmask, edge_width=2)
    assert result.crop_parameters.bounds == EXPECTED_BOUNDS
    assert result.gray_matter.shape == EXPECTED_SHAPE
    outputs = (result.gray_matter, result.white_matter_fill, result.white_matter_mask)
    for output in outputs:
        assert output.shape == EXPECTED_SHAPE
        assert nib.aff2axcodes(output.affine) == ("R", "A", "S")
        assert np.allclose(
            output.affine,
            gm.affine @ np.array(
                [[1, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 1], [0, 0, 0, 1]],
                dtype=float,
            ),
        )
    expected_wmfill = np.asarray(wmfill.dataobj)[1:12, 0:11, 1:10]
    np.testing.assert_array_equal(result.white_matter_fill.dataobj, expected_wmfill)
    np.testing.assert_array_equal(
        result.white_matter_mask.dataobj,
        np.asarray(wmmask.dataobj)[1:12, 0:11, 1:10],
    )


def test_crop_tissue_maps_cli_writes_all_three_ras_products(tmp_path) -> None:
    output_dir = tmp_path / "crop"
    provenance = tmp_path / "crop.json"
    assert main([
        "crop-tissue-maps",
        "--gm", str(FIXTURE / "gm.nii.gz"),
        "--wm-fill", str(FIXTURE / "wmfill.nii.gz"),
        "--wm-mask", str(FIXTURE / "wmmask.nii.gz"),
        "--output-dir", str(output_dir),
        "--edge-width", "2",
        "--provenance", str(provenance),
    ]) == 0
    assert provenance.is_file()
    for filename in (
        "gm_cropped.nii.gz", "wmfill_cropped.nii.gz", "wmmask_cropped.nii.gz"
    ):
        image = nib.load(str(output_dir / filename))
        assert image.shape == EXPECTED_SHAPE
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
