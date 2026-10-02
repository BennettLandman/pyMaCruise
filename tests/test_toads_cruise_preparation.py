from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.cli import main
from pymacruise.toads_cruise_preparation import prepare_toads_cruise_inputs


FIXTURE = Path(__file__).parent / "data" / "synthetic" / "toads_tissue_crop"


def _inputs() -> tuple[nib.spatialimages.SpatialImage, ...]:
    return tuple(
        nib.load(str(FIXTURE / filename))
        for filename in ("gm.nii.gz", "wmfill.nii.gz", "wmmask.nii.gz")
    )


def test_selected_postpass_composes_cleanup_crop_and_topology_in_source_order() -> None:
    gm, wmfill, wmmask = _inputs()
    segmentation = nib.load(str(FIXTURE / "maseg.nii.gz"))
    result = prepare_toads_cruise_inputs(
        gm,
        wmfill,
        wmmask,
        segmentation=segmentation,
        remove_cerebellum=True,
        crop_edge_width=2,
    )
    assert result.stage_order == (
        "removeBrainstemAndCerebellum",
        "cropByPositiveGrayMatter",
        "topologyCorrectWhiteMatterFill",
    )
    assert result.removed_voxel_counts == (8, 8, 8)
    assert result.crop_parameters.bounds == ((1, 11), (0, 10), (1, 9))
    assert result.gray_matter.shape == (11, 11, 9)
    assert result.white_matter_fill.shape == result.gray_matter.shape
    assert result.white_matter_mask.shape == result.gray_matter.shape
    assert result.topology_corrected_white_matter_fill.shape == result.gray_matter.shape
    for image in (
        result.gray_matter,
        result.white_matter_fill,
        result.white_matter_mask,
        result.topology_corrected_white_matter_fill,
    ):
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert result.topology_corrected_white_matter_fill.get_data_dtype() == np.dtype(np.float32)


def test_prepare_cruise_inputs_cli_writes_crop_and_corrected_wm_products(tmp_path) -> None:
    gm, wmfill, wmmask = _inputs()
    segmentation = nib.load(str(FIXTURE / "maseg.nii.gz"))
    expected = prepare_toads_cruise_inputs(
        gm,
        wmfill,
        wmmask,
        segmentation=segmentation,
        remove_cerebellum=True,
        crop_edge_width=2,
    )
    output_dir = tmp_path / "prepared"
    provenance = tmp_path / "run.json"
    assert main([
        "prepare-cruise-inputs",
        "--gm", str(FIXTURE / "gm.nii.gz"),
        "--wm-fill", str(FIXTURE / "wmfill.nii.gz"),
        "--wm-mask", str(FIXTURE / "wmmask.nii.gz"),
        "--segmentation", str(FIXTURE / "maseg.nii.gz"),
        "--remove-cerebellum",
        "--edge-width", "2",
        "--output-dir", str(output_dir),
        "--provenance", str(provenance),
    ]) == 0
    assert provenance.is_file()
    outputs = {
        "gm_cropped.nii.gz": expected.gray_matter,
        "wmfill_cropped.nii.gz": expected.white_matter_fill,
        "wmmask_cropped.nii.gz": expected.white_matter_mask,
        "wmfill_topology_corrected.nii.gz": expected.topology_corrected_white_matter_fill,
    }
    for filename, reference in outputs.items():
        image = nib.load(str(output_dir / filename))
        assert image.shape == reference.shape
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(np.asarray(image.dataobj), np.asarray(reference.dataobj))
