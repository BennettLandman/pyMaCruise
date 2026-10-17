from __future__ import annotations

import json

import nibabel as nib
import numpy as np

from pymacruise.cli import main


def test_atlasem_prelesion_dura_cli_writes_ras_outputs_and_provenance(tmp_path) -> None:
    shape = (17, 17, 17)
    x, y, z = np.indices(shape, dtype=np.float32)
    radius = np.sqrt((x - 8) ** 2 + (y - 8) ** 2 + (z - 8) ** 2)
    mask = (radius < 7.0).astype(np.uint8)
    membership = np.zeros((*shape, 3), dtype=np.float32)
    membership[..., 0] = np.where(mask, 0.05, 0.0)
    membership[..., 1] = np.where((radius >= 3.0) & (radius < 7.0), 0.75, 0.01)
    membership[..., 2] = np.where(radius < 3.0, 0.85, 0.02)
    membership[mask == 0] = 0.0
    t1_values = np.where(mask, 100.0 - 4.0 * radius, 0.0).astype(np.float32)
    # LAS input exercises the command's RAS normalization contract.
    affine = np.diag([-1.0, 1.0, 1.0, 1.0])
    t1_path = tmp_path / "t1.nii.gz"
    membership_path = tmp_path / "memberships.nii.gz"
    mask_path = tmp_path / "brain-mask.nii.gz"
    nib.save(nib.Nifti1Image(t1_values, affine), t1_path)
    nib.save(nib.Nifti1Image(membership, affine), membership_path)
    nib.save(nib.Nifti1Image(mask, affine), mask_path)
    output_dir = tmp_path / "outputs"
    provenance_path = tmp_path / "atlasem-dura.json"

    assert main(
        [
            "atlasem-prelesion-dura",
            "--t1", str(t1_path),
            "--memberships", str(membership_path),
            "--brain-mask", str(mask_path),
            "--class-names", '["Background","CorticalGM","WhiteMatter"]',
            "--output-dir", str(output_dir),
            "--max-iterations", "1",
            "--provenance", str(provenance_path),
        ]
    ) == 0

    expected_names = {
        "atlasem_dura_inputs", "atlasem_segmentation", "dura_stripped_t1",
        "dura_corrected_white_matter", "dura_corrected_gray_matter",
        "dura_brain_mask_level_set", "dura_brain_mask", "dura_wm_topology", "dura_gm_topology",
    }
    assert {path.stem.removesuffix(".nii") for path in output_dir.glob("*.nii.gz")} == expected_names
    for path in output_dir.glob("*.nii.gz"):
        image = nib.load(str(path))
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert np.isfinite(np.asarray(image.dataobj)).all()
    segmentation = nib.load(str(output_dir / "atlasem_segmentation.nii.gz"))
    assert segmentation.get_data_dtype() == np.dtype(np.int8)
    assert segmentation.shape == shape
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert provenance["command"] == "atlasem-prelesion-dura"
    assert provenance["spatial_convention"] == "RAS"
    assert provenance["parameters"] == {
        "class_names": ["Background", "CorticalGM", "WhiteMatter"],
        "max_iterations": 1,
    }
