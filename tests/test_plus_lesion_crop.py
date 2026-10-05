from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.image_io import save_nifti_ras
import pymacruise.plus_lesion_crop as crop_stage
from pymacruise.lesion_preprocessing_crop import crop_lesion_toads_preprocessing_inputs
from pymacruise.lesion_toads_preprocessing import normalize_lesion_toads_images


def test_plus_lesion_crop_normalizes_pads_and_resumes(tmp_path: Path, monkeypatch) -> None:
    shape = (15, 13, 11)
    x, y, z = np.indices(shape, dtype=np.float32)
    mask = ((x - 7.0) / 6.0) ** 2 + ((y - 6.0) / 5.0) ** 2 + ((z - 5.0) / 4.0) ** 2 < 1.0
    t1 = np.zeros(shape, dtype=np.float32)
    flair = np.zeros(shape, dtype=np.float32)
    t1[mask] = (0.2 + 0.01 * x[mask] + 0.005 * z[mask]).astype(np.float32)
    flair[mask] = (0.1 + 0.007 * y[mask] + 0.009 * z[mask]).astype(np.float32)
    affine = np.diag((0.8, 0.9, 1.1, 1.0))
    modalities_path = save_nifti_ras(
        nib.Nifti1Image(np.stack((t1, flair), axis=-1), affine),
        tmp_path / "modalities_t1_flair.nii.gz",
    )
    topology = np.ones(shape, dtype=np.uint8)
    topology[4:11, 3:10, 2:9] = 2
    topology_path = save_nifti_ras(
        nib.Nifti1Image(topology, affine), tmp_path / "topology.nii.gz"
    )

    padding = 2
    output = tmp_path / "lesion-crop"
    result = crop_stage.run_plus_lesion_crop_checkpoint(
        modalities_path, topology_path, output, padding=padding
    )
    assert not result.resumed
    normalized_image = nib.load(str(result.normalized_path))
    cropped_image = nib.load(str(result.cropped_path))
    assert normalized_image.shape == (*shape, 2)
    assert nib.aff2axcodes(normalized_image.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(cropped_image.affine) == ("R", "A", "S")
    normalized_data = np.asarray(normalized_image.dataobj)
    cropped_data = np.asarray(cropped_image.dataobj)
    assert np.isfinite(normalized_data).all() and np.isfinite(cropped_data).all()
    assert normalized_data.min() == 0.0 and normalized_data.max() == 1.0
    assert cropped_data.shape == (
        result.crop_bounds[1] - result.crop_bounds[0] + 1 + 2 * padding,
        result.crop_bounds[3] - result.crop_bounds[2] + 1 + 2 * padding,
        result.crop_bounds[5] - result.crop_bounds[4] + 1 + 2 * padding,
        2,
    )
    offset = np.asarray(
        (
            result.crop_bounds[0] - padding,
            result.crop_bounds[2] - padding,
            result.crop_bounds[4] - padding,
            1.0,
        )
    )
    np.testing.assert_allclose(cropped_image.affine @ (0, 0, 0, 1), affine @ offset)

    direct_normalized = normalize_lesion_toads_images(np.stack((t1, flair)))
    direct_crop = crop_lesion_toads_preprocessing_inputs(
        direct_normalized.normalized_images,
        topology,
        np.zeros(6, dtype=np.float32),
        image_resolution_mm=(0.8, 0.9, 1.1),
        template_resolution_mm=(0.8, 0.9, 1.1),
        padding=padding,
    )
    assert result.crop_bounds == direct_crop.crop.bounds
    np.testing.assert_array_equal(
        cropped_data,
        direct_crop.crop.images.transpose(1, 2, 3, 0),
    )
    assert result.post_crop_transform == tuple(float(v) for v in direct_crop.crop.transform)
    manifest = json.loads(result.manifest_path.read_text())
    assert manifest["status"] == "complete"
    assert manifest["channel_order"] == ["T1", "FLAIR"]
    assert manifest["source_order"] == [
        "normalizeImages",
        "findCroppingBoundaries",
        "cropImages",
        "update crop transform translation",
    ]

    def unexpected(*args, **kwargs):
        raise AssertionError("valid resume should reuse verified normalized and cropped outputs")

    monkeypatch.setattr(crop_stage, "normalize_lesion_toads_images", unexpected)
    monkeypatch.setattr(crop_stage, "crop_lesion_toads_preprocessing_inputs", unexpected)
    resumed = crop_stage.run_plus_lesion_crop_checkpoint(
        modalities_path, topology_path, output, padding=padding, resume=True
    )
    assert resumed.resumed
    assert resumed.crop_bounds == result.crop_bounds
    with pytest.raises(ValueError, match="inputs or settings differ"):
        crop_stage.run_plus_lesion_crop_checkpoint(
            modalities_path, topology_path, output, padding=padding + 1, resume=True
        )


def test_plus_crop_lesion_cli_routes_padding_and_resume(tmp_path: Path, monkeypatch) -> None:
    from pymacruise.cli import main

    captured = {}

    def fake_crop(modalities, topology, output, **kwargs):
        captured.update(modalities=modalities, topology=topology, output=output, **kwargs)
        return crop_stage.PlusLesionCropResult(
            tmp_path / "manifest.json",
            tmp_path / "normalized.nii.gz",
            tmp_path / "cropped.nii.gz",
            (1, 2, 3, 4, 5, 6),
            (0.0,) * 6,
            resumed=True,
        )

    monkeypatch.setattr(crop_stage, "run_plus_lesion_crop_checkpoint", fake_crop)
    assert main(
        [
            "plus-crop-lesion-inputs",
            "--modalities", "modalities.nii.gz",
            "--topology", "topology.nii.gz",
            "--output-dir", "out",
            "--padding", "3",
            "--resume",
        ]
    ) == 0
    assert captured["padding"] == 3
    assert captured["resume"] is True
