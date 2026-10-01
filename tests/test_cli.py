from __future__ import annotations

import json
import subprocess
import sys

import nibabel as nib
import numpy as np
import pytest

from pymacruise.cli import main
from pymacruise.topology import topology_correct
from pymacruise.volume_ops import calculate_volume


def test_cli_version_and_prepare_raw_commands(tmp_path, capsys) -> None:
    with pytest.raises(SystemExit) as exit_info:
        main(["--version"])
    assert exit_info.value.code == 0
    assert "pymacruise" in capsys.readouterr().out

    affine = np.eye(4)
    t1_path = tmp_path / "t1.nii.gz"
    maseg_path = tmp_path / "maseg.nii.gz"
    mask_path = tmp_path / "brain-mask.nii.gz"
    nib.save(nib.Nifti1Image(np.ones((4, 5, 6), dtype=np.float32), affine), t1_path)
    labels = np.zeros((4, 5, 6), dtype=np.uint8)
    labels[1:3, 1:4, 2:5] = 100
    nib.save(nib.Nifti1Image(labels, affine), maseg_path)
    brain_mask = np.zeros((4, 5, 6), dtype=np.uint8)
    brain_mask[1:3, 1:4, 1:5] = 255
    nib.save(nib.Nifti1Image(brain_mask, affine), mask_path)

    assert main(
        [
            "prepare-raw",
            "--t1",
            str(t1_path),
            "--maseg",
            str(maseg_path),
            "--brain-mask",
            str(mask_path),
            "--output-dir",
            str(tmp_path / "prepared"),
        ]
    ) == 0
    assert (tmp_path / "prepared" / "T1.nii.gz").is_file()
    assert (tmp_path / "prepared" / "Seg.nii.gz").is_file()
    assert (tmp_path / "prepared" / "BrainMask.nii.gz").is_file()


def test_calculate_cli_enforces_nii_gz_and_records_ras_operation(tmp_path) -> None:
    affine = np.array(
        [[-2.0, 0.0, 0.0, 8.0], [0.0, 3.0, 0.0, -2.0],
         [0.0, 0.0, 4.0, 5.0], [0.0, 0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    data_a = np.arange(60, dtype=np.float32).reshape(5, 4, 3)
    data_b = np.full((5, 4, 3), 2.5, dtype=np.float32)
    input_a = tmp_path / "a.nii.gz"
    input_b = tmp_path / "b.nii.gz"
    nib.save(nib.Nifti1Image(data_a, affine), input_a)
    nib.save(nib.Nifti1Image(data_b, affine), input_b)
    expected = calculate_volume(nib.load(input_a), nib.load(input_b), "add")
    output = tmp_path / "sum.nii.gz"
    provenance_path = tmp_path / "calculate.json"

    completed = subprocess.run(
        [
            sys.executable, "-m", "pymacruise.cli", "calculate",
            "--input-a", str(input_a), "--input-b", str(input_b),
            "--operation", " add ", "--output", str(output),
            "--provenance", str(provenance_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "reorienting the image to RAS" in completed.stderr
    written = nib.load(output)
    assert nib.aff2axcodes(written.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(np.asarray(written.dataobj), np.asarray(expected.dataobj))
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert provenance["command"] == "calculate"
    assert provenance["spatial_convention"] == "RAS"
    assert provenance["inputs"] == {
        "image_a": str(input_a.resolve()), "image_b": str(input_b.resolve())
    }
    assert provenance["parameters"] == {"operation": "ADD"}
    assert provenance["outputs"] == {"result": str(output.resolve())}

    uncompressed = tmp_path / "a.nii"
    nib.save(nib.Nifti1Image(data_a, affine), uncompressed)
    rejected = subprocess.run(
        [
            sys.executable, "-m", "pymacruise.cli", "calculate",
            "--input-a", str(uncompressed), "--input-b", str(input_b),
            "--operation", "ADD", "--output", str(tmp_path / "bad.nii.gz"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert rejected.returncode != 0
    assert "Image inputs must use .nii.gz" in rejected.stderr


def test_topology_correct_cli_matches_api_with_atlas_refinement(tmp_path) -> None:
    affine = np.diag([-1.0, 1.2, 1.5, 1.0])
    affine[0, 3] = 4.0
    grid = np.indices((5, 5, 5), dtype=np.float32)
    volume_data = (grid[0] * 0.2 + grid[1] * 0.1 + grid[2] * 0.05).astype(np.float32)
    paint_data = np.zeros((5, 5, 5), dtype=np.uint8)
    paint_data[2, 2, 2] = 1
    multi_atlas_wm_data = np.zeros((5, 5, 5), dtype=np.uint8)
    multi_atlas_wm_data[1:4, 1:4, 1:4] = 1
    inputs = {
        "volume": tmp_path / "volume.nii.gz",
        "paint": tmp_path / "paint.nii.gz",
        "wm": tmp_path / "multi-atlas-wm.nii.gz",
    }
    for key, values in (
        ("volume", volume_data),
        ("paint", paint_data),
        ("wm", multi_atlas_wm_data),
    ):
        nib.save(nib.Nifti1Image(values, affine), inputs[key])

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        expected = topology_correct(
            nib.load(inputs["volume"]),
            nib.load(inputs["paint"]),
            multi_atlas_wm=nib.load(inputs["wm"]),
            regularize=0.0,
        )

    corrected_path = tmp_path / "corrected.nii.gz"
    atlas_refined_path = tmp_path / "atlas-refined.nii.gz"
    provenance_path = tmp_path / "topology.json"
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "pymacruise.cli",
            "topology-correct",
            "--volume",
            str(inputs["volume"]),
            "--paint-mask",
            str(inputs["paint"]),
            "--multi-atlas-wm",
            str(inputs["wm"]),
            "--output",
            str(corrected_path),
            "--atlas-refined-output",
            str(atlas_refined_path),
            "--regularize",
            "0",
            "--provenance",
            str(provenance_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "reorienting the image to RAS" in completed.stderr
    corrected = nib.load(corrected_path)
    refined = nib.load(atlas_refined_path)
    assert nib.aff2axcodes(corrected.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(refined.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(corrected.affine, expected.corrected.affine)
    np.testing.assert_array_equal(corrected.dataobj, expected.corrected.dataobj)
    assert expected.atlas_refined is not None
    np.testing.assert_array_equal(refined.affine, expected.atlas_refined.affine)
    np.testing.assert_array_equal(refined.dataobj, expected.atlas_refined.dataobj)

    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert provenance["command"] == "topology-correct"
    assert provenance["spatial_convention"] == "RAS"
    assert provenance["inputs"] == {
        "volume": str(inputs["volume"].resolve()),
        "paint_mask": str(inputs["paint"].resolve()),
        "multi_atlas_wm": str(inputs["wm"].resolve()),
    }
    assert provenance["parameters"] == {
        "paint_threshold": 0.5,
        "lowest_level": 0.0,
        "highest_level": 1.0,
        "regularize": 0.0,
        "atlas_refinement_enabled": True,
    }
    assert provenance["outputs"] == {
        "corrected_volume": str(corrected_path.resolve()),
        "atlas_refined_volume": str(atlas_refined_path.resolve()),
    }


def test_topology_correct_cli_help_and_compressed_output_error(tmp_path) -> None:
    help_result = subprocess.run(
        [sys.executable, "-m", "pymacruise.cli", "topology-correct", "--help"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert help_result.returncode == 0
    assert "--multi-atlas-wm" in help_result.stdout
    assert "--regularize" in help_result.stdout

    values = np.ones((3, 3, 3), dtype=np.float32)
    path = tmp_path / "volume.nii.gz"
    nib.save(nib.Nifti1Image(values, np.eye(4)), path)
    rejected = subprocess.run(
        [
            sys.executable,
            "-m",
            "pymacruise.cli",
            "topology-correct",
            "--volume",
            str(path),
            "--paint-mask",
            str(path),
            "--output",
            str(tmp_path / "corrected.nii"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert rejected.returncode != 0
    assert "Image outputs must use .nii.gz" in rejected.stderr
