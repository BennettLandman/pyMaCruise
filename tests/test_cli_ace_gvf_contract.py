from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import nibabel as nib
import numpy as np
import pytest

from pymacruise.ace import ace_enhance
from pymacruise.cli import main
from pymacruise.gvf import generalized_gvf_image


REPO = Path(__file__).resolve().parents[1]
SYNTHETIC = REPO / "tests" / "data" / "synthetic"


def _run_cli(*arguments: str | Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "pymacruise.cli", *map(str, arguments)],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )


def _save_lps_equivalent(image: nib.spatialimages.SpatialImage, path: Path) -> None:
    """Store the same physical samples in a flipped-x, non-RAS voxel order."""
    data = np.asarray(image.dataobj)
    flip = np.eye(4)
    flip[0, 0] = -1.0
    flip[0, 3] = data.shape[0] - 1
    affine = image.affine @ flip
    nib.save(nib.Nifti1Image(data[::-1, ...].copy(), affine), path)


def test_ace_cli_reorients_inputs_and_matches_api_with_ras_provenance(tmp_path: Path) -> None:
    source = SYNTHETIC / "folded_cortex"
    ras_images = {
        name: nib.load(source / f"{name}.nii.gz")
        for name in ("gray_matter", "white_matter", "multi_atlas_csf")
    }
    lps_paths = {name: tmp_path / f"{name}-lps.nii.gz" for name in ras_images}
    for name, image in ras_images.items():
        _save_lps_equivalent(image, lps_paths[name])

    with pytest.warns(UserWarning, match="reorienting the image to RAS") as api_warnings:
        expected = ace_enhance(
            *(nib.load(lps_paths[name]) for name in ("gray_matter", "white_matter")),
            multi_atlas_csf=nib.load(lps_paths["multi_atlas_csf"]),
            diffusion_iterations=0,
        )
    assert len(api_warnings) == 3

    outputs = {
        "enhanced": tmp_path / "enhanced.nii.gz",
        "skeleton": tmp_path / "skeleton.nii.gz",
        "thinned": tmp_path / "thinned.nii.gz",
    }
    provenance_path = tmp_path / "ace.json"
    completed = _run_cli(
        "ace",
        "--gray", lps_paths["gray_matter"],
        "--white", lps_paths["white_matter"],
        "--atlas-csf", lps_paths["multi_atlas_csf"],
        "--enhanced", outputs["enhanced"],
        "--skeleton", outputs["skeleton"],
        "--thinned-skeleton", outputs["thinned"],
        "--diffusion-iterations", "0",
        "--provenance", provenance_path,
    )
    assert completed.returncode == 0, completed.stderr
    assert completed.stderr.count("reorienting the image to RAS") == 3
    expected_names = {
        "enhanced": "enhanced_gray_matter",
        "skeleton": "skeleton",
        "thinned": "thinned_skeleton",
    }
    for output_name, result_name in expected_names.items():
        actual = nib.load(outputs[output_name])
        expected_image = getattr(expected, result_name)
        assert nib.aff2axcodes(actual.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(actual.affine, expected_image.affine)
        np.testing.assert_array_equal(actual.dataobj, expected_image.dataobj)

    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert provenance["command"] == "ace"
    assert provenance["spatial_convention"] == "RAS"
    assert provenance["inputs"] == {
        "gray": str(lps_paths["gray_matter"].resolve()),
        "white": str(lps_paths["white_matter"].resolve()),
        "atlas_csf": str(lps_paths["multi_atlas_csf"].resolve()),
    }
    assert provenance["parameters"]["diffusion_iterations"] == 0
    assert provenance["parameters"]["merge_method"] == "Merge"
    assert provenance["outputs"] == {
        "enhanced_gray_matter": str(outputs["enhanced"].resolve()),
        "skeleton": str(outputs["skeleton"].resolve()),
        "thinned_skeleton": str(outputs["thinned"].resolve()),
    }


def test_gvf_cli_reorients_to_ras_and_preserves_same_grid_vector_result(tmp_path: Path) -> None:
    source = nib.load(SYNTHETIC / "plus_gvf" / "scalar.nii.gz")
    lps_path = tmp_path / "scalar-lps.nii.gz"
    _save_lps_equivalent(source, lps_path)
    parameters = {"lambda_value": 0.8, "max_iterations": 35, "solver_tolerance": 1e-4}
    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        expected = generalized_gvf_image(nib.load(lps_path), **parameters).vector_field

    output = tmp_path / "vectors.nii.gz"
    provenance_path = tmp_path / "gvf.json"
    completed = _run_cli(
        "gvf", "--input", lps_path, "--output", output,
        "--lambda", str(parameters["lambda_value"]),
        "--max-iterations", str(parameters["max_iterations"]),
        "--solver-tolerance", str(parameters["solver_tolerance"]),
        "--provenance", provenance_path,
    )
    assert completed.returncode == 0, completed.stderr
    assert "reorienting the image to RAS" in completed.stderr
    actual = nib.load(output)
    assert actual.shape == expected.shape
    assert actual.shape[-1] == 3
    assert nib.aff2axcodes(actual.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(actual.affine, expected.affine)
    np.testing.assert_array_equal(actual.dataobj, expected.dataobj)
    assert actual.header.get_intent()[0] == "vector"

    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert provenance["command"] == "gvf"
    assert provenance["spatial_convention"] == "RAS"
    assert provenance["inputs"] == {"scalar_image": str(lps_path.resolve())}
    assert provenance["parameters"] == parameters
    assert provenance["outputs"] == {"vector_field": str(output.resolve())}


def test_ace_and_gvf_cli_help_grid_errors_and_compressed_output_contract(tmp_path: Path) -> None:
    ace_help = _run_cli("ace", "--help")
    gvf_help = _run_cli("gvf", "--help")
    assert ace_help.returncode == gvf_help.returncode == 0
    assert "--atlas-csf" in ace_help.stdout
    assert "--diffusion-iterations" in ace_help.stdout
    assert "--solver-tolerance" in gvf_help.stdout
    assert "--provenance" in ace_help.stdout and "--provenance" in gvf_help.stdout

    image = nib.load(SYNTHETIC / "folded_cortex" / "gray_matter.nii.gz")
    data = np.asarray(image.dataobj)
    wrong_grid_affine = image.affine.copy()
    wrong_grid_affine[0, 3] += 2.0
    wrong_grid = nib.Nifti1Image(data, wrong_grid_affine)
    wrong_grid_path = tmp_path / "wrong-grid.nii.gz"
    nib.save(wrong_grid, wrong_grid_path)
    with pytest.raises(ValueError, match="share one physical grid"):
        main([
            "ace", "--gray", str(SYNTHETIC / "folded_cortex" / "gray_matter.nii.gz"),
            "--white", str(wrong_grid_path), "--enhanced", str(tmp_path / "e.nii.gz"),
            "--skeleton", str(tmp_path / "s.nii.gz"),
            "--thinned-skeleton", str(tmp_path / "t.nii.gz"),
            "--diffusion-iterations", "0",
        ])

    with pytest.raises(ValueError, match="Image outputs must use .nii.gz"):
        main([
            "gvf", "--input", str(SYNTHETIC / "plus_gvf" / "scalar.nii.gz"),
            "--output", str(tmp_path / "vectors.nii"),
        ])
