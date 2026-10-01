from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import nibabel as nib
import numpy as np

from pymacruise.segmentation_refine import refine_segmentation_by_surfaces


REPO = Path(__file__).resolve().parents[1]
FIXTURES = REPO / "tests" / "data" / "synthetic" / "refinement"


def _run_cli(*arguments: str | Path) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        [sys.executable, "-m", "pymacruise.cli", *map(str, arguments)],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, (
        f"CLI failed ({completed.returncode}): {completed.args}\n"
        f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    )
    return completed


def _ras_inputs() -> tuple[nib.spatialimages.SpatialImage, ...]:
    return tuple(
        nib.load(FIXTURES / f"{name}.nii.gz")
        for name in ("segmentation", "inner_levelset", "outer_levelset")
    )


def _save_flipped_nonras_inputs(directory: Path) -> tuple[Path, Path, Path]:
    directory.mkdir(parents=True)
    result = []
    for image, name in zip(
        _ras_inputs(), ("segmentation", "inner", "outer"), strict=True
    ):
        affine = image.affine.copy()
        data = np.asarray(image.dataobj)[::-1]
        affine[:3, 0] *= -1
        affine[0, 3] += (image.shape[0] - 1) * image.affine[0, 0]
        path = directory / f"{name}.nii.gz"
        nib.save(nib.Nifti1Image(data, affine), path)
        result.append(path)
    return tuple(result)


def test_segmentation_refine_cli_preserves_default_and_explicit_corrected_mode(
    tmp_path: Path,
) -> None:
    ras_inputs = _ras_inputs()
    initial_hashes = [
        hashlib.sha256((FIXTURES / f"{name}.nii.gz").read_bytes()).hexdigest()
        for name in ("segmentation", "inner_levelset", "outer_levelset")
    ]
    expected_faithful = refine_segmentation_by_surfaces(*ras_inputs)
    expected_corrected = refine_segmentation_by_surfaces(
        *ras_inputs, source_faithful=False
    )
    faithful_output = tmp_path / "faithful.nii.gz"
    faithful_summary = tmp_path / "faithful-summary.json"
    faithful_provenance = tmp_path / "faithful-provenance.json"
    command = _run_cli(
        "segmentation-refine",
        "--segmentation",
        FIXTURES / "segmentation.nii.gz",
        "--inner-levelset",
        FIXTURES / "inner_levelset.nii.gz",
        "--outer-levelset",
        FIXTURES / "outer_levelset.nii.gz",
        "--output",
        faithful_output,
        "--summary",
        faithful_summary,
        "--provenance",
        faithful_provenance,
    )
    assert "source_faithful=true" in command.stdout
    faithful_written = nib.load(faithful_output)
    assert nib.aff2axcodes(faithful_written.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(
        np.asarray(faithful_written.dataobj), np.asarray(expected_faithful.image.dataobj)
    )
    faithful_report = json.loads(faithful_summary.read_text(encoding="utf-8"))
    assert faithful_report["source_faithful"] is True
    assert faithful_report["rules"]["inner_wm"]["candidate_voxels"] == int(
        expected_faithful.candidate_masks["inner_wm"].sum()
    )
    assert faithful_report["changed_voxels"] == expected_faithful.changed_voxels
    assert faithful_report["unresolved_voxels"] == expected_faithful.unresolved_voxels
    provenance = json.loads(faithful_provenance.read_text(encoding="utf-8"))
    assert provenance["command"] == "segmentation-refine"
    assert provenance["spatial_convention"] == "RAS"
    assert provenance["parameters"]["source_faithful"] is True

    corrected_output = tmp_path / "corrected.nii.gz"
    corrected_summary = tmp_path / "corrected-summary.json"
    corrected_command = _run_cli(
        "segmentation-refine",
        "--segmentation",
        FIXTURES / "segmentation.nii.gz",
        "--inner-levelset",
        FIXTURES / "inner_levelset.nii.gz",
        "--outer-levelset",
        FIXTURES / "outer_levelset.nii.gz",
        "--output",
        corrected_output,
        "--summary",
        corrected_summary,
        "--no-source-faithful",
    )
    assert "source_faithful=false" in corrected_command.stdout
    corrected_written = nib.load(corrected_output)
    np.testing.assert_array_equal(
        np.asarray(corrected_written.dataobj),
        np.asarray(expected_corrected.image.dataobj),
    )
    corrected_report = json.loads(corrected_summary.read_text(encoding="utf-8"))
    assert corrected_report["source_faithful"] is False
    assert np.asarray(faithful_written.dataobj)[3, 3, 3] == 1
    assert np.asarray(corrected_written.dataobj)[3, 3, 3] == 44
    final_hashes = [
        hashlib.sha256((FIXTURES / f"{name}.nii.gz").read_bytes()).hexdigest()
        for name in ("segmentation", "inner_levelset", "outer_levelset")
    ]
    assert initial_hashes == final_hashes


def test_segmentation_refine_cli_reorients_nonras_inputs_with_warning(
    tmp_path: Path,
) -> None:
    paths = _save_flipped_nonras_inputs(tmp_path / "nonras")
    output = tmp_path / "nonras-output.nii.gz"
    summary = tmp_path / "summary.json"
    completed = _run_cli(
        "segmentation-refine",
        "--segmentation",
        paths[0],
        "--inner-levelset",
        paths[1],
        "--outer-levelset",
        paths[2],
        "--output",
        output,
        "--summary",
        summary,
    )
    assert "reorienting the image to RAS" in completed.stderr
    written = nib.load(output)
    expected = refine_segmentation_by_surfaces(*_ras_inputs())
    assert nib.aff2axcodes(written.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(np.asarray(written.dataobj), np.asarray(expected.image.dataobj))


def test_segmentation_refine_cli_rejects_uncompressed_inputs(tmp_path: Path) -> None:
    segmentation, inner, outer = _ras_inputs()
    plain = tmp_path / "segmentation.nii"
    nib.save(segmentation, plain)
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pymacruise.cli",
            "segmentation-refine",
            "--segmentation",
            str(plain),
            "--inner-levelset",
            str(FIXTURES / "inner_levelset.nii.gz"),
            "--outer-levelset",
            str(FIXTURES / "outer_levelset.nii.gz"),
            "--output",
            str(tmp_path / "out.nii.gz"),
            "--summary",
            str(tmp_path / "summary.json"),
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "must use .nii.gz" in result.stderr
