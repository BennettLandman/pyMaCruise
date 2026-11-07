from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

import nibabel as nib
import numpy as np
import pytest

from pymacruise.n4 import N4Config, n4_bias_correct


REPO = Path(__file__).resolve().parents[1]
PHANTOM = REPO / "tests" / "data" / "synthetic" / "n4_cli" / "smooth_bias.nii.gz"


def _run_cli(*arguments: str, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "pymacruise.cli", *map(str, arguments)],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )


def test_n4_cli_matches_api_and_writes_ras_provenance(tmp_path: Path) -> None:
    pytest.importorskip("ants")
    input_path = tmp_path / "smooth-bias.nii.gz"
    output_path = tmp_path / "corrected.nii.gz"
    provenance_path = tmp_path / "n4-run.json"
    source = nib.load(PHANTOM)
    nib.save(source, input_path)

    config = N4Config(
        shrink_factor=2,
        iterations=(3,),
        convergence_tolerance=0.0,
        spline_mesh=(1, 1, 1),
    )
    expected = n4_bias_correct(source, config=config)
    completed = _run_cli(
        "n4",
        "--input", str(input_path),
        "--output", str(output_path),
        "--shrink-factor", "2",
        "--iterations", "3",
        "--provenance", str(provenance_path),
    )
    assert completed.returncode == 0, completed.stderr
    assert "reorienting the image to RAS" in completed.stderr
    written = nib.load(output_path)
    assert nib.aff2axcodes(written.affine) == ("R", "A", "S")
    np.testing.assert_allclose(written.affine, expected.affine, rtol=0, atol=1e-5)
    np.testing.assert_allclose(
        np.asarray(written.dataobj), np.asarray(expected.dataobj), rtol=0, atol=1e-5
    )

    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert provenance["command"] == "n4"
    assert provenance["spatial_convention"] == "RAS"
    assert provenance["inputs"] == {"intensity_image": str(input_path.resolve())}
    assert provenance["parameters"] == {
        "shrink_factor": 2,
        "iterations": [3],
        "convergence_tolerance": 0.0,
        "spline_mesh": [1, 1, 1],
        "rescale_intensities": False,
        "estimation_mask": "all",
    }
    assert provenance["outputs"] == {"corrected_image": str(output_path.resolve())}

    uncompressed_output = tmp_path / "corrected.nii"
    rejected = _run_cli("n4", "--input", str(input_path), "--output", str(uncompressed_output))
    assert rejected.returncode != 0
    assert "Image outputs must use .nii.gz" in rejected.stderr


def test_n4_cli_help_and_missing_ants_guidance_without_eager_import(tmp_path: Path) -> None:
    input_path = tmp_path / "input.nii.gz"
    nib.save(nib.load(PHANTOM), input_path)

    help_result = _run_cli("n4", "--help")
    assert help_result.returncode == 0
    assert "--shrink-factor" in help_result.stdout
    assert "--iterations" in help_result.stdout
    assert "--no-rescale-intensities" in help_result.stdout
    assert "--provenance" in help_result.stdout
    root_help = _run_cli("--help")
    assert root_help.returncode == 0
    assert "n4" in root_help.stdout

    # Mark ANTs as unavailable before importing PyMaCRUISE. The CLI and API
    # modules remain importable; only invoking N4 reports the optional extra.
    env = os.environ.copy()
    unavailable_ants = (
        "import sys; sys.modules['ants'] = None; "
        "from pymacruise.cli import main; "
        "main(['n4', '--input', sys.argv[1], '--output', sys.argv[2]])"
    )
    missing = subprocess.run(
        [sys.executable, "-c", unavailable_ants, str(input_path), str(tmp_path / "out.nii.gz")],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )
    assert missing.returncode == 2
    assert "install PyMaCRUISE with its 'ants' extra" in missing.stderr


def test_base_run_reports_and_records_missing_n4_dependency(tmp_path: Path) -> None:
    fixture = REPO / "tests/data/synthetic/toads_registered_first_iteration_cycle"
    output = tmp_path / "base-run"
    arguments = [
        "base-run", "--t1", str(fixture / "t1.nii.gz"),
        "--brain-mask", str(fixture / "brain-mask.nii.gz"),
        "--maseg", str(REPO / "tests/data/synthetic/base_workflow/maseg.nii.gz"),
        "--atlas-descriptor", str(REPO / "tests/data/synthetic/toads_registered_field_refit/atlas-11.txt"),
        "--output-dir", str(output),
    ]
    unavailable_ants = (
        "import sys; sys.modules['ants'] = None; "
        "from pymacruise.cli import main; sys.exit(main(sys.argv[1:]))"
    )
    for attempt, extra in ((1, []), (2, ["--resume"])):
        result = subprocess.run(
            [sys.executable, "-c", unavailable_ants, *arguments, *extra],
            cwd=REPO, capture_output=True, text=True, check=False,
        )
        assert result.returncode != 0
        assert "ANTsPy is required for N4 correction" in result.stderr
        manifest = json.loads((output / "run.json").read_text())
        assert manifest["status"] == "failed"
        assert manifest["stages"]["n4"]["status"] == "failed"
        assert manifest["stages"]["n4"]["attempt"] == attempt
        assert set(manifest["stages"]) == {"n4"}
        log = output / "stages/n4" / f"attempt-{attempt:02d}/stage.log"
        assert "ANTsPy is required for N4 correction" in log.read_text()
        assert not (log.parent / "t1-n4.nii.gz").exists()

def test_n4_native_estimation_precedes_orientation_only_working_export(tmp_path, monkeypatch):
    from pymacruise import cli
    from pymacruise.image_io import orient_to_ras_without_obliquity
    shape = (6, 7, 8)
    values = np.arange(np.prod(shape), dtype=np.float32).reshape(shape)
    angle = np.deg2rad(9)
    affine = np.eye(4)
    affine[:2, :2] = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    original = nib.Nifti1Image(values, affine)
    native_path = tmp_path / 'native.nii.gz'
    nib.save(original, native_path)
    with pytest.warns(UserWarning):
        working = orient_to_ras_without_obliquity(original)
    working_path = tmp_path / 'working.nii.gz'
    nib.save(working, working_path)
    seen = {}
    def estimate(image, *, config):
        seen['affine'] = image.affine.copy()
        return nib.Nifti1Image(np.asarray(image.dataobj)+10, image.affine)
    monkeypatch.setattr(cli, 'n4_bias_correct', estimate)
    output = tmp_path / 'corrected-working.nii.gz'
    assert cli.main(['n4', '--input', str(native_path), '--output', str(output),
                     '--working-reference', str(working_path), '--orientation-only-working-grid']) == 0
    np.testing.assert_array_equal(seen['affine'], nib.load(native_path).affine)
    saved = nib.load(output)
    np.testing.assert_array_equal(saved.affine, nib.load(working_path).affine)
    np.testing.assert_array_equal(np.asarray(saved.dataobj), values+10)
