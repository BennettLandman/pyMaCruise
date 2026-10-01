from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_dura import run_dura_removal


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "tests" / "data" / "synthetic" / "dura_topology_two_pass"


def test_dura_cli_runs_one_active_outer_pass_and_matches_api(tmp_path: Path) -> None:
    wm_path = FIXTURE / "wm_input.nii.gz"
    gm_path = FIXTURE / "raw_gm.nii.gz"
    mask_path = FIXTURE / "paint_mask.nii.gz"
    wm = nib.load(wm_path)
    gm = nib.load(gm_path)
    mask = nib.load(mask_path)
    t1_path = tmp_path / "synthetic-t1.nii.gz"
    t1_values = np.arange(np.prod(wm.shape), dtype=np.float32).reshape(wm.shape) + 20
    nib.save(nib.Nifti1Image(t1_values, wm.affine), t1_path)

    expected = run_dura_removal(
        nib.load(t1_path), wm, gm, mask, max_iterations=1
    )
    output_dir = tmp_path / "outputs"
    provenance_path = tmp_path / "provenance.json"
    command = [
        sys.executable,
        "-m",
        "pymacruise.cli",
        "dura",
        "--t1",
        str(t1_path),
        "--white",
        str(wm_path),
        "--gray",
        str(gm_path),
        "--white-mask",
        str(mask_path),
        "--output-dir",
        str(output_dir),
        "--max-iterations",
        "1",
        "--provenance",
        str(provenance_path),
    ]
    completed = subprocess.run(
        command, cwd=REPO, capture_output=True, text=True, check=False
    )
    assert completed.returncode == 0, (
        f"CLI failed ({completed.returncode}): {command}\n"
        f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    )

    output_names = (
        "stripped_t1",
        "corrected_white_matter",
        "corrected_gray_matter",
        "brain_mask_level_set",
        "brain_mask",
        "wm_topology",
        "gm_topology",
    )
    for name in output_names:
        output_path = output_dir / f"{name}.nii.gz"
        assert output_path.is_file()
        actual = nib.load(output_path)
        assert nib.aff2axcodes(actual.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(
            np.asarray(actual.dataobj),
            np.asarray(getattr(expected, name).dataobj),
        )
        np.testing.assert_array_equal(actual.affine, wm.affine)

    provenance = json.loads(provenance_path.read_text())
    assert provenance["command"] == "dura"
    assert provenance["spatial_convention"] == "RAS"
    assert provenance["parameters"]["max_iterations"] == 1
    assert provenance["inputs"] == {
        "t1": str(t1_path.resolve()),
        "white_matter": str(wm_path.resolve()),
        "gray_matter": str(gm_path.resolve()),
        "white_matter_mask": str(mask_path.resolve()),
    }
    assert provenance["outputs"] == {
        name: str((output_dir / f"{name}.nii.gz").resolve())
        for name in output_names
    }
