from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.cli import main


def _write_synthetic_descriptor(directory: Path) -> Path:
    shape = (12, 12, 12)
    x, y, z = np.indices(shape, dtype=np.float32)
    radius = np.sqrt((x - 5.5) ** 2 + (y - 5.5) ** 2 + (z - 5.5) ** 2)
    topology = np.where(radius < 2.4, 3, np.where(radius < 4.6, 2, 1)).astype(np.uint8)
    background = np.ones(shape, dtype=np.float32)
    gm = np.exp(-((radius - 3.5) ** 2) / 2.5).astype(np.float32)
    wm = np.exp(-(radius**2) / 18.0).astype(np.float32)
    affine = np.eye(4)
    for name, data in (("topology", topology), ("background", background), ("gm", gm), ("wm", wm)):
        nib.save(nib.Nifti1Image(data, affine), directory / f"{name}.nii.gz")
    descriptor = directory / "atlas.txt"
    descriptor.write_text(
        "\n".join(
            (
                "Structure Atlas File (edit at your own risks)", "", "Structures:\t3",
                "Background\t1\tmask", "Cerebrum-GM\t2\tclosed", "Cerebrum-WM\t3\tclosed",
                "", "Topology Atlas:", "File:\ttopology.nii.gz", "Dimensions:\t12\t12\t12",
                "Resolutions:\t1\t1\t1", "", "Shape Atlas:\t4", "Dimensions:\t12\t12\t12",
                "Resolutions:\t1\t1\t1", "Orientations:\t1\t1\t6\t3",
                "Structure:\tBackground", "File:\tbackground.nii.gz", "Structure:\tCerebrum-GM",
                "File:\tgm.nii.gz", "Structure:\tCerebrum-WM", "File:\twm.nii.gz",
                "Registered Shapes:", "0\t1\t1", "0\t1\t1", "", "Intensity Atlas:\t10",
                "T1_MPRAGE\t0.2\t0.5\t0.8",
            )
        ) + "\n",
        encoding="utf-8",
    )
    return descriptor


def test_atlasem_run_cli_executes_synthetic_descriptor_workflow(tmp_path: Path) -> None:
    shape = (12, 12, 12)
    x, y, z = np.indices(shape, dtype=np.float32)
    radius = np.sqrt((x - 5.5) ** 2 + (y - 5.5) ** 2 + (z - 5.5) ** 2)
    t1_data = np.where(radius < 4.6, 0.35 + 0.025 * x + 0.015 * y, 0.0).astype(np.float32)
    t1_path = tmp_path / "t1.nii.gz"
    nib.save(nib.Nifti1Image(t1_data, np.eye(4)), t1_path)
    descriptor = _write_synthetic_descriptor(tmp_path)
    output = tmp_path / "out"
    provenance = tmp_path / "run.json"

    status = main(
        [
            "atlasem-run", "--t1", str(t1_path), "--descriptor", str(descriptor),
            "--output-dir", str(output), "--crop-border", "2", "--levels", "2",
            "--initial-registration-iterations", "1", "--post-recrop-registration-iterations", "1",
            "--max-em-iterations", "1", "--provenance", str(provenance),
        ]
    )

    assert status == 0
    expected = {
        "atlasem_memberships.nii.gz", "atlasem_segmentation.nii.gz",
        "atlasem_internal_mask.nii.gz", "atlasem_dura_inputs.nii.gz",
    }
    assert {path.name for path in output.iterdir()} == expected
    assert all(nib.aff2axcodes(nib.load(str(output / name)).affine) == ("R", "A", "S") for name in expected)
    assert all(np.isfinite(np.asarray(nib.load(str(output / name)).dataobj)).all() for name in expected)
    manifest = json.loads(provenance.read_text(encoding="utf-8"))
    assert manifest["command"] == "atlasem-run"
    assert manifest["parameters"]["em_iterations_completed"] == 1
