from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.plus_lesion_native_dura_transfer import run_plus_lesion_native_dura_transfer


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_native_dura_transfer_exposes_physical_and_jist_unchanged_modes(tmp_path: Path) -> None:
    source_shape = (4, 5, 6)
    source_affine = np.diag([1.0, 1.5, 2.0, 1.0])
    reference_affine = np.diag([2.0, 2.0, 2.0, 1.0])
    reference_affine[:3, 3] = (-12.0, 4.0, -10.0)
    reference_path = tmp_path / "reference.nii.gz"
    nib.save(nib.Nifti1Image(np.zeros((8, 9, 10), dtype=np.float32), reference_affine), reference_path)

    names = ("csf", "gm", "wm_fill", "wm_mask")
    outputs = {}
    source_data = {}
    for index, name in enumerate(names):
        path = tmp_path / f"{name}.nii.gz"
        values = (np.arange(np.prod(source_shape), dtype=np.float32).reshape(source_shape) + index)
        source_data[name] = values
        nib.save(nib.Nifti1Image(values, source_affine), path)
        outputs[name] = {"path": str(path), "sha256": _sha256(path)}
    manifest = tmp_path / "dura.json"
    manifest.write_text(json.dumps({
        "status": "complete",
        "stage": "plus-lesion-dura-output-checkpoint",
        "spatial_convention": "RAS",
        "outputs": outputs,
    }))

    jist = run_plus_lesion_native_dura_transfer(
        manifest, reference_path, tmp_path / "jist", transfer_mode="jist-unchanged"
    )
    for name, path in jist.outputs.items():
        image = nib.load(str(path))
        assert image.shape == source_shape
        np.testing.assert_array_equal(image.dataobj, source_data[name])
        np.testing.assert_array_equal(image.affine, reference_affine)
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    jist_manifest = json.loads(jist.manifest_path.read_text())
    assert jist_manifest["transfer_mode"] == "jist-unchanged"
    assert jist_manifest["coordinate_mapping"].startswith("preserve RAS source extents")

    physical = run_plus_lesion_native_dura_transfer(
        manifest, reference_path, tmp_path / "physical-affine"
    )
    for path in physical.outputs.values():
        image = nib.load(str(path))
        assert image.shape == (8, 9, 10)
        np.testing.assert_array_equal(image.affine, reference_affine)
    assert json.loads(physical.manifest_path.read_text())["transfer_mode"] == "physical-affine"
