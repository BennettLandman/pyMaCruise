from __future__ import annotations

import gzip
import importlib.util
import struct
from pathlib import Path

import nibabel as nib
import numpy as np

_SCRIPT = Path(__file__).parents[1] / "scripts" / "compare_dura_current_stage.py"
_SPEC = importlib.util.spec_from_file_location("compare_dura_current_stage", _SCRIPT)
assert _SPEC is not None and _SPEC.loader is not None
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)
OUTPUTS = _MODULE.OUTPUTS
compare = _MODULE.compare


def _write_java_volume(path, data: np.ndarray, dtype: np.dtype) -> None:
    with gzip.open(path, "wb") as stream:
        stream.write(struct.pack(">III", *data.shape))
        stream.write(np.asarray(data, dtype=dtype).tobytes(order="C"))


def test_compare_dura_current_stage_uses_saved_filenames_and_bounds_mismatch_examples(tmp_path):
    java_dir = tmp_path / "java"
    python_dir = tmp_path / "python"
    java_dir.mkdir()
    python_dir.mkdir()
    shape = (12, 4, 4)
    affine = np.diag([1.0, 1.0, 1.0, 1.0])

    for name, (filename, dtype) in OUTPUTS.items():
        python_data = np.zeros(shape, dtype=np.float32 if dtype.kind == "f" else np.uint8)
        java_data = python_data.copy()
        if name == "brain_mask":
            python_data.ravel()[0] = 1
            java_data.ravel()[:12] = 1
        nib.save(nib.Nifti1Image(python_data, affine), python_dir / filename)
        _write_java_volume(java_dir / f"{name}.raw.gz", java_data, dtype)

    result = compare(java_dir, python_dir)
    mask = result["outputs"]["brain_mask"]
    assert mask["changed_voxels"] == 11
    assert mask["java_positive_python_zero"]["count"] == 11
    assert len(mask["java_positive_python_zero"]["first_xyz"]) == 10
    assert mask["java_zero_python_positive"]["count"] == 0
    assert result["outputs"]["wm_tc"]["changed_voxels"] == 0
    assert result["outputs"]["gm_tc"]["changed_voxels"] == 0
