from __future__ import annotations

import struct
import json
import hashlib
from pathlib import Path
from types import SimpleNamespace

import nibabel as nib
import numpy as np
import pytest

import scripts.replay_plus_tgdm_saved_arrays as replay
from scripts.replay_plus_tgdm_saved_arrays import read_arrays
from pymacruise.plus_pipeline import PlusPipelineOptions


def test_saved_tgdm_array_reader_checks_shape_and_preserves_scalar_vector_order(tmp_path) -> None:
    shape = (2, 3, 4)
    scalar_count = int(np.prod(shape))
    arrays = [np.full(shape, i, dtype=np.float32) for i in range(5)]
    arrays.append(np.full((*shape, 3), 5, dtype=np.float32))
    path = tmp_path / "arrays.bin"
    with path.open("wb") as stream:
        stream.write(struct.pack(">III", *shape))
        for array in arrays:
            stream.write(np.asarray(array, dtype=">f4").tobytes(order="C"))
    found_shape, found = read_arrays(path)
    assert found_shape == shape
    assert len(found) == 6
    for index, array in enumerate(found[:5]):
        np.testing.assert_array_equal(array, index)
    assert found[5].shape == (*shape, 3)
    np.testing.assert_array_equal(found[5], 5)


def test_saved_tgdm_array_reader_rejects_truncated_inputs(tmp_path) -> None:
    path = tmp_path / "truncated.bin"
    path.write_bytes(struct.pack(">III", 2, 2, 2) + b"short")
    with pytest.raises(ValueError, match="length mismatch"):
        read_arrays(path)


def test_replay_uses_initial_membership_and_omits_optional_mask(tmp_path, monkeypatch) -> None:
    shape = (2, 2, 2)
    arrays = [np.full(shape, i + 0.25, dtype=np.float32) for i in range(5)]
    vector = np.ones((*shape, 3), dtype=np.float32)
    arrays_path = tmp_path / "arrays.bin"
    with arrays_path.open("wb") as stream:
        stream.write(struct.pack(">III", *shape))
        for array in (*arrays, vector):
            stream.write(np.asarray(array, dtype=">f4").tobytes(order="C"))
    digest = hashlib.sha256(arrays_path.read_bytes()).hexdigest()
    provenance_path = tmp_path / "input-provenance.json"
    provenance_path.write_text(json.dumps({
        "held_fixed_input_binary_sha256": digest,
        "shape_xyz": list(shape),
        "affine_ras": np.eye(4).tolist(),
        "settings": PlusPipelineOptions().as_dict(),
        "tgdm_input_mode": "initial_membership",
        "tgdm_optional_brain_mask_level_set": "omitted_per_selected_MaTGDMcrop_wrapper",
    }))
    captured = {}

    def fake_run(*args, **kwargs):
        captured["args"] = args
        captured["kwargs"] = kwargs
        image = nib.Nifti1Image(np.zeros(shape, dtype=np.float32), np.eye(4))
        return SimpleNamespace(
            initial_level_set=image, tissue_sign=image, initial_inner=image,
            inner=image, central=image, outer=image,
            execution=SimpleNamespace(smoothed_initial=SimpleNamespace(level_set=np.zeros(shape))),
        )

    monkeypatch.setattr(replay, "run_plus_tgdm_stage", fake_run)
    output_dir = tmp_path / "replay"
    assert replay.main([
        "--arrays", str(arrays_path), "--input-provenance", str(provenance_path),
        "--output-dir", str(output_dir),
    ]) == 0
    assert len(captured["args"]) == 6
    assert "brain_mask" not in captured["kwargs"]
    np.testing.assert_array_equal(captured["args"][0].dataobj, arrays[0])
    result_provenance = json.loads((output_dir / "input-provenance.json").read_text())
    assert result_provenance["tgdm_input_mode"] == "initial_membership"
    assert result_provenance["python_tgdm_replay"]["brain_mask_argument"] is None


def test_replay_rejects_precomputed_levelset_mode(tmp_path) -> None:
    arrays_path = tmp_path / "arrays.bin"
    arrays_path.write_bytes(b"")
    provenance_path = tmp_path / "input-provenance.json"
    provenance_path.write_text(json.dumps({
        "tgdm_input_mode": "precomputed_initial_level_set",
        "held_fixed_input_binary_sha256": hashlib.sha256(b"").hexdigest(),
    }))
    with pytest.raises(ValueError, match="initial-membership"):
        replay.main([
            "--arrays", str(arrays_path), "--input-provenance", str(provenance_path),
            "--output-dir", str(tmp_path / "replay"),
        ])
