#!/usr/bin/env python3
"""Replay the selected null-bound PLUS TGDM stage on saved held-fixed arrays."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import struct

import nibabel as nib
import numpy as np

from pymacruise.plus_pipeline import PlusPipelineOptions, run_plus_tgdm_stage


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_arrays(path: Path) -> tuple[tuple[int, int, int], tuple[np.ndarray, ...]]:
    payload = path.read_bytes()
    if len(payload) < 12:
        raise ValueError("Held-fixed TGDM input binary is truncated")
    shape = struct.unpack(">III", payload[:12])
    voxel_count = int(np.prod(shape))
    expected_bytes = 12 + voxel_count * 8 * 4
    if len(payload) != expected_bytes:
        raise ValueError(f"Held-fixed TGDM input length mismatch: expected {expected_bytes}, got {len(payload)}")
    raw = np.frombuffer(payload, dtype=">f4", offset=12)
    arrays = tuple(raw[i * voxel_count:(i + 1) * voxel_count].astype(np.float32).reshape(shape)
                   for i in range(5))
    vector = raw[5 * voxel_count:].astype(np.float32).reshape((*shape, 3))
    if any(not np.isfinite(array).all() for array in (*arrays, vector)):
        raise ValueError("Held-fixed TGDM inputs contain non-finite values")
    return shape, (*arrays, vector)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arrays", type=Path, required=True)
    parser.add_argument("--input-provenance", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    arrays_path = args.arrays.resolve()
    provenance_path = args.input_provenance.resolve()
    if not arrays_path.is_file() or not provenance_path.is_file():
        raise FileNotFoundError("Held-fixed arrays and provenance must be regular files")
    provenance = json.loads(provenance_path.read_text())
    if sha256(arrays_path) != provenance.get("held_fixed_input_binary_sha256"):
        raise ValueError("Held-fixed array hash differs from its source provenance")
    if provenance.get("tgdm_input_mode") != "initial_membership":
        raise ValueError("Selected null-bound replay requires the original initial-membership input mode")
    shape, arrays = read_arrays(arrays_path)
    if tuple(provenance.get("shape_xyz", ())) != shape:
        raise ValueError("Held-fixed array dimensions differ from source provenance")
    affine = np.asarray(provenance["affine_ras"], dtype=np.float64)
    if affine.shape != (4, 4) or not np.isfinite(affine).all():
        raise ValueError("Source provenance has no finite 4x4 RAS affine")
    options_record = dict(provenance["settings"])
    options_record["boundary_overlook_labels"] = tuple(options_record["boundary_overlook_labels"])
    options = PlusPipelineOptions(**options_record)
    images = [nib.Nifti1Image(array, affine) for array in arrays[:5]]
    vector = nib.Nifti1Image(arrays[5], affine)
    result = run_plus_tgdm_stage(*images, vector, options=options)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name in ("initial_level_set", "tissue_sign", "initial_inner", "inner", "central", "outer"):
        np.save(args.output_dir / f"python-{name}.npy", np.asarray(getattr(result, name).dataobj))
    np.save(args.output_dir / "python-smoothed_initial.npy",
            np.asarray(result.execution.smoothed_initial.level_set, dtype=np.float32))
    provenance["tgdm_input_mode"] = "initial_membership"
    provenance["tgdm_optional_brain_mask_level_set"] = "omitted_per_selected_MaTGDMcrop_wrapper"
    provenance["python_tgdm_replay"] = {
        "policy": "selected MaTGDMcrop null-bound route",
        "held_fixed_arrays": str(arrays_path),
        "held_fixed_arrays_sha256": sha256(arrays_path),
        "input_provenance_sha256": sha256(provenance_path),
        "python_replay_source": "scripts/replay_plus_tgdm_saved_arrays.py",
        "brain_mask_argument": None,
        "outputs": {
            name: {"path": str(args.output_dir / f"python-{name}.npy"),
                  "sha256": sha256(args.output_dir / f"python-{name}.npy")}
            for name in ("initial_level_set", "tissue_sign", "initial_inner", "inner", "central", "outer")
        },
    }
    (args.output_dir / "input-provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )
    print(f"Python null-bound TGDM replay complete: {args.output_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
