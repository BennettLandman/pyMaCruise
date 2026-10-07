"""Compare selected Java Dura outputs with Python on the brain-resolution base inputs."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import struct
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_dura import dura_initial_level_sets, run_dura_removal
from pymacruise.validation import dice_coefficient


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/data/synthetic/base_workflow_brain_resolution/dura_java_inputs"
SHAPE = (172, 216, 164)
FLOAT_OUTPUTS = (
    "wm_tc", "union", "gm_tc", "inner_distance", "outer_distance", "level_set",
    "stripped_t1", "corrected_wm", "corrected_gm", "brain_mask_level_set", "tgdm_pass_01",
)


def read_java_volume(path: Path, *, byte: bool = False) -> np.ndarray:
    dtype = np.uint8 if byte else np.dtype(">f4")
    with gzip.open(path, "rb") as stream:
        header = stream.read(12)
        if len(header) != 12 or struct.unpack(">III", header) != SHAPE:
            raise ValueError(f"{path} has an invalid Java volume header")
        payload = stream.read()
    expected_bytes = int(np.prod(SHAPE)) * np.dtype(dtype).itemsize
    if len(payload) != expected_bytes:
        raise ValueError(f"{path} contains {len(payload)} bytes, expected {expected_bytes}")
    return np.frombuffer(payload, dtype=dtype).reshape(SHAPE).astype(np.float32 if not byte else np.uint8)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("java_output_dir", type=Path)
    parser.add_argument("--base-stage", type=Path, help="also compare the regenerated Python Dura outputs with a saved workflow stage")
    args = parser.parse_args()

    images = {name: nib.load(FIXTURE / f"{name}.nii.gz") for name in ("t1", "white_matter", "gray_matter", "white_matter_mask")}
    reference = images["t1"]
    for name, image in images.items():
        if image.shape != SHAPE or nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            raise ValueError(f"{name} input is not the expected 172x216x164 RAS grid")
        if not np.array_equal(image.affine, reference.affine):
            raise ValueError(f"{name} affine does not match T1")
    result = run_dura_removal(
        images["t1"], images["white_matter"], images["gray_matter"], images["white_matter_mask"],
        max_iterations=1,
        record_passes=True,
    )
    initial = dura_initial_level_sets(
        np.asarray(result.wm_topology.dataobj, dtype=np.float32),
        np.asarray(result.gm_topology.dataobj, dtype=np.float32),
    )
    gm = np.asarray(images["gray_matter"].dataobj, dtype=np.float32)
    python_outputs = {
        "wm_tc": np.asarray(result.wm_topology.dataobj, dtype=np.float32),
        "union": np.add(np.asarray(result.wm_topology.dataobj, dtype=np.float32), gm, dtype=np.float32),
        "gm_tc": np.asarray(result.gm_topology.dataobj, dtype=np.float32),
        "inner_distance": initial.inner_signed_distance,
        "outer_distance": initial.outer_signed_distance,
        "level_set": np.asarray(result.evolution.level_set, dtype=np.float32),
        "stripped_t1": np.asarray(result.stripped_t1.dataobj, dtype=np.float32),
        "corrected_wm": np.asarray(result.corrected_white_matter.dataobj, dtype=np.float32),
        "corrected_gm": np.asarray(result.corrected_gray_matter.dataobj, dtype=np.float32),
        "brain_mask_level_set": np.asarray(result.brain_mask_level_set.dataobj, dtype=np.float32),
        "brain_mask": np.asarray(result.brain_mask.dataobj, dtype=np.uint8),
        "tgdm_pass_01": np.asarray(result.evolution.pass_level_sets[0], dtype=np.float32),
    }
    output_root = FIXTURE / "java_reference"
    output_root.mkdir(parents=True, exist_ok=True)
    comparisons = {}
    for name, candidate in python_outputs.items():
        java = read_java_volume(args.java_output_dir / f"{name}.raw.gz", byte=name == "brain_mask")
        error = np.abs(candidate.astype(np.float64) - java.astype(np.float64))
        record = {
            "exact": bool(np.array_equal(candidate, java)),
            "changed_voxels": int(np.count_nonzero(candidate != java)),
            "max_absolute_error": float(error.max(initial=0)),
            "mean_absolute_error": float(error.mean()),
        }
        if name in {"inner_distance", "outer_distance", "level_set", "brain_mask_level_set", "tgdm_pass_01"}:
            record["sign_changed_voxels"] = int(np.count_nonzero(np.signbit(candidate) != np.signbit(java)))
        if name == "brain_mask":
            record["dice_foreground"] = dice_coefficient(candidate > 0, java > 0)
        elif name in {"wm_tc", "gm_tc", "level_set", "tgdm_pass_01"}:
            record["dice_at_zero"] = dice_coefficient(candidate <= 0, java <= 0)
        if args.base_stage:
            base_name = {
                "wm_tc": "wm_topology.nii.gz",
                "gm_tc": "gm_topology.nii.gz",
                "level_set": "brain_mask_level_set.nii.gz",
                "brain_mask_level_set": "brain_mask_level_set.nii.gz",
                "brain_mask": "brain_mask.nii.gz",
                "stripped_t1": "stripped_t1.nii.gz",
                "corrected_wm": "corrected_white_matter.nii.gz",
                "corrected_gm": "corrected_gray_matter.nii.gz",
            }.get(name)
            if base_name:
                base = np.asarray(nib.load(args.base_stage / base_name).dataobj)
                record["regenerated_python_matches_base_run"] = bool(np.array_equal(candidate, base))
        comparisons[name] = record
        save_name = "brain_mask.nii.gz" if name == "brain_mask" else f"{name}.nii.gz"
        output = nib.Nifti1Image(java.astype(np.uint8 if name == "brain_mask" else np.float32), reference.affine)
        nib.save(output, output_root / save_name)

    manifest_paths = [*(FIXTURE / f"{name}.nii.gz" for name in images), *(output_root / f"{name}.nii.gz" for name in (*FLOAT_OUTPUTS, "brain_mask"))]
    manifest = {str(path.relative_to(FIXTURE)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(manifest_paths)}
    (FIXTURE / "SHA256SUMS").write_text("".join(f"{digest}  {path}\n" for path, digest in manifest.items()))
    report = {
        "synthetic_only": True,
        "workflow_stage": "first TOADS outputs -> selected Dura wrapper -> second TOADS input",
        "shape": list(SHAPE),
        "coordinate_system": "RAS",
        "parameters": {"membership_level": 0.4, "mask_translation": 2.0, "level_set_translation": 0.0, "tgdm_outer_iterations": 1},
        "outputs": comparisons,
    }
    (FIXTURE / "comparison.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
