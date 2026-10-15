"""Compare selected Java Dura outputs with one saved Python Dura stage."""

from __future__ import annotations

import argparse
import gzip
import json
import struct
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy.ndimage import binary_erosion, distance_transform_edt


OUTPUTS = {
    "wm_tc": ("wm_tc.nii.gz", np.dtype(">f4")),
    "gm_tc": ("gm_tc.nii.gz", np.dtype(">f4")),
    "brain_mask_level_set": ("brain_mask_level_set.nii.gz", np.dtype(">f4")),
    "brain_mask": ("brain_mask.nii.gz", np.dtype("u1")),
    "stripped_t1": ("stripped_t1.nii.gz", np.dtype(">f4")),
    "corrected_wm": ("corrected_wm.nii.gz", np.dtype(">f4")),
    "corrected_gm": ("corrected_gm.nii.gz", np.dtype(">f4")),
}


def read_java_volume(path: Path, dtype: np.dtype) -> np.ndarray:
    with gzip.open(path, "rb") as stream:
        shape = struct.unpack(">III", stream.read(12))
        values = np.frombuffer(stream.read(), dtype=dtype)
    if values.size != int(np.prod(shape)):
        raise ValueError(f"invalid Java volume size in {path}")
    return values.reshape(shape).copy()


def mask_boundary(mask: np.ndarray) -> np.ndarray:
    six_connected = np.zeros((3, 3, 3), dtype=bool)
    six_connected[1, 1, 1] = True
    six_connected[0, 1, 1] = six_connected[2, 1, 1] = True
    six_connected[1, 0, 1] = six_connected[1, 2, 1] = True
    six_connected[1, 1, 0] = six_connected[1, 1, 2] = True
    return mask & ~binary_erosion(mask, structure=six_connected, border_value=0)


def mismatch_summary(mask: np.ndarray, *, limit: int = 10) -> dict[str, object]:
    """Bound JSON size while retaining a few xyz examples for diagnosis."""
    coordinates = np.argwhere(mask)
    return {
        "count": int(coordinates.shape[0]),
        "first_xyz": coordinates[:limit].tolist(),
    }


def compare(java_dir: Path, stage_dir: Path) -> dict[str, object]:
    records: dict[str, object] = {}
    python_mask_image = nib.load(str(stage_dir / "brain_mask.nii.gz"))
    if nib.aff2axcodes(python_mask_image.affine) != ("R", "A", "S"):
        raise ValueError("Python Dura output is not RAS-oriented")
    spacing = tuple(float(value) for value in python_mask_image.header.get_zooms()[:3])
    python_mask = np.asarray(python_mask_image.dataobj) > 0
    for name, (python_name, dtype) in OUTPUTS.items():
        java = read_java_volume(java_dir / f"{name}.raw.gz", dtype)
        image = python_mask_image if name == "brain_mask" else nib.load(str(stage_dir / python_name))
        python = np.asarray(image.dataobj)
        if java.shape != python.shape:
            raise ValueError(f"shape mismatch for {name}: {java.shape} != {python.shape}")
        if name == "brain_mask":
            java_mask = java > 0
            intersection = int(np.count_nonzero(java_mask & python_mask))
            denominator = int(np.count_nonzero(java_mask) + np.count_nonzero(python_mask))
            java_boundary = mask_boundary(java_mask)
            python_boundary = mask_boundary(python_mask)
            to_python = distance_transform_edt(~python_boundary, sampling=spacing)[java_boundary]
            to_java = distance_transform_edt(~java_boundary, sampling=spacing)[python_boundary]
            records[name] = {
                "changed_voxels": int(np.count_nonzero(java_mask != python_mask)),
                "dice": 2.0 * intersection / denominator if denominator else 1.0,
                "symmetric_boundary_voxel_center_mean_distance_mm": float(
                    (to_python.mean() + to_java.mean()) / 2.0
                ),
                "boundary_voxel_center_hausdorff_mm": float(
                    max(to_python.max(initial=0.0), to_java.max(initial=0.0))
                ),
                "java_positive_python_zero": mismatch_summary(java_mask & ~python_mask),
                "java_zero_python_positive": mismatch_summary(~java_mask & python_mask),
            }
            continue

        if not np.isfinite(java).all() or not np.isfinite(python).all():
            raise ValueError(f"non-finite Dura output in {name}")
        delta = np.abs(java.astype(np.float64) - python.astype(np.float64))
        record = {
            "changed_voxels": int(np.count_nonzero(java != python)),
            "max_abs": float(delta.max(initial=0.0)),
            "mean_abs": float(delta.mean()),
        }
        if name in {"wm_tc", "gm_tc", "brain_mask_level_set"}:
            record["zero_level_sign_changes"] = int(
                np.count_nonzero(np.signbit(java) != np.signbit(python))
            )
        records[name] = record
    return {
        "java_output_dir": str(java_dir),
        "python_stage_dir": str(stage_dir),
        "shape": list(python_mask.shape),
        "spacing_mm": list(spacing),
        "coordinate_system": "RAS",
        "mask_distance_definition": "symmetric 6-connected boundary-voxel-center samples",
        "outputs": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("java_output_dir", type=Path)
    parser.add_argument("python_stage_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(compare(args.java_output_dir, args.python_stage_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
