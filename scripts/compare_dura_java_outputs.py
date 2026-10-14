"""Compare a selected Java Dura harness run against the Python Dura API.

The Java output directory must come from
``scripts/run_java_dura_base_reference.sh`` using matching dimensions,
iteration count, and membership level. Inputs are four same-grid RAS NIfTIs:
``t1``, ``white_matter``, ``gray_matter``, and ``white_matter_mask``.
"""

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
from pymacruise.surfaces import extract_levelset_surface
from pymacruise.validation import mesh_surface_distance_metrics


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_java_volume(path: Path, shape: tuple[int, int, int], *, byte: bool = False) -> np.ndarray:
    dtype = np.uint8 if byte else np.dtype(">f4")
    with gzip.open(path, "rb") as stream:
        header = stream.read(12)
        if len(header) != 12 or struct.unpack(">III", header) != shape:
            raise ValueError(f"{path} has an invalid volume header")
        payload = stream.read()
    expected = int(np.prod(shape)) * np.dtype(dtype).itemsize
    if len(payload) != expected:
        raise ValueError(f"{path} has {len(payload)} bytes, expected {expected}")
    return np.frombuffer(payload, dtype=dtype).reshape(shape).astype(np.uint8 if byte else np.float32)


def _comparison(candidate: np.ndarray, reference: np.ndarray, *, label: str) -> dict[str, object]:
    a = np.asarray(candidate)
    b = np.asarray(reference)
    if a.shape != b.shape:
        raise ValueError(f"{label} shapes differ: {a.shape} != {b.shape}")
    error = np.abs(a.astype(np.float64) - b.astype(np.float64))
    result: dict[str, object] = {
        "exact": bool(np.array_equal(a, b)),
        "changed_voxels": int(np.count_nonzero(a != b)),
        "max_absolute_error": float(error.max(initial=0)),
        "mean_absolute_error": float(error.mean()),
    }
    if label in {"inner_distance", "outer_distance", "level_set", "brain_mask_level_set"} or label.startswith("tgdm_pass_"):
        nonpositive_a, nonpositive_b = a <= 0, b <= 0
        denominator = int(nonpositive_a.sum() + nonpositive_b.sum())
        result["sign_changed_voxels"] = int(np.count_nonzero(nonpositive_a != nonpositive_b))
        result["nonpositive_dice"] = float(2 * np.count_nonzero(nonpositive_a & nonpositive_b) / denominator) if denominator else 1.0
    if label == "brain_mask":
        foreground_a, foreground_b = a > 0, b > 0
        denominator = int(foreground_a.sum() + foreground_b.sum())
        result["foreground_dice"] = float(2 * np.count_nonzero(foreground_a & foreground_b) / denominator) if denominator else 1.0
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("t1", "white-matter", "gray-matter", "white-matter-mask"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    parser.add_argument("--java-output-dir", required=True, type=Path)
    parser.add_argument("--python-output-dir", type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--membership-level", type=float, default=0.3)
    parser.add_argument("--iterations", type=int, default=5)
    args = parser.parse_args()

    image_paths = {
        "t1": args.t1,
        "white_matter": args.white_matter,
        "gray_matter": args.gray_matter,
        "white_matter_mask": args.white_matter_mask,
    }
    images = {name: nib.load(str(path)) for name, path in image_paths.items()}
    reference = images["t1"]
    shape = tuple(int(size) for size in reference.shape)
    if len(shape) != 3 or nib.aff2axcodes(reference.affine) != ("R", "A", "S"):
        raise ValueError("Dura comparison inputs must be 3D RAS NIfTIs")
    for name, image in images.items():
        if image.shape != shape or not np.array_equal(image.affine, reference.affine):
            raise ValueError(f"{name} must share the exact T1 grid")

    result = run_dura_removal(
        images["t1"], images["white_matter"], images["gray_matter"], images["white_matter_mask"],
        membership_level=args.membership_level,
        max_iterations=args.iterations,
        record_passes=True,
    )
    initial = dura_initial_level_sets(
        np.asarray(result.wm_topology.dataobj, dtype=np.float32),
        np.asarray(result.gm_topology.dataobj, dtype=np.float32),
        membership_level=args.membership_level,
    )
    python_outputs = {
        "wm_tc": np.asarray(result.wm_topology.dataobj, dtype=np.float32),
        "union": np.add(
            np.asarray(result.wm_topology.dataobj, dtype=np.float32),
            np.asarray(images["gray_matter"].dataobj, dtype=np.float32),
            dtype=np.float32,
        ),
        "gm_tc": np.asarray(result.gm_topology.dataobj, dtype=np.float32),
        "inner_distance": initial.inner_signed_distance,
        "outer_distance": initial.outer_signed_distance,
        "level_set": np.asarray(result.brain_mask_level_set.dataobj, dtype=np.float32),
        "stripped_t1": np.asarray(result.stripped_t1.dataobj, dtype=np.float32),
        "corrected_wm": np.asarray(result.corrected_white_matter.dataobj, dtype=np.float32),
        "corrected_gm": np.asarray(result.corrected_gray_matter.dataobj, dtype=np.float32),
        "brain_mask_level_set": np.asarray(result.brain_mask_level_set.dataobj, dtype=np.float32),
        "brain_mask": np.asarray(result.brain_mask.dataobj, dtype=np.uint8),
    }
    python_outputs.update({
        f"tgdm_pass_{index:02d}": np.asarray(values, dtype=np.float32)
        for index, values in enumerate(result.evolution.pass_level_sets, start=1)
    })

    if args.python_output_dir is not None:
        args.python_output_dir.mkdir(parents=True, exist_ok=True)
        for name, values in python_outputs.items():
            dtype = np.uint8 if name == "brain_mask" else np.float32
            nib.save(
                nib.Nifti1Image(np.asarray(values, dtype=dtype), reference.affine),
                str(args.python_output_dir / f"{name}.nii.gz"),
            )

    comparisons = {}
    java_hashes = {}
    for name, candidate in python_outputs.items():
        path = args.java_output_dir / f"{name}.raw.gz"
        java = _read_java_volume(path, shape, byte=name == "brain_mask")
        comparisons[name] = _comparison(candidate, java, label=name)
        java_hashes[name] = _sha256(path)
    java_levelset = _read_java_volume(args.java_output_dir / "level_set.raw.gz", shape)
    java_image = nib.Nifti1Image(java_levelset, reference.affine)
    python_mesh = extract_levelset_surface(result.brain_mask_level_set, level=0.0)
    java_mesh = extract_levelset_surface(java_image, level=0.0)
    surface_distance = mesh_surface_distance_metrics(
        java_mesh.vertices, java_mesh.faces, python_mesh.vertices, python_mesh.faces
    )
    report = {
        "evidence": "selected original Java Dura wrapper and Python Dura API on identical held-fixed RAS inputs",
        "shape": list(shape),
        "coordinate_system": "RAS",
        "parameters": {
            "membership_level": args.membership_level,
            "outer_iterations": args.iterations,
            "mask_translation_mm": 2.0,
            "level_set_translation_mm": 0.0,
            "alpha": 0.5,
            "alignment_force": 1.0,
            "curvature_force": 0.5,
            "pressure_force": 0.5,
        },
        "input_sha256": {name: _sha256(path) for name, path in image_paths.items()},
        "java_output_sha256": java_hashes,
        "outputs": comparisons,
        "final_dura_mask_surface_comparison": {
            "method": "same Python source-derived isosurface extractor; symmetric vertex-to-triangle surface distance",
            "java_mesh_topology": {"vertices": int(len(java_mesh.vertices)), "faces": int(len(java_mesh.faces))},
            "python_mesh_topology": {"vertices": int(len(python_mesh.vertices)), "faces": int(len(python_mesh.faces))},
            "distances_mm": surface_distance,
        },
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(args.report), "outputs": comparisons}, indent=2))


if __name__ == "__main__":
    main()
