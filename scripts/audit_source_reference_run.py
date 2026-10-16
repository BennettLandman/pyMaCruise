#!/usr/bin/env python3
"""Hash and validate outputs returned from an original MaCRUISE source run."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import nibabel as nib
import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_nifti(path: Path, *, slab_voxels: int) -> dict[str, Any]:
    image = nib.load(str(path))
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError("not a NIfTI-1 or NIfTI-2 image")
    if len(image.shape) not in (3, 4) or any(int(size) < 1 for size in image.shape):
        raise ValueError(f"expected nonempty 3D or 4D image, got shape={image.shape}")
    dtype = image.get_data_dtype()
    if not np.issubdtype(dtype, np.number) or np.issubdtype(dtype, np.complexfloating):
        raise ValueError(f"nonnumeric NIfTI datatype: {image.get_data_dtype()}")

    finite = True
    minimum: float | None = None
    maximum: float | None = None
    spatial_shape = tuple(int(size) for size in image.shape[:3])
    z_chunk = max(1, min(spatial_shape[2], slab_voxels))
    y_chunk = max(1, min(spatial_shape[1], slab_voxels // z_chunk))
    x_chunk = max(1, min(spatial_shape[0], slab_voxels // (y_chunk * z_chunk)))
    frame_count = int(image.shape[3]) if len(image.shape) == 4 else 1
    for frame in range(frame_count):
        for x_start in range(0, spatial_shape[0], x_chunk):
            for y_start in range(0, spatial_shape[1], y_chunk):
                for z_start in range(0, spatial_shape[2], z_chunk):
                    slices = (
                        slice(x_start, x_start + x_chunk),
                        slice(y_start, y_start + y_chunk),
                        slice(z_start, z_start + z_chunk),
                    )
                    if len(image.shape) == 4:
                        slices += (frame,)
                    slab = np.asarray(image.dataobj[slices])
                    finite_values = slab[np.isfinite(slab)]
                    if finite_values.size != slab.size:
                        finite = False
                    if finite_values.size:
                        slab_min = float(np.min(finite_values))
                        slab_max = float(np.max(finite_values))
                        minimum = slab_min if minimum is None else min(minimum, slab_min)
                        maximum = slab_max if maximum is None else max(maximum, slab_max)

    affine = np.asarray(image.affine, dtype=np.float64)
    if not np.isfinite(affine).all():
        raise ValueError("nonfinite NIfTI affine")
    qform, qcode = image.get_qform(coded=True)
    sform, scode = image.get_sform(coded=True)
    axis_codes = tuple(str(code) for code in nib.aff2axcodes(affine))
    return {
        "format": "NIfTI-1" if isinstance(image, nib.Nifti1Image) else "NIfTI-2",
        "shape": list(image.shape),
        "dtype": str(image.get_data_dtype()),
        "zooms": [float(value) for value in image.header.get_zooms()],
        "axis_codes": list(axis_codes),
        "ras_oriented": axis_codes == ("R", "A", "S"),
        "affine": affine.tolist(),
        "qform_code": int(qcode),
        "qform": np.asarray(qform, dtype=np.float64).tolist() if qform is not None else None,
        "sform_code": int(scode),
        "sform": np.asarray(sform, dtype=np.float64).tolist() if sform is not None else None,
        "finite_voxels": finite,
        "minimum": minimum,
        "maximum": maximum,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path, help="returned source-run directory to audit")
    parser.add_argument("--output", type=Path, help="JSON output path (default RUN_DIR/source-run-manifest.json)")
    parser.add_argument("--require-run-info", action="store_true", help="fail unless RUN_INFO.md or RUN_INFO.txt exists")
    parser.add_argument("--require-ras", action="store_true", help="fail on any non-RAS NIfTI output")
    parser.add_argument("--slab-voxels", type=int, default=1_000_000, help="approximate maximum voxels loaded per finite-value check slab")
    args = parser.parse_args()
    root = args.run_dir.expanduser().resolve()
    if not root.is_dir() or root.is_symlink():
        parser.error(f"run directory must be a regular directory: {root}")
    if args.slab_voxels < 1:
        parser.error("--slab-voxels must be positive")
    output = (args.output or (root / "source-run-manifest.json")).expanduser().resolve()
    try:
        output.relative_to(root)
    except ValueError:
        pass

    paths = sorted(path for path in root.rglob("*") if path != output)
    errors: list[str] = []
    warnings: list[str] = []
    files: dict[str, Any] = {}
    info_names = {"run_info.md", "run_info.txt"}
    has_run_info = False
    for path in paths:
        relative = path.relative_to(root).as_posix()
        if path.is_symlink():
            errors.append(f"symlink refused: {relative}")
            continue
        if not path.is_file():
            continue
        record: dict[str, Any] = {"sha256": sha256(path), "size_bytes": path.stat().st_size}
        if path.name.lower() in info_names:
            has_run_info = True
        if path.name.lower().endswith(".nii"):
            errors.append(f"uncompressed NIfTI is not allowed; convert to .nii.gz: {relative}")
        elif path.name.lower().endswith(".nii.gz"):
            try:
                geometry = audit_nifti(path, slab_voxels=args.slab_voxels)
                record["nifti"] = geometry
                if not geometry["finite_voxels"]:
                    errors.append(f"nonfinite NIfTI voxel data: {relative}")
                if not geometry["ras_oriented"]:
                    message = f"non-RAS NIfTI orientation {geometry['axis_codes']}: {relative}"
                    if args.require_ras:
                        errors.append(message)
                    else:
                        warnings.append(message)
            except Exception as exc:
                errors.append(f"could not audit NIfTI {relative}: {type(exc).__name__}: {exc}")
        files[relative] = record

    if args.require_run_info and not has_run_info:
        errors.append("missing RUN_INFO.md or RUN_INFO.txt")
    manifest = {
        "schema_version": 1,
        "run_directory": str(root),
        "require_ras": bool(args.require_ras),
        "require_run_info": bool(args.require_run_info),
        "has_run_info": has_run_info,
        "file_count": len(files),
        "nifti_count": sum("nifti" in record for record in files.values()),
        "files": files,
        "warnings": warnings,
        "errors": errors,
        "passed": not errors,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"manifest": str(output), "files": len(files), "niftis": manifest["nifti_count"], "warnings": warnings, "errors": errors, "passed": manifest["passed"]}, indent=2))
    return 0 if manifest["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
