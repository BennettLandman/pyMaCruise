"""Transfer selected LesionTOADS Dura channels onto native T1 geometry."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.image_io import normalize_image_to_ras, save_nifti_ras
from pymacruise.spatial_transfer import (
    assign_reference_geometry_without_resampling,
    transfer_to_reference_space,
)


@dataclass(frozen=True)
class NativeDuraTransferResult:
    manifest_path: Path
    outputs: dict[str, Path]
    resumed: bool


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def _input_record(path: Path) -> dict[str, str]:
    if not path.is_file() or path.is_symlink():
        raise FileNotFoundError(f"Missing regular transfer input: {path}")
    return {"path": str(path.resolve()), "sha256": _sha256(path)}


def run_plus_lesion_native_dura_transfer(
    dura_manifest_path: str | Path,
    t1_reference_path: str | Path,
    output_dir: str | Path,
    *,
    transfer_mode: str = "physical-affine",
    resume: bool = False,
) -> NativeDuraTransferResult:
    """Transfer GM, WM-fill, and WM-mask Dura images into native T1 space.

    ``physical-affine`` maps RAS physical affines with nearest-neighbor
    sampling and is recommended. ``jist-unchanged`` preserves RAS source
    extents and values, then assigns the reference geometry to reproduce the
    selected JIST behavior for legacy compatibility comparisons.
    """
    if not isinstance(resume, bool):
        raise TypeError("resume must be boolean")
    if transfer_mode not in ("physical-affine", "jist-unchanged"):
        raise ValueError("transfer_mode must be 'physical-affine' or 'jist-unchanged'")
    dura_manifest_file = Path(dura_manifest_path).expanduser().resolve()
    t1_file = Path(t1_reference_path).expanduser().resolve()
    root = Path(output_dir).expanduser().resolve()
    if not dura_manifest_file.is_file() or dura_manifest_file.is_symlink():
        raise FileNotFoundError(f"Missing regular Dura-output manifest: {dura_manifest_file}")
    if not t1_file.is_file() or t1_file.is_symlink() or not t1_file.name.lower().endswith(".nii.gz"):
        raise FileNotFoundError(f"T1 reference must be a regular .nii.gz file: {t1_file}")
    dura_manifest = json.loads(dura_manifest_file.read_text(encoding="utf-8"))
    if (
        dura_manifest.get("status") != "complete"
        or dura_manifest.get("stage") != "plus-lesion-dura-output-checkpoint"
        or dura_manifest.get("spatial_convention") != "RAS"
    ):
        raise ValueError("input manifest is not a completed RAS LesionTOADS Dura checkpoint")
    raw_t1 = nib.load(str(t1_file))
    if not isinstance(raw_t1, (nib.Nifti1Image, nib.Nifti2Image)) or len(raw_t1.shape) != 3:
        raise ValueError("T1 reference must be a 3-D NIfTI image")
    reference = normalize_image_to_ras(raw_t1, source_name="T1 reference for Dura transfer")
    if not np.issubdtype(reference.get_data_dtype(), np.number):
        raise ValueError("T1 reference must contain finite numeric voxels")
    for x_start in range(0, reference.shape[0], 8):
        reference_values = np.asarray(reference.dataobj[x_start : x_start + 8, :, :])
        if not np.isfinite(reference_values).all():
            raise ValueError("T1 reference must contain finite numeric voxels")

    # Preserve CSF as an auxiliary channel for downstream APIs that accept the
    # source's four-channel generator product, even though selected plus only
    # routes GM, WM-fill, and WM-mask.
    channel_names = ("csf", "gm", "wm_fill", "wm_mask")
    source_paths: dict[str, Path] = {}
    for name in channel_names:
        record = dura_manifest.get("outputs", {}).get(name)
        if not isinstance(record, dict):
            raise ValueError(f"Dura checkpoint is missing the {name} output")
        source_path = Path(record.get("path", "")).resolve()
        if (
            not source_path.is_file()
            or source_path.is_symlink()
            or _sha256(source_path) != record.get("sha256")
        ):
            raise ValueError(f"Dura checkpoint {name} output is missing or fails its hash")
        source_paths[name] = source_path
    input_records = {
        "dura_manifest": _input_record(dura_manifest_file),
        "t1_reference": _input_record(t1_file),
        **{f"dura_{name}": _input_record(path) for name, path in source_paths.items()},
    }
    signature_data = {
        "schema_version": 1,
        "stage": "plus-lesion-native-dura-transfer",
        "source_order": ["transfer CSF auxiliary", "transfer GM", "transfer WM-fill", "transfer WM-mask"],
        "interpolation": "nearest-neighbor" if transfer_mode == "physical-affine" else "none",
        "coordinate_mapping": (
            "physical RAS affine" if transfer_mode == "physical-affine"
            else "preserve RAS source extents and values; assign reference geometry"
        ),
        "transfer_mode": transfer_mode,
        "spatial_convention": "RAS",
        "inputs": input_records,
        "source_equivalence_claim": False,
    }
    signature = hashlib.sha256(
        json.dumps(signature_data, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "native-dura-transfer.json"
    if manifest_path.is_symlink():
        raise FileExistsError(f"Refusing to write through symlinked manifest: {manifest_path}")
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"Native Dura transfer exists; use resume=True: {manifest_path}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        # If an upstream stage was repaired during a top-level resume, transfer
        # again from its newly verified Dura products instead of rejecting the
        # updated dependency signature.
        if manifest.get("run_signature") == signature and manifest.get("status") == "complete":
            cached: dict[str, Path] = {}
            for name in channel_names:
                record = manifest.get("outputs", {}).get(name, {})
                path = Path(record.get("path", ""))
                if (
                    not path.is_file()
                    or path.is_symlink()
                    or _sha256(path) != record.get("sha256")
                ):
                    break
                image = nib.load(str(path))
                source_image = nib.load(str(source_paths[name]))
                expected_shape = reference.shape if transfer_mode == "physical-affine" else source_image.shape
                if (
                    not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image))
                    or image.shape != expected_shape
                    or image.get_data_dtype() != np.dtype(np.float32)
                    or nib.aff2axcodes(image.affine) != ("R", "A", "S")
                    or not np.allclose(image.affine, reference.affine, rtol=0, atol=1e-5)
                    or not np.isfinite(np.asarray(image.dataobj)).all()
                ):
                    break
                cached[name] = path
            if len(cached) == len(channel_names):
                return NativeDuraTransferResult(manifest_path, cached, True)
    elif resume:
        raise FileNotFoundError(f"Cannot resume; native Dura manifest does not exist: {manifest_path}")

    manifest = {
        **signature_data,
        "run_signature": signature,
        "manifest_path": str(manifest_path),
        "status": "running",
        "outputs": {},
    }
    _write_json(manifest_path, manifest)
    outputs: dict[str, Path] = {}
    try:
        for name in channel_names:
            source = nib.load(str(source_paths[name]))
            if not isinstance(source, (nib.Nifti1Image, nib.Nifti2Image)) or len(source.shape) != 3:
                raise ValueError(f"Dura {name} input must be a 3-D NIfTI image")
            transferred = (
                transfer_to_reference_space(source, reference, cval=0.0)
                if transfer_mode == "physical-affine"
                else assign_reference_geometry_without_resampling(source, reference)
            )
            final_path = root / f"lesion_toads_{name}_native.nii.gz"
            temporary_path = root / f".{final_path.stem}.partial.nii.gz"
            save_nifti_ras(transferred, temporary_path, dtype=np.float32)
            temporary_path.replace(final_path)
            outputs[name] = final_path
        manifest["outputs"] = {
            name: {"path": str(path.resolve()), "sha256": _sha256(path)}
            for name, path in outputs.items()
        }
        manifest["status"] = "complete"
        _write_json(manifest_path, manifest)
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        _write_json(manifest_path, manifest)
        raise
    return NativeDuraTransferResult(manifest_path, outputs, False)
