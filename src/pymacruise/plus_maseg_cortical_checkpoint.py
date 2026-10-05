"""Restartable MaSeg-to-cortical-mask preparation for plus workflows."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.cortical_separation import cortical_separation
from pymacruise.image_io import normalize_image_to_ras


@dataclass(frozen=True)
class MaSegCorticalCheckpointResult:
    manifest_path: Path
    outputs: dict[str, Path]
    resumed: bool


_OUTPUT_NAMES = ("Atlas_GM", "Atlas_WM", "Atlas_Mask", "Atlas_BrainTissueMask")


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


def run_plus_maseg_cortical_checkpoint(
    maseg_path: str | Path,
    output_dir: str | Path,
    *,
    resume: bool = False,
) -> MaSegCorticalCheckpointResult:
    """Create and persist the four source-selected MaSeg tissue masks."""
    if not isinstance(resume, bool):
        raise TypeError("resume must be boolean")
    source = Path(maseg_path).expanduser().resolve()
    root = Path(output_dir).expanduser().resolve()
    if not source.is_file() or source.is_symlink() or not source.name.lower().endswith(".nii.gz"):
        raise FileNotFoundError(f"MaSeg must be a regular compressed NIfTI input: {source}")
    source_image = nib.load(str(source))
    if not isinstance(source_image, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError("MaSeg must be a NIfTI image")
    normalized_source = normalize_image_to_ras(source_image, source_name="MaSeg")
    input_record = {"path": str(source), "sha256": _sha256(source)}
    signature_data = {
        "schema_version": 1,
        "stage": "plus-maseg-cortical-separation",
        "inputs": {"maseg": input_record},
        "outputs": list(_OUTPUT_NAMES),
        "spatial_convention": "RAS",
    }
    signature = hashlib.sha256(
        json.dumps(signature_data, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "maseg-cortical-separation.json"
    if manifest_path.is_symlink():
        raise FileExistsError(f"Refusing to write through symlinked manifest: {manifest_path}")
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"MaSeg cortical checkpoint exists; use resume=True: {manifest_path}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("run_signature") != signature:
            raise ValueError("Cannot resume MaSeg cortical checkpoint: inputs differ")
        if manifest.get("status") == "complete":
            outputs: dict[str, Path] = {}
            for name in _OUTPUT_NAMES:
                record = manifest.get("outputs", {}).get(name, {})
                path = Path(record.get("path", ""))
                if not path.is_file() or path.is_symlink() or _sha256(path) != record.get("sha256"):
                    break
                image = nib.load(str(path))
                if (
                    not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image))
                    or len(image.shape) != 3
                    or image.get_data_dtype() != np.dtype(np.uint8)
                    or nib.aff2axcodes(image.affine) != ("R", "A", "S")
                    or not np.isfinite(np.asarray(image.dataobj)).all()
                    or image.shape != normalized_source.shape
                    or not np.allclose(image.affine, normalized_source.affine, rtol=0, atol=1e-5)
                ):
                    break
                outputs[name] = path
            if len(outputs) == len(_OUTPUT_NAMES):
                return MaSegCorticalCheckpointResult(manifest_path, outputs, True)
    elif resume:
        raise FileNotFoundError(f"Cannot resume; no MaSeg cortical manifest exists: {manifest_path}")

    manifest = {
        **signature_data,
        "run_signature": signature,
        "manifest_path": str(manifest_path),
        "status": "running",
        "outputs": {},
    }
    _write_json(manifest_path, manifest)
    try:
        images = cortical_separation(source_image, root)
        outputs = {name: root / f"{name}.nii.gz" for name in _OUTPUT_NAMES}
        manifest["source_order"] = ["integer-label validation", "RAS reorientation", "selected MaSeg label-mask generation"]
        manifest["outputs"] = {
            name: {"path": str(path.resolve()), "sha256": _sha256(path)}
            for name, path in outputs.items()
        }
        for name, path in outputs.items():
            image = images[name]
            if (
                image.shape != normalized_source.shape
                or image.get_data_dtype() != np.dtype(np.uint8)
                or nib.aff2axcodes(image.affine) != ("R", "A", "S")
                or not np.allclose(image.affine, normalized_source.affine, rtol=0, atol=1e-5)
            ):
                raise RuntimeError(f"MaSeg cortical output {name} violates the RAS image contract")
        manifest["status"] = "complete"
        _write_json(manifest_path, manifest)
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        _write_json(manifest_path, manifest)
        raise
    return MaSegCorticalCheckpointResult(manifest_path, outputs, False)
