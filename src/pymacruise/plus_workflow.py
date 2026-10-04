"""Restartable upstream registration checkpoint for MaCRUISEplus.

The selected MATLAB runner rigidly registers FLAIR to T1 before entering the
plus pipeline. This module records that one stage around the existing ANTsPy
adapter; it does not reproduce NiftyReg/Aladin's optimizer or assemble the
remaining plus workflow.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.image_io import normalize_image_to_ras, save_nifti_ras
from pymacruise.registration import RegistrationConfig, register_flair_to_t1


@dataclass(frozen=True)
class PlusFlairRegistrationCheckpoint:
    """Persisted result of rigid FLAIR-to-T1 preparation."""

    warped_flair: nib.Nifti1Image
    warped_flair_path: Path
    forward_transforms: tuple[Path, ...]
    inverse_transforms: tuple[Path, ...]
    manifest_path: Path
    resumed: bool


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _signature(t1: Path, flair: Path, config: RegistrationConfig) -> str:
    payload = {
        "schema_version": 1,
        "stage": "plus-flair-to-t1-rigid-registration",
        "inputs": {
            "t1": {"path": str(t1), "sha256": _sha256(t1)},
            "flair": {"path": str(flair), "sha256": _sha256(flair)},
        },
        "registration": asdict(config),
        "spatial_convention": "RAS",
        "backend": "ANTsPy rigid registration; Aladin equivalence not established",
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def _write_manifest(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def _valid_cached_result(
    manifest: dict[str, Any],
    signature: str,
    output_dir: Path,
    fixed: nib.spatialimages.SpatialImage,
) -> PlusFlairRegistrationCheckpoint | None:
    if manifest.get("run_signature") != signature or manifest.get("status") != "complete":
        return None
    output_records = manifest.get("outputs", {})
    warped_path = Path(output_records.get("warped_flair", {}).get("path", ""))
    if not warped_path.is_relative_to(output_dir) or not warped_path.is_file() or warped_path.is_symlink():
        return None
    if _sha256(warped_path) != output_records["warped_flair"].get("sha256"):
        return None
    transforms: dict[str, tuple[Path, ...]] = {}
    for role in ("forward_transforms", "inverse_transforms"):
        records = output_records.get(role, [])
        paths = []
        for record in records:
            path = Path(record.get("path", ""))
            if (
                not path.is_relative_to(output_dir)
                or not path.is_file()
                or path.is_symlink()
                or _sha256(path) != record.get("sha256")
            ):
                return None
            paths.append(path)
        if not paths:
            return None
        transforms[role] = tuple(paths)
    warped = nib.load(str(warped_path))
    if (
        not isinstance(warped, (nib.Nifti1Image, nib.Nifti2Image))
        or warped.shape != fixed.shape
        or not np.allclose(warped.affine, fixed.affine, rtol=0, atol=1e-4)
        or nib.aff2axcodes(warped.affine) != ("R", "A", "S")
        or not np.isfinite(np.asarray(warped.dataobj)).all()
    ):
        return None
    return PlusFlairRegistrationCheckpoint(
        warped_flair=warped,
        warped_flair_path=warped_path,
        forward_transforms=transforms["forward_transforms"],
        inverse_transforms=transforms["inverse_transforms"],
        manifest_path=Path(manifest["manifest_path"]),
        resumed=True,
    )


def run_plus_flair_registration_checkpoint(
    t1_path: str | Path,
    flair_path: str | Path,
    output_dir: str | Path,
    *,
    config: RegistrationConfig = RegistrationConfig(),
    resume: bool = False,
) -> PlusFlairRegistrationCheckpoint:
    """Register FLAIR to the T1 RAS grid and persist a verifiable checkpoint.

    The T1 is the fixed/reference image and FLAIR is moving. Both inputs must
    be compressed NIfTI files. Non-RAS inputs are normalized by the existing
    registration adapter with its standard warning. A valid ``resume=True``
    verifies the manifest, warped image, transforms, and input signature before
    skipping registration. This is one stage of MaCRUISEplus only.
    """
    t1 = Path(t1_path).expanduser().resolve()
    flair = Path(flair_path).expanduser().resolve()
    root = Path(output_dir).expanduser().resolve()
    if not isinstance(config, RegistrationConfig) or config.transform_type != "Rigid":
        raise ValueError("plus FLAIR registration requires a rigid RegistrationConfig")
    for role, path in (("T1", t1), ("FLAIR", flair)):
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"Missing regular {role} input: {path}")
        if not path.name.lower().endswith(".nii.gz"):
            raise ValueError(f"{role} input must be .nii.gz")

    fixed = normalize_image_to_ras(nib.load(str(t1)), source_name="plus T1 registration reference")
    moving = normalize_image_to_ras(nib.load(str(flair)), source_name="plus FLAIR moving image")
    if len(fixed.shape) != 3 or len(moving.shape) != 3:
        raise ValueError("plus T1 and FLAIR inputs must be 3-D scalar images")
    if not np.isfinite(np.asarray(fixed.dataobj)).all() or not np.isfinite(np.asarray(moving.dataobj)).all():
        raise ValueError("plus T1 and FLAIR inputs must contain finite intensities")

    signature = _signature(t1, flair, config)
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "flair-registration.json"
    warped_path = root / "flair_in_t1.nii.gz"
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"Registration checkpoint exists; use resume=True: {manifest_path}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("run_signature") != signature:
            raise ValueError("Cannot resume plus registration: inputs or settings differ")
        cached = _valid_cached_result(manifest, signature, root, fixed)
        if cached is not None:
            return cached
    elif resume:
        raise FileNotFoundError(f"Cannot resume; no plus registration manifest exists: {manifest_path}")

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "stage": "plus-flair-to-t1-rigid-registration",
        "manifest_path": str(manifest_path),
        "run_signature": signature,
        "status": "running",
        "spatial_convention": "RAS",
        "registration_backend": "ANTsPy",
        "equivalence_note": "ANTsPy is the selected Python substitute; NiftyReg/Aladin numerical equivalence is unvalidated.",
        "inputs": {
            "t1": {"path": str(t1), "sha256": _sha256(t1)},
            "flair": {"path": str(flair), "sha256": _sha256(flair)},
        },
        "registration": asdict(config),
    }
    _write_manifest(manifest_path, manifest)
    result = register_flair_to_t1(
        fixed,
        moving,
        output_prefix=root / "ants" / "flair_to_t1_",
        config=config,
    )
    save_nifti_ras(result.warped_moving, warped_path)
    forward = tuple(Path(path).resolve() for path in result.forward_transforms)
    inverse = tuple(Path(path).resolve() for path in result.inverse_transforms)
    if not forward or not inverse or any(not path.is_relative_to(root) for path in (*forward, *inverse)):
        raise RuntimeError("ANTsPy did not persist all transforms inside the checkpoint directory")
    manifest["status"] = "complete"
    manifest["outputs"] = {
        "warped_flair": {"path": str(warped_path), "sha256": _sha256(warped_path)},
        "forward_transforms": [
            {"path": str(path), "sha256": _sha256(path)} for path in forward
        ],
        "inverse_transforms": [
            {"path": str(path), "sha256": _sha256(path)} for path in inverse
        ],
    }
    _write_manifest(manifest_path, manifest)
    return PlusFlairRegistrationCheckpoint(
        result.warped_moving,
        warped_path,
        forward,
        inverse,
        manifest_path,
        resumed=False,
    )
