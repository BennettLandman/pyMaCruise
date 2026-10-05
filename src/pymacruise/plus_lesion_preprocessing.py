"""Prepare selected MaCRUISEplus T1/FLAIR inputs for LesionTOADS.

The selected MATLAB route rigid-registers FLAIR to T1, N4-corrects both
modalities, then strips them with the brain mask. The Python workflow takes
that mask from the user, as required by this project, instead of invoking BET
or SPECTRE. This module implements only that preparation handoff; lesion atlas
registration, LesionTOADS segmentation, and later plus stages remain separate.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras, nifti_ras_on_reference, save_nifti_ras
from pymacruise.n4 import N4Config, n4_bias_correct
from pymacruise.plus_workflow import run_plus_flair_registration_checkpoint
from pymacruise.registration import RegistrationConfig


@dataclass(frozen=True)
class PlusLesionPreprocessingResult:
    """Persisted RAS inputs ready for the selected LesionTOADS handoff."""

    manifest_path: Path
    t1_brain_path: Path
    flair_brain_path: Path
    modalities_path: Path
    resumed: bool


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_hash(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _write_manifest(path: Path, manifest: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def _on_grid(image: nib.spatialimages.SpatialImage, reference: nib.spatialimages.SpatialImage) -> bool:
    return len(image.shape) in (3, 4) and image.shape[:3] == reference.shape and np.allclose(
        image.affine, reference.affine, rtol=0, atol=1e-4
    )


def _valid_result(
    manifest: dict[str, Any], signature: str, root: Path, reference: nib.spatialimages.SpatialImage
) -> PlusLesionPreprocessingResult | None:
    if manifest.get("status") != "complete" or manifest.get("run_signature") != signature:
        return None
    registration_record = manifest.get("registration_manifest", {})
    registration_path = Path(registration_record.get("path", ""))
    if (
        not registration_path.is_relative_to(root)
        or not registration_path.is_file()
        or registration_path.is_symlink()
        or _sha256(registration_path) != registration_record.get("sha256")
    ):
        return None
    outputs = manifest.get("outputs")
    if not isinstance(outputs, dict):
        return None
    paths: dict[str, Path] = {}
    for name in ("t1_brain", "flair_brain", "modalities_t1_flair"):
        record = outputs.get(name)
        if not isinstance(record, dict):
            return None
        path = Path(record.get("path", ""))
        if not path.is_relative_to(root) or not path.is_file() or path.is_symlink():
            return None
        if _sha256(path) != record.get("sha256"):
            return None
        image = nib.load(str(path))
        if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
            return None
        if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            return None
        if not _on_grid(image, reference):
            return None
        expected_ndim = 4 if name == "modalities_t1_flair" else 3
        if len(image.shape) != expected_ndim or (expected_ndim == 4 and image.shape[-1] != 2):
            return None
        if not np.isfinite(np.asarray(image.dataobj)).all():
            return None
        paths[name] = path
    return PlusLesionPreprocessingResult(
        manifest_path=Path(manifest["manifest_path"]),
        t1_brain_path=paths["t1_brain"],
        flair_brain_path=paths["flair_brain"],
        modalities_path=paths["modalities_t1_flair"],
        resumed=True,
    )


def run_plus_lesion_preprocessing_checkpoint(
    t1_path: str | Path,
    flair_path: str | Path,
    brain_mask_path: str | Path,
    output_dir: str | Path,
    *,
    registration_config: RegistrationConfig = RegistrationConfig(),
    n4_config: N4Config = N4Config(),
    resume: bool = False,
) -> PlusLesionPreprocessingResult:
    """Register FLAIR to T1, correct both images, and apply the supplied mask.

    Input images must be compressed 3-D NIfTI. The brain mask must match the
    T1 grid; FLAIR may arrive on a different grid because rigid registration
    brings it onto T1. N4 follows registration in the source route. Negative
    corrected values are set to zero when the supplied brain mask is applied,
    matching the selected ``removeNegative`` cleanup. The 4-D handoff stores
    channel-last ``(T1, FLAIR)`` modalities for downstream LesionTOADS code.
    """
    if not isinstance(registration_config, RegistrationConfig):
        raise TypeError("registration_config must be RegistrationConfig")
    if not isinstance(n4_config, N4Config):
        raise TypeError("n4_config must be N4Config")
    if not isinstance(resume, bool):
        raise TypeError("resume must be boolean")

    t1_file = Path(t1_path).expanduser().resolve()
    flair_file = Path(flair_path).expanduser().resolve()
    mask_file = Path(brain_mask_path).expanduser().resolve()
    root = Path(output_dir).expanduser().resolve()
    for role, path in (("T1", t1_file), ("FLAIR", flair_file), ("brain mask", mask_file)):
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"Missing regular {role} input: {path}")
        if not path.name.lower().endswith(".nii.gz"):
            raise ValueError(f"{role} input must be .nii.gz")

    t1 = load_nifti_ras(t1_file)
    mask = load_nifti_ras(mask_file)
    if len(t1.shape) != 3 or len(mask.shape) != 3:
        raise ValueError("T1 and brain mask must be 3-D scalar images")
    if not _on_grid(mask, t1):
        raise ValueError("brain mask must match the T1 reference grid")
    for role, image in (("T1", t1), ("brain mask", mask)):
        for start in range(0, image.shape[0], 8):
            if not np.isfinite(np.asarray(image.dataobj[start : start + 8])).all():
                raise ValueError(f"{role} contains nonfinite values")
    brain = np.asarray(mask.dataobj) > 0
    if not brain.any():
        raise ValueError("brain mask must contain at least one nonzero voxel")

    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "lesion-preprocessing.json"
    signature = _json_hash(
        {
            "stage": "plus-lesion-input-preparation",
            "inputs": {
                "t1": {"path": str(t1_file), "sha256": _sha256(t1_file)},
                "flair": {"path": str(flair_file), "sha256": _sha256(flair_file)},
                "brain_mask": {"path": str(mask_file), "sha256": _sha256(mask_file)},
            },
            "registration": asdict(registration_config),
            "n4": asdict(n4_config),
            "channel_order": ["T1", "FLAIR"],
            "spatial_convention": "RAS",
        }
    )
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"Lesion preprocessing checkpoint exists; use resume=True: {manifest_path}")
        cached_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if cached_manifest.get("run_signature") != signature:
            raise ValueError("Cannot resume plus lesion preprocessing: inputs or settings differ")
        cached = _valid_result(cached_manifest, signature, root, t1)
        if cached is not None:
            return cached
    elif resume:
        raise FileNotFoundError(f"Cannot resume; no lesion preprocessing manifest exists: {manifest_path}")

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "stage": "plus-lesion-input-preparation",
        "manifest_path": str(manifest_path),
        "run_signature": signature,
        "status": "running",
        "spatial_convention": "RAS",
        "source_order": ["rigid FLAIR-to-T1 registration", "N4 T1", "N4 registered FLAIR", "brain-mask strip and negative cleanup"],
        "registration_backend": "ANTsPy",
        "inputs": {
            "t1": {"path": str(t1_file), "sha256": _sha256(t1_file)},
            "flair": {"path": str(flair_file), "sha256": _sha256(flair_file)},
            "brain_mask": {"path": str(mask_file), "sha256": _sha256(mask_file)},
        },
        "settings": {"registration": asdict(registration_config), "n4": asdict(n4_config)},
        "channel_order": ["T1", "FLAIR"],
        "equivalence_note": "Synthetic/input-preparation composition; ANTsPy versus NiftyReg/Aladin and the complete LesionTOADS workflow remain unvalidated.",
    }
    _write_manifest(manifest_path, manifest)

    registration = run_plus_flair_registration_checkpoint(
        t1_file,
        flair_file,
        root / "flair-registration",
        config=registration_config,
        resume=resume and (root / "flair-registration" / "flair-registration.json").exists(),
    )
    manifest["registration_manifest"] = {
        "path": str(registration.manifest_path),
        "sha256": _sha256(registration.manifest_path),
    }
    _write_manifest(manifest_path, manifest)

    corrected_t1 = n4_bias_correct(t1, config=n4_config)
    corrected_flair = n4_bias_correct(registration.warped_flair, config=n4_config)
    t1_values = np.asarray(corrected_t1.dataobj, dtype=np.float32)
    flair_values = np.asarray(corrected_flair.dataobj, dtype=np.float32)
    if t1_values.shape != t1.shape or flair_values.shape != t1.shape:
        raise RuntimeError("N4 output did not preserve the T1 reference grid")
    t1_values = np.where(brain, np.maximum(t1_values, np.float32(0.0)), np.float32(0.0))
    flair_values = np.where(brain, np.maximum(flair_values, np.float32(0.0)), np.float32(0.0))
    t1_output = save_nifti_ras(nifti_ras_on_reference(t1, t1_values, dtype=np.float32), root / "t1_brain.nii.gz")
    flair_output = save_nifti_ras(
        nifti_ras_on_reference(t1, flair_values, dtype=np.float32), root / "flair_brain.nii.gz"
    )
    modalities = np.stack((t1_values, flair_values), axis=-1)
    modalities_output = save_nifti_ras(
        nifti_ras_on_reference(t1, modalities, dtype=np.float32), root / "modalities_t1_flair.nii.gz"
    )
    manifest["status"] = "complete"
    manifest["outputs"] = {
        "t1_brain": {"path": str(t1_output), "sha256": _sha256(t1_output)},
        "flair_brain": {"path": str(flair_output), "sha256": _sha256(flair_output)},
        "modalities_t1_flair": {"path": str(modalities_output), "sha256": _sha256(modalities_output)},
    }
    manifest["registration_output"] = {
        "path": str(registration.warped_flair_path),
        "sha256": _sha256(registration.warped_flair_path),
    }
    _write_manifest(manifest_path, manifest)
    return PlusLesionPreprocessingResult(
        manifest_path, t1_output, flair_output, modalities_output, resumed=False
    )
