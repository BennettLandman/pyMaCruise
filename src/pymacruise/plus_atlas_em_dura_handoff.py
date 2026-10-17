"""Restartable pre-LesionTOADS AtlasEM/Dura input handoff for MaCRUISEplus."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.atlas_em_descriptor import parse_atlas_em_descriptor
from pymacruise.atlas_em_execution import run_atlas_em_descriptor_workflow
from pymacruise.atlas_em_registration import (
    ATLAS_EM_PLUGIN_DEFAULT_INITIAL_REGISTRATION_ITERATIONS,
    ATLAS_EM_PLUGIN_DEFAULT_LEVELS,
    ATLAS_EM_PLUGIN_DEFAULT_POST_RECROP_REGISTRATION_ITERATIONS,
)
from pymacruise.image_io import load_nifti_ras, save_nifti_ras


@dataclass(frozen=True)
class PlusAtlasEMDuraHandoffResult:
    """AtlasEM/Dura outputs and Dura-conditioned T1/FLAIR handoff."""

    manifest_path: Path
    modalities_path: Path
    dura_stripped_t1_path: Path
    dura_masked_flair_path: Path
    output_paths: dict[str, Path]
    resumed: bool


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _record(path: Path) -> dict[str, str]:
    if not path.is_file() or path.is_symlink():
        raise FileNotFoundError(f"Missing regular AtlasEM/Dura handoff file: {path}")
    return {"path": str(path.resolve()), "sha256": _sha256(path)}


def _write_manifest(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def _cached_result(
    manifest_path: Path,
    manifest: dict[str, Any],
    signature: str,
    root: Path,
    reference: nib.spatialimages.SpatialImage,
) -> PlusAtlasEMDuraHandoffResult | None:
    if manifest.get("status") != "complete" or manifest.get("run_signature") != signature:
        return None
    outputs = manifest.get("outputs")
    if not isinstance(outputs, dict):
        return None
    paths: dict[str, Path] = {}
    for name, record in outputs.items():
        if not isinstance(record, dict):
            return None
        path = Path(record.get("path", ""))
        if not path.is_relative_to(root) or not path.is_file() or path.is_symlink():
            return None
        if _sha256(path) != record.get("sha256"):
            return None
        if path.suffixes[-2:] != [".nii", ".gz"]:
            return None
        image = nib.load(str(path))
        if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
            return None
        if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            return None
        if image.shape[:3] != reference.shape[:3] or not np.allclose(
            image.affine, reference.affine, rtol=0, atol=1e-4
        ):
            return None
        if not np.isfinite(np.asarray(image.dataobj)).all():
            return None
        paths[name] = path
    required = {
        "lesion_modalities", "dura_stripped_t1", "dura_masked_flair",
    }
    if not required.issubset(paths):
        return None
    modalities = nib.load(str(paths["lesion_modalities"]))
    if len(modalities.shape) != 4 or modalities.shape[-1] != 2:
        return None
    return PlusAtlasEMDuraHandoffResult(
        manifest_path=manifest_path,
        modalities_path=paths["lesion_modalities"],
        dura_stripped_t1_path=paths["dura_stripped_t1"],
        dura_masked_flair_path=paths["dura_masked_flair"],
        output_paths=paths,
        resumed=True,
    )


def run_plus_atlas_em_dura_handoff_checkpoint(
    t1_brain_path: str | Path,
    flair_brain_path: str | Path,
    brain_mask_path: str | Path,
    atlasem_descriptor_path: str | Path,
    output_dir: str | Path,
    *,
    levels: int = ATLAS_EM_PLUGIN_DEFAULT_LEVELS,
    initial_registration_iterations: int = ATLAS_EM_PLUGIN_DEFAULT_INITIAL_REGISTRATION_ITERATIONS,
    post_recrop_registration_iterations: int = ATLAS_EM_PLUGIN_DEFAULT_POST_RECROP_REGISTRATION_ITERATIONS,
    max_em_iterations: int = 99,
    max_membership_change: float = 0.01,
    smoothing: float = 10.0,
    max_dura_iterations: int = 5,
    resume: bool = False,
) -> PlusAtlasEMDuraHandoffResult:
    """Run AtlasEM and Dura before LesionTOADS and mask registered FLAIR.

    T1 and registered FLAIR are the brain-masked outputs of the preceding
    N4/input-preparation stage. The AtlasEM descriptor is distinct from the
    LesionTOADS atlas descriptor. The returned channel-last compressed NIfTI
    contains Dura-stripped T1 and Dura-masked registered FLAIR, on the same
    native RAS T1 grid, ready for the existing LesionTOADS crop stage.
    """
    if not isinstance(resume, bool):
        raise TypeError("resume must be boolean")
    root = Path(output_dir).expanduser().resolve()
    manifest_path = root / "atlasem-dura-handoff.json"
    t1_path = Path(t1_brain_path).expanduser().resolve()
    flair_path = Path(flair_brain_path).expanduser().resolve()
    mask_path = Path(brain_mask_path).expanduser().resolve()
    descriptor_path = Path(atlasem_descriptor_path).expanduser().resolve()
    for name, path in (
        ("T1", t1_path), ("registered FLAIR", flair_path),
        ("brain mask", mask_path), ("AtlasEM descriptor", descriptor_path),
    ):
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"Missing regular {name} input: {path}")
        if name != "AtlasEM descriptor" and not path.name.lower().endswith(".nii.gz"):
            raise ValueError(f"{name} input must be .nii.gz")
    descriptor = parse_atlas_em_descriptor(descriptor_path)
    if descriptor.missing_files:
        raise FileNotFoundError(
            "AtlasEM descriptor references missing assets: "
            + ", ".join(str(path) for path in descriptor.missing_files)
        )
    t1 = load_nifti_ras(t1_path)
    flair = load_nifti_ras(flair_path)
    mask = load_nifti_ras(mask_path)
    if len(t1.shape) != 3 or len(flair.shape) != 3 or len(mask.shape) != 3:
        raise ValueError("AtlasEM/Dura handoff T1, FLAIR, and brain mask must be 3-D")
    for name, image in (("registered FLAIR", flair), ("brain mask", mask)):
        if image.shape != t1.shape or not np.allclose(
            image.affine, t1.affine, rtol=0, atol=1e-4
        ):
            raise ValueError(f"{name} must match the T1 RAS grid")
    mask_values = np.asarray(mask.dataobj)
    if mask_values.dtype != np.bool_ and not np.isin(mask_values, (0, 1)).all():
        raise ValueError("brain mask must be binary")

    input_records = {
        "t1_brain": _record(t1_path),
        "flair_brain": _record(flair_path),
        "brain_mask": _record(mask_path),
        "atlasem_descriptor": _record(descriptor_path),
        "atlasem_assets": {str(path): _record(path)["sha256"] for path in descriptor.referenced_files},
    }
    settings: dict[str, Any] = {
        "modality": "T1_MPRAGE",
        "levels": int(levels),
        "initial_registration_iterations": int(initial_registration_iterations),
        "post_recrop_registration_iterations": int(post_recrop_registration_iterations),
        "max_em_iterations": int(max_em_iterations),
        "max_membership_change": float(max_membership_change),
        "smoothing": float(smoothing),
        "max_dura_iterations": int(max_dura_iterations),
        "dura_flair_policy": "Dura brain mask, clamp negative values to zero",
    }
    signature = hashlib.sha256(
        json.dumps(
            {"stage": "plus-prelesion-atlasem-dura-handoff", "inputs": input_records,
             "settings": settings, "spatial_convention": "RAS"},
            sort_keys=True, separators=(",", ":"), allow_nan=False,
        ).encode()
    ).hexdigest()
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"AtlasEM/Dura handoff exists; use resume=True: {manifest_path}")
        cached_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if cached_manifest.get("run_signature") != signature:
            raise ValueError("Cannot resume AtlasEM/Dura handoff: inputs or settings differ")
        cached = _cached_result(manifest_path, cached_manifest, signature, root, t1)
        if cached is not None:
            return cached
    elif resume:
        raise FileNotFoundError(f"Cannot resume; no AtlasEM/Dura handoff exists: {manifest_path}")

    root.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "stage": "plus-prelesion-atlasem-dura-handoff",
        "manifest_path": str(manifest_path),
        "run_signature": signature,
        "status": "running",
        "spatial_convention": "RAS",
        "inputs": input_records,
        "settings": settings,
        "source_order": [
            "AtlasEM on prepared T1", "selected AtlasEM-output Dura removal",
            "Dura-mask registered FLAIR", "form Dura-conditioned T1/FLAIR handoff",
        ],
        "equivalence_claim": "synthetic stage integration only; selected atlas assets and complete PLUS source equivalence remain unverified",
    }
    _write_manifest(manifest_path, manifest)
    atlasem = run_atlas_em_descriptor_workflow(
        t1, mask, descriptor_path, root / "atlasem-dura",
        modality=settings["modality"], levels=levels,
        initial_registration_iterations=initial_registration_iterations,
        post_recrop_registration_iterations=post_recrop_registration_iterations,
        max_em_iterations=max_em_iterations,
        max_membership_change=max_membership_change,
        smoothing=smoothing,
        max_dura_iterations=max_dura_iterations,
    )
    dura_t1_path = atlasem.output_paths["dura_stripped_t1"]
    dura_mask_path = atlasem.output_paths["dura_brain_mask"]
    dura_t1 = load_nifti_ras(dura_t1_path)
    dura_mask_image = load_nifti_ras(dura_mask_path)
    if dura_t1.shape != t1.shape or dura_mask_image.shape != t1.shape:
        raise RuntimeError("AtlasEM/Dura returned products outside the T1 grid")
    if not np.allclose(dura_t1.affine, t1.affine, rtol=0, atol=1e-4) or not np.allclose(
        dura_mask_image.affine, t1.affine, rtol=0, atol=1e-4
    ):
        raise RuntimeError("AtlasEM/Dura returned a product with a mismatched T1 affine")
    dura_mask = np.asarray(dura_mask_image.dataobj) > 0
    flair_values = np.asarray(flair.dataobj, dtype=np.float32)
    dura_flair = np.where(dura_mask, np.maximum(flair_values, 0.0), 0.0).astype(np.float32)
    dura_t1_values = np.asarray(dura_t1.dataobj, dtype=np.float32)
    modalities_values = np.stack((dura_t1_values, dura_flair), axis=-1)
    dura_flair_path = save_nifti_ras(
        nib.Nifti1Image(dura_flair, t1.affine, header=t1.header.copy()),
        root / "dura_masked_flair.nii.gz", dtype=np.float32,
    )
    modalities_path = save_nifti_ras(
        nib.Nifti1Image(modalities_values, t1.affine, header=t1.header.copy()),
        root / "dura_conditioned_t1_flair.nii.gz", dtype=np.float32,
    )
    outputs = {name: Path(path) for name, path in atlasem.output_paths.items()}
    outputs["dura_stripped_t1"] = Path(dura_t1_path)
    outputs["dura_masked_flair"] = dura_flair_path
    outputs["lesion_modalities"] = modalities_path
    manifest["status"] = "complete"
    manifest["outputs"] = {name: _record(path) for name, path in outputs.items()}
    _write_manifest(manifest_path, manifest)
    return PlusAtlasEMDuraHandoffResult(
        manifest_path, modalities_path, Path(dura_t1_path), dura_flair_path, outputs, False
    )
