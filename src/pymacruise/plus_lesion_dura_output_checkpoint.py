"""Persist source-selected LesionTOADS Dura-removal channels at P41m."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.image_io import nifti_ras_on_reference, save_nifti_ras
from pymacruise.lesion_dura_outputs import generate_lesion_dura_removal_outputs


@dataclass(frozen=True)
class LesionDuraOutputCheckpointResult:
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


def _record(path: Path) -> dict[str, str]:
    if not path.is_file() or path.is_symlink():
        raise FileNotFoundError(f"Missing regular checkpoint input: {path}")
    return {"path": str(path.resolve()), "sha256": _sha256(path)}


def _checked_output(manifest: dict[str, Any], name: str) -> Path:
    record = manifest.get("outputs", {}).get(name)
    if not isinstance(record, dict):
        raise ValueError(f"Completed input manifest has no {name!r} output")
    path = Path(record.get("path", "")).resolve()
    if not path.is_file() or path.is_symlink() or _sha256(path) != record.get("sha256"):
        raise ValueError(f"Input output {name!r} is missing or fails its manifest hash")
    return path


def _complete_manifest(path: Path, expected_stage: str) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink():
        raise FileNotFoundError(f"Missing regular {expected_stage} manifest: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("status") != "complete" or value.get("stage") != expected_stage:
        raise ValueError(f"Input manifest is not a completed {expected_stage}: {path}")
    return value


def _verified_output(
    record: dict[str, Any], name: str, reference: nib.spatialimages.SpatialImage
) -> Path | None:
    path = Path(record.get("path", ""))
    if (
        not path.is_file()
        or path.is_symlink()
        or _sha256(path) != record.get("sha256")
    ):
        return None
    image = nib.load(str(path))
    if (
        not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image))
        or len(image.shape) != 3
        or image.get_data_dtype() != np.dtype(np.float32)
        or nib.aff2axcodes(image.affine) != ("R", "A", "S")
        or image.shape != reference.shape
        or not np.allclose(image.affine, reference.affine, rtol=0, atol=1e-5)
        or not np.isfinite(np.asarray(image.dataobj)).all()
    ):
        return None
    return path


def run_plus_lesion_dura_output_checkpoint(
    registered_recrop_manifest_path: str | Path,
    registered_loop_manifest_path: str | Path,
    output_dir: str | Path,
    *,
    resume: bool = False,
) -> LesionDuraOutputCheckpointResult:
    """Apply `generateDuraRemovalOutputs` to persisted P41m state.

    The four channel images are saved on the registered recrop grid. Native
    T1-grid transfer and the later MaCRUISEplus surface workflow are separate
    stages and are not inferred here.
    """
    if not isinstance(resume, bool):
        raise TypeError("resume must be boolean")
    recrop_path = Path(registered_recrop_manifest_path).expanduser().resolve()
    loop_path = Path(registered_loop_manifest_path).expanduser().resolve()
    root = Path(output_dir).expanduser().resolve()
    recrop = _complete_manifest(recrop_path, "plus-lesion-registered-recrops")
    loop = _complete_manifest(loop_path, "plus-lesion-registered-loop")

    if recrop.get("spatial_convention") != "RAS" or loop.get("spatial_convention") != "RAS":
        raise ValueError("registered recrop and LesionTOADS loop checkpoints must be RAS")
    modalities_path = _checked_output(recrop, "modalities")
    summary_path = _checked_output(recrop, "summary")
    memberships_path = _checked_output(loop, "memberships")
    segmentation_path = _checked_output(loop, "segmentation")
    input_files = {
        "recropped_modalities": _record(modalities_path),
        "registered_recrop_summary": _record(summary_path),
        "registered_memberships": _record(memberships_path),
        "registered_segmentation": _record(segmentation_path),
    }
    upstream_manifests = {
        "registered_recrop_manifest": _record(recrop_path),
        "registered_loop_manifest": _record(loop_path),
    }
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    class_names = tuple(summary.get("class_names", ()))
    if len(class_names) != 13:
        raise ValueError("registered recrop must describe the selected 13-class atlas order")
    reference4d = nib.load(str(modalities_path))
    memberships_image = nib.load(str(memberships_path))
    segmentation_image = nib.load(str(segmentation_path))
    if not all(isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)) for image in (
        reference4d, memberships_image, segmentation_image
    )):
        raise ValueError("P41m state inputs must be NIfTI images")
    if (
        len(reference4d.shape) != 4
        or reference4d.shape[-1] != 2
        or memberships_image.shape != (*reference4d.shape[:3], 17)
        or segmentation_image.shape != reference4d.shape[:3]
        or any(nib.aff2axcodes(image.affine) != ("R", "A", "S") for image in (
            reference4d, memberships_image, segmentation_image
        ))
        or any(not np.allclose(image.affine, reference4d.affine, rtol=0, atol=1e-5) for image in (
            memberships_image, segmentation_image
        ))
    ):
        raise ValueError("P41m memberships and segmentation must share the registered recrop RAS grid")
    reference = nib.Nifti1Image(
        np.zeros(reference4d.shape[:3], dtype=np.uint8),
        reference4d.affine,
        header=reference4d.header,
    )
    signature_data = {
        "schema_version": 1,
        "stage": "plus-lesion-dura-output-checkpoint",
        "source_order": ["generateDuraRemovalOutputs"],
        "channel_order": ["csf", "gm", "wm_fill", "wm_mask"],
        "class_names": list(class_names),
        "inputs": input_files,
        "upstream_manifests": upstream_manifests,
        "spatial_convention": "RAS",
        "source_equivalence_claim": False,
    }
    signature_payload = {
        key: value for key, value in signature_data.items() if key != "upstream_manifests"
    }
    signature = hashlib.sha256(
        json.dumps(signature_payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "lesion-dura-output-checkpoint.json"
    if manifest_path.is_symlink():
        raise FileExistsError(f"Refusing to write through symlinked manifest: {manifest_path}")
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"Dura output checkpoint exists; use resume=True: {manifest_path}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        # A completed upstream loop may have been repaired during resume, so
        # its manifest/hash can legitimately change. Recompute this downstream
        # checkpoint from the newly verified inputs instead of rejecting the
        # parent workflow's dependency invalidation.
        if manifest.get("run_signature") == signature and manifest.get("status") == "complete":
            verified = {
                name: _verified_output(record, name, reference)
                for name, record in manifest.get("outputs", {}).items()
            }
            if set(verified) == {"csf", "gm", "wm_fill", "wm_mask"} and all(verified.values()):
                return LesionDuraOutputCheckpointResult(
                    manifest_path, {name: path for name, path in verified.items() if path is not None}, True
                )
    elif resume:
        raise FileNotFoundError(f"Cannot resume; no Dura output manifest exists: {manifest_path}")

    manifest = {
        **signature_data,
        "run_signature": signature,
        "manifest_path": str(manifest_path),
        "status": "running",
        "outputs": {},
    }
    _write_json(manifest_path, manifest)
    try:
        packed = np.asarray(memberships_image.dataobj, dtype=np.float32)
        segmentation = np.asarray(segmentation_image.dataobj)
        result = generate_lesion_dura_removal_outputs(
            packed[..., :13],
            segmentation,
            class_names,
            packed[..., 13],
            packed[..., 14],
        )
        channel_values = {
            "csf": result.csf,
            "gm": result.gray_matter,
            "wm_fill": result.white_matter_fill,
            "wm_mask": result.white_matter_mask,
        }
        outputs: dict[str, Path] = {}
        records: dict[str, dict[str, str]] = {}
        for name, values in channel_values.items():
            final_path = root / f"lesion_toads_{name}.nii.gz"
            temporary_path = root / f".{final_path.stem}.partial.nii.gz"
            save_nifti_ras(
                nifti_ras_on_reference(reference, values, dtype=np.float32),
                temporary_path,
                dtype=np.float32,
            )
            temporary_path.replace(final_path)
            outputs[name] = final_path
            records[name] = {"path": str(final_path.resolve()), "sha256": _sha256(final_path)}
        manifest["outputs"] = records
        manifest["status"] = "complete"
        _write_json(manifest_path, manifest)
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        _write_json(manifest_path, manifest)
        raise
    return LesionDuraOutputCheckpointResult(manifest_path, outputs, False)
