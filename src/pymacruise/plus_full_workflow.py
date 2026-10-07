"""Restartable plus input-to-terminal LesionTOADS workflow orchestration."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.image_io import normalize_image_to_ras
from pymacruise.lesion_atlas_priors import load_lesion_atlas_topology
from pymacruise.n4 import N4Config
from pymacruise.plus_maseg_cortical_checkpoint import run_plus_maseg_cortical_checkpoint
from pymacruise.plus_file_workflow import run_plus_file_workflow
from pymacruise.plus_pipeline import PlusPipelineOptions
from pymacruise.plus_lesion_atlas_init import run_plus_lesion_initial_atlas_alignment_checkpoint
from pymacruise.plus_lesion_crop import run_plus_lesion_crop_checkpoint
from pymacruise.plus_lesion_dura_output_checkpoint import run_plus_lesion_dura_output_checkpoint
from pymacruise.plus_lesion_preprocessing import run_plus_lesion_preprocessing_checkpoint
from pymacruise.plus_lesion_registered_recrop import run_plus_lesion_registered_recrop_checkpoint
from pymacruise.plus_lesion_loop_workflow import run_plus_lesion_loop_workflow
from pymacruise.plus_lesion_native_dura_transfer import run_plus_lesion_native_dura_transfer
from pymacruise.registration import RegistrationConfig


@dataclass(frozen=True)
class PlusFullWorkflowResult:
    """Top-level files and optional terminal outputs for a plus run."""

    manifest_path: Path
    preprocessing_manifest: Path
    crop_manifest: Path
    initial_alignment_manifest: Path
    registered_recrop_manifest: Path
    lesion_workflow_manifest: Path
    classification_path: Path | None
    lesion_subtype_path: Path | None
    dura_output_paths: dict[str, Path]
    native_dura_output_paths: dict[str, Path]
    surface_output_paths: dict[str, Path]
    resumed: bool


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def _file_record(path: Path) -> dict[str, str]:
    if not path.is_file() or path.is_symlink():
        raise FileNotFoundError(f"Missing regular workflow file: {path}")
    return {"path": str(path.resolve()), "sha256": _sha256(path)}


def _verify_stage(path: Path, stage: str | None) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink():
        raise FileNotFoundError(f"Missing {stage} manifest: {path}")
    record = json.loads(path.read_text(encoding="utf-8"))
    if record.get("status") != "complete" or (stage is not None and record.get("stage") != stage):
        raise RuntimeError(f"Stage {stage!r} did not complete: {path}")
    for name, output in record.get("outputs", {}).items():
        output_path = Path(output.get("path", ""))
        if (
            not output_path.is_file()
            or output_path.is_symlink()
            or _sha256(output_path) != output.get("sha256")
        ):
            raise RuntimeError(f"Stage {stage!r} output {name!r} failed hash verification")
    return record


def run_plus_full_workflow(
    t1_path: str | Path,
    flair_path: str | Path,
    brain_mask_path: str | Path,
    maseg_path: str | Path,
    atlas_descriptor_path: str | Path,
    output_dir: str | Path,
    *,
    registration_config: RegistrationConfig = RegistrationConfig(),
    n4_config: N4Config = N4Config(),
    crop_padding: int = 10,
    init_align_iterations: int = 10,
    pyramid_levels: int = 2,
    lesion_iterations: int = 99,
    first_thinning_limit: float = 3.0,
    first_thinning_volume_ratio: float = 0.5,
    source_zero_variance_nan: bool = True,
    max_distance: float = 1e-4,
    loop_convergence_threshold: float | None = None,
    n_iter_max: int = 99,
    loop_thinning_limit: float = 3.0,
    max_passes: int = 99,
    pipeline_options: PlusPipelineOptions = PlusPipelineOptions(),
    resume: bool = False,
) -> PlusFullWorkflowResult:
    """Run raw inputs through LesionTOADS and terminal CRUISE surface outputs.

    The required brain mask is always a user input on the T1 grid. BET and
    SPECTRE are not invoked. The supplied MaSeg and brain mask must match T1.
    """
    if not isinstance(registration_config, RegistrationConfig):
        raise TypeError("registration_config must be RegistrationConfig")
    if not isinstance(n4_config, N4Config):
        raise TypeError("n4_config must be N4Config")
    if not isinstance(pipeline_options, PlusPipelineOptions):
        raise TypeError("pipeline_options must be PlusPipelineOptions")
    if not isinstance(resume, bool):
        raise TypeError("resume must be boolean")
    t1 = Path(t1_path).expanduser().resolve()
    flair = Path(flair_path).expanduser().resolve()
    mask = Path(brain_mask_path).expanduser().resolve()
    maseg = Path(maseg_path).expanduser().resolve()
    descriptor = Path(atlas_descriptor_path).expanduser().resolve()
    root = Path(output_dir).expanduser().resolve()
    # Reject grid mismatches before costly preprocessing/registration. Compare
    # physical RAS grids after axis reorientation; no interpolation is done.
    t1_image = nib.load(str(t1))
    if not isinstance(t1_image, (nib.Nifti1Image, nib.Nifti2Image)) or len(t1_image.shape) != 3:
        raise ValueError("T1 must be a 3-D NIfTI reference image")
    t1_ras = normalize_image_to_ras(t1_image, source_name="T1")
    for path, name in ((mask, "brain mask"), (maseg, "MaSeg")):
        image = nib.load(str(path))
        if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)) or len(image.shape) != 3:
            raise ValueError(f"{name} must be a 3-D NIfTI image")
        image_ras = normalize_image_to_ras(image, source_name=name)
        if image_ras.shape != t1_ras.shape or not np.allclose(
            image_ras.affine, t1_ras.affine, rtol=0, atol=1e-5
        ):
            raise ValueError(f"{name} must match the T1 reference grid after RAS normalization")
    inputs = {
        "t1": _file_record(t1),
        "flair": _file_record(flair),
        "brain_mask": _file_record(mask),
        "maseg": _file_record(maseg),
        "atlas_descriptor": _file_record(descriptor),
    }
    # Resolve and validate the descriptor-owned topology file before starting
    # child stages; it is required by the selected crop computation.
    topology = load_lesion_atlas_topology(descriptor)
    inputs["topology"] = _file_record(topology.path)
    settings: dict[str, Any] = {
        "registration": asdict(registration_config),
        "n4": asdict(n4_config),
        "pipeline_options": pipeline_options.as_dict(),
        "crop_padding": int(crop_padding),
        "init_align_iterations": int(init_align_iterations),
        "pyramid_levels": int(pyramid_levels),
        "lesion_iterations": int(lesion_iterations),
        "first_thinning_limit": float(first_thinning_limit),
        "first_thinning_volume_ratio": float(first_thinning_volume_ratio),
        "source_zero_variance_nan": bool(source_zero_variance_nan),
        "max_distance": float(max_distance),
        "loop_convergence_threshold": (
            None if loop_convergence_threshold is None else float(loop_convergence_threshold)
        ),
        "n_iter_max": int(n_iter_max),
        "loop_thinning_limit": float(loop_thinning_limit),
        "max_passes": int(max_passes),
    }
    signature_data = {
        "schema_version": 1,
        "stage": "plus-full-workflow",
        "stage_order": [
            "preprocessing", "crop", "initial_alignment", "registered_recrop",
            "lesion_workflow", "dura_outputs", "native_dura_transfer",
            "maseg_cortical_separation", "plus_surface_workflow",
        ],
        "spatial_convention": "RAS",
        "inputs": inputs,
        "settings": settings,
        "source_equivalence_claim": False,
    }
    signature = hashlib.sha256(
        json.dumps(signature_data, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "plus-full-workflow.json"
    if manifest_path.is_symlink():
        raise FileExistsError(f"Refusing to write through symlinked manifest: {manifest_path}")
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"Plus workflow manifest exists; use resume=True: {manifest_path}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("run_signature") != signature:
            raise ValueError("Cannot resume plus workflow: inputs or settings differ")
    else:
        if resume:
            raise FileNotFoundError(f"Cannot resume; no plus workflow manifest exists: {manifest_path}")
        manifest = {
            **signature_data,
            "run_signature": signature,
            "status": "running",
            "active_stage": None,
            "stages": {},
        }
        _write_json(manifest_path, manifest)

    stage_paths = {
        "preprocessing": root / "preprocessing" / "lesion-preprocessing.json",
        "crop": root / "crop" / "lesion-crop.json",
        "initial_alignment": root / "initial_alignment" / "lesion-initial-atlas-alignment.json",
        "registered_recrop": root / "registered_recrop" / "lesion-registered-recrops.json",
        "lesion_workflow": root / "lesion_workflow" / "plus-lesion-loop-workflow.json",
        "dura_outputs": root / "dura-outputs" / "lesion-dura-output-checkpoint.json",
        "native_dura_transfer": root / "native-dura" / "native-dura-transfer.json",
        "maseg_cortical_separation": root / "maseg-cortical-separation" / "maseg-cortical-separation.json",
        "plus_surface_workflow": root / "plus-surface" / "run.json",
    }
    definitions = (
        ("preprocessing", "plus-lesion-input-preparation"),
        ("crop", "plus-lesion-normalization-crop"),
        ("initial_alignment", "plus-lesion-initial-atlas-alignment"),
        ("registered_recrop", "plus-lesion-registered-recrops"),
        ("lesion_workflow", "plus-lesion-loop-workflow"),
        ("dura_outputs", "plus-lesion-dura-output-checkpoint"),
        ("native_dura_transfer", "plus-lesion-native-dura-transfer"),
        ("maseg_cortical_separation", "plus-maseg-cortical-separation"),
        ("plus_surface_workflow", None),
    )
    results: dict[str, Any] = {}
    for name, child_stage in definitions:
        # Resume completed child manifests, but let a downstream stage that
        # never started run fresh after an interruption in its parent.
        child_resume = bool(resume and stage_paths[name].is_file())
        manifest["status"] = "running"
        manifest["active_stage"] = name
        manifest["stages"][name] = {"status": "running"}
        _write_json(manifest_path, manifest)
        try:
            if name == "preprocessing":
                result = run_plus_lesion_preprocessing_checkpoint(
                    t1, flair, mask, root / name,
                    registration_config=registration_config,
                    n4_config=n4_config,
                    resume=child_resume,
                )
            elif name == "crop":
                prior = results["preprocessing"]
                result = run_plus_lesion_crop_checkpoint(
                    prior.modalities_path, topology.path, root / name,
                    padding=crop_padding, resume=child_resume,
                )
            elif name == "initial_alignment":
                prior = results["crop"]
                result = run_plus_lesion_initial_atlas_alignment_checkpoint(
                    prior.cropped_path, prior.manifest_path, descriptor, root / name,
                    init_align_iterations=init_align_iterations,
                    pyramid_levels=pyramid_levels,
                    resume=child_resume,
                )
            elif name == "registered_recrop":
                crop = results["crop"]
                alignment = results["initial_alignment"]
                result = run_plus_lesion_registered_recrop_checkpoint(
                    crop.normalized_path, crop.manifest_path, alignment.manifest_path,
                    descriptor, root / name, resume=child_resume,
                )
            elif name == "lesion_workflow":
                recrop = results["registered_recrop"]
                result = run_plus_lesion_loop_workflow(
                    recrop.manifest_path, descriptor, root / name,
                    iterations=lesion_iterations,
                    first_thinning_limit=first_thinning_limit,
                    first_thinning_volume_ratio=first_thinning_volume_ratio,
                    source_zero_variance_nan=source_zero_variance_nan,
                    max_distance=max_distance,
                    loop_convergence_threshold=loop_convergence_threshold,
                    n_iter_max=n_iter_max,
                    thinning_limit=loop_thinning_limit,
                    max_passes=max_passes,
                    resume=child_resume,
                )
            elif name == "dura_outputs":
                recrop = results["registered_recrop"]
                loop = results["lesion_workflow"]
                result = run_plus_lesion_dura_output_checkpoint(
                    recrop.manifest_path,
                    loop.registered_loop_manifest,
                    root / "dura-outputs",
                    resume=child_resume,
                )
            elif name == "native_dura_transfer":
                dura = results["dura_outputs"]
                result = run_plus_lesion_native_dura_transfer(
                    dura.manifest_path, t1, root / "native-dura", resume=child_resume
                )
            elif name == "maseg_cortical_separation":
                result = run_plus_maseg_cortical_checkpoint(
                    maseg, root / "maseg-cortical-separation", resume=child_resume
                )
            else:
                native = results["native_dura_transfer"]
                cortical = results["maseg_cortical_separation"]
                plus = run_plus_file_workflow(
                    t1_reference_path=t1,
                    segmentation_path=maseg,
                    dura_csf_path=native.outputs["csf"],
                    dura_gm_path=native.outputs["gm"],
                    dura_wm_fill_path=native.outputs["wm_fill"],
                    dura_wm_mask_path=native.outputs["wm_mask"],
                    multi_atlas_gm_path=cortical.outputs["Atlas_GM"],
                    multi_atlas_wm_path=cortical.outputs["Atlas_WM"],
                    multi_atlas_mask_path=cortical.outputs["Atlas_Mask"],
                    brain_mask_path=mask,
                    output_dir=root / "plus-surface",
                    options=pipeline_options,
                    resume=child_resume,
                )
                result = plus
        except Exception as exc:
            manifest["status"] = "failed"
            manifest["active_stage"] = name
            manifest["stages"][name] = {
                "status": "failed",
                "error": f"{type(exc).__name__}: {exc}",
            }
            _write_json(manifest_path, manifest)
            raise
        child_path = Path(result.manifest_path).resolve()
        child = _verify_stage(child_path, child_stage)
        stage_paths[name] = child_path
        results[name] = result
        manifest["stages"][name] = {
            "status": "complete",
            "manifest": _file_record(child_path),
            "output_count": len(child.get("outputs", {})),
            "resumed": bool(
                getattr(result, "resumed", False)
                or getattr(result, "resumed_stages", ())
            ),
        }
        _write_json(manifest_path, manifest)

    loop_result = results["lesion_workflow"]
    dura_result = results["dura_outputs"]
    native_dura_result = results["native_dura_transfer"]
    plus_surface_result = results["plus_surface_workflow"]
    manifest["status"] = "complete"
    manifest["active_stage"] = None
    manifest["terminal"] = loop_result.classification_path is not None
    manifest["final_outputs"] = {
        "dura_checkpoint": _file_record(dura_result.manifest_path),
        **{
            name: _file_record(path)
            for name, path in dura_result.outputs.items()
        },
        **{
            f"native_{name}": _file_record(path)
            for name, path in native_dura_result.outputs.items()
        },
        **{
            f"plus_{name}": _file_record(path)
            for name, path in plus_surface_result.outputs.items()
        },
    }
    _write_json(manifest_path, manifest)
    return PlusFullWorkflowResult(
        manifest_path,
        stage_paths["preprocessing"],
        stage_paths["crop"],
        stage_paths["initial_alignment"],
        stage_paths["registered_recrop"],
        stage_paths["lesion_workflow"],
        loop_result.classification_path,
        loop_result.lesion_subtype_path,
        dura_result.outputs,
        native_dura_result.outputs,
        plus_surface_result.outputs,
        resume,
    )
