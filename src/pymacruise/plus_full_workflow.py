"""Restartable plus input-to-terminal LesionTOADS workflow orchestration."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

from pymacruise.lesion_atlas_priors import load_lesion_atlas_topology
from pymacruise.n4 import N4Config
from pymacruise.plus_lesion_atlas_init import run_plus_lesion_initial_atlas_alignment_checkpoint
from pymacruise.plus_lesion_crop import run_plus_lesion_crop_checkpoint
from pymacruise.plus_lesion_preprocessing import run_plus_lesion_preprocessing_checkpoint
from pymacruise.plus_lesion_registered_recrop import run_plus_lesion_registered_recrop_checkpoint
from pymacruise.plus_lesion_loop_workflow import run_plus_lesion_loop_workflow
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


def _verify_stage(path: Path, stage: str) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink():
        raise FileNotFoundError(f"Missing {stage} manifest: {path}")
    record = json.loads(path.read_text(encoding="utf-8"))
    if record.get("status") != "complete" or record.get("stage") != stage:
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
    max_distance: float = 0.0,
    loop_convergence_threshold: float | None = None,
    n_iter_max: int = 99,
    loop_thinning_limit: float = 3.0,
    max_passes: int = 1,
    resume: bool = False,
) -> PlusFullWorkflowResult:
    """Run preprocessing through a terminal or continuing LesionTOADS result.

    The required brain mask is always a user input on the T1 grid. BET and
    SPECTRE are not invoked. Dura-to-CRUISE surface processing remains a
    downstream workflow; this runner ends at the selected P41m checkpoint.
    """
    if not isinstance(registration_config, RegistrationConfig):
        raise TypeError("registration_config must be RegistrationConfig")
    if not isinstance(n4_config, N4Config):
        raise TypeError("n4_config must be N4Config")
    if not isinstance(resume, bool):
        raise TypeError("resume must be boolean")
    t1 = Path(t1_path).expanduser().resolve()
    flair = Path(flair_path).expanduser().resolve()
    mask = Path(brain_mask_path).expanduser().resolve()
    descriptor = Path(atlas_descriptor_path).expanduser().resolve()
    root = Path(output_dir).expanduser().resolve()
    inputs = {
        "t1": _file_record(t1),
        "flair": _file_record(flair),
        "brain_mask": _file_record(mask),
        "atlas_descriptor": _file_record(descriptor),
    }
    # Resolve and validate the descriptor-owned topology file before starting
    # child stages; it is required by the selected crop computation.
    topology = load_lesion_atlas_topology(descriptor)
    inputs["topology"] = _file_record(topology.path)
    settings: dict[str, Any] = {
        "registration": asdict(registration_config),
        "n4": asdict(n4_config),
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
            "preprocessing", "crop", "initial_alignment", "registered_recrop", "lesion_workflow"
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
        "initial_alignment": root / "initial-alignment" / "lesion-initial-atlas-alignment.json",
        "registered_recrop": root / "registered-recrop" / "lesion-registered-recrop.json",
        "lesion_workflow": root / "lesion-workflow" / "plus-lesion-loop-workflow.json",
    }
    definitions = (
        ("preprocessing", "plus-lesion-input-preparation"),
        ("crop", "plus-lesion-normalization-crop"),
        ("initial_alignment", "plus-lesion-initial-atlas-alignment"),
        ("registered_recrop", "plus-lesion-registered-recrops"),
        ("lesion_workflow", "plus-lesion-loop-workflow"),
    )
    results: dict[str, Any] = {}
    for name, child_stage in definitions:
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
                    resume=resume,
                )
            elif name == "crop":
                prior = results["preprocessing"]
                result = run_plus_lesion_crop_checkpoint(
                    prior.modalities_path, topology.path, root / name,
                    padding=crop_padding, resume=resume,
                )
            elif name == "initial_alignment":
                prior = results["crop"]
                result = run_plus_lesion_initial_atlas_alignment_checkpoint(
                    prior.cropped_path, prior.manifest_path, descriptor, root / name,
                    init_align_iterations=init_align_iterations,
                    pyramid_levels=pyramid_levels,
                    resume=resume,
                )
            elif name == "registered_recrop":
                crop = results["crop"]
                alignment = results["initial_alignment"]
                result = run_plus_lesion_registered_recrop_checkpoint(
                    crop.normalized_path, crop.manifest_path, alignment.manifest_path,
                    descriptor, root / name, resume=resume,
                )
            else:
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
                    resume=resume,
                )
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
            "resumed": bool(getattr(result, "resumed", False)),
        }
        _write_json(manifest_path, manifest)

    loop_result = results["lesion_workflow"]
    manifest["status"] = "complete"
    manifest["active_stage"] = None
    manifest["terminal"] = loop_result.classification_path is not None
    manifest["final_outputs"] = _file_record(loop_result.manifest_path)
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
        resume,
    )
