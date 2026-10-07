"""Compose persisted plus LesionTOADS checkpoints into one resumable command.

This workflow connects existing P41i, P41k, and P41m file stages. It adds no
segmentation algorithm and does not claim terminal LesionTOADS classification
or source-composite equivalence.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .plus_lesion_first_iteration import run_plus_lesion_first_iteration_checkpoint
from .plus_lesion_registered_branch import run_plus_lesion_registered_branch_checkpoint
from .plus_lesion_registered_loop import persist_plus_lesion_registered_loop_checkpoint


@dataclass(frozen=True)
class PlusLesionLoopWorkflowResult:
    manifest_path: Path
    first_iteration_manifest: Path
    registered_branch_manifest: Path
    registered_loop_manifest: Path
    resumed: bool
    classification_path: Path | None = None
    lesion_subtype_path: Path | None = None


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


def _file_record(path: Path) -> dict[str, str]:
    if not path.is_file() or path.is_symlink():
        raise FileNotFoundError(f"Expected a regular stage manifest: {path}")
    return {"path": str(path.resolve()), "sha256": _sha256(path)}


def _verify_child(path: Path, *, stage: str) -> dict[str, Any]:
    record = json.loads(path.read_text(encoding="utf-8"))
    if record.get("status") != "complete" or record.get("stage") != stage:
        raise RuntimeError(f"Child stage {stage!r} did not complete: {path}")
    for name, output in record.get("outputs", {}).items():
        output_path = Path(output.get("path", ""))
        if (
            not output_path.is_file()
            or output_path.is_symlink()
            or _sha256(output_path) != output.get("sha256")
        ):
            raise RuntimeError(f"Child stage {stage!r} output {name!r} failed hash verification")
    return record


def run_plus_lesion_loop_workflow(
    registered_recrop_manifest_path: str | Path,
    atlas_descriptor_path: str | Path,
    output_dir: str | Path,
    *,
    iterations: int = 99,
    first_thinning_limit: float = 3.0,
    first_thinning_volume_ratio: float = 0.5,
    source_zero_variance_nan: bool = True,
    max_distance: float = 1e-4,
    n_iter_max: int = 99,
    thinning_limit: float = 3.0,
    max_passes: int = 99,
    loop_convergence_threshold: float | None = None,
    correction_fields_path: str | Path | None = None,
    correction_polynomial_degree: int | None = None,
    initial_correction_degree: int | None = None,
    resume: bool = False,
) -> PlusLesionLoopWorkflowResult:
    """Run P41i -> P41k -> P41m with a hash-checked top-level manifest.

    Completed child stages are reused only through their own validated
    checkpoint APIs. The top-level record marks each stage running before it
    starts and records the child manifest hash after completion. If a process
    stops during a child stage, resume delegates recovery to that child's
    checkpoint implementation.
    """
    root = Path(output_dir).expanduser().resolve()
    parent_path = Path(registered_recrop_manifest_path).expanduser().resolve()
    descriptor = Path(atlas_descriptor_path).expanduser().resolve()
    if not parent_path.is_file() or parent_path.is_symlink():
        raise FileNotFoundError(f"Missing regular registered-recrop manifest: {parent_path}")
    if not descriptor.is_file() or descriptor.is_symlink():
        raise FileNotFoundError(f"Missing regular atlas descriptor: {descriptor}")
    settings = {
        "iterations": int(iterations),
        "first_thinning_limit": float(first_thinning_limit),
        "first_thinning_volume_ratio": float(first_thinning_volume_ratio),
        "source_zero_variance_nan": source_zero_variance_nan,
        "max_distance": float(max_distance),
        "n_iter_max": int(n_iter_max),
        "thinning_limit": float(thinning_limit),
        "max_passes": int(max_passes),
        "loop_convergence_threshold": (
            None if loop_convergence_threshold is None else float(loop_convergence_threshold)
        ),
        "correction_polynomial_degree": correction_polynomial_degree,
        "initial_correction_degree": initial_correction_degree,
    }
    correction_fields_record = None
    if correction_fields_path is not None:
        fields_path = Path(correction_fields_path).expanduser().resolve()
        if not fields_path.name.lower().endswith(".nii.gz"):
            raise ValueError("correction fields must use compressed NIfTI format (.nii.gz)")
        correction_fields_record = _file_record(fields_path)
    if loop_convergence_threshold is not None and (
        not np.isfinite(loop_convergence_threshold) or loop_convergence_threshold < 0
    ):
        raise ValueError("loop_convergence_threshold must be finite and nonnegative")
    signature_data = {
        "stage": "plus-lesion-loop-workflow",
        "registered_recrop_manifest": _file_record(parent_path),
        "atlas_descriptor": _file_record(descriptor),
        "correction_fields": correction_fields_record,
        "settings": settings,
        "spatial_convention": "RAS",
    }
    signature = hashlib.sha256(
        json.dumps(signature_data, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "plus-lesion-loop-workflow.json"
    if manifest_path.is_symlink():
        raise FileExistsError(f"Refusing to write through symlinked workflow manifest: {manifest_path}")
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"Workflow manifest exists; use resume=True: {manifest_path}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("run_signature") != signature:
            raise ValueError("Cannot resume plus LesionTOADS workflow: inputs or settings differ")
    else:
        if resume:
            raise FileNotFoundError(f"Cannot resume; workflow manifest does not exist: {manifest_path}")
        manifest = {
            **signature_data,
            "run_signature": signature,
            "status": "running",
            "source_equivalence_claim": False,
            "stages": {},
        }
        _write_json(manifest_path, manifest)

    child_paths = {
        "first_iteration": root / "first-iteration" / "lesion-first-iteration.json",
        "registered_branch": root / "registered-branch" / "lesion-registered-branch.json",
        "registered_loop": root / "registered-loop" / "lesion-registered-loop.json",
    }
    child_definitions = (
        ("first_iteration", "plus-lesion-first-iteration"),
        ("registered_branch", "plus-lesion-registered-branch"),
        ("registered_loop", "plus-lesion-registered-loop"),
    )
    for name, stage in child_definitions:
        child_resume = bool(resume and child_paths[name].is_file())
        manifest["status"] = "running"
        manifest["active_stage"] = name
        manifest["stages"][name] = {"status": "running"}
        _write_json(manifest_path, manifest)
        try:
            if name == "first_iteration":
                result = run_plus_lesion_first_iteration_checkpoint(
                    parent_path,
                    descriptor,
                    root / "first-iteration",
                    iterations=iterations,
                    first_thinning_limit=first_thinning_limit,
                    first_thinning_volume_ratio=first_thinning_volume_ratio,
                    source_zero_variance_nan=source_zero_variance_nan,
                    resume=child_resume,
                )
            elif name == "registered_branch":
                result = run_plus_lesion_registered_branch_checkpoint(
                    parent_path,
                    descriptor,
                    root / "registered-branch",
                    iterations=iterations,
                    first_thinning_limit=first_thinning_limit,
                    first_thinning_volume_ratio=first_thinning_volume_ratio,
                    source_zero_variance_nan=source_zero_variance_nan,
                    max_distance=max_distance,
                    first_iteration_manifest_path=child_paths["first_iteration"],
                    resume=child_resume,
                )
            else:
                result = persist_plus_lesion_registered_loop_checkpoint(
                    child_paths["registered_branch"],
                    descriptor,
                    root / "registered-loop",
                    n_iter_max=n_iter_max,
                    thinning_limit=thinning_limit,
                    source_zero_variance_nan=source_zero_variance_nan,
                    max_passes=max_passes,
                    convergence_threshold=loop_convergence_threshold,
                    correction_fields_path=correction_fields_path,
                    correction_polynomial_degree=correction_polynomial_degree,
                    initial_correction_degree=initial_correction_degree,
                    resume=child_resume,
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
        if child_path != child_paths[name].resolve():
            # The stage APIs own their stable manifest names; the final path
            # is recorded so a future naming change remains provenance-safe.
            child_paths[name] = child_path
        child = _verify_child(child_path, stage=stage)
        manifest["stages"][name] = {
            "status": "complete",
            "manifest": _file_record(child_path),
            "output_count": len(child.get("outputs", {})),
            "resumed": bool(getattr(result, "resumed", False)),
        }
        _write_json(manifest_path, manifest)

    manifest["status"] = "complete"
    manifest["active_stage"] = None
    manifest["workflow_executed"] = True
    manifest["final_outputs"] = manifest["stages"]["registered_loop"]["manifest"]
    final_loop = _verify_child(child_paths["registered_loop"], stage="plus-lesion-registered-loop")
    loop_summary = json.loads(Path(final_loop["outputs"]["summary"]["path"]).read_text(encoding="utf-8"))
    terminal = bool(loop_summary.get("stop"))
    manifest["terminal"] = terminal
    manifest["terminal_stop_reason"] = loop_summary.get("stop_reason")
    manifest["result_scope"] = (
        "synthetic selected LesionTOADS terminal classification and lesion subtype"
        if terminal else
        "synthetic selected LesionTOADS continuation through a bounded registered-loop segment"
    )
    _write_json(manifest_path, manifest)
    return PlusLesionLoopWorkflowResult(
        manifest_path=manifest_path,
        first_iteration_manifest=child_paths["first_iteration"],
        registered_branch_manifest=child_paths["registered_branch"],
        registered_loop_manifest=child_paths["registered_loop"],
        resumed=resume,
        classification_path=(
            Path(final_loop["outputs"]["classification"]["path"])
            if "classification" in final_loop["outputs"] else None
        ),
        lesion_subtype_path=(
            Path(final_loop["outputs"]["lesion_subtype"]["path"])
            if "lesion_subtype" in final_loop["outputs"] else None
        ),
    )
