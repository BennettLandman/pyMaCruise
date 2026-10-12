"""Restartable plus input-to-terminal LesionTOADS workflow orchestration."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.image_io import (
    is_axis_aligned_ras_grid,
    normalize_image_to_ras,
    regrid_to_orthogonal_ras,
    save_nifti_ras,
)
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


def _implementation_signature(source_root: Path | None = None) -> str:
    """Fingerprint package code so a resume cannot accept stale stage outputs.

    The plus runner has child workflows whose input and option signatures are
    stable across implementation edits. A source change therefore needs a
    fresh output generation while preserving the prior generation for review.
    """
    root = Path(source_root) if source_root is not None else Path(__file__).resolve().parent
    sources = sorted(root.rglob("*.py"))
    if not sources:
        raise RuntimeError(f"No Python implementation files found beneath {root}")
    return hashlib.sha256(json.dumps([
        {"path": path.relative_to(root).as_posix(), "sha256": _sha256(path)}
        for path in sources
    ], sort_keys=True, separators=(",", ":")).encode()).hexdigest()


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
    regrid_inputs: bool = True,
    dura_transfer_mode: str = "physical-affine",
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
    if not isinstance(regrid_inputs, bool):
        raise TypeError("regrid_inputs must be boolean")
    if dura_transfer_mode not in ("physical-affine", "jist-unchanged"):
        raise ValueError("dura_transfer_mode must be 'physical-affine' or 'jist-unchanged'")
    if not isinstance(resume, bool):
        raise TypeError("resume must be boolean")
    t1 = Path(t1_path).expanduser().resolve()
    flair = Path(flair_path).expanduser().resolve()
    mask = Path(brain_mask_path).expanduser().resolve()
    maseg = Path(maseg_path).expanduser().resolve()
    descriptor = Path(atlas_descriptor_path).expanduser().resolve()
    root = Path(output_dir).expanduser().resolve()
    # Reject source-grid mismatches before costly work. Regridding changes only
    # an oblique T1 grid; its paired mask and MaSeg use nearest-neighbor sampling.
    t1_image = nib.load(str(t1))
    if not isinstance(t1_image, (nib.Nifti1Image, nib.Nifti2Image)) or len(t1_image.shape) != 3:
        raise ValueError("T1 must be a 3-D NIfTI reference image")
    t1_ras = normalize_image_to_ras(t1_image, source_name="T1")
    if not regrid_inputs and not is_axis_aligned_ras_grid(t1_ras):
        raise ValueError(
            "TOADS requires an axis-aligned RAS grid; regridding is disabled, "
            "so provide an orthogonal RAS T1 or enable --regrid-inputs"
        )
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
        "regrid_inputs": regrid_inputs,
        "dura_transfer_mode": dura_transfer_mode,
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
    implementation_signature = _implementation_signature()
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
        previous_implementation = manifest.get("implementation_signature")
        if previous_implementation == implementation_signature:
            stage_root = Path(manifest.get("stage_root", root)).resolve()
            if stage_root != root and root not in stage_root.parents:
                raise ValueError("Cannot resume plus workflow: stage root escapes output directory")
        else:
            # Keep the completed prior generation intact. Child workflows use
            # fixed manifests and input-only signatures, so a code revision
            # must run in a new directory instead of silently reusing them.
            stage_root = root / "code-revisions" / implementation_signature[:16]
            manifest["previous_implementation_signature"] = previous_implementation
            manifest["implementation_signature"] = implementation_signature
            manifest["stage_root"] = str(stage_root)
            manifest["stages"] = {}
            manifest["active_stage"] = None
            manifest["status"] = "running"
            _write_json(manifest_path, manifest)
    else:
        if resume:
            raise FileNotFoundError(f"Cannot resume; no plus workflow manifest exists: {manifest_path}")
        manifest = {
            **signature_data,
            "run_signature": signature,
            "implementation_signature": implementation_signature,
            "stage_root": str(root),
            "status": "running",
            "active_stage": None,
            "stages": {},
        }
        stage_root = root
        _write_json(manifest_path, manifest)

    regridded_resources: list[dict[str, str]] = []
    if regrid_inputs and not is_axis_aligned_ras_grid(t1_ras):
        print(
            "Plus workflow: regridding oblique RAS T1, brain mask, and MaSeg "
            "to an orthogonal RAS grid (T1 linear; mask/labels nearest neighbor)."
        )
        orthogonal_t1 = regrid_to_orthogonal_ras(t1_ras, interpolation_order=1)
        orthogonal_mask = regrid_to_orthogonal_ras(
            normalize_image_to_ras(nib.load(str(mask)), source_name="brain mask"),
            reference=orthogonal_t1,
            interpolation_order=0,
        )
        orthogonal_maseg = regrid_to_orthogonal_ras(
            normalize_image_to_ras(nib.load(str(maseg)), source_name="MaSeg"),
            reference=orthogonal_t1,
            interpolation_order=0,
        )
        regridded_dir = root / "resources" / "regridded-inputs"
        regridded_dir.mkdir(parents=True, exist_ok=True)
        derived = {
            "t1": (regridded_dir / "T1.nii.gz", orthogonal_t1),
            "brain_mask": (regridded_dir / "BrainMask.nii.gz", orthogonal_mask),
            "maseg": (regridded_dir / "MaSeg.nii.gz", orthogonal_maseg),
        }
        for role, (path, image) in derived.items():
            save_nifti_ras(image, path)
            regridded_resources.append({"role": role, "path": str(path), "sha256": _sha256(path)})
        t1, mask, maseg = derived["t1"][0], derived["brain_mask"][0], derived["maseg"][0]

    stage_root.mkdir(parents=True, exist_ok=True)
    manifest["implementation_signature"] = implementation_signature
    manifest["stage_root"] = str(stage_root)
    manifest["derived_resources"] = {"regridded_inputs": regridded_resources}
    _write_json(manifest_path, manifest)

    stage_paths = {
        "preprocessing": stage_root / "preprocessing" / "lesion-preprocessing.json",
        "crop": stage_root / "crop" / "lesion-crop.json",
        "initial_alignment": stage_root / "initial_alignment" / "lesion-initial-atlas-alignment.json",
        "registered_recrop": stage_root / "registered_recrop" / "lesion-registered-recrops.json",
        "lesion_workflow": stage_root / "lesion_workflow" / "plus-lesion-loop-workflow.json",
        "dura_outputs": stage_root / "dura-outputs" / "lesion-dura-output-checkpoint.json",
        "native_dura_transfer": stage_root / "native-dura" / "native-dura-transfer.json",
        "maseg_cortical_separation": stage_root / "maseg-cortical-separation" / "maseg-cortical-separation.json",
        "plus_surface_workflow": stage_root / "plus-surface" / "run.json",
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
                    t1, flair, mask, stage_root / name,
                    registration_config=registration_config,
                    n4_config=n4_config,
                    resume=child_resume,
                )
            elif name == "crop":
                prior = results["preprocessing"]
                result = run_plus_lesion_crop_checkpoint(
                    prior.modalities_path, topology.path, stage_root / name,
                    padding=crop_padding, resume=child_resume,
                )
            elif name == "initial_alignment":
                prior = results["crop"]
                result = run_plus_lesion_initial_atlas_alignment_checkpoint(
                    prior.cropped_path, prior.manifest_path, descriptor, stage_root / name,
                    init_align_iterations=init_align_iterations,
                    pyramid_levels=pyramid_levels,
                    resume=child_resume,
                )
            elif name == "registered_recrop":
                crop = results["crop"]
                alignment = results["initial_alignment"]
                result = run_plus_lesion_registered_recrop_checkpoint(
                    crop.normalized_path, crop.manifest_path, alignment.manifest_path,
                    descriptor, stage_root / name, resume=child_resume,
                )
            elif name == "lesion_workflow":
                recrop = results["registered_recrop"]
                result = run_plus_lesion_loop_workflow(
                    recrop.manifest_path, descriptor, stage_root / name,
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
                    stage_root / "dura-outputs",
                    resume=child_resume,
                )
            elif name == "native_dura_transfer":
                dura = results["dura_outputs"]
                result = run_plus_lesion_native_dura_transfer(
                    dura.manifest_path, t1, stage_root / "native-dura",
                    transfer_mode=dura_transfer_mode,
                    resume=child_resume,
                )
            elif name == "maseg_cortical_separation":
                result = run_plus_maseg_cortical_checkpoint(
                    maseg, stage_root / "maseg-cortical-separation", resume=child_resume
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
                    multi_atlas_mask_path=cortical.outputs["Atlas_BrainTissueMask"],
                    brain_mask_path=mask,
                    output_dir=stage_root / "plus-surface",
                    dura_transfer_mode=dura_transfer_mode,
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
