"""Restartable file-based MaCRUISEplus continuation from P39 Dura products.

This workflow consumes already aligned/supplied segmentation and atlas maps plus
persisted P39 Dura channels. It does not execute atlas registration, TOADS, or
LesionTOADS. Stage records make transfer/crop, topology-through-TGDM, and
surface/thickness restartable independently.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.cropping import (
    CropParameters,
    crop_image,
    crop_image_with_reference_bounds,
    find_crop_parameters,
)
from pymacruise.image_io import normalize_image_to_ras, save_nifti_ras
from pymacruise.lesion_post_loop_outputs import (
    LesionToadsOutputProduct,
    LesionToadsPostLoopOutputs,
)
from pymacruise.plus_lesion_inputs import prepare_plus_lesion_inputs
from pymacruise.plus_pipeline import PlusPipelineOptions, run_plus_post_tgdm_stages
from pymacruise.plus_topology_composition import run_plus_topology_and_tgdm
from pymacruise.surfaces import write_legacy_vtk
from pymacruise.tgdm import TGDMMaCRUISEImageResult


@dataclass(frozen=True)
class PlusFileWorkflowResult:
    manifest_path: Path
    outputs: dict[str, Path]
    resumed_stages: tuple[str, ...]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def _write_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def _load(path: Path, name: str) -> nib.spatialimages.SpatialImage:
    if not path.is_file() or path.is_symlink():
        raise FileNotFoundError(f"Missing regular {name} input: {path}")
    if not path.name.lower().endswith(".nii.gz"):
        raise ValueError(f"{name} input must be .nii.gz")
    image = nib.load(str(path))
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)) or len(image.shape) != 3:
        raise ValueError(f"{name} must be a 3-D NIfTI image")
    image = normalize_image_to_ras(image, source_name=name)
    # Validate in slabs: a brain-resolution volume must not be duplicated in
    # RAM just to check finiteness before the stage algorithms allocate their
    # working arrays.
    if not np.issubdtype(image.get_data_dtype(), np.number):
        raise ValueError(f"{name} must contain numeric values")
    for x_start in range(0, image.shape[0], 8):
        values = np.asarray(image.dataobj[x_start : x_start + 8, :, :])
        if not np.isfinite(values).all():
            raise ValueError(f"{name} must contain finite values")
    return image


def _record_valid(
    record: dict[str, Any] | None,
    dependency_signature: str,
    root: Path,
) -> bool:
    if not record or record.get("status") != "complete":
        return False
    if record.get("dependency_signature") != dependency_signature:
        return False
    outputs = record.get("outputs")
    if not isinstance(outputs, dict) or not outputs:
        return False
    for item in outputs.values():
        try:
            path = Path(item["path"])
            if (
                not path.is_relative_to(root)
                or not path.is_file()
                or path.is_symlink()
                or _sha256(path) != item["sha256"]
            ):
                return False
            if path.name.lower().endswith(".nii.gz"):
                image = nib.load(str(path))
                if (
                    not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image))
                    or len(image.shape) not in (3, 4)
                    or nib.aff2axcodes(image.affine) != ("R", "A", "S")
                    or not np.isfinite(np.asarray(image.dataobj)).all()
                ):
                    return False
        except (KeyError, OSError, TypeError, ValueError):
            return False
    return True


def _output_records(paths: dict[str, Path]) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    for name, path in paths.items():
        if not path.is_file() or path.is_symlink():
            raise RuntimeError(f"Stage did not create regular output {name}: {path}")
        if path.name.lower().endswith(".nii.gz"):
            image = nib.load(str(path))
            if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
                raise RuntimeError(f"Stage output is not NIfTI: {path}")
            if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
                raise RuntimeError(f"Stage output is not RAS: {path}")
            if not np.isfinite(np.asarray(image.dataobj)).all():
                raise RuntimeError(f"Stage output contains nonfinite values: {path}")
        records[name] = {"path": str(path.resolve()), "sha256": _sha256(path)}
    return records


def _crop_from_manifest(record: dict[str, Any]) -> CropParameters:
    crop = record["crop_parameters"]
    return CropParameters(
        original_shape=tuple(crop["original_shape"]),
        bounds=tuple(tuple(pair) for pair in crop["bounds"]),
        background=float(crop["background"]),
        original_affine=np.asarray(crop["original_affine"], dtype=np.float64),
    )


def _new_attempt(stage_root: Path) -> Path:
    stage_root.mkdir(parents=True, exist_ok=True)
    existing = [
        int(path.name.split("-")[-1])
        for path in stage_root.glob("attempt-*")
        if path.name.split("-")[-1].isdigit()
    ]
    destination = stage_root / f"attempt-{max(existing, default=0) + 1:02d}"
    destination.mkdir()
    return destination


def run_plus_file_workflow(
    *,
    t1_reference_path: str | Path,
    segmentation_path: str | Path,
    dura_csf_path: str | Path,
    dura_gm_path: str | Path,
    dura_wm_fill_path: str | Path,
    dura_wm_mask_path: str | Path,
    multi_atlas_gm_path: str | Path,
    multi_atlas_wm_path: str | Path,
    multi_atlas_mask_path: str | Path,
    brain_mask_path: str | Path,
    output_dir: str | Path,
    crop_edge_width: int = 5,
    dura_transfer_mode: str = "physical-affine",
    options: PlusPipelineOptions = PlusPipelineOptions(),
    resume: bool = False,
) -> PlusFileWorkflowResult:
    """Run the file-based plus continuation and resume verified stages.

    Inputs are a T1 reference, segmentation, four persisted P39 Dura channels,
    atlas GM/WM/mask, and brain mask. All maps are normalized to RAS and must
    match the segmentation reference after the P39 channels are transferred.
    Atlas/TOADS registration and the full TOADS/LesionTOADS algorithms remain
    upstream inputs; this is not a complete MaCRUISEplus workflow.
    """
    if not isinstance(options, PlusPipelineOptions):
        raise TypeError("options must be PlusPipelineOptions")
    if isinstance(crop_edge_width, bool) or not isinstance(crop_edge_width, int) or crop_edge_width < 0:
        raise ValueError("crop_edge_width must be a nonnegative integer")
    if dura_transfer_mode not in ("physical-affine", "jist-unchanged"):
        raise ValueError("dura_transfer_mode must be 'physical-affine' or 'jist-unchanged'")
    named_paths = {
        "t1_reference": Path(t1_reference_path).expanduser().resolve(),
        "segmentation": Path(segmentation_path).expanduser().resolve(),
        "dura_csf": Path(dura_csf_path).expanduser().resolve(),
        "dura_gm": Path(dura_gm_path).expanduser().resolve(),
        "dura_wm_fill": Path(dura_wm_fill_path).expanduser().resolve(),
        "dura_wm_mask": Path(dura_wm_mask_path).expanduser().resolve(),
        "multi_atlas_gm": Path(multi_atlas_gm_path).expanduser().resolve(),
        "multi_atlas_wm": Path(multi_atlas_wm_path).expanduser().resolve(),
        "multi_atlas_mask": Path(multi_atlas_mask_path).expanduser().resolve(),
        "brain_mask": Path(brain_mask_path).expanduser().resolve(),
    }
    images = {name: _load(path, name.replace("_", " ")) for name, path in named_paths.items()}
    reference = images["t1_reference"]
    for name in ("segmentation", "brain_mask"):
        image = images[name]
        if image.shape != reference.shape or not np.allclose(
            image.affine, reference.affine, rtol=0, atol=1e-4
        ):
            raise ValueError(f"{name} must match the T1 reference grid")
    dura_images = [images[name] for name in ("dura_csf", "dura_gm", "dura_wm_fill", "dura_wm_mask")]
    if dura_transfer_mode == "physical-affine":
        for name in ("dura_csf", "dura_gm", "dura_wm_fill", "dura_wm_mask"):
            image = images[name]
            if image.shape != reference.shape or not np.allclose(image.affine, reference.affine, rtol=0, atol=1e-4):
                raise ValueError(f"{name} must match the T1 grid for physical-affine transfer")
    else:
        source_shape = dura_images[0].shape
        for name, image in zip(("dura_csf", "dura_gm", "dura_wm_fill", "dura_wm_mask"), dura_images, strict=True):
            if image.shape != source_shape or not np.allclose(image.affine, reference.affine, rtol=0, atol=1e-4):
                raise ValueError(
                    f"{name} must share the JIST source extents and carry the T1 reference header"
                )
    settings = {
        "crop_edge_width": crop_edge_width,
        "dura_transfer_mode": dura_transfer_mode,
        "plus_pipeline": options.as_dict(),
    }
    input_records = {
        name: {"path": str(path), "sha256": _sha256(path)}
        for name, path in named_paths.items()
    }
    run_signature = _json_hash({"schema_version": 1, "spatial_convention": "RAS", "inputs": input_records, "settings": settings})
    root = Path(output_dir).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "run.json"
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"Plus run manifest exists; use resume=True: {manifest_path}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("run_signature") != run_signature:
            manifest["previous_run_signature"] = manifest.get("run_signature")
            manifest["run_signature"] = run_signature
            manifest["inputs"] = input_records
            manifest["settings"] = settings
    else:
        if resume:
            raise FileNotFoundError(f"Cannot resume; no plus run manifest exists: {manifest_path}")
        manifest = {
            "schema_version": 1,
            "workflow": "plus downstream from P39 Dura products",
            "status": "running",
            "spatial_convention": "RAS",
            "run_signature": run_signature,
            "inputs": input_records,
            "settings": settings,
            "equivalence_note": "Synthetic file-based composition; full source workflow equivalence is unvalidated.",
            "stages": {},
            "history": [],
        }
        _write_json(manifest_path, manifest)

    # A process may have died during any stage of a previously complete run.
    # Persist the active state before validating/recomputing checkpoints so a
    # failed resume cannot leave a misleading top-level `complete` status.
    manifest["status"] = "running"
    _write_json(manifest_path, manifest)

    resumed: list[str] = []
    stages: dict[str, dict[str, Any]] = manifest["stages"]

    # Stage 1: transfer P39 channels and crop all common-grid companions.
    crop_dependency = _json_hash({"inputs": input_records, "settings": settings})
    crop_record = stages.get("transfer_crop")
    if _record_valid(crop_record, crop_dependency, root):
        resumed.append("transfer_crop")
        crop_outputs = {name: Path(item["path"]) for name, item in crop_record["outputs"].items()}
        crop_parameters = _crop_from_manifest(crop_record)
    else:
        root_attempt = _new_attempt(root / "stages" / "transfer_crop")
        stages["transfer_crop"] = {
            "status": "running",
            "attempt": root_attempt.name,
            "dependency_signature": crop_dependency,
        }
        manifest["history"].append({"stage": "transfer_crop", "status": "running", "attempt": root_attempt.name})
        _write_json(manifest_path, manifest)
        segmentation = images["segmentation"]
        if dura_transfer_mode == "physical-affine":
            channels = np.stack(
                [np.asarray(images[name].dataobj, dtype=np.float32) for name in ("dura_csf", "dura_gm", "dura_wm_fill", "dura_wm_mask")]
            )
            labels = np.asarray(segmentation.dataobj)
            dura_output = LesionToadsPostLoopOutputs(
                classification=np.zeros(labels.shape, dtype=np.int8),
                membership_classification=None,
                lesion_classification=np.zeros(labels.shape, dtype=np.int16),
                result_mode="dura removal inputs",
                products=(LesionToadsOutputProduct("dura_removal_inputs", channels),),
                source_order=("dura_removal_inputs",),
            )
            prepared = prepare_plus_lesion_inputs(
                dura_output,
                reference,
                segmentation,
                images["multi_atlas_gm"],
                images["multi_atlas_wm"],
                images["multi_atlas_mask"],
                images["brain_mask"],
                crop_edge_width=crop_edge_width,
            )
            crop_parameters = prepared.crop_parameters
            crop_images = {
                "multi_atlas_gm": prepared.multi_atlas_gm,
                "multi_atlas_wm": prepared.multi_atlas_wm,
                "multi_atlas_mask": prepared.multi_atlas_mask,
                "brain_mask": prepared.brain_mask,
                "segmentation": prepared.segmentation,
                "toads_gm": prepared.toads_gm,
                "toads_wm_fill": prepared.toads_wm_fill,
                "toads_wm_mask": prepared.toads_wm_mask,
            }
        else:
            crop_parameters = find_crop_parameters(
                images["multi_atlas_gm"], threshold=0.0, edge_width=crop_edge_width
            )
            crop_images = {
                name: crop_image(image, crop_parameters)
                for name, image in (
                    ("multi_atlas_gm", images["multi_atlas_gm"]),
                    ("multi_atlas_wm", images["multi_atlas_wm"]),
                    ("multi_atlas_mask", images["multi_atlas_mask"]),
                    ("brain_mask", images["brain_mask"]),
                    ("segmentation", segmentation),
                )
            }
            crop_images.update({
                "toads_gm": crop_image_with_reference_bounds(images["dura_gm"], crop_parameters),
                "toads_wm_fill": crop_image_with_reference_bounds(images["dura_wm_fill"], crop_parameters),
                "toads_wm_mask": crop_image_with_reference_bounds(images["dura_wm_mask"], crop_parameters),
            })
        crop_outputs = {}
        for name, image in crop_images.items():
            path = root_attempt / f"{name}.nii.gz"
            save_nifti_ras(image, path)
            crop_outputs[name] = path
        crop_json = root_attempt / "crop.json"
        _write_json(crop_json, {
            "original_shape": list(crop_parameters.original_shape),
            "bounds": [list(pair) for pair in crop_parameters.bounds],
            "background": crop_parameters.background,
            "original_affine": crop_parameters.original_affine.tolist(),
        })
        crop_outputs["crop_parameters"] = crop_json
        crop_record = {
            "status": "complete",
            "attempt": root_attempt.name,
            "dependency_signature": crop_dependency,
            "crop_parameters": {
                "original_shape": list(crop_parameters.original_shape),
                "bounds": [list(pair) for pair in crop_parameters.bounds],
                "background": crop_parameters.background,
                "original_affine": crop_parameters.original_affine.tolist(),
            },
            "outputs": _output_records(crop_outputs),
        }
        stages["transfer_crop"] = crop_record
        manifest["history"].append({"stage": "transfer_crop", "status": "complete", "attempt": root_attempt.name})
        _write_json(manifest_path, manifest)

    # Stage 2: selected plus topology, ACE, GVF, and nested TGDM.
    crop_sig = _json_hash(crop_record["outputs"])
    tgdm_dependency = _json_hash({"crop": crop_sig, "settings": settings})
    tgdm_record = stages.get("topology_ace_gvf_tgdm")
    if _record_valid(tgdm_record, tgdm_dependency, root):
        resumed.append("topology_ace_gvf_tgdm")
        tgdm_outputs = {name: Path(item["path"]) for name, item in tgdm_record["outputs"].items()}
        tgdm_images = {name: _load(path, name) for name, path in tgdm_outputs.items()}
    else:
        root_attempt = _new_attempt(root / "stages" / "topology_ace_gvf_tgdm")
        stages["topology_ace_gvf_tgdm"] = {
            "status": "running",
            "attempt": root_attempt.name,
            "dependency_signature": tgdm_dependency,
        }
        manifest["history"].append({"stage": "topology_ace_gvf_tgdm", "status": "running", "attempt": root_attempt.name})
        _write_json(manifest_path, manifest)
        cropped = {name: _load(Path(crop_record["outputs"][name]["path"]), name) for name in (
            "segmentation", "toads_wm_fill", "multi_atlas_gm", "multi_atlas_wm", "multi_atlas_mask", "toads_gm", "toads_wm_mask"
        )}
        composed = run_plus_topology_and_tgdm(
            cropped["segmentation"],
            cropped["toads_wm_mask"],
            cropped["multi_atlas_gm"],
            cropped["multi_atlas_wm"],
            cropped["multi_atlas_mask"],
            cropped["toads_gm"],
            cropped["toads_wm_fill"],
            cropped["toads_wm_mask"],
            options=options,
        )
        result = composed.tgdm
        tgdm_images = {name: getattr(result, name) for name in (
            "initial_level_set", "tissue_sign", "initial_inner", "inner", "central", "outer"
        )}
        tgdm_outputs = {}
        for name, image in tgdm_images.items():
            path = root_attempt / f"{name}.nii.gz"
            save_nifti_ras(image, path)
            tgdm_outputs[name] = path
        tgdm_record = {
            "status": "complete",
            "attempt": root_attempt.name,
            "dependency_signature": tgdm_dependency,
            "outputs": _output_records(tgdm_outputs),
        }
        stages["topology_ace_gvf_tgdm"] = tgdm_record
        manifest["history"].append({"stage": "topology_ace_gvf_tgdm", "status": "complete", "attempt": root_attempt.name})
        _write_json(manifest_path, manifest)

    # Stage 3: restore full grid, extract meshes, and calculate thickness.
    tail_dependency = _json_hash({
        "tgdm": _json_hash(tgdm_record["outputs"]),
        "crop": crop_sig,
        # MaTGDMcrop uses a 20-valued exterior when restoring level sets.
        # Include this contract so old zero-background meshes are regenerated.
        "uncrop_background": 20.0,
    })
    tail_record = stages.get("surfaces_thickness")
    if _record_valid(tail_record, tail_dependency, root):
        resumed.append("surfaces_thickness")
        tail_outputs = {name: Path(item["path"]) for name, item in tail_record["outputs"].items()}
    else:
        root_attempt = _new_attempt(root / "stages" / "surfaces_thickness")
        stages["surfaces_thickness"] = {
            "status": "running",
            "attempt": root_attempt.name,
            "dependency_signature": tail_dependency,
        }
        manifest["history"].append({"stage": "surfaces_thickness", "status": "running", "attempt": root_attempt.name})
        _write_json(manifest_path, manifest)
        images_for_tail = {name: _load(Path(tgdm_record["outputs"][name]["path"]), name) for name in (
            "initial_level_set", "tissue_sign", "initial_inner", "inner", "central", "outer"
        )}
        tgdm = TGDMMaCRUISEImageResult(**images_for_tail, execution=None)  # execution is unused by the post-TGDM adapters
        post = run_plus_post_tgdm_stages(tgdm, crop_parameters)
        tail_images = {
            "inner_level_set": post.level_sets.inner,
            "central_level_set": post.level_sets.central,
            "outer_level_set": post.level_sets.outer,
            "thickness_volume": post.thickness.volume.image,
        }
        tail_outputs = {}
        for name, image in tail_images.items():
            path = root_attempt / f"{name}.nii.gz"
            save_nifti_ras(image, path)
            tail_outputs[name] = path
        for name, mesh in (
            ("inner_surface", post.surfaces.inner),
            ("central_surface", post.surfaces.central),
            ("outer_surface", post.surfaces.outer),
        ):
            path = root_attempt / f"{name}.vtk"
            write_legacy_vtk(mesh, path)
            tail_outputs[name] = path
        metrics = root_attempt / "thickness.json"
        def finite_metric(value: float) -> float | None:
            """Represent unavailable empty-region statistics as JSON null."""
            numeric = float(value)
            return numeric if np.isfinite(numeric) else None

        crossing_voxels = int(post.thickness.volume.crossing_voxel_count)
        _write_json(metrics, {
            "spatial_convention": "RAS",
            "surface_units": "mm",
            "volume": {
                "mean_mm": finite_metric(post.thickness.volume.mean_mm),
                "standard_deviation_mm": finite_metric(post.thickness.volume.standard_deviation_mm),
                "minimum_mm": finite_metric(post.thickness.volume.minimum_mm),
                "maximum_mm": finite_metric(post.thickness.volume.maximum_mm),
                "crossing_voxel_count": crossing_voxels,
                "available": crossing_voxels > 0,
            },
            "surface": {
                "mean_mm": finite_metric(post.thickness.surface.mean_mm),
                "standard_deviation_mm": finite_metric(post.thickness.surface.standard_deviation_mm),
                "nonzero_vertex_count": post.thickness.surface.nonzero_vertex_count,
            },
        })
        tail_outputs["thickness_metrics"] = metrics
        tail_record = {
            "status": "complete",
            "attempt": root_attempt.name,
            "dependency_signature": tail_dependency,
            "outputs": _output_records(tail_outputs),
        }
        stages["surfaces_thickness"] = tail_record
        manifest["history"].append({"stage": "surfaces_thickness", "status": "complete", "attempt": root_attempt.name})
        _write_json(manifest_path, manifest)

    manifest["status"] = "complete"
    manifest["workflow_scope"] = "synthetic or supplied prealigned P39-to-plus continuation; upstream TOADS and LesionTOADS are external inputs"
    _write_json(manifest_path, manifest)
    return PlusFileWorkflowResult(manifest_path, tail_outputs, tuple(resumed))
