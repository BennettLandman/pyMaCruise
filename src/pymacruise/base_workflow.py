"""Restartable synthetic/native base MaCRUISE workflow orchestration.

This module composes existing stage CLIs. It adds checkpointing, provenance,
and stage contracts; it does not replace any image-processing algorithm.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from contextlib import redirect_stderr, redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
from typing import Any
import uuid

import nibabel as nib
import numpy as np

from pymacruise.cli import main as cli_main
from pymacruise.image_io import normalize_image_to_ras


@dataclass(frozen=True)
class BaseRunOptions:
    registration_levels: int = 1
    registration_iterations: int = 1
    registration_max_iterations: int = 1
    toads_iterations: int = 2
    dura_iterations: int = 1
    crop_edge_width: int = 2
    ace_diffusion_iterations: int = 1
    gvf_iterations: int = 200
    tgdm_iterations: int = 1


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


def _write_manifest(path: Path, manifest: dict[str, Any]) -> None:
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n")
    tmp.replace(path)


def _output_records(
    paths: list[Path],
    reference: nib.spatialimages.SpatialImage,
    *,
    require_reference_grid: bool = True,
) -> list[dict[str, Any]]:
    records = []
    for path in paths:
        if not path.is_file() or path.is_symlink():
            raise RuntimeError(f"Stage did not create a regular output file: {path}")
        if path.name.lower().endswith(".nii.gz"):
            image = nib.load(str(path))
            if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
                raise RuntimeError(f"Stage output is not NIfTI: {path}")
            if len(image.shape) not in (3, 4):
                raise RuntimeError(f"Stage output is not 3D or component-last 4D NIfTI: {path}")
            if require_reference_grid and (
                image.shape != reference.shape
                or not np.allclose(image.affine, reference.affine, rtol=0, atol=1e-4)
            ):
                raise RuntimeError(f"Stage output is not on the native grid: {path}")
            if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
                raise RuntimeError(f"Stage output is not RAS: {path}")
            values = np.asanyarray(image.dataobj)
            if not np.issubdtype(values.dtype, np.number) or not np.isfinite(values).all():
                raise RuntimeError(f"Stage output has nonfinite or nonnumeric voxels: {path}")
        record = {"path": str(path.resolve()), "sha256": _sha256(path), "bytes": path.stat().st_size}
        if path.name.lower().endswith(".nii.gz"):
            record["shape"] = [int(size) for size in image.shape]
            record["axis_codes"] = list(nib.aff2axcodes(image.affine))
        records.append(record)
    return records


def _stage_is_valid(record: dict[str, Any] | None, dependency_signature: str | None = None) -> bool:
    if not record or record.get("status") != "complete" or not record.get("outputs"):
        return False
    if dependency_signature is not None and record.get("dependency_signature") != dependency_signature:
        return False
    for output in record["outputs"]:
        path = Path(output["path"])
        if not path.is_file() or path.is_symlink() or _sha256(path) != output["sha256"]:
            return False
    return True


def _toads_telemetry(text: str) -> dict[str, Any] | None:
    match = re.search(
        r"iterations=(\d+), .*?terminal_energy_ratio=([-+0-9.eE]+), "
        r"energy_threshold_reached=(true|false), iteration_limit_reached=(true|false)",
        text,
    )
    if match is None:
        return None
    return {
        "iterations": int(match.group(1)),
        "terminal_energy_ratio": float(match.group(2)),
        "energy_threshold_reached": match.group(3) == "true",
        "iteration_limit_reached": match.group(4) == "true",
    }


def _materialize_dura_atlas(source_descriptor: Path, destination: Path) -> Path:
    """Copy a descriptor and its image assets, applying wrapper tissue names.

    The TOADS source itself does not need this rename; the Dura/CRUISE output
    adapter classifies tissue products by these names. Copying keeps the run
    self-contained and leaves the supplied atlas untouched.
    """
    source_root = source_descriptor.parent.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    text = source_descriptor.read_text(encoding="utf-8")
    mapping = {
        "Csf1": "SulcalCSF", "Csf2": "Ventricles",
        "Gm1": "CerebellarGM", "Gm2": "CerebralGM", "Gm3": "Caudate",
        "Gm4": "Thalamus", "Gm5": "Putamen",
        "Wm1": "Brainstem", "Wm2": "CerebellarWM", "Wm3": "CerebralWM",
    }
    source_names = set(mapping)
    output_names = set(mapping.values())
    present = {
        line.partition(":")[2].strip()
        for line in text.splitlines()
        if line.startswith("Structure:")
    }
    if source_names <= present:
        for old, new in mapping.items():
            text = text.replace(old, new)
    elif not output_names <= present:
        missing = sorted(output_names - present)
        raise ValueError(
            "Base atlas must provide the selected ten Dura tissue structures; "
            "missing " + ", ".join(missing)
        )
    for line in source_descriptor.read_text(encoding="utf-8").splitlines():
        if not line.startswith("File:"):
            continue
        relative = Path(line.partition(":")[2].strip())
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe atlas asset path in descriptor: {relative}")
        source = (source_root / relative).resolve()
        if not source.is_relative_to(source_root) or not source.is_file():
            raise FileNotFoundError(f"Missing or unsafe atlas asset: {source}")
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            if _sha256(target) != _sha256(source):
                raise ValueError(f"Existing run atlas resource differs from input: {target}")
        else:
            shutil.copy2(source, target)
    target_descriptor = destination / source_descriptor.name
    desired_bytes = text.encode("utf-8")
    if target_descriptor.exists():
        if target_descriptor.read_bytes() != desired_bytes:
            raise ValueError(f"Existing derived atlas descriptor differs from expected content: {target_descriptor}")
    else:
        tmp = target_descriptor.with_suffix(target_descriptor.suffix + ".tmp")
        tmp.write_bytes(desired_bytes)
        tmp.replace(target_descriptor)
    return target_descriptor


def run_base_workflow(
    *,
    t1_path: str | Path,
    brain_mask_path: str | Path,
    atlas_descriptor: str | Path,
    maseg_path: str | Path,
    output_dir: str | Path,
    options: BaseRunOptions = BaseRunOptions(),
    resume: bool = False,
) -> Path:
    """Run N4 → TOADS → Dura → TOADS → optional MaSeg cleanup → CRUISE.

    A stage is skipped only when its complete output set still matches the
    recorded SHA-256 values and the run input/configuration signature matches.
    A failed/incomplete attempt is preserved; retry uses a fresh attempt
    directory, making interruption recovery non-destructive.
    """
    t1_path = Path(t1_path).expanduser().resolve()
    brain_mask_path = Path(brain_mask_path).expanduser().resolve()
    atlas_descriptor = Path(atlas_descriptor).expanduser().resolve()
    output_dir = Path(output_dir).expanduser().resolve()
    maseg_path = Path(maseg_path).expanduser().resolve()
    inputs = {"t1": t1_path, "brain_mask": brain_mask_path, "atlas_descriptor": atlas_descriptor}
    inputs["maseg"] = maseg_path
    for role, path in inputs.items():
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"Missing regular {role} input: {path}")
    if not t1_path.name.lower().endswith(".nii.gz") or not brain_mask_path.name.lower().endswith(".nii.gz"):
        raise ValueError("T1 and brain-mask inputs must be .nii.gz")
    if not maseg_path.name.lower().endswith(".nii.gz"):
        raise ValueError("MaSeg input must be .nii.gz")
    reference = normalize_image_to_ras(nib.load(str(t1_path)), source_name="base workflow T1")
    mask = normalize_image_to_ras(nib.load(str(brain_mask_path)), source_name="base workflow brain mask")
    if reference.shape != mask.shape or not np.allclose(reference.affine, mask.affine, rtol=0, atol=1e-4):
        raise ValueError("Brain mask must match the input T1 grid")
    maseg_image = normalize_image_to_ras(nib.load(str(maseg_path)), source_name="base workflow MaSeg")
    maseg_values = np.asanyarray(maseg_image.dataobj)
    if maseg_image.shape != reference.shape or not np.allclose(maseg_image.affine, reference.affine, rtol=0, atol=1e-4):
        raise ValueError("MaSeg input must match the input T1 grid")
    if not np.issubdtype(maseg_values.dtype, np.integer) or not np.isfinite(maseg_values).all():
        raise ValueError("MaSeg input must contain finite integer labels")
    if not reference.shape or len(reference.shape) != 3:
        raise ValueError("Input T1 must be a 3D image")
    atlas_assets: list[Path] = []
    atlas_root = atlas_descriptor.parent.resolve()
    for line in atlas_descriptor.read_text(encoding="utf-8").splitlines():
        if line.startswith("File:"):
            relative = Path(line.partition(":")[2].strip())
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError(f"Unsafe atlas asset path in descriptor: {relative}")
            asset = (atlas_root / relative).resolve()
            if not asset.is_relative_to(atlas_root) or not asset.is_file():
                raise FileNotFoundError(f"Missing or unsafe atlas asset: {asset}")
            atlas_assets.append(asset)
    if not atlas_assets:
        raise ValueError("Atlas descriptor does not reference any asset files")

    settings = asdict(options)
    if any(isinstance(value, bool) or not isinstance(value, int) or value < 1 for value in settings.values()):
        raise ValueError("All synthetic base workflow iteration/count options must be positive integers")
    signature_payload = {
        "schema_version": 1,
        "spatial_convention": "RAS",
        "inputs": {name: {"path": str(path), "sha256": _sha256(path)} for name, path in inputs.items()},
        "atlas_assets": [{"path": str(path), "sha256": _sha256(path)} for path in atlas_assets],
        "options": settings,
    }
    signature = _json_hash(signature_payload)
    output_dir.mkdir(parents=True, exist_ok=True)
    run_atlas = _materialize_dura_atlas(atlas_descriptor, output_dir / "resources" / "dura-atlas")
    manifest_path = output_dir / "run.json"
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"Run manifest exists; use --resume to continue: {manifest_path}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("run_signature") != signature:
            raise ValueError("Cannot resume: inputs or run settings differ from recorded manifest")
    else:
        if resume:
            raise FileNotFoundError(f"Cannot resume; no run manifest exists: {manifest_path}")
        manifest = {
            **signature_payload,
            "run_signature": signature,
            "status": "running",
            "stages": {},
            "history": [],
            "workflow_scope": "synthetic/native-grid base route; not a full Java-equivalence claim",
            "derived_resources": [{"path": str(run_atlas), "sha256": _sha256(run_atlas)}],
        }
        _write_manifest(manifest_path, manifest)

    # Use a new immutable attempt directory whenever a prior attempt did not
    # finish cleanly. This preserves partial products for diagnosis.
    def run_stage(
        name: str,
        outputs: list[Path],
        command: list[str],
        dependencies: list[Path],
    ) -> None:
        dependency_signature = _json_hash(
            [{"path": str(path.resolve()), "sha256": _sha256(path)} for path in dependencies]
        )
        previous = manifest["stages"].get(name)
        if _stage_is_valid(previous, dependency_signature):
            print(f"Resume: verified {name}; skipping")
            return
        attempts = 1 + sum(item.get("stage") == name for item in manifest["history"])
        stage_root = output_dir / "stages" / name
        planned_stage_dir = stage_root / f"attempt-{attempts:02d}"
        stage_dir = planned_stage_dir
        while stage_dir.exists():
            # A process may have stopped after creating its attempt directory
            # but before recording a completed stage in run.json. Preserve that
            # partial attempt and choose a fresh directory on restart.
            attempts += 1
            stage_dir = output_dir / "stages" / name / f"attempt-{attempts:02d}"
        stage_dir.mkdir(parents=True, exist_ok=False)
        attempt_outputs = [
            stage_dir / path.relative_to(planned_stage_dir)
            if path.is_relative_to(planned_stage_dir)
            else path
            for path in outputs
        ]
        stdout_buffer, stderr_buffer = io.StringIO(), io.StringIO()
        log_path = stage_dir / "stage.log"
        try:
            with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
                cli_main(command_for_dir(command, stage_dir))
            log_path.write_text(
                stdout_buffer.getvalue()
                + ("\n[stderr]\n" if stderr_buffer.getvalue() else "")
                + stderr_buffer.getvalue(),
                encoding="utf-8",
            )
            print(stdout_buffer.getvalue(), end="")
            records = _output_records(
                [*attempt_outputs, log_path], reference,
                require_reference_grid=name != "cruise-core",
            )
            if name == "cruise-core":
                by_name = {Path(item["path"]).name: Path(item["path"]) for item in records}
                for filename in ("inner.nii.gz", "central.nii.gz", "outer.nii.gz"):
                    image = nib.load(str(by_name[filename]))
                    if image.shape != reference.shape or not np.allclose(image.affine, reference.affine, rtol=0, atol=1e-4):
                        raise RuntimeError(f"CRUISE native output is not on the T1 grid: {filename}")
                crop_names = (
                    "gm_cropped.nii.gz", "wmfill_cropped.nii.gz", "wmmask_cropped.nii.gz",
                    "wmfill_topology_corrected.nii.gz", "ace_enhanced_gm.nii.gz",
                    "ace_skeleton.nii.gz", "ace_thinned_skeleton.nii.gz",
                )
                crop_images = [nib.load(str(by_name[filename])) for filename in crop_names]
                if any(image.shape != crop_images[0].shape or not np.allclose(image.affine, crop_images[0].affine, rtol=0, atol=1e-4) for image in crop_images[1:]):
                    raise RuntimeError("CRUISE cropped tissue/ACE outputs do not share one grid")
                gvf = nib.load(str(by_name["gvf.nii.gz"]))
                if gvf.shape != (*crop_images[0].shape, 3) or not np.allclose(gvf.affine, crop_images[0].affine, rtol=0, atol=1e-4):
                    raise RuntimeError("CRUISE GVF output does not share the cropped grid")
            stage_record: dict[str, Any] = {
                "status": "complete", "attempt": attempts,
                "outputs": records,
                "dependency_signature": dependency_signature,
                "command": command_for_dir(command, stage_dir),
            }
            if name in {"toads-first", "toads-second"}:
                telemetry = _toads_telemetry(stdout_buffer.getvalue())
                if telemetry is None:
                    raise RuntimeError(f"Could not record selected TOADS stop telemetry for {name}")
                stage_record["toads_stop"] = telemetry
                stage_record["converged"] = telemetry["energy_threshold_reached"]
            manifest["stages"][name] = stage_record
            manifest["history"].append({"stage": name, "attempt": attempts, "status": "complete"})
            _write_manifest(manifest_path, manifest)
            print(f"Stage check passed: {name} ({len(records)} outputs)")
        except Exception as exc:
            log_path.write_text(
                stdout_buffer.getvalue()
                + ("\n[stderr]\n" if stderr_buffer.getvalue() else "")
                + stderr_buffer.getvalue()
                + f"\n[exception] {type(exc).__name__}: {exc}\n",
                encoding="utf-8",
            )
            print(stdout_buffer.getvalue(), end="")
            manifest["stages"][name] = {"status": "failed", "attempt": attempts, "error": f"{type(exc).__name__}: {exc}"}
            manifest["history"].append({"stage": name, "attempt": attempts, "status": "failed"})
            _write_manifest(manifest_path, manifest)
            raise

    # Command output roots are fixed to each attempt directory. Resolve the
    # latest valid attempt by pointing the stable aliases at its products.
    def command_for_dir(command: list[str], stage_dir: Path) -> list[str]:
        return [arg.replace("{STAGE}", str(stage_dir)) for arg in command]

    def active_output(
        stage_name: str, relative: str, dependencies: list[Path] | None = None
    ) -> Path:
        record = manifest["stages"].get(stage_name)
        expected = None
        if dependencies is not None:
            expected = _json_hash(
                [{"path": str(path.resolve()), "sha256": _sha256(path)} for path in dependencies]
            )
        if _stage_is_valid(record, expected):
            match = next((Path(item["path"]) for item in record["outputs"] if Path(item["path"]).name == Path(relative).name), None)
            if match is not None:
                return match
        # current attempt path is predictable before it runs
        count = 1 + sum(item.get("stage") == stage_name for item in manifest["history"])
        return output_dir / "stages" / stage_name / f"attempt-{count:02d}" / relative

    n4 = active_output("n4", "t1-n4.nii.gz", [t1_path])
    run_stage("n4", [n4], ["n4", "--input", str(t1_path), "--output", "{STAGE}/t1-n4.nii.gz"], [t1_path])
    n4 = active_output("n4", "t1-n4.nii.gz", [t1_path])

    toads_common = [
        "--atlas-descriptor", str(run_atlas),
        "--registration-levels", str(options.registration_levels),
        "--registration-iterations", str(options.registration_iterations),
        "--registration-max-iterations", str(options.registration_max_iterations),
        "--maximum-iterations", str(options.toads_iterations),
    ]
    atlas_dependencies = [run_atlas, *(run_atlas.parent / asset.relative_to(atlas_root) for asset in atlas_assets)]
    first_dependencies = [n4, brain_mask_path, *atlas_dependencies]
    first_dura_files = [
        active_output("toads-first", "classification.nii.gz", first_dependencies),
        active_output("toads-first", "dura-inputs/t1-n4_wmfill.nii.gz", first_dependencies),
        active_output("toads-first", "dura-inputs/t1-n4_gm.nii.gz", first_dependencies),
        active_output("toads-first", "dura-inputs/t1-n4_wmmask.nii.gz", first_dependencies),
        active_output("toads-first", "dura-inputs/t1-n4_csf.nii.gz", first_dependencies),
    ]
    run_stage("toads-first", first_dura_files, [
        "toads-loop", "--t1-n4", str(n4), "--brain-mask", str(brain_mask_path),
        "--output", "{STAGE}/classification.nii.gz", "--dura-output-dir", "{STAGE}/dura-inputs",
        *toads_common,
    ], [n4, brain_mask_path, *atlas_dependencies])
    first_label = active_output("toads-first", "classification.nii.gz")
    first_outputs = manifest["stages"]["toads-first"]["outputs"]
    first_files = {Path(item["path"]).name: Path(item["path"]) for item in first_outputs}
    first_wm = next(path for name, path in first_files.items() if name.endswith("_wmfill.nii.gz"))
    first_gm = next(path for name, path in first_files.items() if name.endswith("_gm.nii.gz"))
    first_wmmask = next(path for name, path in first_files.items() if name.endswith("_wmmask.nii.gz"))
    dura_dependencies = [n4, first_wm, first_gm, first_wmmask]
    dura_products = [active_output("dura", f"{name}.nii.gz", dura_dependencies) for name in (
        "stripped_t1", "corrected_white_matter", "corrected_gray_matter", "brain_mask_level_set", "brain_mask", "wm_topology", "gm_topology"
    )]
    run_stage("dura", dura_products, [
        "dura", "--t1", str(n4), "--white", str(first_wm), "--gray", str(first_gm),
        "--white-mask", str(first_wmmask), "--output-dir", "{STAGE}",
        "--max-iterations", str(options.dura_iterations),
    ], [n4, first_wm, first_gm, first_wmmask])
    stripped = active_output("dura", "stripped_t1.nii.gz", dura_dependencies)
    brain = active_output("dura", "brain_mask.nii.gz", dura_dependencies)
    second_dependencies = [stripped, brain, *atlas_dependencies]
    second_outputs = [
        active_output("toads-second", "classification.nii.gz", second_dependencies),
        active_output("toads-second", "dura-inputs/stripped_t1_wmfill.nii.gz", second_dependencies),
        active_output("toads-second", "dura-inputs/stripped_t1_gm.nii.gz", second_dependencies),
        active_output("toads-second", "dura-inputs/stripped_t1_wmmask.nii.gz", second_dependencies),
        active_output("toads-second", "dura-inputs/stripped_t1_csf.nii.gz", second_dependencies),
    ]
    run_stage("toads-second", second_outputs, [
        "toads-loop", "--t1-n4", str(stripped), "--brain-mask", str(brain),
        "--output", "{STAGE}/classification.nii.gz", "--dura-output-dir", "{STAGE}/dura-inputs",
        *toads_common,
    ], [stripped, brain, *atlas_dependencies])
    second_outputs = manifest["stages"]["toads-second"]["outputs"]
    second_files = {Path(item["path"]).name: Path(item["path"]) for item in second_outputs}
    second_gm = next(path for name, path in second_files.items() if name.endswith("_gm.nii.gz"))
    second_wm = next(path for name, path in second_files.items() if name.endswith("_wmfill.nii.gz"))
    second_wmmask = next(path for name, path in second_files.items() if name.endswith("_wmmask.nii.gz"))
    maseg_args = ["--segmentation", str(maseg_path), "--remove-cerebellum"]
    core_dependencies = [second_gm, second_wm, second_wmmask]
    core_dependencies.append(maseg_path)
    core_outputs = [active_output("cruise-core", f"{name}.nii.gz", core_dependencies) for name in (
        "gm_cropped", "wmfill_cropped", "wmmask_cropped", "wmfill_topology_corrected",
        "ace_enhanced_gm", "ace_skeleton", "ace_thinned_skeleton", "gvf", "inner", "central", "outer"
    )]
    core_outputs.append(active_output("cruise-core", "provenance.json", core_dependencies))
    run_stage("cruise-core", core_outputs, [
        "toads-cruise-core", "--gm", str(second_gm), "--wm-fill", str(second_wm), "--wm-mask", str(second_wmmask),
        *maseg_args, "--edge-width", str(options.crop_edge_width),
        "--ace-diffusion-iterations", str(options.ace_diffusion_iterations),
        "--gvf-iterations", str(options.gvf_iterations),
        "--tgdm-initial-iterations", str(options.tgdm_iterations),
        "--tgdm-inner-iterations", str(options.tgdm_iterations),
        "--tgdm-central-iterations", str(options.tgdm_iterations),
        "--output-dir", "{STAGE}",
        "--provenance", "{STAGE}/provenance.json",
    ], core_dependencies)
    core_files = {Path(item["path"]).name: Path(item["path"]) for item in manifest["stages"]["cruise-core"]["outputs"]}
    refine_dir = output_dir / "stages" / "maseg-refinement" / "attempt-01"
    refined_output = active_output("maseg-refinement", "SEG_refine.nii.gz", [
        maseg_path, core_files["inner.nii.gz"], core_files["outer.nii.gz"]
    ])
    refine_summary = refined_output.with_name("SEG_refine-summary.json")
    run_stage("maseg-refinement", [refined_output, refine_summary], [
        "segmentation-refine", "--segmentation", str(maseg_path),
        "--inner-levelset", str(core_files["inner.nii.gz"]),
        "--outer-levelset", str(core_files["outer.nii.gz"]),
        "--output", "{STAGE}/SEG_refine.nii.gz",
        "--summary", "{STAGE}/SEG_refine-summary.json",
    ], [maseg_path, core_files["inner.nii.gz"], core_files["outer.nii.gz"]])
    manifest["status"] = "complete"
    manifest["converged"] = all(
        manifest["stages"][stage]["toads_stop"]["energy_threshold_reached"]
        for stage in ("toads-first", "toads-second")
    )
    manifest["result_interpretation"] = (
        "base-stage-run-converged"
        if manifest["converged"]
        else "stage-run-complete-but-toads-did-not-meet-energy-threshold"
    )
    manifest["source_equivalence_claim"] = False
    manifest["final_outputs"] = (
        manifest["stages"]["cruise-core"]["outputs"]
        + manifest["stages"]["maseg-refinement"]["outputs"]
    )
    _write_manifest(manifest_path, manifest)
    return manifest_path
