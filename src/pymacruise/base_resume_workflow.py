"""Restartable MaCRUISE base tail starting from a completed TOADS second pass.

The source run is read only. Its completed ``toads-second`` products are
validated and recorded as immutable upstream inputs in a distinct manifest;
MaSeg preparation, CRUISE core, source-faithful refinement, and surfaces then
run as independently hash-checked stages.
"""

from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import re
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.base_workflow import (
    _implementation_signature,
    _json_hash,
    _output_records,
    _sha256,
    _stage_is_valid,
    _write_manifest,
)
from pymacruise.cli import main as cli_main
from pymacruise.contracts import has_finite_integer_values
from pymacruise.image_io import normalize_image_to_ras
from pymacruise.surfaces import read_legacy_vtk


_V3_CROP_EDGE_WIDTH = 5


def _selected_v3_profile_record(
    manifest: dict[str, Any],
    stages: dict[str, Any],
    source_options: dict[str, Any],
    toads_command: list[str],
    dura_command: list[str],
    dura_membership_level: float,
) -> dict[str, Any]:
    """Describe whether recorded upstream settings match the selected V3 profile.

    This is deliberately informational: importing a complete, hash-verified
    checkpoint with shorter settings remains supported. Missing wrapper
    defaults are called out as unrecorded instead of being inferred from a
    command line that omitted them.
    """
    expected_options: dict[str, int] = {
        "registration_levels": 4,
        "registration_iterations": 50,
        "post_refresh_registration_iterations": 5,
        "registration_max_iterations": 5,
        "toads_iterations": 50,
        "dura_iterations": 5,
        "crop_edge_width": _V3_CROP_EDGE_WIDTH,
    }
    checked: dict[str, dict[str, Any]] = {}
    mismatches: list[dict[str, Any]] = []

    def check(name: str, expected: Any, actual: Any, evidence: str, *, numeric: bool = False) -> None:
        matches = bool(np.isclose(float(actual), float(expected), rtol=0.0, atol=1e-12)) if numeric else actual == expected
        checked[name] = {
            "expected": expected,
            "actual": actual,
            "matches": matches,
            "evidence": evidence,
        }
        if not matches:
            mismatches.append({
                "setting": name,
                "expected": expected,
                "actual": actual,
                "evidence": evidence,
            })

    for name, expected in expected_options.items():
        check(name, expected, source_options.get(name), f"run.json options.{name}")
    check(
        "dura_command_max_iterations", 5,
        int(dura_command[dura_command.index("--max-iterations") + 1]),
        "Dura command --max-iterations",
    )
    check(
        "dura_membership_level", 0.3, dura_membership_level,
        "Dura command --membership-level (effective value)", numeric=True,
    )

    # The generated TOADS command makes these wrapper settings explicit. The
    # helper is only called after _validate_source_run has verified presence
    # and consistency with run.json options.
    command_option_names = {
        "--registration-levels": "registration_levels",
        "--registration-iterations": "registration_iterations",
        "--post-refresh-registration-iterations": "post_refresh_registration_iterations",
        "--registration-max-iterations": "registration_max_iterations",
        "--maximum-iterations": "toads_iterations",
    }
    for flag, option_name in command_option_names.items():
        check(
            f"toads_{option_name}", expected_options[option_name],
            int(toads_command[toads_command.index(flag) + 1]),
            f"toads-second command {flag}",
        )
    if "--maximum-distance" in toads_command:
        check(
            "toads_maximum_distance", 1e-4,
            float(toads_command[toads_command.index("--maximum-distance") + 1]),
            "toads-second command --maximum-distance", numeric=True,
        )
    else:
        checked["toads_maximum_distance"] = {
            "expected": 1e-4,
            "actual": None,
            "matches": None,
            "evidence": "omitted from command; CLI default is not persisted in the source run manifest",
        }
    if "--convergence-threshold" in toads_command:
        check(
            "toads_convergence_threshold", 1e-4,
            float(toads_command[toads_command.index("--convergence-threshold") + 1]),
            "toads-second command --convergence-threshold", numeric=True,
        )
    else:
        checked["toads_convergence_threshold"] = {
            "expected": 1e-4,
            "actual": None,
            "matches": None,
            "evidence": "omitted from command; CLI default is not persisted in the source run manifest",
        }
    if "--modality" in toads_command:
        check(
            "toads_modality", "T1_MPRAGE",
            toads_command[toads_command.index("--modality") + 1],
            "toads-second command --modality",
        )
    else:
        checked["toads_modality"] = {
            "expected": "T1_MPRAGE",
            "actual": None,
            "matches": None,
            "evidence": "omitted from command; CLI default is not persisted in the source run manifest",
        }

    n4_stage = stages.get("n4", {})
    n4_command = n4_stage.get("command") if isinstance(n4_stage, dict) else None
    if isinstance(n4_command, list) and "--input" in n4_command:
        n4_input = Path(n4_command[n4_command.index("--input") + 1]).resolve()
        declared_t1 = Path(manifest["inputs"]["t1"]["path"]).resolve()
        check(
            "n4_input_matches_declared_t1", str(declared_t1), str(n4_input),
            "N4 command --input compared with run.json inputs.t1.path",
        )
    else:
        checked["n4_input_matches_declared_t1"] = {
            "expected": "N4 runs from the T1 declared in run.json",
            "actual": None,
            "matches": None,
            "evidence": "N4 input command is not recorded in the source run manifest",
        }

    expected_stage_order = ["n4", "toads-first", "dura", "toads-second"]
    history = manifest.get("history")
    if isinstance(history, list):
        latest_events: dict[str, tuple[int, dict[str, Any]]] = {}
        for index, item in enumerate(history):
            if (
                isinstance(item, dict)
                and item.get("status") == "complete"
                and item.get("stage") in expected_stage_order
            ):
                latest_events[item["stage"]] = (index, item)
        recorded_events: list[tuple[int, str]] = []
        for name in expected_stage_order:
            event = latest_events.get(name)
            stage_record = stages.get(name)
            if event is None or not isinstance(stage_record, dict):
                continue
            recorded_attempt = stage_record.get("attempt")
            if isinstance(recorded_attempt, int) and not isinstance(recorded_attempt, bool):
                event = next((
                    (index, item) for index, item in enumerate(history)
                    if isinstance(item, dict) and item.get("stage") == name
                    and item.get("status") == "complete"
                    and item.get("attempt") == recorded_attempt
                ), None)
                if event is None:
                    continue
            recorded_events.append((event[0], name))
        if len(recorded_events) == len(expected_stage_order):
            recorded_order = [name for _, name in sorted(recorded_events)]
            check(
                "upstream_stage_order", expected_stage_order, recorded_order,
                "latest complete history event matching each stage's recorded attempt",
            )
        else:
            checked["upstream_stage_order"] = {
                "expected": expected_stage_order,
                "actual": None,
                "matches": None,
                "evidence": "history does not contain a completed event for every recorded upstream attempt",
            }
    else:
        checked["upstream_stage_order"] = {
            "expected": expected_stage_order,
            "actual": None,
            "matches": None,
            "evidence": "ordered stage history is not recorded",
        }

    # The source wrapper settings below are known from source inspection but
    # are not serialized in this run manifest. Their absence prevents an
    # affirmative full-profile match even if every recorded setting agrees.
    regularization_flags = {
        "wm_topology_regularization": "--wm-topology-regularization",
        "gm_topology_regularization": "--gm-topology-regularization",
    }
    actual_dura_regularization = {}
    for name, flag in regularization_flags.items():
        actual_dura_regularization[name] = (
            float(dura_command[dura_command.index(flag) + 1])
            if flag in dura_command else None
        )
    expected_dura_regularization = {
        "wm_topology_regularization": 0.05,
        "gm_topology_regularization": 0.02,
    }
    regularization_observed = all(value is not None for value in actual_dura_regularization.values())
    checked["dura_source_wrapper_parameter_set"] = {
        "expected": expected_dura_regularization,
        "actual": actual_dura_regularization,
        "matches": (
            all(np.isclose(actual_dura_regularization[name], expected) for name, expected in expected_dura_regularization.items())
            if regularization_observed else None
        ),
        "evidence": "explicit Dura command arguments for selected V3 WM/GM topology regularization",
    }
    checked["toads_source_constructor_parameter_set"] = {
        "expected": "selected run_Toads_jist.m constructor profile",
        "actual": None,
        "matches": None,
        "evidence": "the command/run options omit several constructor parameters",
    }

    scope = manifest.get("workflow_scope")
    modality = (
        toads_command[toads_command.index("--modality") + 1]
        if "--modality" in toads_command else None
    )
    return {
        "profile": "selected-v3-base-upstream-recorded-settings",
        "matches": False if mismatches else (True if all(value["matches"] is not None for value in checked.values()) else None),
        "recorded_settings_match": not mismatches,
        "fully_observed": all(value["matches"] is not None for value in checked.values()),
        "mismatches": mismatches,
        "checked_settings": checked,
        "indicators": {
            "workflow_scope": scope,
            "workflow_scope_mentions_v3": isinstance(scope, str) and "v3" in scope.lower(),
            "spatial_convention": manifest.get("spatial_convention"),
            "source_stage_sequence": [name for name in ("n4", "toads-first", "dura", "toads-second") if name in stages],
            "toads_command_modality_explicit": modality,
            "source_profile_defaults_not_recorded": [
                name for name, value in checked.items() if value["matches"] is None
            ],
            "dura_wrapper_parameters_not_recorded": [
                "curvature=0.5", "level=1.0", "alignment=0.5", "pressure=0.5",
                "target=0.3", "mask=2.0", "level2=0.0",
            ],
            "toads_constructor_parameters_not_recorded": [
                "atlas coefficient=2.0", "smoothing=0.2", "rigid alignment",
                "levels=4", "firstIter=50", "topIter=5", "lastIter=0",
                "maximum distance=1e-4 unless explicit command flag is present",
            ],
            "scope_note": "matches recorded upstream settings only; this does not establish full V3 workflow equivalence",
        },
    }


def _record_file(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"Expected a regular upstream input file: {path}")
    return {"path": str(path.resolve()), "sha256": _sha256(path), "bytes": path.stat().st_size}


def _validate_source_run(source_manifest_path: Path) -> tuple[dict[str, Any], dict[str, Path], dict[str, Any]]:
    """Validate source run identity, completed stage records, hashes and grid."""
    if not source_manifest_path.is_file() or source_manifest_path.is_symlink():
        raise FileNotFoundError(f"Missing source run manifest: {source_manifest_path}")
    try:
        manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot read source run manifest {source_manifest_path}: {exc}") from exc
    if manifest.get("schema_version") != 1 or manifest.get("status") != "complete":
        raise ValueError("Source run must be a completed schema-version-1 workflow manifest")
    if manifest.get("spatial_convention") != "RAS":
        raise ValueError("Source run must declare RAS spatial convention")
    scope = manifest.get("workflow_scope")
    if not isinstance(scope, str) or "base route" not in scope.lower() or "plus" in scope.lower():
        raise ValueError("Source manifest is not the supported completed base-workflow profile")
    run_signature = manifest.get("run_signature")
    if not isinstance(run_signature, str) or re.fullmatch(r"[0-9a-f]{64}", run_signature) is None:
        raise ValueError("Source run must record a SHA-256 run signature")
    implementation = manifest.get("implementation_signature")
    if not isinstance(implementation, str) or re.fullmatch(r"[0-9a-f]{64}", implementation) is None:
        raise ValueError("Source run must record a SHA-256 implementation signature")
    source_inputs = manifest.get("inputs")
    required_inputs = ("t1", "brain_mask", "atlas_descriptor", "maseg")
    if not isinstance(source_inputs, dict) or any(name not in source_inputs for name in required_inputs):
        raise ValueError("Source base run must identify T1, brain mask, atlas descriptor, and MaSeg inputs")
    for role in required_inputs:
        record = source_inputs[role]
        if not isinstance(record, dict) or not isinstance(record.get("path"), str):
            raise ValueError(f"Source run has a malformed {role} input record")
        path = Path(record["path"])
        if not path.is_file() or path.is_symlink() or _sha256(path) != record.get("sha256"):
            raise ValueError(f"Source run {role} input is missing or its recorded SHA-256 changed")
    atlas_assets = manifest.get("atlas_assets")
    if not isinstance(atlas_assets, list) or not atlas_assets:
        raise ValueError("Source base run must record its selected atlas assets")
    for record in atlas_assets:
        if not isinstance(record, dict) or not isinstance(record.get("path"), str):
            raise ValueError("Source run has a malformed atlas asset record")
        path = Path(record["path"])
        if not path.is_file() or path.is_symlink() or _sha256(path) != record.get("sha256"):
            raise ValueError(f"Source atlas asset is missing or its recorded SHA-256 changed: {path}")
    descriptor_path = Path(source_inputs["atlas_descriptor"]["path"])
    declared_assets: set[Path] = set()
    try:
        descriptor_lines = descriptor_path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise ValueError(f"Cannot read source atlas descriptor: {exc}") from exc
    for line in descriptor_lines:
        if not line.startswith("File:"):
            continue
        relative = Path(line.partition(":")[2].strip())
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe source atlas asset path: {relative}")
        asset_path = (descriptor_path.parent / relative).resolve()
        if not asset_path.is_relative_to(descriptor_path.parent.resolve()) or not asset_path.is_file():
            raise ValueError(f"Source atlas descriptor references a missing/unsafe asset: {relative}")
        declared_assets.add(asset_path)
    recorded_assets = {Path(record["path"]).resolve() for record in atlas_assets}
    if declared_assets != recorded_assets:
        raise ValueError("Source atlas descriptor assets disagree with the run manifest atlas asset set")
    audit_record = manifest.get("atlas_conversion_audit")
    if audit_record is not None:
        if not isinstance(audit_record, dict) or not isinstance(audit_record.get("path"), str):
            raise ValueError("Source run has a malformed atlas conversion audit record")
        audit_path = Path(audit_record["path"])
        if not audit_path.is_file() or audit_path.is_symlink() or _sha256(audit_path) != audit_record.get("sha256"):
            raise ValueError("Source atlas conversion audit is missing or its recorded SHA-256 changed")
    t1_record = source_inputs.get("t1") if isinstance(source_inputs, dict) else None
    if not isinstance(t1_record, dict) or not isinstance(t1_record.get("path"), str):
        raise ValueError("Source run does not record its T1 input")
    source_t1 = Path(t1_record["path"])
    if not source_t1.is_file() or source_t1.is_symlink() or _sha256(source_t1) != t1_record.get("sha256"):
        raise ValueError("Source run T1 input is missing or no longer matches its recorded SHA-256")
    t1_image = nib.load(str(source_t1))
    if len(t1_image.shape) != 3 or nib.aff2axcodes(t1_image.affine) != ("R", "A", "S"):
        raise ValueError("Source run T1 must be a 3D RAS NIfTI")
    t1_values = np.asanyarray(t1_image.dataobj)
    if not np.issubdtype(t1_values.dtype, np.number) or not np.isfinite(t1_values).all():
        raise ValueError("Source run T1 contains nonfinite or nonnumeric data")
    brain_image = nib.load(source_inputs["brain_mask"]["path"])
    if len(brain_image.shape) != 3 or brain_image.shape != t1_image.shape or not np.allclose(brain_image.affine, t1_image.affine, rtol=0, atol=1e-4):
        raise ValueError("Source base-run brain mask does not match the T1 grid")
    if nib.aff2axcodes(brain_image.affine) != ("R", "A", "S"):
        raise ValueError("Source base-run brain mask is not RAS")
    brain_values = np.asanyarray(brain_image.dataobj)
    if not np.isfinite(brain_values).all():
        raise ValueError("Source base-run brain mask contains nonfinite voxels")

    stages = manifest.get("stages", {})
    if not isinstance(stages, dict) or any(
        not isinstance(stages.get(name), dict) or stages[name].get("status") != "complete"
        for name in ("n4", "toads-first", "dura", "toads-second")
    ):
        raise ValueError("Source manifest is not a completed N4 → TOADS → Dura → TOADS base workflow")
    for stage_name in ("n4", "toads-first", "dura"):
        records = stages[stage_name].get("outputs")
        if not isinstance(records, list) or not records:
            raise ValueError(f"Source {stage_name} stage has no output hash records")
        stage_dependency = stages[stage_name].get("dependency_signature")
        if not isinstance(stage_dependency, str) or re.fullmatch(r"[0-9a-f]{64}", stage_dependency) is None:
            raise ValueError(f"Source {stage_name} stage has no valid dependency signature")
        for record in records:
            if not isinstance(record, dict) or not isinstance(record.get("path"), str):
                raise ValueError(f"Source {stage_name} stage has a malformed output record")
            path = Path(record["path"])
            if not path.is_file() or path.is_symlink() or path.stat().st_size != record.get("bytes") or _sha256(path) != record.get("sha256"):
                raise ValueError(f"Source {stage_name} stage output hash/size changed: {path}")
    stage = stages.get("toads-second")
    if not isinstance(stage, dict) or stage.get("status") != "complete":
        raise ValueError("Source run has no completed toads-second stage")
    if isinstance(stage.get("attempt"), bool) or not isinstance(stage.get("attempt"), int) or stage["attempt"] < 1:
        raise ValueError("Source toads-second stage has an invalid attempt number")
    dependency_signature = stage.get("dependency_signature")
    if not isinstance(dependency_signature, str) or re.fullmatch(r"[0-9a-f]{64}", dependency_signature) is None:
        raise ValueError("Source toads-second stage has no valid dependency signature")
    if not isinstance(stage.get("toads_stop"), dict) or not isinstance(stage.get("converged"), bool):
        raise ValueError("Source toads-second stage is missing its recorded stop settings/telemetry")
    command = stage.get("command")
    if not isinstance(command, list) or not command or command[0] != "toads-loop":
        raise ValueError("Source toads-second stage is not the selected toads-loop command")
    # These checks protect against importing an unrelated or modified stage
    # while allowing the source wrapper's actual iteration counts to vary.
    checked_flags = (
        "--t1-n4", "--brain-mask", "--atlas-descriptor", "--output", "--dura-output-dir",
        "--registration-levels", "--registration-iterations", "--post-refresh-registration-iterations",
        "--registration-max-iterations", "--maximum-iterations",
    )
    for required_flag in checked_flags:
        if required_flag not in command:
            raise ValueError(f"Source toads-second command omits {required_flag}")
    source_options = manifest.get("options")
    if not isinstance(source_options, dict):
        raise ValueError("Source run does not record workflow settings")
    required_settings = (
        "registration_levels", "registration_iterations", "post_refresh_registration_iterations",
        "registration_max_iterations", "toads_iterations", "dura_iterations",
        "crop_edge_width", "ace_diffusion_iterations", "gvf_iterations", "tgdm_iterations",
    )
    if any(name not in source_options for name in required_settings):
        raise ValueError("Source base-run settings are incomplete for the selected V3 workflow")
    for name in required_settings:
        value = source_options[name]
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"Source base-run setting {name} must be a positive integer")
    dura_command = stages["dura"].get("command")
    if not isinstance(dura_command, list) or not dura_command or dura_command[0] != "dura":
        raise ValueError("Source Dura stage does not record the selected Dura command")
    if "--max-iterations" not in dura_command:
        raise ValueError("Source Dura command omits selected setting --max-iterations")
    if dura_command[dura_command.index("--max-iterations") + 1] != str(source_options["dura_iterations"]):
        raise ValueError("Source Dura iteration setting disagrees with run settings")
    if "--membership-level" in dura_command:
        try:
            dura_membership_level = float(dura_command[dura_command.index("--membership-level") + 1])
        except (ValueError, IndexError) as exc:
            raise ValueError("Source Dura membership setting is malformed") from exc
    else:
        # Older saved base-run manifests omit the parser's selected .3 default
        # from the command and run-options object; preserve that effective value.
        dura_membership_level = 0.3
    if not np.isfinite(dura_membership_level) or not 0.0 <= dura_membership_level <= 1.0:
        raise ValueError("Source Dura membership setting is outside [0, 1]")
    option_names = {
        "--registration-levels": "registration_levels",
        "--registration-iterations": "registration_iterations",
        "--post-refresh-registration-iterations": "post_refresh_registration_iterations",
        "--registration-max-iterations": "registration_max_iterations",
        "--maximum-iterations": "toads_iterations",
    }
    for flag, option_name in option_names.items():
        recorded = command[command.index(flag) + 1]
        if recorded != str(source_options.get(option_name)):
            raise ValueError(f"Source toads-second {flag} setting disagrees with the run settings")
    derived_resources = manifest.get("derived_resources")
    if not isinstance(derived_resources, list) or not derived_resources:
        raise ValueError("Source base run does not record its materialized atlas/config resources")
    verified_derived: set[Path] = set()
    for record in derived_resources:
        if not isinstance(record, dict) or not isinstance(record.get("path"), str):
            raise ValueError("Source run has a malformed derived resource record")
        path = Path(record["path"])
        if not path.is_file() or path.is_symlink() or _sha256(path) != record.get("sha256"):
            raise ValueError(f"Source derived configuration resource changed: {path}")
        verified_derived.add(path.resolve())
    atlas_command_path = Path(command[command.index("--atlas-descriptor") + 1]).resolve()
    if atlas_command_path not in verified_derived:
        raise ValueError("Source toads-second atlas descriptor is not the verified materialized base atlas")
    source_assets_by_relative = {
        Path(record["path"]).resolve().relative_to(descriptor_path.parent.resolve()): Path(record["path"]).resolve()
        for record in atlas_assets
    }
    try:
        materialized_lines = atlas_command_path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise ValueError(f"Cannot read materialized source atlas descriptor: {exc}") from exc
    materialized_assets: set[Path] = set()
    for line in materialized_lines:
        if not line.startswith("File:"):
            continue
        relative = Path(line.partition(":")[2].strip())
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe materialized atlas asset path: {relative}")
        materialized = (atlas_command_path.parent / relative).resolve()
        source = source_assets_by_relative.get(relative)
        if source is None or not materialized.is_file() or _sha256(materialized) != _sha256(source):
            raise ValueError(f"Materialized atlas asset does not match the source run: {relative}")
        materialized_assets.add(materialized)
    if not materialized_assets:
        raise ValueError("Verified materialized source atlas descriptor contains no assets")
    dura_products = {
        Path(record["path"]).name: Path(record["path"])
        for record in stages["dura"].get("outputs", [])
        if isinstance(record, dict) and isinstance(record.get("path"), str)
    }
    for filename in ("stripped_t1.nii.gz", "brain_mask.nii.gz"):
        path = dura_products.get(filename)
        if path is None or not path.is_file() or path.is_symlink():
            raise ValueError(f"Source completed Dura stage does not provide {filename}")
        matching_flag = "--t1-n4" if filename == "stripped_t1.nii.gz" else "--brain-mask"
        if Path(command[command.index(matching_flag) + 1]).resolve() != path.resolve():
            raise ValueError(f"Source toads-second {matching_flag} is not the completed Dura output")
    maximum_iterations = command[command.index("--maximum-iterations") + 1]
    stop = stage["toads_stop"]
    if not isinstance(stop.get("iterations"), int) or stop["iterations"] < 1:
        raise ValueError("Source toads-second stop telemetry has no valid iteration count")
    if not isinstance(stop.get("energy_threshold_reached"), bool) or not isinstance(stop.get("iteration_limit_reached"), bool):
        raise ValueError("Source toads-second stop telemetry is incomplete")
    if stage["converged"] != stop["energy_threshold_reached"]:
        raise ValueError("Source toads-second convergence marker disagrees with stop telemetry")
    if stop["iterations"] > int(maximum_iterations):
        raise ValueError("Source toads-second stop telemetry exceeds the configured iteration limit")

    outputs = stage.get("outputs")
    if not isinstance(outputs, list) or not outputs:
        raise ValueError("Source toads-second stage has no recorded outputs")
    by_name: dict[str, Path] = {}
    for record in outputs:
        if not isinstance(record, dict) or not isinstance(record.get("path"), str):
            raise ValueError("Malformed source toads-second output record")
        path = Path(record["path"])
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"Recorded source stage output is missing or unsafe: {path}")
        if _sha256(path) != record.get("sha256") or path.stat().st_size != record.get("bytes"):
            raise ValueError(f"Recorded source stage output hash/size changed: {path}")
        if path.name.lower().endswith(".nii.gz"):
            image = nib.load(str(path))
            if len(image.shape) != 3 or image.shape != t1_image.shape:
                raise ValueError(f"Source TOADS output does not match source T1 dimensions: {path}")
            if not np.allclose(image.affine, t1_image.affine, rtol=0, atol=1e-4):
                raise ValueError(f"Source TOADS output does not match source T1 affine: {path}")
            if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
                raise ValueError(f"Source TOADS output is not RAS: {path}")
            values = np.asanyarray(image.dataobj)
            if not np.issubdtype(values.dtype, np.number) or not np.isfinite(values).all():
                raise ValueError(f"Source TOADS output is not finite numeric data: {path}")
            if record.get("shape") != list(image.shape) or record.get("axis_codes") != ["R", "A", "S"]:
                raise ValueError(f"Source output geometry does not match its manifest record: {path}")
        if path.name in by_name:
            raise ValueError(f"Duplicate source TOADS output basename: {path.name}")
        by_name[path.name] = path.resolve()
    required = {
        "toads_gm": "stripped_t1_gm.nii.gz",
        "toads_wmfill": "stripped_t1_wmfill.nii.gz",
        "toads_wmmask": "stripped_t1_wmmask.nii.gz",
    }
    missing = [filename for filename in required.values() if filename not in by_name]
    if missing:
        raise ValueError("Source toads-second output set is missing: " + ", ".join(missing))
    selected = {role: by_name[filename] for role, filename in required.items()}
    profile_record = _selected_v3_profile_record(
        manifest, stages, source_options, command, dura_command, dura_membership_level
    )
    upstream = {
        "source_manifest": _record_file(source_manifest_path),
        "source_t1": _record_file(source_t1),
        "source_run_signature": run_signature,
        "source_implementation_signature": implementation,
        "source_stage": "toads-second",
        "source_stage_attempt": stage.get("attempt"),
        "source_stage_settings": {
            "workflow_options": source_options,
            "dura_command": dura_command,
            "dura_membership_level": dura_membership_level,
            "command": command,
            "toads_stop": stop,
            "converged": stage["converged"],
        },
        "selected_v3_profile_match": profile_record,
        "source_stage_outputs": [dict(record) for record in outputs],
    }
    return manifest, selected, upstream


def run_base_resume_from_toads(
    *, source_run: str | Path, maseg_path: str | Path, output_dir: str | Path,
    resume: bool = False,
) -> Path:
    """Run the base-workflow tail from a verified ``toads-second`` checkpoint."""
    source_run = Path(source_run).expanduser().resolve()
    maseg_path = Path(maseg_path).expanduser().resolve()
    output_dir = Path(output_dir).expanduser().resolve()
    if output_dir == source_run.parent or output_dir.is_relative_to(source_run.parent):
        raise ValueError("Tail output directory must be separate from the source run directory")
    source_manifest, toads_maps, upstream = _validate_source_run(source_run)
    if not maseg_path.is_file() or maseg_path.is_symlink() or not maseg_path.name.lower().endswith(".nii.gz"):
        raise ValueError("MaSeg input must be a regular .nii.gz file")
    source_t1 = Path(source_manifest["inputs"]["t1"]["path"])
    reference = normalize_image_to_ras(nib.load(str(source_t1)), source_name="TOADS source T1")
    maseg = normalize_image_to_ras(nib.load(str(maseg_path)), source_name="base tail MaSeg")
    maseg_values = np.asanyarray(maseg.dataobj)
    if len(maseg.shape) != 3 or maseg.shape != reference.shape or not np.allclose(maseg.affine, reference.affine, rtol=0, atol=1e-4):
        raise ValueError("MaSeg must match the source run T1 grid")
    if not has_finite_integer_values(maseg_values):
        raise ValueError("MaSeg input must contain finite integer labels")
    recorded_maseg = source_manifest["inputs"]["maseg"]
    if _sha256(maseg_path) != recorded_maseg.get("sha256"):
        raise ValueError("MaSeg must match the source base-run input by SHA-256")
    source_options = source_manifest["options"]
    settings = {
        "crop_edge_width": _V3_CROP_EDGE_WIDTH,
        "crop_edge_width_source_run": source_options.get("crop_edge_width"),
        "ace_diffusion_iterations": 25,
        "ace_use_updated_csf_in_ace": 2,
        "gvf_lambda": 1.0,
        "gvf_iterations": 200,
        "gvf_solver_tolerance": 1e-6,
        "tgdm_initial_iterations": 6,
        "tgdm_inner_iterations": 4,
        "tgdm_central_iterations": 7,
        "wm_iso_level": 0.5,
        "gm_confident_level": 0.001,
        "remove_cerebellum": False,
        "refinement_source_faithful": True,
    }
    payload = {
        "schema_version": 1,
        "workflow_profile": "base-resume-from-toads",
        "spatial_convention": "RAS",
        "source_run_manifest": upstream["source_manifest"],
        "source_run_signature": upstream["source_run_signature"],
        "source_toads_second": upstream,
        "inputs": {"maseg": _record_file(maseg_path)},
        "options": settings,
    }
    signature = _json_hash(payload)
    implementation = _implementation_signature()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "run.json"
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"Tail run manifest exists; use --resume to continue: {manifest_path}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("run_signature") != signature:
            raise ValueError("Cannot resume: imported source run, MaSeg, or tail settings changed")
    else:
        if resume:
            raise FileNotFoundError(f"Cannot resume; no tail manifest exists: {manifest_path}")
        manifest = {
            **payload,
            "run_signature": signature,
            "implementation_signature": implementation,
            "status": "running",
            "stages": {},
            "history": [],
            "source_equivalence_claim": False,
        }
        _write_manifest(manifest_path, manifest)
    if manifest.get("implementation_signature") != implementation:
        manifest["previous_implementation_signature"] = manifest.get("implementation_signature")
        manifest["implementation_signature"] = implementation
    manifest["status"] = "running"
    _write_manifest(manifest_path, manifest)

    expected_stage_signatures: dict[str, str] = {}

    def execute_stage(name: str, relative_outputs: tuple[str, ...], command_template: list[str], dependencies: list[Path]) -> dict[str, Path]:
        dependency_signature = _json_hash({
            "implementation_signature": implementation,
            "dependencies": [{"path": str(p.resolve()), "sha256": _sha256(p)} for p in dependencies],
        })
        expected_stage_signatures[name] = dependency_signature
        previous = manifest["stages"].get(name)
        if _stage_is_valid(previous, dependency_signature):
            return {Path(item["path"]).name: Path(item["path"]) for item in previous["outputs"]}
        attempt = 1 + sum(item.get("stage") == name for item in manifest["history"])
        stage_root = output_dir / "stages" / name
        while (stage_root / f"attempt-{attempt:02d}").exists():
            attempt += 1
        stage_dir = stage_root / f"attempt-{attempt:02d}"
        stage_dir.mkdir(parents=True)
        stdout, stderr = io.StringIO(), io.StringIO()
        log_path = stage_dir / "stage.log"
        command = [arg.replace("{STAGE}", str(stage_dir)) for arg in command_template]
        try:
            with log_path.open("w", encoding="utf-8") as log:
                class Tee:
                    def write(self, value: str) -> int:
                        stdout.write(value); n = log.write(value); log.flush(); return n
                    def flush(self) -> None:
                        stdout.flush(); log.flush()
                with redirect_stdout(Tee()), redirect_stderr(stderr):
                    cli_main(command)
                if stderr.getvalue():
                    log.write("\n[stderr]\n" + stderr.getvalue()); log.flush()
            outputs = [stage_dir / rel for rel in relative_outputs]
            records = _output_records([*outputs, log_path], reference, require_reference_grid=name not in {"base-tissue-preparation", "cruise-core"})
            if name == "base-tissue-preparation":
                full_names = [n for n in relative_outputs if n.endswith(".nii.gz") and "_cropped" not in n and not n.startswith("boundary_csf") and n not in {"fused_gm.nii.gz", "fused_wmfill.nii.gz", "fused_mask.nii.gz"}]
                full_images = [nib.load(str(stage_dir / n)) for n in full_names]
                if any(im.shape != reference.shape or not np.allclose(im.affine, reference.affine, rtol=0, atol=1e-4) for im in full_images):
                    raise RuntimeError("Full-grid MaSeg/TOADS preparation maps do not match source T1")
                crop_names = [n for n in relative_outputs if n.endswith("_cropped.nii.gz")]
                crop_names.extend(n for n in ("boundary_csf_raw.nii.gz", "boundary_csf_toads_cleaned.nii.gz", "boundary_csf.nii.gz") if n in relative_outputs)
                crop_names.extend(n for n in ("fused_gm.nii.gz", "fused_wmfill.nii.gz", "fused_mask.nii.gz") if n in relative_outputs)
                images = [nib.load(str(stage_dir / n)) for n in crop_names]
                if not images or any(im.shape != images[0].shape or not np.allclose(im.affine, images[0].affine, rtol=0, atol=1e-4) for im in images[1:]):
                    raise RuntimeError("MaSeg/TOADS preparation crops do not share a grid")
                if nib.aff2axcodes(images[0].affine) != ("R", "A", "S"):
                    raise RuntimeError("Preparation crop is not RAS")
            if name == "cruise-core":
                for filename in ("inner.nii.gz", "central.nii.gz", "outer.nii.gz"):
                    image = nib.load(str(stage_dir / filename))
                    if image.shape != reference.shape or not np.allclose(image.affine, reference.affine, rtol=0, atol=1e-4):
                        raise RuntimeError(f"CRUISE output is not on source T1 grid: {filename}")
                cropped_names = ("gm_cropped.nii.gz", "wmfill_cropped.nii.gz", "wmmask_cropped.nii.gz",
                                 "wmfill_topology_corrected.nii.gz", "ace_enhanced_gm.nii.gz",
                                 "ace_skeleton.nii.gz", "ace_thinned_skeleton.nii.gz")
                cropped = [nib.load(str(stage_dir / n)) for n in cropped_names]
                if any(im.shape != cropped[0].shape or not np.allclose(im.affine, cropped[0].affine, rtol=0, atol=1e-4) for im in cropped[1:]):
                    raise RuntimeError("CRUISE cropped scalar products do not share one grid")
                gvf = nib.load(str(stage_dir / "gvf.nii.gz"))
                if gvf.shape != (*cropped[0].shape, 3) or not np.allclose(gvf.affine, cropped[0].affine, rtol=0, atol=1e-4):
                    raise RuntimeError("CRUISE GVF product does not share the cropped grid")
            if name == "surface-thickness":
                for filename in ("inner_surface.vtk", "central_surface.vtk", "outer_surface.vtk", "inner_surface_thickness.vtk"):
                    mesh = read_legacy_vtk(stage_dir / filename)
                    if mesh.coordinate_frame != "RAS" or mesh.units != "mm":
                        raise RuntimeError(f"Surface product does not declare RAS millimetres: {filename}")
                    if mesh.vertices.ndim != 2 or mesh.vertices.shape[1] != 3 or not np.isfinite(mesh.vertices).all():
                        raise RuntimeError(f"Surface product contains invalid RAS coordinates: {filename}")
            manifest["stages"][name] = {
                "status": "complete", "attempt": attempt, "command": command,
                "dependency_signature": dependency_signature, "outputs": records,
            }
            manifest["history"].append({"stage": name, "attempt": attempt, "status": "complete"})
            _write_manifest(manifest_path, manifest)
            return {Path(item["path"]).name: Path(item["path"]) for item in records}
        except Exception as exc:
            manifest["stages"][name] = {"status": "failed", "attempt": attempt, "error": f"{type(exc).__name__}: {exc}"}
            manifest["history"].append({"stage": name, "attempt": attempt, "status": "failed"})
            _write_manifest(manifest_path, manifest)
            raise

    prep_names = (
        "atlas_gm.nii.gz", "atlas_wm.nii.gz", "atlas_mask.nii.gz", "atlas_brain_tissue_mask.nii.gz",
        "toads_gm_transferred.nii.gz", "toads_wmfill_transferred.nii.gz", "toads_mask_transferred.nii.gz",
        "atlas_gm_cropped.nii.gz", "atlas_wm_cropped.nii.gz", "atlas_mask_cropped.nii.gz",
        "atlas_brain_tissue_mask_cropped.nii.gz", "toads_gm_cropped.nii.gz", "toads_wmfill_cropped.nii.gz",
        "toads_mask_cropped.nii.gz", "boundary_csf_raw.nii.gz", "boundary_csf_toads_cleaned.nii.gz",
        "boundary_csf.nii.gz", "fused_gm.nii.gz", "fused_wmfill.nii.gz", "fused_mask.nii.gz", "provenance.json",
    )
    prep_deps = [maseg_path, *toads_maps.values()]
    prep = execute_stage("base-tissue-preparation", prep_names, [
        "base-tissue-preparation", "--maseg", str(maseg_path), "--gm", str(toads_maps["toads_gm"]),
        "--wm-fill", str(toads_maps["toads_wmfill"]), "--wm-mask", str(toads_maps["toads_wmmask"]),
        "--edge-width", str(_V3_CROP_EDGE_WIDTH), "--gm-confident-level", "0.001",
        "--output-dir", "{STAGE}", "--provenance", "{STAGE}/provenance.json",
    ], prep_deps)
    core_names = (
        "gm_cropped.nii.gz", "wmfill_cropped.nii.gz", "wmmask_cropped.nii.gz", "wmfill_topology_corrected.nii.gz",
        "ace_enhanced_gm.nii.gz", "ace_skeleton.nii.gz", "ace_thinned_skeleton.nii.gz", "gvf.nii.gz",
        "inner.nii.gz", "central.nii.gz", "outer.nii.gz", "provenance.json",
    )
    prep_prov = prep["provenance.json"]
    core_deps = [prep[n] for n in ("fused_gm.nii.gz", "fused_wmfill.nii.gz", "fused_mask.nii.gz", "boundary_csf.nii.gz")] + [prep_prov]
    core = execute_stage("cruise-core", core_names, [
        "toads-cruise-core", "--gm", str(core_deps[0]), "--wm-fill", str(core_deps[1]), "--wm-mask", str(core_deps[2]),
        "--boundary-csf", str(core_deps[3]),
        "--base-preparation-provenance", str(prep_prov), "--edge-width", str(_V3_CROP_EDGE_WIDTH),
        "--ace-diffusion-iterations", "25", "--use-updated-csf-in-ace", str(settings["ace_use_updated_csf_in_ace"]),
        "--gvf-iterations", "200", "--wm-iso-level", str(settings["wm_iso_level"]),
        "--tgdm-initial-iterations", "6", "--tgdm-inner-iterations", "4",
        "--tgdm-central-iterations", "7", "--output-dir", "{STAGE}",
        "--provenance", "{STAGE}/provenance.json",
    ], core_deps)
    refine_names = ("SEG_refine.nii.gz", "SEG_refine-summary.json")
    refine = execute_stage("maseg-refinement", refine_names, [
        "segmentation-refine", "--segmentation", str(maseg_path), "--inner-levelset", str(core["inner.nii.gz"]),
        "--outer-levelset", str(core["outer.nii.gz"]), "--source-faithful", "--output", "{STAGE}/SEG_refine.nii.gz",
        "--summary", "{STAGE}/SEG_refine-summary.json",
    ], [maseg_path, core["inner.nii.gz"], core["outer.nii.gz"]])
    surface_names = ("inner_surface.vtk", "central_surface.vtk", "outer_surface.vtk", "inner_surface_thickness.vtk", "thickness_levelset.nii.gz", "thickness_summary.json", "provenance.json")
    surfaces = execute_stage("surface-thickness", surface_names, [
        "base-surface-thickness", "--inner-levelset", str(core["inner.nii.gz"]), "--central-levelset", str(core["central.nii.gz"]),
        "--outer-levelset", str(core["outer.nii.gz"]), "--output-dir", "{STAGE}", "--provenance", "{STAGE}/provenance.json",
    ], [core[n] for n in ("inner.nii.gz", "central.nii.gz", "outer.nii.gz")])
    # Recheck both the immutable source checkpoint and every tail-stage record
    # before sealing the manifest. This catches upstream edits during execution
    # and prevents an interrupted or stale stage from appearing as complete.
    _, _, verified_upstream = _validate_source_run(source_run)
    if verified_upstream != upstream:
        raise RuntimeError("Imported source run or TOADS outputs changed while the tail was running")
    for name, dependency_signature in expected_stage_signatures.items():
        if not _stage_is_valid(manifest["stages"].get(name), dependency_signature):
            raise RuntimeError(f"Cannot complete tail manifest; stage {name} failed final hash/signature verification")
    manifest["status"] = "complete"
    manifest["source_equivalence_claim"] = False
    manifest["final_outputs"] = [
        *manifest["stages"]["cruise-core"]["outputs"],
        *manifest["stages"]["maseg-refinement"]["outputs"],
        *manifest["stages"]["surface-thickness"]["outputs"],
    ]
    _write_manifest(manifest_path, manifest)
    return manifest_path
