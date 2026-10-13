#!/usr/bin/env python3
"""Run source-selected plus surface stages from saved loop checkpoints.

This composes the persisted LesionTOADS loop result with the existing Dura,
native-space transfer, MaSeg separation, and plus surface APIs. It does not
claim the stages were executed by one uninterrupted ``plus-full-run`` call.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pymacruise.plus_file_workflow import run_plus_file_workflow  # noqa: E402
from pymacruise.plus_lesion_dura_output_checkpoint import (  # noqa: E402
    run_plus_lesion_dura_output_checkpoint,
)
from pymacruise.plus_lesion_native_dura_transfer import (  # noqa: E402
    run_plus_lesion_native_dura_transfer,
)
from pymacruise.plus_maseg_cortical_checkpoint import (  # noqa: E402
    run_plus_maseg_cortical_checkpoint,
)
from pymacruise.plus_pipeline import PlusPipelineOptions  # noqa: E402
from pymacruise.tgdm import TGDMBaselineParameters  # noqa: E402


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _record(path: Path) -> dict[str, str]:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink():
        raise FileNotFoundError(f"Expected a regular file: {resolved}")
    return {"path": str(resolved), "sha256": _sha256(resolved)}


def _write_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--loop-manifest", required=True, type=Path)
    parser.add_argument("--t1", required=True, type=Path)
    parser.add_argument("--brain-mask", required=True, type=Path)
    parser.add_argument("--maseg", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--dura-transfer-mode",
        choices=("physical-affine", "jist-unchanged"),
        default="physical-affine",
        help="native Dura transfer behavior (default: physical-affine)",
    )
    parser.add_argument(
        "--allow-checkpoint",
        action="store_true",
        help="compose surfaces from a nonterminal registered-loop checkpoint",
    )
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    loop_path = args.loop_manifest.expanduser().resolve()
    loop = json.loads(loop_path.read_text(encoding="utf-8"))
    if (
        loop.get("status") != "complete"
        or loop.get("stage") != "plus-lesion-registered-loop"
        or loop.get("spatial_convention") != "RAS"
    ):
        raise ValueError("input must be a complete RAS registered-loop checkpoint")
    summary_record = loop.get("outputs", {}).get("summary", {})
    summary_path = Path(summary_record.get("path", "")).resolve()
    if not summary_path.is_file() or _sha256(summary_path) != summary_record.get("sha256"):
        raise ValueError("registered-loop summary is missing or fails its checkpoint hash")
    loop_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    terminal = loop_summary.get("stop") is True
    if not terminal and not args.allow_checkpoint:
        raise ValueError(
            "nonterminal loop state requires --allow-checkpoint; outputs will remain marked nonterminal"
        )
    recrop_record = loop.get("inputs", {}).get("registered_recrop_manifest", {})
    recrop_path = Path(recrop_record.get("path", "")).resolve()
    if not recrop_path.is_file() or _sha256(recrop_path) != recrop_record.get("sha256"):
        raise ValueError("registered recrop manifest is missing or fails its loop-checkpoint hash")

    output_root = args.output_dir.expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    dura = run_plus_lesion_dura_output_checkpoint(
        recrop_path,
        loop_path,
        output_root / "lesion-dura",
        resume=args.resume,
    )
    native = run_plus_lesion_native_dura_transfer(
        dura.manifest_path,
        args.t1,
        output_root / "native-dura",
        transfer_mode=args.dura_transfer_mode,
        resume=args.resume,
    )
    cortical = run_plus_maseg_cortical_checkpoint(
        args.maseg,
        output_root / "maseg-cortical-separation",
        resume=args.resume,
    )
    plus = run_plus_file_workflow(
        t1_reference_path=args.t1,
        segmentation_path=args.maseg,
        dura_csf_path=native.outputs["csf"],
        dura_gm_path=native.outputs["gm"],
        dura_wm_fill_path=native.outputs["wm_fill"],
        dura_wm_mask_path=native.outputs["wm_mask"],
        multi_atlas_gm_path=cortical.outputs["Atlas_GM"],
        multi_atlas_wm_path=cortical.outputs["Atlas_WM"],
        multi_atlas_mask_path=cortical.outputs["Atlas_BrainTissueMask"],
        brain_mask_path=args.brain_mask,
        output_dir=output_root / "plus-surfaces",
        dura_transfer_mode=args.dura_transfer_mode,
        options=PlusPipelineOptions(),
        resume=args.resume,
    )

    stages = {
        "registered_recrop": recrop_path,
        "registered_loop": loop_path,
        "lesion_dura": dura.manifest_path,
        "native_dura_transfer": native.manifest_path,
        "maseg_cortical_separation": cortical.manifest_path,
        "plus_surface_workflow": plus.manifest_path,
    }
    inputs = {
        "t1": _record(args.t1),
        "brain_mask": _record(args.brain_mask),
        "maseg": _record(args.maseg),
        "registered_loop": _record(loop_path),
        "registered_recrop": _record(recrop_path),
    }
    surface_names = ("inner_level_set", "central_level_set", "outer_level_set")
    final_outputs = {
        f"plus_{name}": _record(plus.outputs[name])
        for name in surface_names
    }
    stage_records = {name: _record(path) for name, path in stages.items()}
    settings = {
        "plus_pipeline_options": PlusPipelineOptions().as_dict(),
        "tgdm_surface_evolution": asdict(TGDMBaselineParameters()),
        "dura_transfer_mode": args.dura_transfer_mode,
        "file_workflow_run_manifest": _record(plus.manifest_path),
    }
    composition_stage = (
        "synthetic-plus-terminal-surface-composition"
        if terminal else "synthetic-plus-checkpoint-surface-composition"
    )
    signature_data = {
        "schema_version": 1,
        "stage": composition_stage,
        "spatial_convention": "RAS",
        "stop_reason": loop_summary.get("stop_reason"),
        "completed_iteration": loop_summary.get("completed_iteration"),
        "inputs": inputs,
        "stages": stage_records,
        "settings": settings,
        "final_outputs": final_outputs,
        "source_equivalence_claim": False,
    }
    run_signature = hashlib.sha256(
        json.dumps(signature_data, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()
    manifest = {
        **signature_data,
        "run_signature": run_signature,
        "status": "complete",
        "terminal": bool(terminal),
        "workflow_scope": "source-selected plus surface composition from hash-checked registered-loop state; "
        + ("terminal; stages were composed from checkpoints, not one uninterrupted plus-full-run invocation"
           if terminal else "nonterminal checkpoint audit; stages were composed from saved loop outputs"),
    }
    manifest_path = output_root / f"{composition_stage}.json"
    _write_json(manifest_path, manifest)
    label = "Terminal" if terminal else "Nonterminal checkpoint"
    print(f"{label} synthetic plus surface composition: {manifest_path}")
    print(f"Stop reason: {manifest['stop_reason']}")
    print(f"Completed loop iteration: {manifest['completed_iteration']}")


if __name__ == "__main__":
    main()
