"""Restartable selected-source MaSeg refinement for the PLUS workflow."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.image_io import save_nifti_ras
from pymacruise.segmentation_refine import refine_segmentation_by_surfaces


@dataclass(frozen=True)
class PlusSegmentationRefinementResult:
    manifest_path: Path
    output_path: Path
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


def run_plus_segmentation_refinement_checkpoint(
    segmentation_path: str | Path,
    inner_levelset_path: str | Path,
    outer_levelset_path: str | Path,
    output_dir: str | Path,
    *,
    resume: bool = False,
) -> PlusSegmentationRefinementResult:
    """Apply the selected V3 `run_SegmentationRefine` rules with hash-checked resume."""
    inputs = {
        "segmentation": Path(segmentation_path).resolve(),
        "inner_levelset": Path(inner_levelset_path).resolve(),
        "outer_levelset": Path(outer_levelset_path).resolve(),
    }
    for path in inputs.values():
        # Match the other input readers: resolving an input symlink is allowed.
        if not path.is_file():
            raise FileNotFoundError(f"Expected regular PLUS refinement input: {path}")
    root = Path(output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "plus-segmentation-refinement.json"
    output_path = root / "SEG_refine.nii.gz"
    signature_data = {
        "schema_version": 1,
        "stage": "plus-segmentation-refinement",
        "inputs": {name: {"path": str(path), "sha256": _sha256(path)} for name, path in inputs.items()},
        "settings": {
            "source": "MaCRUISE SVN r9360 run_SegmentationRefine / MAsegRefineBySurf",
            "inner_coefficient_mm": 0.2,
            "outer_coefficient_mm": 0.2,
            "inner_central_coefficient_mm": 0.2,
            "outer_central_coefficient_mm": 0.2,
            "deal_with_central": True,
            "neighborhood_radius": 2,
            "max_sweeps": 50,
            "source_faithful": True,
        },
        "spatial_convention": "RAS",
        "implementation_signature": {
            filename: _sha256(Path(__file__).with_name(filename))
            for filename in ("plus_segmentation_refinement.py", "segmentation_refine.py", "image_io.py")
        },
    }
    signature = hashlib.sha256(
        json.dumps(signature_data, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if manifest_path.is_symlink():
        raise FileExistsError(f"Refusing to write through symlinked manifest: {manifest_path}")
    if manifest_path.exists():
        if not resume:
            raise FileExistsError(f"PLUS refinement manifest exists; use resume=True: {manifest_path}")
        previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        if previous.get("run_signature") != signature:
            raise ValueError("Cannot resume PLUS refinement: inputs or settings differ")
        record = previous.get("outputs", {}).get("SEG_refine.nii.gz", {})
        if (previous.get("status") == "complete" and output_path.is_file()
                and not output_path.is_symlink() and record.get("sha256") == _sha256(output_path)):
            return PlusSegmentationRefinementResult(manifest_path, output_path, True)

    manifest = {
        **signature_data,
        "run_signature": signature,
        "status": "running",
        "outputs": {},
    }
    _write_json(manifest_path, manifest)
    try:
        result = refine_segmentation_by_surfaces(
            nib.load(str(inputs["segmentation"])),
            nib.load(str(inputs["inner_levelset"])),
            nib.load(str(inputs["outer_levelset"])),
            inner_coefficient_mm=0.2,
            outer_coefficient_mm=0.2,
            inner_central_coefficient_mm=0.2,
            outer_central_coefficient_mm=0.2,
            deal_with_central=True,
            neighborhood_radius=2,
            max_sweeps=50,
            source_faithful=True,
        )
        # Preserve segmentation geometry on the normalized full native grid.
        save_nifti_ras(result.image, output_path)
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        _write_json(manifest_path, manifest)
        raise
    manifest["status"] = "complete"
    central_changes = result.changed_masks.get("central_fill")
    manifest["rule_metrics"] = {
        name: {
            "candidate_voxels": int(result.candidate_masks[name].sum()),
            "changed_voxels": int(
                result.changed_voxels[name]
                if name in result.changed_voxels
                else np.count_nonzero(result.candidate_masks[name] & central_changes)
            ),
            "unresolved_voxels": int(
                result.unresolved_voxels[name]
                if name in result.unresolved_voxels
                else result.candidate_masks[name].sum()
                - np.count_nonzero(result.candidate_masks[name] & central_changes)
            ),
            "sweeps": int(result.sweeps.get(name, result.sweeps.get("central_fill", 0))),
        }
        for name in result.candidate_masks
    }
    manifest["outputs"] = {
        output_path.name: {"path": str(output_path), "sha256": _sha256(output_path)}
    }
    _write_json(manifest_path, manifest)
    return PlusSegmentationRefinementResult(manifest_path, output_path, False)
