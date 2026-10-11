"""Opt-in replay and lineage diagnostics for saved nested-TGDM crop inputs.

This module is intentionally separate from production workflow code. It replays
only the nested-TGDM stage from already-saved crop dependencies and writes its
intermediate level sets to a new diagnostic directory.
"""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import argparse
import json
from pathlib import Path
import platform
import resource
import sys
import time
from typing import Any

import nibabel as nib
import numpy as np
from scipy import ndimage

from pymacruise.cropping import CropParameters, crop_image
from pymacruise.image_io import save_nifti_ras
from pymacruise.tgdm import TGDMBaselineParameters, tgdm_nested_arrays


_CROP_INPUTS = {
    "initial_membership": "wmfill_topology_corrected.nii.gz",
    "white_matter": "wmfill_cropped.nii.gz",
    "gray_matter": "gm_cropped.nii.gz",
    "ace_gray_matter": "ace_enhanced_gm.nii.gz",
    "ace_skeleton": "ace_thinned_skeleton.nii.gz",
    "gvf": "gvf.nii.gz",
}
_OUTPUT_FIELDS = {
    "01_initial_membership_levelset": lambda result: result.initial_level_set,
    "02_smoothed_levelset": lambda result: result.smoothed_initial.level_set,
    "03_initial_inner_levelset": lambda result: result.initial_inner.level_set,
    "04_central_levelset": lambda result: result.central.level_set,
    "05_central_shrink_levelset": lambda result: result.central_shrink.level_set,
    "06_delivered_inner_levelset": lambda result: result.inner.level_set,
}


def levelset_component_summary(level_set: np.ndarray) -> dict[str, Any]:
    """Summarize nonpositive components without retaining label volumes."""
    values = np.asarray(level_set)
    if values.ndim != 3 or not np.isfinite(values).all():
        raise ValueError("Level set must be a finite 3D array")
    negative = values <= 0
    summary: dict[str, Any] = {"negative_voxels": int(np.count_nonzero(negative))}
    for connectivity, name in ((1, "6"), (2, "18"), (3, "26")):
        labels, count = ndimage.label(
            negative,
            structure=ndimage.generate_binary_structure(3, connectivity),
        )
        sizes = np.bincount(labels.ravel())[1:]
        summary[f"connectivity_{name}"] = {
            "components": int(count),
            "singleton_components": int(np.count_nonzero(sizes == 1)),
            "components_under_10_voxels": int(np.count_nonzero(sizes < 10)),
            "largest_component_voxels": int(sizes.max(initial=0)),
            "smallest_component_voxels": int(sizes.min()) if count else None,
        }
        del labels, sizes
    return summary


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _max_rss_gib() -> float:
    peak = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    # macOS reports bytes; Linux and the other supported POSIX targets report KiB.
    return peak / (1024**3 if sys.platform == "darwin" else 1024**2)


def _array_comparison(actual: np.ndarray, expected: np.ndarray) -> dict[str, Any]:
    if actual.shape != expected.shape:
        return {
            "same_shape": False,
            "actual_shape": list(actual.shape),
            "expected_shape": list(expected.shape),
        }
    written = np.ones(actual.shape, dtype=bool)
    written[-1, :, :] = False
    written[:, -1, :] = False
    written[:, :, -1] = False

    def metrics(mask: np.ndarray) -> dict[str, Any]:
        left = actual[mask]
        right = expected[mask]
        difference = np.abs(left.astype(np.float64) - right.astype(np.float64))
        return {
            "voxels": int(left.size),
            "bitwise_equal": bool(np.array_equal(left, right)),
            "different_voxels": int(np.count_nonzero(left != right)),
            "max_absolute_error": float(difference.max(initial=0.0)),
            "mean_absolute_error": float(difference.mean()),
        }

    return {
        "same_shape": True,
        "comparison_note": (
            "The selected JIST uncrop uses exclusive high bounds; final crop planes "
            "are unwritten and therefore compared separately."
        ),
        "legacy_uncrop_written_region": metrics(written),
        "exclusive_high_boundary_planes": metrics(~written),
        "all_crop_voxels": metrics(np.ones(actual.shape, dtype=bool)),
    }


def _singleton_negative_coordinates(level_set: np.ndarray) -> np.ndarray:
    labels, _count = ndimage.label(
        np.asarray(level_set) <= 0,
        structure=ndimage.generate_binary_structure(3, 1),
    )
    sizes = np.bincount(labels.ravel())
    singleton_ids = np.flatnonzero(sizes == 1)
    singleton_ids = singleton_ids[singleton_ids != 0]
    return np.argwhere(np.isin(labels, singleton_ids))


def run_saved_tgdm_diagnostic(
    core_attempt_dir: str | Path,
    output_dir: str | Path,
    *,
    overwrite: bool = False,
) -> Path:
    """Replay TGDM from one saved base-core attempt and save six crop fields.

    No TOADS, ACE, GVF, or full MaCRUISE stage is rerun. The destination must
    be new unless ``overwrite`` is explicitly requested.
    """
    core_dir = Path(core_attempt_dir).resolve()
    destination = Path(output_dir).resolve()
    if not core_dir.is_dir():
        raise FileNotFoundError(f"Core attempt directory not found: {core_dir}")
    if destination.exists() and any(destination.iterdir()) and not overwrite:
        raise FileExistsError(
            f"Diagnostic output directory is not empty: {destination}; use --overwrite"
        )

    started = time.monotonic()
    provenance_path = core_dir / "provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    bounds_json = provenance["parameters"]["crop_bounds"]
    bounds = tuple(tuple(int(v) for v in axis) for axis in bounds_json)
    reference = nib.load(str(core_dir / "inner.nii.gz"))
    crop = CropParameters(
        original_shape=tuple(int(v) for v in reference.shape),
        bounds=bounds,  # type: ignore[arg-type]
        background=20.0,
        original_affine=np.asarray(reference.affine, dtype=np.float64),
    )
    stage_root = core_dir.parent.parent
    inputs: dict[str, tuple[Path, nib.spatialimages.SpatialImage, np.ndarray]] = {}
    crop_reference_affine: np.ndarray | None = None
    for key, filename in _CROP_INPUTS.items():
        path = core_dir / filename
        image = nib.load(str(path))
        if len(image.shape) != 3 + (1 if key == "gvf" else 0):
            # `gvf` is a 4D vector field, all remaining fields are scalar.
            raise ValueError(f"Unexpected dimensionality for {key}: {image.shape}")
        expected_grid = crop.cropped_shape
        if tuple(image.shape[:3]) != expected_grid:
            raise ValueError(
                f"{key} grid {image.shape[:3]} differs from recorded crop {expected_grid}"
            )
        if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            raise ValueError(f"{key} is not RAS-oriented")
        if key == "gvf" and image.shape[3] != 3:
            raise ValueError(f"GVF must have exactly three vector components; got {image.shape}")
        if crop_reference_affine is None:
            crop_reference_affine = np.asarray(image.affine, dtype=np.float64)
        elif not np.allclose(image.affine, crop_reference_affine, rtol=0, atol=1e-5):
            raise ValueError(f"{key} affine differs from the TGDM crop grid")
        array = np.asarray(image.dataobj, dtype=np.float32)
        if not np.isfinite(array).all():
            raise ValueError(f"{key} contains NaN or infinite values")
        inputs[key] = (path, image, array)

    params = TGDMBaselineParameters(**provenance["parameters"]["tgdm"])
    execution = tgdm_nested_arrays(
        inputs["initial_membership"][2],
        inputs["white_matter"][2],
        inputs["gray_matter"][2],
        inputs["ace_gray_matter"][2],
        ace_skeleton=inputs["ace_skeleton"][2],
        gvf=inputs["gvf"][2],
        parameters=params,
    )
    input_metadata = {
        key: {
            "path": str(path),
            "sha256": _sha256(path),
            "shape": list(image.shape),
            "axis_codes": list(nib.aff2axcodes(image.affine)),
        }
        for key, (path, image, _array) in inputs.items()
    }
    # The evolution result is self-contained; release input array buffers before
    # writing and analyzing the saved crop fields.
    del inputs

    destination.mkdir(parents=True, exist_ok=True)
    reference_crop = nib.load(str(core_dir / _CROP_INPUTS["initial_membership"]))
    outputs: dict[str, Any] = {}
    for stem, select in _OUTPUT_FIELDS.items():
        values = np.asarray(select(execution), dtype=np.float32)
        path = destination / f"{stem}.nii.gz"
        image = nib.Nifti1Image(values, reference_crop.affine, header=reference_crop.header.copy())
        save_nifti_ras(image, path, dtype=np.float32)
        outputs[stem] = {
            "path": str(path),
            "sha256": _sha256(path),
            "shape": list(values.shape),
            "levelset_components": levelset_component_summary(values),
        }

    comparisons: dict[str, Any] = {}
    for result_name in ("inner", "central", "outer"):
        saved = nib.load(str(core_dir / f"{result_name}.nii.gz"))
        saved_crop = crop_image(saved, crop)
        actual = (
            execution.inner.level_set
            if result_name == "inner"
            else execution.central.level_set
            if result_name == "central"
            else execution.outer.level_set
        )
        comparisons[result_name] = _array_comparison(
            actual, np.asarray(saved_crop.dataobj, dtype=np.float32)
        )
        if result_name in {"inner", "central"}:
            replay_singletons = _singleton_negative_coordinates(actual)
            saved_singletons = _singleton_negative_coordinates(
                np.asarray(saved_crop.dataobj, dtype=np.float32)
            )
            comparisons[result_name]["negative_singletons_6_connected"] = {
                "replay_count": int(len(replay_singletons)),
                "saved_count": int(len(saved_singletons)),
                "coordinates_match_exactly": bool(
                    np.array_equal(replay_singletons, saved_singletons)
                ),
                "matching_coordinates": int(
                    len(set(map(tuple, replay_singletons)) & set(map(tuple, saved_singletons)))
                ),
            }

    report = {
        "evidence_category": "diagnostic replay of saved human core crop; not a full-workflow rerun or source-equivalence claim",
        "core_attempt_dir": str(core_dir),
        "output_dir": str(destination),
        "source_provenance_sha256": _sha256(provenance_path),
        "parameters": asdict(params),
        "crop_bounds_inclusive": [list(axis) for axis in crop.bounds],
        "crop_shape": list(crop.cropped_shape),
        "inputs": input_metadata,
        "outputs": outputs,
        "replay_vs_saved_final_fields": comparisons,
        "evolution_accounting": {
            "smoothed_initial_completed_steps": execution.smoothed_initial.completed_steps,
            "initial_inner_completed_steps": execution.initial_inner.completed_steps,
            "central_completed_steps": execution.central.completed_steps,
            "central_shrink_completed_steps": execution.central_shrink.completed_steps,
            "delivered_inner_completed_steps": execution.inner.completed_steps,
            "outer_completed_steps": execution.outer.completed_steps,
        },
        "peak_process_rss_gib": _max_rss_gib(),
        "runtime_seconds": time.monotonic() - started,
        "python": platform.python_version(),
    }
    report_path = destination / "lineage.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-attempt-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--overwrite", action="store_true", help="replace files in the diagnostic output directory"
    )
    args = parser.parse_args()
    print(run_saved_tgdm_diagnostic(args.core_attempt_dir, args.output_dir, overwrite=args.overwrite))


if __name__ == "__main__":
    main()
