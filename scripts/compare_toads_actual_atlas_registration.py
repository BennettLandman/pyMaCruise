#!/usr/bin/env python3
"""Compare bundled-atlas TOADS startup registration with selected Java source."""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_alignment_startup import (
    prepare_toads_alignment_startup,
    register_toads_alignment_startup,
    register_toads_from_atlas_assets,
)
from pymacruise.toads_atlas import (
    load_toads_atlas,
    load_toads_atlas_data,
    load_toads_intensity_priors,
    toads_atlas_shapes_in_structure_order,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FIXTURE = ROOT / "tests/data/synthetic/base_workflow_brain_resolution"
DEFAULT_ATLAS = ROOT / "data/atlases/cruise-atlas-10obj-toads2008.txt"
JAVA_HARNESS = ROOT / "tests/java_reference/run_toads_pyramid_actual_registration_reference.sh"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t1", type=Path, default=DEFAULT_FIXTURE / "t1.nii.gz")
    parser.add_argument("--brain-mask", type=Path, default=DEFAULT_FIXTURE / "brain-mask.nii.gz")
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "build/toads-actual-atlas-registration",
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=1,
        help="registration iterations per level; 1 gives a bounded source-comparison slice",
    )
    parser.add_argument(
        "--verify-full-startup",
        action="store_true",
        help="also run the bundled-atlas wrapper through post-registration constructor state",
    )
    args = parser.parse_args()
    if args.max_iterations < 1:
        parser.error("--max-iterations must be positive")

    for path in (args.t1, args.brain_mask, args.atlas):
        if not path.is_file():
            parser.error(f"input does not exist: {path}")
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    atlas = load_toads_atlas(args.atlas)
    atlas_data = load_toads_atlas_data(atlas)
    intensity_priors = load_toads_intensity_priors(
        atlas.descriptor, ["T1_MPRAGE"]
    )
    mask_indices = [
        index
        for index, structure in enumerate(atlas.structures)
        if structure.topology.lower() == "mask"
    ]
    if len(mask_indices) != 1:
        raise ValueError("atlas descriptor must identify one mask/background structure")
    background_label = atlas.structures[mask_indices[0]].label
    t1 = nib.load(args.t1)
    brain_mask = nib.load(args.brain_mask)
    startup = prepare_toads_alignment_startup(
        t1,
        brain_mask,
        atlas_data.topology,
        intensity_priors,
        atlas.shape_voxel_sizes_mm,
        padding=10,
        background_label=background_label,
    )
    shapes = toads_atlas_shapes_in_structure_order(atlas, atlas_data)
    membership_path = output_dir / "memberships-f32le.bin"
    shapes_path = output_dir / "atlas-shapes-f32le.bin"
    startup.initial_memberships.tofile(membership_path)
    shapes.tofile(shapes_path)

    python_result = register_toads_alignment_startup(
        startup,
        shapes,
        levels=4,
        max_iterations=args.max_iterations,
        allow_nonfinite_trials=False,
    )
    image_spacing = startup.geometry.signed_image_resolutions_mm
    template_spacing = atlas.shape_voxel_sizes_mm
    harness_command = [
        str(JAVA_HARNESS),
        str(membership_path),
        str(shapes_path),
        *(str(value) for value in startup.initial_memberships.shape[:3]),
        *(str(value) for value in shapes.shape[:3]),
        str(shapes.shape[-1]),
        *(str(float(value)) for value in image_spacing),
        *(str(float(value)) for value in template_spacing),
        "4",
        str(args.max_iterations),
        *(str(float(value)) for value in startup.transform),
    ]
    java = subprocess.run(
        harness_command, check=True, text=True, capture_output=True
    )
    (output_dir / "java-harness.log").write_text(java.stdout, encoding="utf-8")
    java_result = json.loads(
        next(line for line in java.stdout.splitlines() if line.startswith("{"))
    )
    python_bits = python_result.registration.transform.view(np.int32).tolist()
    transform_match = python_bits == java_result["transform_bits"]
    if not transform_match:
        raise AssertionError(
            "selected Java/Python final transform mismatch: "
            f"python={python_bits}, java={java_result['transform_bits']}"
        )

    source_level_numbers = (4, 3, 2)
    level_reports = []
    for python_level, source_level in zip(
        python_result.registration.levels, source_level_numbers, strict=True
    ):
        start = java.stdout.index(f"level {source_level}\n")
        end = java.stdout.find(f"level {source_level - 1}\n", start + 1)
        section = java.stdout[start:] if end < 0 else java.stdout[start:end]
        python_trials = sum(
            item.rejected_trials + int(item.accepted)
            for item in python_level.result.iterations
        )
        python_accepted = sum(
            item.accepted for item in python_level.result.iterations
        )
        java_trials = section.count("trial transform:")
        java_accepted = section.count(" (E:")
        if (python_trials, python_accepted) != (java_trials, java_accepted):
            raise AssertionError(
                f"scale {python_level.scale} line-search trace mismatch: "
                f"python={(python_trials, python_accepted)}, "
                f"java={(java_trials, java_accepted)}"
            )
        level_reports.append(
            {
                "block_scale": python_level.scale,
                "finite_transform": bool(
                    np.isfinite(python_level.result.transform).all()
                ),
                "iterations": len(python_level.result.iterations),
                "accepted_updates": python_accepted,
                "line_search_trials": python_trials,
                "java_accepted_updates": java_accepted,
                "java_line_search_trials": java_trials,
                "terminal_relative_change": (
                    float(python_level.result.iterations[-1].relative_change)
                    if python_level.result.iterations
                    else None
                ),
                "stop_reason": python_level.result.stop_reason,
                "iteration_cap_reached": (
                    python_level.result.stop_reason == "iteration_limit"
                ),
            }
        )

    report = {
        "scope": "selected actual VandyCRUISE rigid pyramid on bundled atlas assets and synthetic T1/mask",
        "source": java_result["source"],
        "t1_sha256": _sha256(args.t1),
        "brain_mask_sha256": _sha256(args.brain_mask),
        "atlas_descriptor_sha256": _sha256(args.atlas),
        "membership_shape": list(startup.initial_memberships.shape),
        "shape_prior_shape": list(shapes.shape),
        "class_count": int(shapes.shape[-1]),
        "levels": level_reports,
        "python_transform_float32_bits": python_bits,
        "java_transform_float32_bits": java_result["transform_bits"],
        "transform_bit_exact": transform_match,
        "max_iterations_per_level": args.max_iterations,
        "python_source_faithful_nonfinite_trials": False,
    }
    if args.verify_full_startup:
        expected_transform_bits = python_bits
        del startup, atlas_data, intensity_priors, shapes, python_result
        gc.collect()
        stages: list[str] = []
        full_startup = register_toads_from_atlas_assets(
            t1,
            brain_mask,
            atlas.descriptor,
            levels=4,
            max_iterations=args.max_iterations,
            padding=10,
            background_ratio=0.0,
            refinement_iterations=1,
            relation_limit=3,
            progress_callback=stages.append,
        )
        full_transform_bits = (
            full_startup.alignment.registration.transform.view(np.int32).tolist()
        )
        if full_transform_bits != expected_transform_bits:
            raise AssertionError(
                "full wrapper startup changed the initial registration transform"
            )
        scientific_arrays = {
            "initial_memberships": full_startup.segmentation_initialization.memberships,
            "transformed_shape_priors": full_startup.transformed_shape_priors,
            "class_standard_deviations": full_startup.initial_variances.standard_deviations,
            "image_variance_weights": full_startup.initial_image_variance.image_variance_weights,
            "centroid_similarity": full_startup.initial_centroid_similarity.similarity,
            "relation_map": full_startup.relation_state.relation.relation_map,
            "boundary_distances": full_startup.relation_state.boundary_distances,
            "closest_labels": full_startup.relation_state.closest_labels,
        }
        nonfinite_arrays = [
            name
            for name, values in scientific_arrays.items()
            if not np.isfinite(np.asarray(values)).all()
        ]
        if nonfinite_arrays:
            raise AssertionError(
                f"nonfinite selected startup products: {nonfinite_arrays}"
            )
        report["full_startup"] = {
            "stage_callbacks": stages,
            "initial_registration_transform_bit_exact": True,
            "post_refresh_transform_finite": bool(
                np.isfinite(full_startup.post_refresh_registration.final_transform).all()
            ),
            "constructor_segmentation_shape": list(
                full_startup.segmentation_constructor.segmentation.shape
            ),
            "finite_startup_products": sorted(scientific_arrays),
        }
    report_path = output_dir / "comparison.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"report={report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
