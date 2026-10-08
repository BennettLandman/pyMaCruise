#!/usr/bin/env python3
"""Compare bundled-atlas TOADS startup registration with selected Java source."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_alignment_startup import (
    prepare_toads_alignment_startup,
    register_toads_alignment_startup,
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

    report = {
        "scope": "selected actual VandyCRUISE rigid pyramid on bundled atlas assets and synthetic T1/mask",
        "source": java_result["source"],
        "t1_sha256": _sha256(args.t1),
        "brain_mask_sha256": _sha256(args.brain_mask),
        "atlas_descriptor_sha256": _sha256(args.atlas),
        "membership_shape": list(startup.initial_memberships.shape),
        "shape_prior_shape": list(shapes.shape),
        "class_count": int(shapes.shape[-1]),
        "levels": [
            {
                "block_scale": level.scale,
                "finite_transform": bool(np.isfinite(level.result.transform).all()),
                "accepted_updates": sum(item.accepted for item in level.result.iterations),
                "line_search_trials": sum(
                    item.rejected_trials + int(item.accepted)
                    for item in level.result.iterations
                ),
            }
            for level in python_result.registration.levels
        ],
        "python_transform_float32_bits": python_bits,
        "java_transform_float32_bits": java_result["transform_bits"],
        "transform_bit_exact": transform_match,
        "max_iterations_per_level": args.max_iterations,
        "python_source_faithful_nonfinite_trials": False,
    }
    report_path = output_dir / "comparison.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"report={report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
