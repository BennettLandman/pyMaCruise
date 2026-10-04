from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_alignment_startup import register_toads_from_atlas_assets
from pymacruise.toads_atlas import (
    load_toads_atlas,
    toads_atlas_shapes_in_structure_order,
)


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "data/synthetic/toads_downsampled_vandy_atlas"
PROJECT = ROOT.parent
ATLAS_FILE = ROOT.parent / "data/atlases/cruise-atlas-10obj-toads2008.txt"
OUTPUT = FIXTURE / "post_refresh_registration_java_reference.json"
BUILD = PROJECT / "build/java-reference/toads-vandy-post-refresh-inputs"


def _write_f32le(path: Path, array: np.ndarray) -> None:
    np.asarray(array, dtype="<f4").tofile(path)


def main() -> None:
    atlas = load_toads_atlas(ATLAS_FILE)
    alignment = register_toads_from_atlas_assets(
        nib.load(str(FIXTURE / "t1.nii.gz")),
        nib.load(str(FIXTURE / "brain-mask.nii.gz")),
        atlas,
        levels=1,
        max_iterations=1,
        refinement_iterations=1,
        modality="T1_MPRAGE",
        padding=10,
        relation_limit=3,
        energy_tolerance=1e-6,
        minimum_lambda=1e-3,
    )
    with np.load(FIXTURE / "two_iteration_java_reference.npz", allow_pickle=False) as reference:
        memberships = np.asarray(reference["iteration1_memberships"], dtype=np.float32)
    shapes = toads_atlas_shapes_in_structure_order(atlas, alignment.atlas_data)
    BUILD.mkdir(parents=True, exist_ok=True)
    membership_path = BUILD / "memberships-f32le.bin"
    shapes_path = BUILD / "shapes-f32le.bin"
    _write_f32le(membership_path, memberships)
    _write_f32le(shapes_path, shapes)
    dimensions = memberships.shape[:3]
    template_dimensions = shapes.shape[:3]
    image_res = alignment.alignment.setup.signed_image_resolutions_mm
    template_res = alignment.alignment.setup.template_resolutions_mm
    transform = alignment.template_refresh.recropped_transform
    args = [
        str(membership_path), str(shapes_path),
        *(str(int(v)) for v in dimensions),
        *(str(int(v)) for v in template_dimensions),
        *(str(float(v)) for v in image_res),
        *(str(float(v)) for v in template_res),
        *(str(float(v)) for v in transform),
    ]
    reference_root = ROOT / "java_reference"
    result = subprocess.run(
        [str(reference_root / "generate_toads_vandy_post_refresh_reference.sh"), *args],
        check=True,
        capture_output=True,
        text=True,
        env=os.environ.copy(),
    )
    output_line = next(line for line in reversed(result.stdout.splitlines()) if line.startswith("{"))
    reference = json.loads(output_line)
    harness = reference_root / "ToadsVandyPostRefreshRegistrationHarness.java"
    reference["harness_sha256"] = hashlib.sha256(harness.read_bytes()).hexdigest()
    OUTPUT.write_text(json.dumps(reference, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
