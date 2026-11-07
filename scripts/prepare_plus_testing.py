#!/usr/bin/env python3
"""Prepare a self-contained PLUS tester directory from the bundled atlas."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import sys
import nibabel as nib
import numpy as np

from pymacruise.atlas_resources import get_bundled_atlas_asset, read_bundled_atlas_source
from pymacruise.image_io import is_axis_aligned_ras_grid, normalize_image_to_ras
from adapt_ibsr11objv3_atlas import adapt, _source_assets


def prepare(output: Path, inputs: Path) -> Path:
    output = output.expanduser().resolve()
    inputs = inputs.expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Use a new or empty testing directory: {output}")
    input_records = {}
    for role, filename in (("t1", "T1.nii.gz"), ("flair", "FLAIR.nii.gz"),
                           ("brain_mask", "brain_mask.nii.gz"), ("maseg", "MaSeg.nii.gz")):
        path = inputs / filename
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"Missing regular {role} input: {path}")
        input_records[role] = {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    images = {name: normalize_image_to_ras(nib.load(record["path"]), source_name=name)
              for name, record in input_records.items()}
    if not is_axis_aligned_ras_grid(images["t1"]):
        raise ValueError("This native-T1 comparison kit requires an axis-aligned RAS T1. Prepare an agreed common grid first for oblique T1 cases.")
    for role in ("brain_mask", "maseg"):
        if images[role].shape != images["t1"].shape or not np.allclose(images[role].affine, images["t1"].affine, atol=1e-4, rtol=0):
            raise ValueError(f"{role} must match the T1 physical grid")
    output.mkdir(parents=True, exist_ok=True)
    native = output / "native-atlas"
    native.mkdir()
    inventory = json.loads(get_bundled_atlas_asset("ibsr11objV3/SOURCE_MANIFEST.json").read_text())
    for name in inventory["files"]:
        # The original descriptor references ibsr11objV3/*.raw relative to
        # its parent. Keep the original descriptor bytes and that exact layout.
        destination = native / (name if name.endswith('.txt') else f"ibsr11objV3/{name}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(read_bundled_atlas_source("ibsr11objV3", name))
    descriptor = native / inventory["descriptor"]
    parsed = _source_assets(descriptor)
    # The selected DemonToadDeformableAtlas loader expects consecutive
    # Structure/File pairs. The preserved older descriptor inserts per-map
    # Dimensions/Resolutions/Center rows, causing that loader to stop after
    # Background. Remove only those redundant rows; retain source voxels,
    # global geometry, profiles, structure order and native file references.
    rows = []
    in_shape_record = False
    for line in parsed["lines"]:
        if line.startswith("Intensity Atlas:"):
            in_shape_record = False
        if line.startswith("Structure:"):
            in_shape_record = True
        if in_shape_record and line.startswith(("Dimensions:", "Resolutions:", "Center:")):
            continue
        rows.append(line.rstrip() if line.strip() else " ")
    java_descriptor = native / "ibsr11objV3-atlasem-java.txt"
    java_descriptor.write_text("\n".join(rows) + "\n")
    adapted = adapt(descriptor, output / "python-atlas")
    repo = Path(__file__).resolve().parents[1]
    lesion = get_bundled_atlas_asset("lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt")
    argv = [sys.executable, "-m", "pymacruise.cli", "plus-full-run"]
    for role, flag in (("t1", "--t1"), ("flair", "--flair"), ("brain_mask", "--brain-mask"), ("maseg", "--maseg")):
        argv += [flag, input_records[role]["path"]]
    argv += ["--atlas-descriptor", str(lesion), "--atlasem-descriptor", str(adapted),
             "--output-dir", str(output / "python-run"), "--regrid-inputs",
             "--init-align-iterations", "50", "--pyramid-levels", "4",
             "--lesion-iterations", "99", "--n-iter-max", "99", "--max-passes", "99",
             "--first-thinning-limit", "3", "--loop-thinning-limit", "3",
             "--loop-convergence-threshold", "0.0001", "--dura-transfer-mode", "physical-affine"]
    script = output / "run-python.sh"
    script.write_text("#!/usr/bin/env bash\nset -euo pipefail\n"
        "export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1\n"
        "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1\n"
        f"cd {shlex.quote(str(repo))}\nexec {shlex.join(argv)} \"$@\"\n")
    script.chmod(0o755)
    settings = {
        "schema_version": 1, "inputs": input_records,
        "original_atlasem_descriptor": str(descriptor),
        "source_atlasem_descriptor": str(java_descriptor),
        "source_atlasem_descriptor_sha256": hashlib.sha256(java_descriptor.read_bytes()).hexdigest(),
        "python_atlasem_descriptor": str(adapted),
        "atlasem_adaptation_manifest": str(adapted.parent / "adaptation-manifest.json"),
        "python_command": argv,
        "assumptions": [
            "Selected-mask scientific PLUS profile: user mask replaces BET/SPECTRE; visualization/MNI/FreeSurfer exports excluded.",
            "ibsr11objV3 replaces historical ibsr8obj only at AtlasEM; LesionTOADS uses its separate 2012 atlas.",
            "Preserve FLAIR physical geometry during registration; regrid oblique T1/mask/MaSeg only when necessary.",
            "ANTsPy rigid registration substitutes for source Aladin; optimizer equivalence is unverified.",
            "LesionTOADS initial atlas alignment uses 50 iterations and four pyramid levels, overriding CLI 10/2 convenience defaults.",
            "Physical-affine Dura transfer is the testing default; source JIST transfer differences must be reported.",
            "AtlasEM covariance floor 1e-6 is a Python regularization difference from singular source behavior.",
            "Native descriptor geometry governs conversion; legacy XML headers have inconsistent spacing/orientation.",
            "Java descriptor removes redundant per-map grid/center rows for the selected loader; native raw bytes and global geometry are unchanged. Full runtime remains unexecuted here.",
        ],
    }
    (output / "testing-kit.json").write_text(json.dumps(settings, indent=2) + "\n")
    # Optional final argument of the adapted original driver.
    (output / "source-config.json").write_text(json.dumps({
        "input_dir": str(inputs), "output_dir": str(output / "source-run"),
        "atlasem_descriptor": str(java_descriptor),
    }, indent=2) + "\n")
    return script


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, default=Path("source-reference-run/inputs/100260"))
    args = parser.parse_args()
    try:
        script = prepare(args.output_dir, args.input_dir)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(2, f"PLUS kit preparation failed: {exc}\n")
    print(f"Prepared {script.parent}\nRun: bash {shlex.quote(str(script))}")


if __name__ == "__main__":
    main()
