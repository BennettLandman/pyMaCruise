"""Generate selected MaCRUISE ACE outputs from the copied Java sources."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import struct
import subprocess

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parents[1]
INPUTS = ROOT / "tests/data/synthetic/ace"
OUTPUT = INPUTS / "java_reference"
BUILD = ROOT / "build/java-reference/ace"
SOURCE_ROOT = ROOT / "tests/java_reference"


def write_inputs(path: Path, images: list[nib.spatialimages.SpatialImage]) -> tuple[int, int, int]:
    arrays = [np.asarray(image.dataobj, dtype=np.float32) for image in images]
    if any(array.shape != arrays[0].shape for array in arrays):
        raise ValueError("ACE fixture images must share a shape")
    with path.open("wb") as stream:
        stream.write(struct.pack(">III", *arrays[0].shape))
        for array in arrays:
            stream.write(array.astype(">f4", copy=False).tobytes(order="C"))
    return arrays[0].shape


def read_output(path: Path, shape: tuple[int, int, int]) -> tuple[np.ndarray, ...]:
    payload = path.read_bytes()
    outputs = []
    offset = 0
    for _ in range(3):
        found = struct.unpack(">III", payload[offset : offset + 12])
        offset += 12
        if found != shape:
            raise ValueError(f"ACE output shape {found} does not match {shape}")
        count = int(np.prod(shape))
        values = np.frombuffer(payload, dtype=">f4", count=count, offset=offset).astype(np.float32)
        if values.size != count:
            raise ValueError("ACE Java output is truncated")
        outputs.append(values.reshape(shape))
        offset += count * 4
    if offset != len(payload):
        raise ValueError(f"ACE Java output has {len(payload) - offset} trailing bytes")
    return tuple(outputs)


def read_intermediate(path: Path, shape: tuple[int, int, int]) -> np.ndarray:
    payload = path.read_bytes()
    found = struct.unpack(">III", payload[:12])
    if found != shape:
        raise ValueError(f"ACE intermediate {path.name} has shape {found}, expected {shape}")
    values = np.frombuffer(payload, dtype=">f4", offset=12)
    if values.size != int(np.prod(shape)):
        raise ValueError(f"ACE intermediate {path.name} is truncated")
    return values.astype(np.float32).reshape(shape)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--place-to-add-csf", type=int, choices=(0, 1, 2, 3), default=2)
    parser.add_argument(
        "--diffusion-iterations", type=int, default=25,
        help="selected Java AnisotropicDiffusion iterations (default: source value 25)",
    )
    parser.add_argument("--input-root", type=Path, default=INPUTS,
                        help="directory containing gray_matter, white_matter, and multi_atlas_csf NIfTIs")
    parser.add_argument("--output-dir", type=Path,
                        help="reference output directory (defaults to INPUT_ROOT/java_reference)")
    args = parser.parse_args()
    if args.diffusion_iterations < 0:
        parser.error("--diffusion-iterations must be non-negative")
    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        raise RuntimeError("A JDK with javac and java is required to generate ACE references")
    build = BUILD / f"place-{args.place_to_add_csf}"
    classes = build / "classes"
    intermediate = build / "intermediate"
    input_root = args.input_root.resolve()
    output = (args.output_dir or (input_root / "java_reference")).resolve()
    build.mkdir(parents=True, exist_ok=True)
    classes.mkdir(parents=True, exist_ok=True)
    intermediate.mkdir(parents=True, exist_ok=True)
    source_root = SOURCE_ROOT / "upstream"
    source_override = None
    original_ace_source = (
        source_root
        / "edu/vanderbilt/masi/algorithms/CRUISE/ace/MaAnatomicallyConsistentEnhance.java"
    )
    if args.diffusion_iterations != 25:
        source_text = original_ace_source.read_text(encoding="utf-8")
        original_setting = "private static final int max_iter = 25;"
        if source_text.count(original_setting) != 1:
            raise RuntimeError("selected Java ACE max_iter declaration changed; review source mapping")
        source_override = build / "source_override"
        override_path = (
            source_override
            / "edu/vanderbilt/masi/algorithms/CRUISE/ace/MaAnatomicallyConsistentEnhance.java"
        )
        override_path.parent.mkdir(parents=True, exist_ok=True)
        override_path.write_text(
            source_text.replace(
                original_setting,
                f"private static final int max_iter = {args.diffusion_iterations};",
            ),
            encoding="utf-8",
        )
    output.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            javac,
            "-d",
            str(classes),
            "-sourcepath",
            f"{SOURCE_ROOT / 'stubs'}:{source_override or source_root}:{source_root}",
            str(ROOT / "tests/java_reference/ACEReferenceHarness.java"),
        ],
        check=True,
        cwd=ROOT,
    )
    names = ("gray_matter.nii.gz", "white_matter.nii.gz", "multi_atlas_csf.nii.gz")
    images = [nib.load(input_root / name) for name in names]
    shape = write_inputs(build / "ace_input.bin", images)
    subprocess.run(
        [
            java,
            "-Xmx2g",
            "-cp",
            str(classes),
            "ACEReferenceHarness",
            str(build / "ace_input.bin"),
            str(build / "ace_output.bin"),
            str(intermediate),
            str(args.place_to_add_csf),
        ],
        check=True,
        cwd=ROOT,
    )
    enhanced, skeleton, thinned = read_output(build / "ace_output.bin", shape)
    reference_paths = {
        "enhanced_gray_matter": output / "enhanced_gm_java.nii.gz",
        "cleaned_skeleton": output / "csf_skeleton_java.nii.gz",
        "thinned_skeleton": output / "thinned_skeleton_java.nii.gz",
    }
    arrays = {"enhanced_gray_matter": enhanced, "cleaned_skeleton": skeleton, "thinned_skeleton": thinned}
    for key, path in reference_paths.items():
        dtype = np.int32 if key != "enhanced_gray_matter" else np.float32
        nib.save(nib.Nifti1Image(arrays[key].astype(dtype), images[0].affine), path)

    stage_output = output / "stages"
    stage_output.mkdir(parents=True, exist_ok=True)
    stage_names = []
    for binary_path in sorted(intermediate.glob("*.bin")):
        stage = read_intermediate(binary_path, shape)
        stage_name = binary_path.stem
        stage_path = stage_output / f"{stage_name}_java.nii.gz"
        nib.save(nib.Nifti1Image(stage, images[0].affine), stage_path)
        stage_names.append(stage_path.relative_to(output).as_posix())

    metadata = {
        "purpose": (
            "selected MaACE Java wrapper after-Gaussian synthetic reference"
            if args.place_to_add_csf == 3
            else "selected MaCRUISE ACE wrapper Java synthetic reference"
        ),
        "provenance": "synthetic; committed analytic fixture; no subject data",
        "frame": "RAS",
        "units": "mm",
        "shape": list(shape),
        "input_files": list(names),
        "parameters": {
            "intensity_threshold": 126.9,
            "gaussian_alpha": 1.8,
            "tensor_alpha": 1.6,
            "diffusion_iterations": args.diffusion_iterations,
            "multi_atlas_csf_place": {
                0: "None",
                1: "BeforeDiffuse",
                2: "AfterDiffuse",
                3: "AfterGaussian",
            }[args.place_to_add_csf],
            "multi_atlas_csf_merge": "Merge",
            "narrow_band": {"delta_t": 0.1, "outer_iterations": 1, "inner_updates": 10, "curvature_threshold": 5.0},
        },
        "outputs": {key: path.name for key, path in reference_paths.items()},
        "java_intermediate_stages": stage_names,
    }
    if source_override is not None:
        metadata["java_reference_build_note"] = (
            "Build-only copy of MaAnatomicallyConsistentEnhance.java changes only "
            f"max_iter from 25 to {args.diffusion_iterations}; original upstream source is unchanged."
        )
    (output / "spec.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    files = sorted(path for path in output.rglob("*") if path.is_file() and path.name != "SHA256SUMS")
    (output / "SHA256SUMS").write_text(
        "\n".join(
            f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(output).as_posix()}"
            for path in files
        )
        + "\n",
        encoding="ascii",
    )


if __name__ == "__main__":
    main()
