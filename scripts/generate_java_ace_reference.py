"""Generate selected MaCRUISE ACE outputs from the copied Java sources."""

from __future__ import annotations

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
CLASSES = BUILD / "classes"
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
    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        raise RuntimeError("A JDK with javac and java is required to generate ACE references")
    BUILD.mkdir(parents=True, exist_ok=True)
    CLASSES.mkdir(parents=True, exist_ok=True)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            javac,
            "-d",
            str(CLASSES),
            "-sourcepath",
            f"{SOURCE_ROOT / 'stubs'}:{SOURCE_ROOT / 'upstream'}",
            str(ROOT / "tests/java_reference/ACEReferenceHarness.java"),
        ],
        check=True,
        cwd=ROOT,
    )
    names = ("gray_matter.nii.gz", "white_matter.nii.gz", "multi_atlas_csf.nii.gz")
    images = [nib.load(INPUTS / name) for name in names]
    shape = write_inputs(BUILD / "ace_input.bin", images)
    subprocess.run(
        [
            java,
            "-cp",
            str(CLASSES),
            "ACEReferenceHarness",
            str(BUILD / "ace_input.bin"),
            str(BUILD / "ace_output.bin"),
            str(BUILD / "intermediate"),
        ],
        check=True,
        cwd=ROOT,
    )
    enhanced, skeleton, thinned = read_output(BUILD / "ace_output.bin", shape)
    reference_paths = {
        "enhanced_gray_matter": OUTPUT / "enhanced_gm_java.nii.gz",
        "cleaned_skeleton": OUTPUT / "csf_skeleton_java.nii.gz",
        "thinned_skeleton": OUTPUT / "thinned_skeleton_java.nii.gz",
    }
    arrays = {"enhanced_gray_matter": enhanced, "cleaned_skeleton": skeleton, "thinned_skeleton": thinned}
    for key, path in reference_paths.items():
        dtype = np.int32 if key != "enhanced_gray_matter" else np.float32
        nib.save(nib.Nifti1Image(arrays[key].astype(dtype), images[0].affine), path)

    stage_output = OUTPUT / "stages"
    stage_output.mkdir(parents=True, exist_ok=True)
    stage_names = []
    for binary_path in sorted((BUILD / "intermediate").glob("*.bin")):
        stage = read_intermediate(binary_path, shape)
        stage_name = binary_path.stem
        stage_path = stage_output / f"{stage_name}_java.nii.gz"
        nib.save(nib.Nifti1Image(stage, images[0].affine), stage_path)
        stage_names.append(stage_path.relative_to(OUTPUT).as_posix())

    metadata = {
        "purpose": "selected MaCRUISE ACE wrapper Java synthetic reference",
        "provenance": "synthetic; committed concentric tissue phantom",
        "frame": "RAS",
        "units": "mm",
        "shape": list(shape),
        "input_files": list(names),
        "parameters": {
            "intensity_threshold": 126.9,
            "gaussian_alpha": 1.8,
            "tensor_alpha": 1.6,
            "diffusion_iterations": 25,
            "multi_atlas_csf_place": "AfterDiffuse",
            "multi_atlas_csf_merge": "Merge",
            "narrow_band": {"delta_t": 0.1, "outer_iterations": 1, "inner_updates": 10, "curvature_threshold": 5.0},
        },
        "outputs": {key: path.name for key, path in reference_paths.items()},
        "java_intermediate_stages": stage_names,
    }
    (OUTPUT / "spec.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    files = sorted(path for path in OUTPUT.rglob("*") if path.is_file() and path.name != "SHA256SUMS")
    (OUTPUT / "SHA256SUMS").write_text(
        "\n".join(
            f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(OUTPUT).as_posix()}"
            for path in files
        )
        + "\n",
        encoding="ascii",
    )


if __name__ == "__main__":
    main()
