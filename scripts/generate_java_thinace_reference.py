"""Generate VandyCRUISE ThinACE references from deterministic phantoms.

Requires a JDK. Java binary interchange files are transient under build/;
committed image inputs and reference outputs are compressed NIfTI only.
"""

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
OUTPUT = ROOT / "tests/data/synthetic/skeleton"
BUILD = ROOT / "build/java-reference/thinace"
CLASSES = BUILD / "classes"
HARNESS = ROOT / "tests/java_reference/ThinACEReferenceHarness.java"
SOURCE_ROOT = ROOT / "tests/java_reference"


def phantoms() -> dict[str, np.ndarray]:
    line = np.zeros((15, 15, 15), dtype=np.int32)
    line[2:13, 7, 7] = 255

    cube = np.zeros_like(line)
    cube[2:13, 2:13, 2:13] = 255

    grid25 = np.indices((25, 25, 25), dtype=np.float32)
    radius25 = sum((grid25[axis] - 12.0) ** 2 for axis in range(3))
    spheres = {
        f"sphere_r{radius}": np.where(
            sum((grid25[axis] - 12.0) ** 2 for axis in range(3)) <= radius**2,
            255,
            0,
        ).astype(np.int32)
        for radius in (7, 8, 9, 10, 11)
    }

    grid35 = np.indices((35, 35, 35), dtype=np.float32)
    spheres.update(
        {
            f"sphere_r{radius}": np.where(
                sum((grid35[axis] - 17.0) ** 2 for axis in range(3)) <= radius**2,
                255,
                0,
            ).astype(np.int32)
            for radius in (12, 14, 16)
        }
    )

    grid33 = np.indices((33, 33, 33), dtype=np.float32)
    center33 = (16.0, 16.0, 16.0)
    radii33 = (11.0, 9.0, 7.0)
    radius33 = sum(
        ((grid33[axis] - center33[axis]) / radii33[axis]) ** 2
        for axis in range(3)
    )
    nested = np.zeros((33, 33, 33), dtype=np.int32)
    nested[(radius33 <= 1.0) & (radius33 > 0.4)] = 128  # gray matter shell
    nested[radius33 <= 0.4] = 255  # white matter core
    return {"line": line, "cube": cube, **spheres, "nested_gm_wm": nested}


def write_binary(path: Path, values: np.ndarray) -> None:
    with path.open("wb") as stream:
        stream.write(struct.pack(">III", *values.shape))
        stream.write(values.astype(">i4", copy=False).tobytes(order="C"))


def read_binary(path: Path, shape: tuple[int, int, int]) -> np.ndarray:
    raw = path.read_bytes()
    found = struct.unpack(">III", raw[:12])
    if found != shape:
        raise ValueError(f"Java output shape {found} does not match input shape {shape}")
    values = np.frombuffer(raw, dtype=">i4", offset=12)
    if values.size != int(np.prod(shape)):
        raise ValueError(f"Java output has {values.size} voxels; expected {np.prod(shape)}")
    return values.astype(np.int32).reshape(shape)


def main() -> None:
    javac = shutil.which("javac")
    java = shutil.which("java")
    if not javac or not java:
        raise RuntimeError("A JDK with javac and java is required to regenerate Java references")
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
            str(HARNESS),
        ],
        check=True,
        cwd=ROOT,
    )

    affine = np.eye(4)
    metadata: dict[str, object] = {
        "purpose": "VandyCRUISE ThinACE synthetic parity reference",
        "provenance": "synthetic; generated from analytic masks",
        "frame": "RAS",
        "units": "mm",
        "voxel_spacing_mm": [1, 1, 1],
        "cases": {},
    }
    for name, values in phantoms().items():
        input_image = nib.Nifti1Image(values.astype(np.uint8), affine)
        input_path = OUTPUT / f"{name}_input.nii.gz"
        reference_path = OUTPUT / f"{name}_thinace_java.nii.gz"
        nib.save(input_image, input_path)
        binary_input = BUILD / f"{name}_input.bin"
        binary_output = BUILD / f"{name}_output.bin"
        binary_distance = BUILD / f"{name}_distance.bin"
        write_binary(binary_input, values)
        subprocess.run(
            [
                java,
                "-cp",
                str(CLASSES),
                "ThinACEReferenceHarness",
                str(binary_input),
                str(binary_output),
                str(binary_distance),
            ],
            check=True,
            cwd=ROOT,
        )
        reference = read_binary(binary_output, values.shape)
        nib.save(nib.Nifti1Image(reference, affine), reference_path)
        metadata["cases"][name] = {
            "shape": list(values.shape),
            "input_values": [int(value) for value in np.unique(values)],
            "reference_values": [int(value) for value in np.unique(reference)],
            "input": input_path.name,
            "java_reference": reference_path.name,
        }

    (OUTPUT / "spec.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    filenames = sorted(path.name for path in OUTPUT.iterdir() if path.is_file() and path.name != "SHA256SUMS")
    lines = [
        f"{hashlib.sha256((OUTPUT / filename).read_bytes()).hexdigest()}  {filename}"
        for filename in filenames
    ]
    (OUTPUT / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="ascii")


if __name__ == "__main__":
    main()
