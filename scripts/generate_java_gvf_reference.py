"""Generate small synthetic references with the copied VandyCRUISE Java GVF."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import struct
import subprocess

import nibabel as nib
import numpy as np
from pymacruise.gvf import generalized_gvf


ROOT = Path(__file__).parents[1]
OUTPUT = ROOT / "tests/data/synthetic/gvf/java_reference"
BUILD = ROOT / "build/java-reference/gvf"
CLASSES = BUILD / "classes"
SOURCE_ROOT = ROOT / "tests/java_reference"


def phantoms() -> dict[str, np.ndarray]:
    axes = np.indices((16, 12, 10), dtype=np.float32)
    x, y, z = axes
    ramp = 0.2 * x + 0.7 * y + 0.4 * z
    step = ((x >= 7) & (y + z >= 9)).astype(np.float32)
    gaussian = np.exp(-(((x - 7.0) / 3.4) ** 2 + ((y - 5.0) / 2.7) ** 2 + ((z - 4.0) / 2.1) ** 2))
    asymmetric = 0.6 * gaussian + 0.4 * np.exp(-(((x - 11.0) / 2.0) ** 2 + ((y - 3.0) / 3.5) ** 2 + ((z - 6.0) / 1.7) ** 2))
    return {
        "ramp_xyz": ramp.astype(np.float32),
        "oblique_step": step.astype(np.float32),
        "asymmetric_blobs": asymmetric.astype(np.float32),
    }


def write_binary(path: Path, values: np.ndarray) -> None:
    with path.open("wb") as stream:
        stream.write(struct.pack(">III", *values.shape))
        stream.write(values.astype(">f4", copy=False).tobytes(order="C"))


def read_binary(path: Path, shape: tuple[int, int, int]) -> np.ndarray:
    data = path.read_bytes()
    found = struct.unpack(">III", data[:12])
    if found != shape:
        raise ValueError(f"Java output shape {found} does not match input {shape}")
    values = np.frombuffer(data, dtype=">f4", offset=12)
    if values.size != int(np.prod(shape)) * 3:
        raise ValueError("Java vector output has an unexpected length")
    return values.astype(np.float32).reshape((*shape, 3))


def main() -> None:
    javac, java = shutil.which("javac"), shutil.which("java")
    if not javac or not java:
        raise RuntimeError("A JDK with javac and java is required to generate Java references")
    BUILD.mkdir(parents=True, exist_ok=True)
    CLASSES.mkdir(parents=True, exist_ok=True)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    sources = [
        str(path)
        for path in sorted((SOURCE_ROOT / "upstream/edu/jhu/ece/iacl/algorithms/gvf").glob("*.java"))
    ]
    sources.append(str(ROOT / "tests/java_reference/GVFReferenceHarness.java"))
    subprocess.run(
        [javac, "-d", str(CLASSES), "-sourcepath", f"{SOURCE_ROOT / 'stubs'}:{SOURCE_ROOT / 'upstream'}", *sources],
        check=True,
        cwd=ROOT,
    )

    affine = np.eye(4)
    metadata: dict[str, object] = {
        "purpose": "VandyCRUISE generalized GVF synthetic parity references",
        "provenance": "synthetic; analytic scalar fields",
        "frame": "RAS",
        "units": "voxel-axis component values; 1 mm isotropic fixtures",
        "lambda": 1.0,
        "iterations": 200,
        "normalize_vectors": False,
        "cases": {},
    }
    for name, values in phantoms().items():
        input_path = OUTPUT / f"{name}_input.nii.gz"
        reference_path = OUTPUT / f"{name}_gvf_java.nii.gz"
        nib.save(nib.Nifti1Image(values, affine), input_path)
        binary_input, binary_output = BUILD / f"{name}_input.bin", BUILD / f"{name}_output.bin"
        write_binary(binary_input, values)
        subprocess.run(
            [java, "-cp", str(CLASSES), "GVFReferenceHarness", str(binary_input), str(binary_output), "1.0", "200"],
            check=True,
            cwd=ROOT,
        )
        reference = read_binary(binary_output, values.shape)
        candidate = generalized_gvf(values, lambda_value=1.0, max_iterations=200)
        absolute_error = np.abs(candidate - reference)
        header = nib.Nifti1Header()
        header.set_data_dtype(np.float32)
        header.set_intent("vector")
        nib.save(nib.Nifti1Image(reference, affine, header), reference_path)
        metadata["cases"][name] = {
            "shape": list(values.shape),
            "input": input_path.name,
            "java_reference": reference_path.name,
            "python_vs_java": {
                "mean_absolute_error": float(absolute_error.mean()),
                "maximum_absolute_error": float(absolute_error.max()),
                "root_mean_square_error": float(np.sqrt(np.mean(absolute_error**2))),
            },
        }

    (OUTPUT / "spec.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    files = sorted(path for path in OUTPUT.iterdir() if path.is_file() and path.name != "SHA256SUMS")
    (OUTPUT / "SHA256SUMS").write_text(
        "\n".join(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}" for path in files) + "\n",
        encoding="ascii",
    )


if __name__ == "__main__":
    main()
