"""Compare Python CG and selected VandyCRUISE GVF on a downsampled NIfTI.

This is a bounded diagnostic for larger scalar volumes, not a substitute for
the committed analytic GVF references or a full-resolution Java run.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import struct
import subprocess
import warnings

import nibabel as nib
import numpy as np

from pymacruise.gvf import generalized_gvf
from pymacruise.image_io import normalize_image_to_ras


ROOT = Path(__file__).resolve().parents[1]
JAVA_SOURCE = ROOT / "tests/java_reference"
UPSTREAM = JAVA_SOURCE / "upstream"


def _write_java_input(path: Path, values: np.ndarray) -> None:
    with path.open("wb") as stream:
        stream.write(struct.pack(">III", *values.shape))
        stream.write(values.astype(">f4", copy=False).tobytes(order="C"))


def _read_java_output(path: Path, shape: tuple[int, int, int]) -> np.ndarray:
    data = path.read_bytes()
    found = struct.unpack(">III", data[:12])
    if found != shape:
        raise ValueError(f"Java output shape {found} does not match input {shape}")
    values = np.frombuffer(data, dtype=">f4", offset=12)
    if values.size != int(np.prod(shape)) * 3:
        raise ValueError("Java vector output has an unexpected length")
    return values.astype(np.float32).reshape((*shape, 3))


def _compile_harness(classes: Path) -> None:
    javac = shutil.which("javac")
    if javac is None:
        raise RuntimeError("A JDK with javac and java is required for this comparison")
    classes.mkdir(parents=True, exist_ok=True)
    sources = sorted(
        path
        for path in (UPSTREAM / "edu/jhu/ece/iacl/algorithms/gvf").glob("*.java")
        if path.name != "FastMarchingGradient.java"
    )
    sources.append(JAVA_SOURCE / "GVFReferenceHarness.java")
    subprocess.run(
        [
            javac,
            "-d",
            str(classes),
            "-sourcepath",
            f"{JAVA_SOURCE / 'stubs'}:{UPSTREAM}",
            *(str(path) for path in sources),
        ],
        check=True,
        cwd=ROOT,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="RAS scalar NIfTI, for example ACE enhanced GM")
    parser.add_argument("--output-dir", type=Path, required=True, help="Ignored build directory for temporary harness files")
    parser.add_argument("--factor", type=int, default=4, help="Stride for each voxel axis (default: 4)")
    parser.add_argument("--iterations", type=int, default=200, help="Fixed Java V-cycle and Python CG cap")
    parser.add_argument("--lambda-value", type=float, default=1.0)
    parser.add_argument("--solver-tolerance", type=float, default=1e-6)
    args = parser.parse_args()
    if args.factor < 1 or args.iterations < 1:
        parser.error("--factor and --iterations must be positive")

    image = nib.load(str(args.input))
    ras = normalize_image_to_ras(image, source_name="GVF probe input")
    if len(ras.shape) != 3:
        raise ValueError("GVF probe input must be a 3D scalar NIfTI")
    values = np.asarray(ras.dataobj, dtype=np.float32)[:: args.factor, :: args.factor, :: args.factor].copy()
    if any(length < 2 for length in values.shape):
        raise ValueError(f"Downsampled dimensions must be at least two, found {values.shape}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    classes = args.output_dir / "classes"
    _compile_harness(classes)
    binary_input = args.output_dir / "gvf-input.bin"
    binary_output = args.output_dir / "gvf-java.bin"
    _write_java_input(binary_input, values)
    java = shutil.which("java")
    if java is None:
        raise RuntimeError("A JDK with java is required for this comparison")
    subprocess.run(
        [
            java,
            "-cp",
            str(classes),
            "GVFReferenceHarness",
            str(binary_input.resolve()),
            str(binary_output.resolve()),
            str(args.lambda_value),
            str(args.iterations),
        ],
        check=True,
        cwd=ROOT,
    )
    java_result = _read_java_output(binary_output, values.shape)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        python_result = generalized_gvf(
            values,
            lambda_value=args.lambda_value,
            max_iterations=args.iterations,
            solver_tolerance=args.solver_tolerance,
        )

    metrics = []
    for component in range(3):
        difference = np.abs(python_result[..., component] - java_result[..., component])
        metrics.append(
            {
                "component": component,
                "mae": float(difference.mean()),
                "max_abs_error": float(difference.max()),
                "rmse": float(np.sqrt(np.mean(difference.astype(np.float64) ** 2))),
            }
        )
    result = {
        "input": str(args.input),
        "ras_shape": list(ras.shape),
        "factor": args.factor,
        "probe_shape": list(values.shape),
        "lambda": args.lambda_value,
        "iterations": args.iterations,
        "python_solver": "VandyCRUISE multigrid V-cycle (Numba port)",
        "python_residual_threshold": args.solver_tolerance,
        "java_solver": "VandyCRUISE fixed-count multigrid V-cycles",
        "python_solver_warnings": [str(item.message) for item in caught],
        "finite": bool(np.isfinite(java_result).all() and np.isfinite(python_result).all()),
        "component_metrics": metrics,
    }
    report = args.output_dir / "comparison.json"
    report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"report={report}")


if __name__ == "__main__":
    main()
