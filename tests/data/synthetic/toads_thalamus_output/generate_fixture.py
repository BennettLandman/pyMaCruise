"""Regenerate the synthetic NIfTI inputs and pinned excerpt output fixture."""

from __future__ import annotations

import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).resolve().parent
HARNESS = Path(__file__).resolve().parents[3] / "java_reference" / "ToadsThalamusOutputReferenceHarness.java"
SOURCE_SHA256 = "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
WRAPPER_SHA256 = "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
AFFINE_RAS = np.eye(4, dtype=np.float64)


def _membership_values() -> np.ndarray:
    nx, ny, nz, classes = 4, 3, 5, 6
    values = np.empty((nx, ny, nz, classes), dtype=np.float32)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                for k in range(classes):
                    if (x + 2 * y + z + k) % 11 == 0:
                        values[x, y, z, k] = np.float32(0.0)
                    else:
                        numerator = np.float32(x * 71 + y * 19 + z * 7 + k * 11 + 1)
                        values[x, y, z, k] = np.float32(numerator / np.float32(100.0))
    return values


def _best_values() -> np.ndarray:
    patterns = np.asarray(
        [[0, 4, 1, 4], [4, 4, 3, 2], [1, 2, 3, 1], [0, 2, 3, 1]],
        dtype=np.int16,
    )
    result = np.empty((4, 3, 5, 4), dtype=np.int16)
    for x in range(4):
        for y in range(3):
            for z in range(5):
                voxel = (x * 3 + y) * 5 + z
                result[x, y, z] = patterns[voxel % len(patterns)]
    return result


def main() -> None:
    memberships = _membership_values()
    best = _best_values()
    for name, array in (("memberships", memberships), ("best_classes", best)):
        nib.save(nib.Nifti1Image(array, AFFINE_RAS), HERE / f"{name}.nii.gz")

    with tempfile.TemporaryDirectory(prefix="java-build-", dir=HERE) as build_dir:
        subprocess.run(["javac", "-d", build_dir, str(HARNESS)], check=True)
        completed = subprocess.run(
            ["java", "-cp", build_dir, "ToadsThalamusOutputReferenceHarness"],
            check=True,
            capture_output=True,
            text=True,
        )
    reference = json.loads(completed.stdout)
    reference.update(
        source_method="MultipleToadSegmentation.generateThalamusOutput",
        source_sha256=SOURCE_SHA256,
        wrapper_sha256=WRAPPER_SHA256,
        harness_sha256=hashlib.sha256(HARNESS.read_bytes()).hexdigest(),
        fixture_sha256={
            filename: hashlib.sha256((HERE / filename).read_bytes()).hexdigest()
            for filename in ("memberships.nii.gz", "best_classes.nii.gz")
        },
        affine_ras=AFFINE_RAS.tolist(),
    )
    for mode in ("exact", "approx"):
        expected = np.asarray(reference["java_reference"][mode], dtype=np.float32).reshape(
            tuple(reference["shape_xyz"])
        )
        nib.save(nib.Nifti1Image(expected, AFFINE_RAS), HERE / f"expected_{mode}.nii.gz")
        reference["fixture_sha256"][f"expected_{mode}.nii.gz"] = hashlib.sha256(
            (HERE / f"expected_{mode}.nii.gz").read_bytes()
        ).hexdigest()
    (HERE / "java_reference.json").write_text(json.dumps(reference, indent=2) + "\n")


if __name__ == "__main__":
    main()
