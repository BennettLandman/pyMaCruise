"""Regenerate RAS input NIfTIs and Java harness references for this fixture."""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[4]
HARNESS = ROOT / "tests/java_reference/ToadsCruiseOutputsReferenceHarness.java"
HERE = Path(__file__).resolve().parent
NAMES = (
    "SulcalCSF", "Sulcal-CSF", "CorticalGM", "Cerebrum-GM", "CerebralGM",
    "WhiteMatter", "CerebellarGM", "CerebellumXYZ", "Brainstem", "Background",
    "Sulcalcsf", "Cerebrum-gm", "DeepWM", "Cerebellum",
)
SHAPE = (3, 4, 5)
APPROX = 3


def main() -> None:
    x, y, z = np.indices(SHAPE)
    memberships = np.empty((*SHAPE, len(NAMES)), dtype=np.float32)
    for k in range(len(NAMES)):
        memberships[..., k] = (x * 101 + y * 23 + z * 7 + k * 13 + 1) / 1000.0
    segmentation = ((x + 2 * y + 2 * z) % len(NAMES)).astype(np.int16)
    first = (x + 2 * y + z) % len(NAMES)
    best = np.stack([first, (first + 5) % len(NAMES), (first + 9) % len(NAMES)], axis=-1).astype(np.int16)
    affine = np.asarray(
        [[1.3, 0, 0, 21.0], [0, 1.7, 0, -9.0], [0, 0, 2.2, 13.0], [0, 0, 0, 1]],
        dtype=np.float64,
    )
    nib.save(nib.Nifti1Image(memberships, affine), HERE / "memberships.nii.gz")
    nib.save(nib.Nifti1Image(segmentation, affine), HERE / "segmentation.nii.gz")
    nib.save(nib.Nifti1Image(best, affine), HERE / "best_classes.nii.gz")

    with tempfile.TemporaryDirectory(prefix="toads-cruise-reference-") as tmp:
        subprocess.run(["javac", "-d", tmp, str(HARNESS)], check=True)
        lines = subprocess.check_output(
            ["java", "-cp", tmp, "ToadsCruiseOutputsReferenceHarness"], text=True
        ).splitlines()
    result: dict[str, dict[str, list[float]]] = {"exact": {}, "approx": {}}
    for line in lines:
        fields = line.split("\t")
        result[fields[0]][fields[1]] = [float(v) for v in fields[2:]]
    reference = {
        "source": "VandyCRUISE MultipleToadSegmentation.generateCruiseOutputs excerpt",
        "source_sha256": "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884",
        "wrapper_sha256": "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d",
        "shape_xyz": list(SHAPE),
        "atlas_names": list(NAMES),
        "channels": ["csf", "gray_matter", "white_matter_fill", "white_matter_mask"],
        "java_reference": result,
    }
    (HERE / "java_reference.json").write_text(json.dumps(reference, indent=2) + "\n")


if __name__ == "__main__":
    main()
