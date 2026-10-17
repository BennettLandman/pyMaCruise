"""Generate a compact two-class AtlasEM pyramid atlas for PLUS integration tests."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parent
SHAPE = (24, 24, 24)
AFFINE = np.diag([1.0, 1.0, 1.0, 1.0])


def main() -> None:
    x, y, z = np.indices(SHAPE)
    center = (np.asarray(SHAPE, dtype=np.float32) - np.float32(1.0)) / np.float32(2.0)
    radius = np.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2 + (z - center[2]) ** 2)
    cortex = radius < np.float32(8.0)
    topology = np.where(cortex, 2, 1).astype(np.uint8)
    background = np.where(cortex, 0.05, 0.95).astype(np.float32)
    cortical_prior = np.where(cortex, 0.95, 0.05).astype(np.float32)
    for filename, values in (
        ("topology.nii.gz", topology),
        ("background.nii.gz", background),
        ("cortex.nii.gz", cortical_prior),
    ):
        nib.save(nib.Nifti1Image(values, AFFINE), ROOT / filename)

    descriptor = "\n".join((
        "Structure Atlas File (edit at your own risks)", "", "Structures:\t2",
        "Background\t1\tmask", "Cortex\t2\tclosed", "", "Topology Atlas:",
        "File:\ttopology.nii.gz", "Dimensions:\t24\t24\t24", "Resolutions:\t1\t1\t1",
        "", "Shape Atlas:\t4", "Dimensions:\t24\t24\t24", "Resolutions:\t1\t1\t1",
        "Orientations:\t1\t1\t6\t3", "Structure:\tBackground", "File:\tbackground.nii.gz",
        "Dimensions:\t24\t24\t24", "Resolutions:\t1\t1\t1", "Center:\t0\t0\t0",
        "Structure:\tCortex", "File:\tcortex.nii.gz", "Dimensions:\t24\t24\t24",
        "Resolutions:\t1\t1\t1", "Center:\t0\t0\t0", "",
        "Registered Shapes (0/1: true/false)", "0\t1", "0\t1", "Intensity Atlas:\t12",
        "T1_MPRAGE\t0.2\t0.8", "Intensity Variance:", "T1_MPRAGE\t1\t1",
        "Lesions:", "T1_MPRAGE\t0.2\tFLAIR\t0.2", "Black Holes:", "T1_MPRAGE\t0.05",
        "Slope:\t3.5", "Modality weights:", "T1_MPRAGE\t0.4\t0.6",
        "Optimized factors:", "T1_MPRAGE\t1\t1\t1\t1", "",
    ))
    (ROOT / "atlasem.txt").write_text(descriptor, encoding="utf-8")
    hashes = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(ROOT.iterdir())
        if path.name != "SHA256SUMS" and path.is_file()
    }
    (ROOT / "SHA256SUMS").write_text(
        "".join(f"{digest}  {name}\n" for name, digest in hashes.items()),
        encoding="utf-8",
    )
    (ROOT / "provenance.json").write_text(json.dumps({
        "synthetic_only": True,
        "purpose": "end-to-end AtlasEM/Dura-enabled PLUS integration and restart test",
        "spatial_convention": "RAS",
        "shape_xyz": list(SHAPE),
        "spacing_mm": [1.0, 1.0, 1.0],
        "classes": ["Background", "Cortex"],
        "registration_levels_supported": 4,
        "selected_atlas_equivalence": False,
    }, indent=2) + "\n", encoding="utf-8")
    hashes = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(ROOT.iterdir())
        if path.name not in {"SHA256SUMS", "provenance.json"} and path.is_file()
    }
    (ROOT / "SHA256SUMS").write_text(
        "".join(f"{digest}  {name}\n" for name, digest in hashes.items()),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
