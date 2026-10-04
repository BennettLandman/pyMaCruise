from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parent
SHAPE = (9, 9, 9)
AFFINE = np.eye(4, dtype=np.float64)


def main() -> None:
    x, y, z = np.indices(SHAPE, dtype=np.int32)
    images = np.empty((*SHAPE, 2), dtype=np.float32)
    images[..., 0] = np.float32(0.25) + np.float32(0.011) * ((3 * x + 5 * y + 7 * z) % 31)
    images[..., 1] = np.float32(0.35) + np.float32(0.013) * ((7 * x + 2 * y + 11 * z) % 29)
    segmentation = ((x + 2 * y + 3 * z) % 13).astype(np.int16)
    shapes = np.empty((*SHAPE, 13), dtype=np.float32)
    for class_index in range(13):
        offset = (x + 2 * y + 3 * z + class_index) % 25
        shapes[..., class_index] = np.float32(0.25) + np.float32(0.01) * offset.astype(np.float32)

    images_path = ROOT / "images_ras.nii.gz"
    segmentation_path = ROOT / "constructor_labels_ras.nii.gz"
    shapes_path = ROOT / "transformed_shapes_ras.nii.gz"
    nib.save(nib.Nifti1Image(images, AFFINE), images_path)
    nib.save(nib.Nifti1Image(segmentation, AFFINE), segmentation_path)
    nib.save(nib.Nifti1Image(shapes, AFFINE), shapes_path)
    spec = {
        "synthetic": True,
        "dimensions_xyz": list(SHAPE),
        "channels_last": ["T1_MPRAGE", "FLAIR"],
        "class_count": 13,
        "affine_ras_mm": AFFINE.tolist(),
        "image_sha256_float32le": hashlib.sha256(images.astype("<f4", copy=False).tobytes()).hexdigest(),
        "labels_sha256_int16le": hashlib.sha256(segmentation.astype("<i2", copy=False).tobytes()).hexdigest(),
        "shapes_sha256_float32le": hashlib.sha256(shapes.astype("<f4", copy=False).tobytes()).hexdigest(),
        "provenance": "Deterministic synthetic voxel ramps, interleaved 13-class constructor labels, and positive transformed shape probabilities.",
        "java_boundary": "Tests validate each selected Java component and exact composition; this generated volume is not a composite Java segmentation reference.",
    }
    (ROOT / "spec.json").write_text(json.dumps(spec, indent=2) + "\n")


if __name__ == "__main__":
    main()
