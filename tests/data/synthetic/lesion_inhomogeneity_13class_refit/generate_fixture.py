from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / "toads_exact_polynomial_field_degrees23"


def main() -> None:
    source_images = {
        name: nib.load(str(SOURCE / (name + ".nii.gz")))
        for name in ("image", "memberships", "segmentation", "initial_field")
    }
    image = np.asarray(source_images["image"].dataobj, dtype=np.float32)
    memberships = np.asarray(source_images["memberships"].dataobj, dtype=np.float32)
    segmentation = np.asarray(source_images["segmentation"].dataobj, dtype=np.uint8)

    expanded_memberships = np.zeros(image.shape + (13,), dtype=np.float32)
    expanded_memberships[..., :11] = memberships
    expanded_segmentation = segmentation.copy()
    candidates = np.flatnonzero(expanded_segmentation.ravel() != 0)
    lesion_voxels = candidates[::19]
    added_object_voxels = candidates[7::23]
    expanded_segmentation.ravel()[lesion_voxels] = 11
    expanded_segmentation.ravel()[added_object_voxels] = 12
    expanded_memberships[..., 11].ravel()[lesion_voxels] = np.float32(0.35)
    expanded_memberships[..., 12].ravel()[added_object_voxels] = np.float32(0.65)

    outputs = {
        "image": image,
        "memberships": expanded_memberships,
        "segmentation": expanded_segmentation,
        "initial_field": np.asarray(source_images["initial_field"].dataobj, dtype=np.float32),
    }
    for name, data in outputs.items():
        reference = source_images[name]
        nib.save(nib.Nifti1Image(data, reference.affine, reference.header), str(ROOT / (name + ".nii.gz")))

    centers = np.asarray(json.loads((SOURCE / "centroids.json").read_text()), dtype=np.float32)
    centers = np.concatenate((centers, np.asarray([91.0, 117.0], dtype=np.float32)))
    (ROOT / "centroids.json").write_text(json.dumps(centers.tolist(), indent=2) + "\n")
    class_types = ("mask", "optimized_csf", "optimized_csf", *("obj",) * 7, "lesion", "obj", "obj")
    (ROOT / "class_types.json").write_text(json.dumps(class_types, indent=2) + "\n")


if __name__ == "__main__":
    main()
