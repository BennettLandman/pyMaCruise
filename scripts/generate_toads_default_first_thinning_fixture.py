#!/usr/bin/env python3
"""Generate the analytic folded 11-class first-thinning RAS fixture."""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "toads_default_first_thinning"
N = 36
CLASSES = 11


def _label(x: int, y: int, z: int) -> int:
    return (x + 2 * y + 3 * z + (y * z) % 3 + (x * y) % 2) % CLASSES


def _intensity(channel: int, x: int, y: int, z: int, label: int) -> np.float32:
    if channel == 0:
        delta = (x * 13 + y * 7 + z * 3) % 17 - 8
        return np.float32(np.float32(label + 1) / np.float32(12.0) + np.float32(delta) / np.float32(1700.0))
    delta = (x * 5 + y * 11 + z * 7) % 19 - 9
    return np.float32(np.float32(11 - label) / np.float32(12.0) + np.float32(delta) / np.float32(1900.0))


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    labels = np.empty((N, N, N), dtype=np.int16)
    images = np.empty((N, N, N, 2), dtype=np.float32)
    field = np.empty((N, N, N), dtype=np.float32)
    transformed_shapes = np.empty((N, N, N, CLASSES), dtype=np.float32)
    closest_labels = np.empty((N, N, N, 4), dtype=np.int16)
    for x in range(N):
        for y in range(N):
            for z in range(N):
                label = _label(x, y, z)
                labels[x, y, z] = label
                images[x, y, z, 0] = _intensity(0, x, y, z, label)
                images[x, y, z, 1] = _intensity(1, x, y, z, label)
                field[x, y, z] = np.float32(1.0) + np.float32((x + 3 * y + 5 * z) % 9) / np.float32(100.0)
                for k in range(CLASSES):
                    transformed_shapes[x, y, z, k] = (
                        np.float32(0.20)
                        + np.float32((7 * x + 5 * y + 3 * z + 2 * k) % 70)
                        / np.float32(100.0)
                    )
                closest_labels[x, y, z] = (
                    (x + y + z) % CLASSES,
                    (x + 2 * y + 3 * z + 1) % CLASSES,
                    (x + 3 * y + 2 * z + 4) % CLASSES,
                    (2 * x + y + z + 7) % CLASSES,
                )
    affine = np.asarray(
        [[1.1, 0.0, 0.0, -18.0], [0.0, 1.3, 0.0, 4.0], [0.0, 0.0, 1.7, 22.0], [0.0, 0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    image = nib.Nifti1Image(images, affine)
    image.header.set_xyzt_units("mm", "sec")
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    nib.save(image, ROOT / "subject_images.nii.gz")
    segmentation = nib.Nifti1Image(labels, affine)
    segmentation.header.set_xyzt_units("mm")
    segmentation.set_qform(affine, code=1)
    segmentation.set_sform(affine, code=1)
    nib.save(segmentation, ROOT / "initial_segmentation.nii.gz")
    correction_field = nib.Nifti1Image(field, affine)
    correction_field.header.set_xyzt_units("mm")
    correction_field.set_qform(affine, code=1)
    correction_field.set_sform(affine, code=1)
    nib.save(correction_field, ROOT / "fixed_inhomogeneity_field.nii.gz")
    shapes = nib.Nifti1Image(transformed_shapes, affine)
    shapes.header.set_xyzt_units("mm")
    shapes.set_qform(affine, code=1)
    shapes.set_sform(affine, code=1)
    nib.save(shapes, ROOT / "membership_transformed_shapes.nii.gz")
    relations = nib.Nifti1Image(closest_labels, affine)
    relations.header.set_xyzt_units("mm")
    relations.set_qform(affine, code=1)
    relations.set_sform(affine, code=1)
    nib.save(relations, ROOT / "membership_closest_labels.nii.gz")
    centroids = np.empty((2, CLASSES), dtype=np.float32)
    for k in range(CLASSES):
        centroids[0, k] = np.float32(k + 1) / np.float32(12.0)
        centroids[1, k] = np.float32(11 - k) / np.float32(12.0)
    (ROOT / "spec.json").write_text(
        json.dumps(
            {
                "shape_xyz_channels": [N, N, N, 2],
                "classes": CLASSES,
                "crop_border": 10,
                "volume_ratio": 0.5,
                "thinning_limit": 3.0,
                "centroids": centroids.tolist(),
                "affine_ras_mm": affine.tolist(),
                "label_formula": "(x + 2*y + 3*z + (y*z)%3 + (x*y)%2) % 11",
            },
            indent=2,
        )
        + "\n"
    )
    (ROOT / "membership_stage_spec.json").write_text(
        json.dumps(
            {
                "image_variance_weights": [0.61, 0.34],
                "modality_weights": [0.73, 0.46],
                "smoothing": 0.043,
                "prior_factor": 0.64,
                "prior_coefficient": 0.16,
                "prior_scale": 0.05,
                "centroid_smoothness": 0.1,
                "category_exponents": [[0.8, 1.0, 1.0], [0.8, 1.0, 1.0]],
                "class_variance_weights": [
                    [0.81, 0.93, 1.17, 0.76, 0.88, 1.02, 1.14, 0.69, 0.96, 1.21, 0.73],
                    [1.19, 0.72, 0.87, 1.08, 0.79, 1.25, 0.91, 1.13, 0.68, 0.98, 1.16],
                ],
                "smoothing_priors": [0.65, 0.71, 0.84, 0.91, 1.02, 1.13, 0.76, 0.88, 1.05, 1.17, 0.82],
                "atlas_shape_priors": [0.10, 0.08, 0.09, 0.11, 0.12, 0.09, 0.10, 0.11, 0.07, 0.06, 0.07],
                "transformed_shape_formula": "0.20 + ((7*x + 5*y + 3*z + 2*k) % 70) / 100",
                "closest_label_formula": [
                    "(x+y+z)%11", "(x+2*y+3*z+1)%11",
                    "(x+3*y+2*z+4)%11", "(2*x+y+z+7)%11"
                ],
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
