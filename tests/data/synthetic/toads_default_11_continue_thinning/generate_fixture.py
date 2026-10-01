"""Write deterministic true-RAS inputs for the selected 11-class chain."""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np


def main() -> None:
    root = Path(__file__).parent
    size, classes = 14, 11
    x, y, z = np.indices((size, size, size), dtype=np.int16)
    labels = ((x * 3 + y * 5 + z * 7 + (x * y) % classes) % classes).astype(np.int16)
    skeleton = np.full(labels.shape, -1, dtype=np.int16)
    for k in range(classes):
        i, j, l = 1 + k, 1 + (k * 3) % 12, 1 + (k * 5) % 12
        labels[i, j, l] = k
        skeleton[i, j, l] = k

    memberships = np.empty((*labels.shape, classes), dtype=np.float32)
    for k in range(classes):
        memberships[..., k] = np.where(labels == k, np.float32(0.91), np.float32(0.009))
    ordering = -((x * 13 + y * 7 + z * 3) % 19).astype(np.float32) / np.float32(10.0)
    shapes = np.zeros((*labels.shape, classes), dtype=np.float32)
    delta = (x * 13 + y * 7 + z * 3) % 17 - 8
    images = (labels.astype(np.float32) + np.float32(1.0)) / np.float32(12.0)
    images += delta.astype(np.float32) / np.float32(1700.0)
    field = np.float32(1.0) + ((x + 3 * y + 5 * z) % 9).astype(np.float32) / np.float32(100.0)
    second_delta = (x * 5 + y * 11 + z * 7) % 19 - 9
    second_image = (np.float32(11.0) - labels.astype(np.float32)) / np.float32(12.0)
    second_image += second_delta.astype(np.float32) / np.float32(1900.0)
    second_field = np.float32(1.0) + ((2 * x + 5 * y + 7 * z) % 11).astype(np.float32) / np.float32(100.0)
    closest = np.empty((*labels.shape, 4), dtype=np.int8)
    closest[..., 0] = (x + y + z) % classes
    closest[..., 1] = (x + 2 * y + 3 * z + 1) % classes
    closest[..., 2] = (x + 3 * y + 2 * z + 4) % classes
    closest[..., 3] = (2 * x + y + z + 7) % classes
    affine = np.diag([0.8, 1.2, 1.7, 1.0]).astype(np.float64)
    affine[:3, 3] = [-12.0, 3.5, 8.0]
    for name, values in (
        ("initial_segmentation.nii.gz", labels),
        ("initial_skeleton.nii.gz", skeleton),
        ("exact_memberships.nii.gz", memberships),
        ("growing_ordering.nii.gz", ordering),
        ("transformed_shapes.nii.gz", shapes),
        ("image_values.nii.gz", images),
        ("intensity_field.nii.gz", field),
        ("image_values_channel2.nii.gz", second_image),
        ("intensity_field_channel2.nii.gz", second_field),
        ("membership_closest_labels.nii.gz", closest),
    ):
        nib.save(nib.Nifti1Image(values, affine), str(root / name))
    stats = {
        "initial_centroids": [(k + 1) / 12.0 for k in range(classes)],
        "category_exponents": [0.8, 1.0, 1.0],
        "centroid_smoothness": 0.1,
        "intensity_field_is_fixed": True,
    }
    (root / "statistics_inputs.json").write_text(
        json.dumps(stats, indent=2) + "\n", encoding="utf-8"
    )
    membership_parameters = {
        "modality_weights": [0.73, 0.46],
        "image_variances": [0.61, 0.34],
        "class_variance_weights": [
            [0.81, 0.93, 1.17, 0.76, 0.88, 1.02, 1.14, 0.69, 0.96, 1.21, 0.73],
            [1.19, 0.72, 0.87, 1.08, 0.79, 1.25, 0.91, 1.13, 0.68, 0.98, 1.16],
        ],
        "smoothing_priors": [0.65, 0.71, 0.84, 0.91, 1.02, 1.13, 0.76, 0.88, 1.05, 1.17, 0.82],
        "atlas_shape_priors": [0.10, 0.08, 0.09, 0.11, 0.12, 0.09, 0.10, 0.11, 0.07, 0.06, 0.07],
        "smoothing": 0.043,
        "prior_factor": 0.64,
        "prior_coefficient": 0.16,
        "prior_scale": 0.05,
        "category_exponents": [0.8, 1.0, 1.0],
    }
    (root / "membership_parameters.json").write_text(
        json.dumps(membership_parameters, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
