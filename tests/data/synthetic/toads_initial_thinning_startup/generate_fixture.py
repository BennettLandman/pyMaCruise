"""Generate the deterministic folded RAS fixture for P18cyn."""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent
SHAPE = (25, 25, 25)
CLASSES = 11
AFFINE = np.array(
    [[0.8, 0.0, 0.0, 4.0], [0.0, 1.1, 0.0, -7.0], [0.0, 0.0, 1.4, 12.0], [0, 0, 0, 1]],
    dtype=np.float64,
)


def _save(name: str, values: np.ndarray) -> None:
    image = nib.Nifti1Image(values, AFFINE)
    image.header.set_xyzt_units("mm")
    image.set_qform(AFFINE, code=1)
    image.set_sform(AFFINE, code=1)
    nib.save(image, str(ROOT / name))


def main() -> None:
    x, y, z = np.indices(SHAPE, dtype=np.float32)
    folded_coordinate = (
        x + 1.6 * np.sin((y + 0.7 * z) * np.float32(0.38))
        + 0.75 * np.sin((z - 0.35 * y) * np.float32(0.57))
    )
    folded_coordinate = np.clip(folded_coordinate, 0.0, SHAPE[0] - 1e-4)
    labels = np.floor(folded_coordinate * np.float32(CLASSES / SHAPE[0])).astype(np.uint8)

    initial = np.empty((*SHAPE, CLASSES), dtype=np.float32)
    for k in range(CLASSES):
        cyclic_distance = np.minimum((labels.astype(np.int16) - k) % CLASSES,
                                     (k - labels.astype(np.int16)) % CLASSES).astype(np.float32)
        raw = np.exp(-(cyclic_distance * cyclic_distance) / np.float32(0.5)).astype(np.float32)
        initial[..., k] = raw
    initial /= initial.sum(axis=-1, keepdims=True, dtype=np.float32)

    images = np.empty((*SHAPE, 2), dtype=np.float32)
    field = np.empty_like(images)
    bias = np.empty_like(images)
    for channel in range(2):
        base = (labels.astype(np.float32) + np.float32(1.0)) / np.float32(12.0)
        if channel == 1:
            base = np.float32(1.0) - base
        ripple = (np.sin((x * np.float32(0.31) + y * np.float32(0.17) + z * np.float32(0.23)))
                  * np.float32(0.012))
        images[..., channel] = np.asarray(base + ripple, dtype=np.float32)
        field[..., channel] = np.asarray(
            np.float32(1.0) + np.float32(0.035) * np.sin(y * np.float32(0.21) + z * np.float32(0.13) + channel),
            dtype=np.float32,
        )
        bias[..., channel] = np.asarray(
            np.float32(0.004) * np.cos(x * np.float32(0.19) - z * np.float32(0.16) + channel),
            dtype=np.float32,
        )

    shapes = np.empty((*SHAPE, CLASSES), dtype=np.float32)
    for k in range(CLASSES):
        wave = (
            np.sin((y + z * np.float32(0.8) + k * np.float32(0.3)) * np.float32(0.31))
            + np.cos((x - z * np.float32(0.4) + k) * np.float32(0.27))
        )
        shapes[..., k] = np.asarray(
            np.float32(0.52) + np.float32(0.19) * wave
            + np.float32(0.08) * np.sin((folded_coordinate + k) * np.float32(0.45)),
            dtype=np.float32,
        )
    shapes = np.clip(shapes, np.float32(0.05), np.float32(0.95)).astype(np.float32)

    _save("segmentation_ras.nii.gz", labels)
    _save("initial_memberships_ras.nii.gz", initial)
    _save("images_ras.nii.gz", images)
    _save("fields_ras.nii.gz", field)
    _save("psf_bias_ras.nii.gz", bias)
    _save("transformed_shapes_ras.nii.gz", shapes)

    spec = {
        "shape_xyz": list(SHAPE),
        "class_count": CLASSES,
        "class_order": ["MASK"] + ["OPTIMIZED_CSF"] * 2 + ["OPTIMIZED_GM"] * 5 + ["OPTIMIZED_WM"] * 3,
        "relation_limit": 2,
        "thinning_volume_ratio": 0.5,
        "thinning_limit": 3.0,
        "membership_smoothing": 0.043,
        "prior_factor": 0.64,
        "prior_coefficient": 0.16,
        "parameters": {
            "centroids": [[.28,.31,.34,.37,.40,.43,.46,.49,.52,.55,.58], [.24,.27,.30,.33,.36,.39,.42,.45,.48,.51,.54]],
            "exponent_by_channel_class": [[1,.78,.78,1.16,1.16,1.16,1.16,1.16,1.55,1.55,1.55], [1,.91,.91,1.29,1.29,1.29,1.29,1.29,1.73,1.73,1.73]],
            "modality_weights": [.73,.46],
            "class_variance_weights": [[.81,.93,1.17,.76,.88,1.02,1.14,.69,.96,1.21,.73], [1.19,.72,.87,1.08,.79,1.25,.91,1.13,.68,.98,1.16]],
            "image_variances": [.61,.34],
            "smoothing_priors": [.65,.71,.84,.91,1.02,1.13,.76,.88,1.05,1.17,.82],
            "atlas_shape_priors": [.10,.08,.09,.11,.12,.09,.10,.11,.07,.06,.07],
            "centroid_similarity": [
                [0,.11,.14,.18,.13,.16,.12,.15,.09,.10,.08], [.12,0,.13,.17,.11,.16,.14,.10,.09,.12,.08],
                [.15,.12,0,.13,.18,.12,.16,.11,.10,.07,.09], [.14,.16,.12,0,.13,.18,.11,.14,.08,.10,.12],
                [.11,.13,.17,.12,0,.11,.19,.14,.10,.09,.07], [.16,.12,.14,.17,.11,0,.12,.18,.09,.08,.10],
                [.13,.15,.11,.14,.18,.12,0,.13,.08,.09,.11], [.12,.10,.14,.16,.13,.18,.11,0,.07,.09,.08],
                [.10,.09,.12,.08,.11,.09,.10,.07,0,.13,.16], [.11,.12,.08,.10,.09,.07,.11,.08,.14,0,.12],
                [.08,.10,.09,.12,.07,.11,.10,.09,.17,.13,0],
            ],
        },
        "fixed_upstream": [
            "centroids", "class_variance_weights", "image_variances", "centroid_similarity",
            "modality_weights", "category exponents", "transformed atlas shape samples", "atlas shape priors",
        ],
        "description": "Deterministic asymmetric folded 11-class label field, smooth nonuniform initial memberships, two synthetic channels, distinct channel fields, PSF bias, and strictly positive spatially folded shape samples.",
    }
    (ROOT / "spec.json").write_text(json.dumps(spec, indent=2) + "\n")


if __name__ == "__main__":
    main()
