"""Generate the committed folded-cortex exact-field degree 2/3 phantom."""

from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent


def main() -> None:
    n = 21
    coordinates = np.linspace(-1.0, 1.0, n, dtype=np.float32)
    x, y, z = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij")
    radius = np.sqrt(x * x + y * y + z * z)
    azimuth = np.arctan2(y, x)
    polar = np.arctan2(np.sqrt(x * x + y * y), z)
    fold = np.float32(0.06) * np.sin(np.float32(5.0) * azimuth + np.float32(2.0) * polar)
    segmentation = ((radius >= np.float32(0.32) + fold) & (radius <= np.float32(0.76) + fold)).astype(np.uint8)
    image = np.asarray(0.38 + 0.19 * x + 0.12 * y * z + 0.07 * x * x - 0.035 * y * y * z, dtype=np.float32)
    centroids = np.linspace(0.22, 0.82, 11, dtype=np.float32)
    memberships = np.exp(-((image[..., None] - centroids) / np.float32(0.18)) ** 2).astype(np.float32)
    memberships /= memberships.sum(axis=-1, keepdims=True)
    initial_field = np.ones(image.shape, dtype=np.float32)
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    for name, data in (
        ("image", image),
        ("memberships", memberships),
        ("segmentation", segmentation),
        ("initial_field", initial_field),
    ):
        image_file = nib.Nifti1Image(data, affine)
        image_file.set_qform(affine, code=1)
        image_file.set_sform(affine, code=1)
        nib.save(image_file, str(ROOT / f"{name}.nii.gz"))
    (ROOT / "centroids.json").write_text(json.dumps([float(value) for value in centroids], indent=2) + "\n")
    (ROOT / "README.md").write_text(
        "# TOADS exact polynomial field degrees 2 and 3\n\n"
        "A deterministic 21³ synthetic RAS phantom with an angularly folded cortical shell, "
        "a smooth positive intensity field, eleven soft class-membership maps, and an initial "
        "unity correction field. All images are stored as `.nii.gz`. The fixture provides "
        "enough stride-three samples to identify the 10-term degree-2 and 20-term degree-3 "
        "Chebyshev models. Regenerate with `python generate_fixture.py`.\n"
    )


if __name__ == "__main__":
    main()
