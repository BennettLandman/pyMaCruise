"""Generate a folded-cortex case for the selected TOADS rigid pyramid."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parent


def _sigmoid(values: np.ndarray) -> np.ndarray:
    return (1.0 / (1.0 + np.exp(-np.clip(values, -40.0, 40.0)))).astype(
        np.float32
    )


def _folded_tissue_fields(center_x: float) -> tuple[np.ndarray, ...]:
    shape = (48, 48, 48)
    x, y, z = np.indices(shape, dtype=np.float32)
    x = x - np.float32(center_x)
    y = y - np.float32(23.5)
    z = z - np.float32(23.5)
    radius = np.sqrt(x * x + y * y + z * z)
    phi = np.arctan2(y, x)
    theta = np.arctan2(z, np.sqrt(x * x + y * y) + np.float32(1e-6))
    fold = np.float32(1.25) * np.sin(
        np.float32(5.0) * phi
        + np.float32(0.7) * np.sin(np.float32(3.0) * theta)
    ) * np.sin(
        np.float32(4.0) * theta
        + np.float32(0.2) * np.cos(np.float32(2.0) * phi)
    )
    inner = np.float32(10.5) + np.float32(0.35) * fold
    outer = np.float32(17.0) + fold
    white_matter = _sigmoid((inner - radius) / np.float32(1.0))
    gray_matter = _sigmoid((radius - inner) / np.float32(0.8)) * _sigmoid(
        (outer - radius) / np.float32(0.8)
    )
    csf = _sigmoid((radius - outer) / np.float32(1.0))
    return white_matter, gray_matter, csf


def main() -> None:
    # The input cortex is offset 1.5 voxels from the atlas to force nontrivial
    # registration updates at all three source pyramid scales.
    wm, gm, csf = _folded_tissue_fields(23.5)
    memberships = np.stack((np.ones_like(wm), wm, gm, csf), axis=-1)
    tissue = memberships[..., 1:]
    tissue /= np.maximum(tissue.sum(axis=-1, keepdims=True), np.float32(1e-6))

    wm, gm, csf = _folded_tissue_fields(25.0)
    shapes = np.stack((np.ones_like(wm), wm, gm, csf), axis=-1)
    shapes[..., 1:] /= np.maximum(
        shapes[..., 1:].sum(axis=-1, keepdims=True), np.float32(1e-6)
    )

    affine = np.eye(4, dtype=np.float64)
    ROOT.mkdir(parents=True, exist_ok=True)
    nib.save(
        nib.Nifti1Image(memberships.astype(np.float32), affine),
        ROOT / "memberships.nii.gz",
    )
    nib.save(
        nib.Nifti1Image(shapes.astype(np.float32), affine),
        ROOT / "atlas-shape-priors.nii.gz",
    )


if __name__ == "__main__":
    main()
