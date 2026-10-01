"""Generate the synthetic initial-membership image for plus TGDM composition."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent


def main() -> None:
    size = 25
    coords = np.indices((size, size, size), dtype=np.float32)
    x, y, z = (coords[axis] - 12.0 for axis in range(3))
    radius = np.sqrt(x * x + y * y + z * z)
    theta = np.arctan2(y, x)
    phi = np.arctan2(z, np.sqrt(x * x + y * y))
    folds = 0.85 * np.sin(3.0 * theta + 0.4) * np.cos(2.0 * phi)
    folds += 0.35 * np.sin(5.0 * theta - 0.7) * np.cos(phi)
    membership = (radius + folds < 9.7).astype(np.float32) * 200.0
    affine = np.diag((0.8, 0.8, 1.2, 1.0))
    output = ROOT / "topology_corrected_initial_membership.nii.gz"
    nib.save(nib.Nifti1Image(membership, affine), str(output))
    provenance = {
        "kind": "project-generated folded radial membership phantom",
        "shape": [size, size, size],
        "affine_ras_mm": affine.tolist(),
        "foreground_membership": 200.0,
        "background_membership": 0.0,
        "sha256": {output.name: hashlib.sha256(output.read_bytes()).hexdigest()},
    }
    (ROOT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
