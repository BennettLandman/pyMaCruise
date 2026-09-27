"""Create committed inputs for the source-mapped segmentation refinement."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np


OUTPUT = Path(__file__).parent
SHAPE = (11, 11, 11)
AFFINE = np.diag([0.5, 0.75, 1.0, 1.0])


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    labels = np.zeros(SHAPE, dtype=np.uint16)
    inner = np.full(SHAPE, 1.0, dtype=np.float32)
    outer = np.full(SHAPE, 1.0, dtype=np.float32)

    # A clean cortical region for central-mask voting.
    labels[4:9, 4:9, 4:9] = 120
    outer[4:9, 4:9, 4:9] = -1.0
    # Inner-side erroneous cortical label, with local WM mode 44.
    labels[1:6, 1:6, 1:6] = 44
    inner[1:6, 1:6, 1:6] = -1.0
    labels[3, 3, 3] = 101
    inner[3, 3, 3] = -0.5
    # A tie between 44 and 45 exercises the lowest-label tie rule.
    labels[2, 8, 2] = 109
    labels[1, 8, 1] = 44
    labels[1, 8, 2] = 44
    labels[3, 8, 1] = 45
    labels[3, 8, 2] = 45
    inner[1, 8, 1] = -1.0
    inner[1, 8, 2] = -1.0
    inner[3, 8, 1] = -1.0
    inner[3, 8, 2] = -1.0
    inner[2, 8, 2] = -0.5
    # Beyond-outer-surface cortical label is removed.
    labels[8, 3, 3] = 105
    inner[8, 3, 3] = 1.0
    outer[8, 3, 3] = 0.5
    # Central WM and missing central-GM voxels receive a local cortical mode.
    labels[5, 5, 5] = 45
    inner[5, 5, 5] = 0.5
    outer[5, 5, 5] = -0.5
    labels[7, 7, 7] = 0
    inner[7, 7, 7] = 0.5
    outer[7, 7, 7] = -0.5
    # Border voxel has no valid neighbors of the required tissue class.
    labels[0, 0, 10] = 102
    inner[0, 0, 10] = -0.5
    outer[0, 0, 10] = -1.0
    labels[10, 10, 10] = 7

    arrays = {"segmentation": labels, "inner_levelset": inner, "outer_levelset": outer}
    for name, values in arrays.items():
        nib.save(nib.Nifti1Image(values, AFFINE), str(OUTPUT / f"{name}.nii.gz"))
    spec = {
        "purpose": "Synthetic P32 MaCRUISE surface-guided segmentation refinement",
        "shape_voxels": list(SHAPE),
        "spacing_mm": [0.5, 0.75, 1.0],
        "orientation": "RAS",
        "source_coefficients_mm": 0.2,
        "source_max_spacing_threshold_voxels": 0.2,
        "target_voxels_ijk": {
            "inner_fill_44": [3, 3, 3],
            "inner_tie_low_label_44": [2, 8, 2],
            "outer_clear_zero": [8, 3, 3],
            "central_wm_fill_120": [5, 5, 5],
            "central_gm_fill_120": [7, 7, 7],
            "unresolved_border_label_102": [0, 0, 10],
            "unrelated_label_7": [10, 10, 10],
        },
        "seed": None,
    }
    (OUTPUT / "spec.json").write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")
    names = [*arrays, "spec"]
    entries = []
    for name in names:
        path = OUTPUT / f"{name}.nii.gz" if name != "spec" else OUTPUT / "spec.json"
        entries.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}")
    (OUTPUT / "SHA256SUMS").write_text("\n".join(entries) + "\n", encoding="ascii")


if __name__ == "__main__":
    main()
