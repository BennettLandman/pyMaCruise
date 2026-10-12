"""Create an exact 2 mm decimation of the committed 1 mm folded cortex case."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

SOURCE = Path(__file__).resolve().parents[1] / "labeled_folded_cortex"
OUTPUT = Path(__file__).resolve().parent
DISCRETE = {"ace_skeleton", "brain_mask", "initial_membership", "maseg_labels"}


def main() -> None:
    records: dict[str, object] = {}
    for source_path in sorted(SOURCE.glob("*.nii.gz")):
        image = nib.load(str(source_path))
        if image.shape != (65, 65, 65):
            raise ValueError(f"expected a 65^3 source, got {image.shape}")
        values = np.asarray(image.dataobj)[::2, ::2, ::2]
        if source_path.name.removesuffix(".nii.gz") not in DISCRETE:
            values = values.astype(np.float32)
        affine = image.affine.copy()
        affine[:3, :3] *= 2.0
        output = nib.Nifti1Image(values, affine)
        output.set_qform(affine, code=1)
        output.set_sform(affine, code=1)
        output.set_data_dtype(values.dtype)
        destination = OUTPUT / source_path.name
        nib.save(output, str(destination))
        records[destination.name] = {
            "sha256": hashlib.sha256(destination.read_bytes()).hexdigest(),
            "shape": list(values.shape),
            "dtype": str(values.dtype),
        }

    provenance = {
        "description": "2 mm exact decimation of the 1 mm labeled folded cortex TGDM case",
        "source_case": "../labeled_folded_cortex",
        "sampling": "every second voxel in each axis; same physical field of view",
        "orientation": "RAS",
        "affine_ras_mm": nib.load(str(OUTPUT / "t1.nii.gz")).affine.tolist(),
        "images": records,
    }
    (OUTPUT / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
