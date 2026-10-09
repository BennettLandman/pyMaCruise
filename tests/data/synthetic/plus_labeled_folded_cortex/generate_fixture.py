"""Create same-grid folded WM/GM/CSF inputs for complete plus validation."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).parent
BASE = ROOT.parent / "labeled_folded_cortex"


def main() -> None:
    reference = nib.load(str(BASE / "t1.nii.gz"))
    labels = np.asarray(nib.load(str(BASE / "maseg_labels.nii.gz")).dataobj, dtype=np.uint16)
    if not {0, 5, 44, 100}.issubset(set(int(v) for v in np.unique(labels))):
        raise ValueError("base fixture must provide labels 0/5/44/100")
    # Tissue intensities are deterministic and modality ordered. They are
    # synthetic contrasts, not a learned or clinical MRI intensity model.
    t1 = np.zeros(labels.shape, dtype=np.float32)
    flair = np.full(labels.shape, 80.0, dtype=np.float32)
    t1[labels == 5], flair[labels == 5] = 20.0, 90.0
    t1[labels == 44], flair[labels == 44] = 130.0, 180.0
    t1[labels == 100], flair[labels == 100] = 90.0, 150.0
    mask = np.asarray(nib.load(str(BASE / "brain_mask.nii.gz")).dataobj, dtype=np.uint8)
    arrays = {
        "t1.nii.gz": t1,
        "flair.nii.gz": flair,
        "brain-mask.nii.gz": mask,
        "maseg.nii.gz": labels,
    }
    records = {}
    for name, values in arrays.items():
        image = nib.Nifti1Image(values, reference.affine, reference.header.copy())
        image.set_qform(reference.affine, code=1)
        image.set_sform(reference.affine, code=1)
        image.set_data_dtype(values.dtype)
        path = ROOT / name
        nib.save(image, str(path))
        records[name] = {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "shape": list(values.shape),
            "dtype": str(values.dtype),
            "axis_codes": list(nib.aff2axcodes(image.affine)),
            "affine_ras_mm": image.affine.tolist(),
        }
    (ROOT / "provenance.json").write_text(
        json.dumps({
            "description": "65-cubed labeled folded cortex for restartable synthetic MaCRUISEplus validation",
            "source_fixture": "../labeled_folded_cortex",
            "labels": {"0": "background", "5": "sulcal CSF", "44": "white matter", "100": "cortical gray matter"},
            "intensity_assumption": "deterministic tissue-specific T1 and FLAIR values; not clinical MRI",
            "images": records,
        }, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
