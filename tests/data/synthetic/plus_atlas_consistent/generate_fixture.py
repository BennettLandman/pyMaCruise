"""Generate a deterministic LesionTOADS-atlas-aligned synthetic PLUS case.

The case uses the bundled topology labels and descriptor intensity means. It is
an integration fixture, not independent ground truth for atlas methods.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_atlas_priors import load_lesion_atlas_priors, load_lesion_atlas_topology

ROOT = Path(__file__).parent
REPO = Path(__file__).resolve().parents[4]
ATLAS = REPO / "data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"
TOPOLOGY_ID_TO_MASEG = {
    1: 0,
    5: 5,
    6: 5,
    14: 100,
    15: 100,
    16: 100,
    17: 100,
    18: 100,
    22: 100,
    24: 44,
    25: 44,
    26: 44,
    27: 44,
}


def main() -> None:
    topology = load_lesion_atlas_topology(ATLAS)
    priors = load_lesion_atlas_priors(ATLAS, ("T1_MPRAGE", "FLAIR"))
    topology_image = topology.image
    class_ids = np.asarray(topology_image.dataobj, dtype=np.uint8)
    observed = set(int(v) for v in np.unique(class_ids))
    if observed != set(TOPOLOGY_ID_TO_MASEG):
        raise ValueError(f"unexpected topology IDs: {sorted(observed)}")

    class_to_index = {int(item.label): index for index, item in enumerate(priors.structures)}
    t1_means = priors.intensity_priors[0] * np.float32(255.0)
    flair_means = priors.intensity_priors[1] * np.float32(255.0)
    t1 = np.zeros(class_ids.shape, dtype=np.float32)
    flair = np.full(class_ids.shape, 80.0, dtype=np.float32)
    maseg = np.zeros(class_ids.shape, dtype=np.uint16)
    for topology_id, maseg_id in TOPOLOGY_ID_TO_MASEG.items():
        mask = class_ids == topology_id
        class_index = class_to_index[topology_id]
        t1[mask] = t1_means[class_index]
        flair[mask] = flair_means[class_index]
        maseg[mask] = maseg_id
    # A positive exterior baseline keeps cubic ANTs interpolation nonnegative
    # at the brain boundary, as in the existing folded PLUS fixture.
    flair[class_ids == 1] = 80.0

    # Add smooth, deterministic within-class intensity variation. This avoids a
    # perfectly flat modality while keeping the labels and atlas anatomy fixed.
    i, j, k = np.ogrid[:class_ids.shape[0], :class_ids.shape[1], :class_ids.shape[2]]
    texture = (
        0.012 * np.sin(i * 0.071 + j * 0.023)
        + 0.009 * np.cos(j * 0.053 - k * 0.031)
        + 0.006 * np.sin(k * 0.047 + i * 0.019)
    ).astype(np.float32)
    brain = class_ids != 1
    t1[brain] *= 1.0 + texture[brain]
    flair[brain] *= 1.0 + texture[brain]
    mask = brain.astype(np.uint8)

    arrays = {
        "t1.nii.gz": (t1, np.float32),
        "flair.nii.gz": (flair, np.float32),
        "brain-mask.nii.gz": (mask, np.uint8),
        "maseg.nii.gz": (maseg, np.uint16),
    }
    records = {}
    for name, (values, dtype) in arrays.items():
        image = nib.Nifti1Image(values.astype(dtype, copy=False), topology_image.affine, topology_image.header.copy())
        image.set_qform(topology_image.affine, code=1)
        image.set_sform(topology_image.affine, code=1)
        image.set_data_dtype(dtype)
        path = ROOT / name
        nib.save(image, str(path))
        records[name] = {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "shape": list(values.shape),
            "dtype": str(np.dtype(dtype)),
            "axis_codes": list(nib.aff2axcodes(image.affine)),
            "affine_ras_mm": image.affine.tolist(),
        }
    (ROOT / "provenance.json").write_text(
        json.dumps({
            "description": "Atlas-aligned synthetic MaCRUISEplus integration fixture",
            "evidence_category": "synthetic integration; not independent atlas or source-method ground truth",
            "source_atlas_descriptor": str(ATLAS.relative_to(REPO)),
            "topology_id_to_maseg_label": {str(k): v for k, v in TOPOLOGY_ID_TO_MASEG.items()},
            "intensity_model": "selected T1_MPRAGE and FLAIR descriptor means scaled by 255 with deterministic smooth within-class modulation; FLAIR exterior baseline 80 keeps cubic interpolation nonnegative",
            "spatial_assumptions": "RAS NIfTIs on bundled topology-atlas grid; atlas classes supply anatomy and synthetic MaSeg tissue grouping",
            "images": records,
        }, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
