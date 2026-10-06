"""Generate an asymmetric 36×45×34 plus-workflow RAS integration phantom."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy.ndimage import zoom


ROOT = Path(__file__).parent
SOURCE = ROOT.parent / "plus_pre_tgdm_composition"
SHAPE = (36, 45, 34)


def main() -> None:
    source = nib.load(str(SOURCE / "multi_atlas_gm.nii.gz"))
    source_shape = np.asarray(source.shape, dtype=np.int64)
    shape = np.asarray(SHAPE, dtype=np.int64)
    factors = shape.astype(np.float64) / source_shape
    gm = zoom(
        (np.asarray(source.dataobj) > 0).astype(np.uint8),
        factors,
        order=0,
        prefilter=False,
    ).astype(bool)
    wm_source = nib.load(str(SOURCE / "multi_atlas_wm.nii.gz"))
    wm = zoom(
        (np.asarray(wm_source.dataobj) > 0).astype(np.uint8),
        factors,
        order=0,
        prefilter=False,
    ).astype(bool)
    mask_source = nib.load(str(SOURCE / "multi_atlas_brain_mask.nii.gz"))
    brain = zoom(
        (np.asarray(mask_source.dataobj) > 0).astype(np.uint8),
        factors,
        order=0,
        prefilter=False,
    ).astype(bool)
    if any(values.shape != SHAPE for values in (gm, wm, brain)):
        raise ValueError("resampling did not produce the requested synthetic grid")

    affine = source.affine.copy()
    center_world = nib.affines.apply_affine(source.affine, (source_shape - 1) / 2)
    affine[:3, :3] = source.affine[:3, :3] * (
        (source_shape - 1) / (shape - 1)
    )
    affine[:3, 3] = center_world - affine[:3, :3] @ ((shape - 1) / 2)

    x, y, _z = np.indices(SHAPE, dtype=np.float32)
    t1 = np.zeros(SHAPE, dtype=np.float32)
    t1[gm] = 90.0 + 0.25 * x[gm] + 0.1 * y[gm]
    t1[wm] = 130.0
    t1[~brain] = 0.0
    maseg = np.zeros(SHAPE, dtype=np.uint16)
    maseg[gm] = 100
    maseg[wm] = 4
    flair_affine = affine.copy()
    flair_affine[0, 3] += 0.25

    arrays = {
        "t1.nii.gz": (t1, affine),
        "flair.nii.gz": ((t1 * 0.7 + 10).astype(np.float32), flair_affine),
        "brain-mask.nii.gz": (brain.astype(np.uint8), affine),
        "maseg.nii.gz": (maseg, affine),
    }
    ROOT.mkdir(parents=True, exist_ok=True)
    records: dict[str, dict[str, object]] = {}
    for name, (values, grid) in arrays.items():
        image = nib.Nifti1Image(values, grid)
        image.set_qform(grid, code=1)
        image.set_sform(grid, code=1)
        path = ROOT / name
        nib.save(image, str(path))
        records[name] = {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "shape": list(values.shape),
            "dtype": str(values.dtype),
            "axis_codes": list(nib.aff2axcodes(grid)),
            "affine": grid.tolist(),
        }
    (ROOT / "provenance.json").write_text(
        json.dumps(
            {
                "description": "Project-generated folded-cortex plus integration fixture",
                "source_fixture": "../plus_pre_tgdm_composition",
                "generation": "nearest-neighbor label/mask scaling; deterministic synthetic modality values",
                "images": records,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
