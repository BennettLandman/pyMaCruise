"""Generate a deterministic asymmetric 64×80×64 folded plus phantom."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy.ndimage import zoom


ROOT = Path(__file__).parent
SOURCE = ROOT.parent / "plus_pre_tgdm_composition"
SHAPE = (64, 80, 64)


def _nearest_mask(path: Path, shape: tuple[int, int, int]) -> np.ndarray:
    source = nib.load(str(path))
    factors = np.asarray(shape, dtype=np.float64) / np.asarray(source.shape, dtype=np.float64)
    result = zoom(
        (np.asarray(source.dataobj) > 0).astype(np.uint8),
        factors,
        order=0,
        prefilter=False,
    ).astype(bool)
    if result.shape != shape:
        raise ValueError(f"resampling {path.name} produced unexpected shape {result.shape}")
    return result


def main() -> None:
    shape = np.asarray(SHAPE, dtype=np.int64)
    source_t1 = nib.load(str(SOURCE / "multi_atlas_gm.nii.gz"))
    source_shape = np.asarray(source_t1.shape, dtype=np.int64)
    gm = _nearest_mask(SOURCE / "multi_atlas_gm.nii.gz", SHAPE)
    wm = _nearest_mask(SOURCE / "multi_atlas_wm.nii.gz", SHAPE)
    brain = _nearest_mask(SOURCE / "multi_atlas_brain_mask.nii.gz", SHAPE)

    affine = source_t1.affine.copy()
    center_world = nib.affines.apply_affine(source_t1.affine, (source_shape - 1) / 2)
    affine[:3, :3] = source_t1.affine[:3, :3] * ((source_shape - 1) / (shape - 1))
    affine[:3, 3] = center_world - affine[:3, :3] @ ((shape - 1) / 2)

    x = np.arange(SHAPE[0], dtype=np.float32)[:, None, None]
    y = np.arange(SHAPE[1], dtype=np.float32)[None, :, None]
    t1 = np.zeros(SHAPE, dtype=np.float32)
    t1[gm] = 90.0 + 0.25 * np.broadcast_to(x, SHAPE)[gm] + 0.1 * np.broadcast_to(y, SHAPE)[gm]
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
    records: dict[str, dict[str, object]] = {}
    ROOT.mkdir(parents=True, exist_ok=True)
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
                "description": "Synthetic folded-cortex plus workflow resource fixture",
                "source_fixture": "../plus_pre_tgdm_composition",
                "generation": "deterministic nearest-neighbor support scaling with synthetic T1/FLAIR values",
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
