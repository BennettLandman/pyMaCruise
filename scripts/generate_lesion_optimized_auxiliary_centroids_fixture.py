"""Add optimized-factor cutoff and out-of-range cases to the special-class phantom."""

from pathlib import Path

import nibabel as nib
import numpy as np


BASE = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic"
SOURCE = BASE / "lesion_special_centroids" / "special_centroid_inputs.nii.gz"
ROOT = BASE / "lesion_optimized_auxiliary_centroids"


def main() -> None:
    source = nib.load(str(SOURCE))
    frames = np.asarray(source.dataobj, dtype=np.float32).copy()
    # Tiny optimized class winner: dist = membership² × intensity <= ZERO.
    frames[2, 0, 0, 3] = np.float32(1e-6)  # hard class 2 membership
    frames[2, 0, 0, 10:12] = np.float32(0.0)
    # Large optimized LESION-class winner drives the Java table out of range.
    frames[3, 1, 0, 2] = np.float32(1.0)  # hard class 1 membership
    frames[3, 1, 0, 10:12] = np.float32(0.0)
    frames[3, 1, 0, 6] = np.float32(5.0)
    ROOT.mkdir(parents=True, exist_ok=True)
    result = nib.Nifti1Image(frames, source.affine, source.header.copy())
    result.set_qform(source.get_qform(), code=int(source.header["qform_code"]))
    result.set_sform(source.get_sform(), code=int(source.header["sform_code"]))
    nib.save(result, str(ROOT / "optimized_auxiliary_centroid_inputs.nii.gz"))


if __name__ == "__main__":
    main()
