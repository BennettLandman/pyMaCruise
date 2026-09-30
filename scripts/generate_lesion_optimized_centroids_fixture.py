"""Generate the RAS synthetic input for optimized LesionTOADS centroids."""

from pathlib import Path

import nibabel as nib
import numpy as np


BASE = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic"
SOURCE = BASE / "lesion_variances" / "variance_inputs.nii.gz"
ROOT = BASE / "lesion_optimized_centroids"


def main() -> None:
    source = nib.load(str(SOURCE))
    frames = np.asarray(source.dataobj, dtype=np.float32)
    x, y, z = np.indices(frames.shape[:3])
    frames[..., 0] = ((x + 2 * y + z) % 4).astype(np.float32)
    extra_membership = np.float32(0.1) + (
        ((x * 5 + y * 7 + z * 3 + 2) % 11).astype(np.float32) * np.float32(0.04)
    )
    packed = np.concatenate(
        (frames[..., :5], extra_membership[..., None], frames[..., 5:]), axis=3
    )
    # One exact-zero T1 residual (also with the field applied) and one tiny
    # uncorrected residual whose square remains below the source ZERO cutoff.
    packed[1, 0, 0, 6] = np.float32(0.55)
    packed[0, 0, 2, 6] = np.float32(1.050001)
    ROOT.mkdir(parents=True, exist_ok=True)
    result = nib.Nifti1Image(packed, source.affine, source.header.copy())
    result.set_qform(source.get_qform(), code=int(source.header["qform_code"]))
    result.set_sform(source.get_sform(), code=int(source.header["sform_code"]))
    nib.save(result, str(ROOT / "optimized_centroid_inputs.nii.gz"))


if __name__ == "__main__":
    main()
