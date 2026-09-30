"""Append a deterministic OUTLIER membership map to the centroid RAS phantom."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "lesion_outlier_centroids"
SOURCE = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "lesion_variances" / "variance_inputs.nii.gz"


def main() -> None:
    image = nib.load(str(SOURCE))
    packed = np.asarray(image.dataobj, dtype=np.float32)
    x, y, z = np.indices(packed.shape[:3])
    levels = np.asarray([0.0, 0.25, 0.5, 0.75, 1.0], dtype=np.float32)
    outlier = levels[(x + 2 * y + 3 * z) % len(levels)]
    output = np.concatenate((packed, outlier[..., None]), axis=3)
    ROOT.mkdir(parents=True, exist_ok=True)
    result = nib.Nifti1Image(output, image.affine, image.header.copy())
    result.set_qform(image.get_qform(), code=int(image.header["qform_code"]))
    result.set_sform(image.get_sform(), code=int(image.header["sform_code"]))
    nib.save(result, str(ROOT / "outlier_centroid_inputs.nii.gz"))


if __name__ == "__main__":
    main()
