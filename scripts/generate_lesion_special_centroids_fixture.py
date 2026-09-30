"""Generate a compact RAS phantom for LESION/BLACKHOLE centroid routing."""

from pathlib import Path

import nibabel as nib
import numpy as np


BASE = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic"
SOURCE = BASE / "lesion_variances" / "variance_inputs.nii.gz"
ROOT = BASE / "lesion_special_centroids"


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
    # Force source branch/tie cases on the LESION and BLACKHOLE hard classes.
    packed[0, 0, 1, 2] = np.float32(0.8)  # class 1 membership ties lesion map
    packed[0, 0, 1, 10] = np.float32(0.8)
    packed[0, 0, 1, 11] = np.float32(0.2)
    packed[0, 0, 2, 3] = np.float32(0.3)  # lesion and black-hole tie above class
    packed[0, 0, 2, 10] = np.float32(0.8)
    packed[0, 0, 2, 11] = np.float32(0.8)
    packed[1, 0, 1, 3] = np.float32(0.6)  # all three tie; class wins by precedence
    packed[1, 0, 1, 10] = np.float32(0.6)
    packed[1, 0, 1, 11] = np.float32(0.6)
    packed[0, 1, 0, 3] = np.float32(0.3)  # black-hole wins
    packed[0, 1, 0, 10] = np.float32(0.2)
    packed[0, 1, 0, 11] = np.float32(0.9)
    packed[1, 0, 0, 2] = np.float32(0.2)  # lesion wins with exact zero distance
    packed[1, 0, 0, 10] = np.float32(0.9)
    packed[1, 0, 0, 11] = np.float32(0.1)
    packed[1, 0, 0, 6] = np.float32(0.85)
    ROOT.mkdir(parents=True, exist_ok=True)
    result = nib.Nifti1Image(packed, source.affine, source.header.copy())
    result.set_qform(source.get_qform(), code=int(source.header["qform_code"]))
    result.set_sform(source.get_sform(), code=int(source.header["sform_code"]))
    nib.save(result, str(ROOT / "special_centroid_inputs.nii.gz"))


if __name__ == "__main__":
    main()
