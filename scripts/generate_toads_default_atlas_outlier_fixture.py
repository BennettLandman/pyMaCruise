"""Build the committed P18czg fixture from the P18czf synthetic atlas case."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parents[1] / "tests/data/synthetic"
SOURCE = ROOT / "toads_default_atlas_dispatch"
OUTPUT = ROOT / "toads_default_atlas_outlier"


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name in (
        "images", "field", "psf_bias", "folded_ribbon", "initial_memberships",
        "transformed_shapes", "closest_labels",
    ):
        image = nib.load(str(SOURCE / f"{name}.nii.gz"))
        nib.save(image, str(OUTPUT / f"{name}.nii.gz"))

    ribbon = np.asarray(nib.load(str(SOURCE / "folded_ribbon.nii.gz")).dataobj)
    nx, ny, nz = ribbon.shape
    x, y, z = np.ogrid[:nx, :ny, :nz]
    index = (2 * x + 3 * y + z) % 9
    outlier_map = np.asarray(np.float32(0.12) + index.astype(np.float32) * np.float32(0.047), dtype=np.float32)
    reference = nib.load(str(SOURCE / "images.nii.gz"))
    outlier_image = nib.Nifti1Image(outlier_map, reference.affine, header=reference.header.copy())
    outlier_image.set_qform(reference.affine, code=1)
    outlier_image.set_sform(reference.affine, code=1)
    outlier_image.header.set_data_dtype(np.float32)
    outlier_image.header.set_xyzt_units("mm")
    nib.save(outlier_image, str(OUTPUT / "initial_outlier_membership.nii.gz"))
    print(f"wrote P18czg synthetic inputs to {OUTPUT}")


if __name__ == "__main__":
    main()
