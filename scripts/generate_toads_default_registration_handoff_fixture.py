"""Copy the folded default-atlas shape prior as the registration template."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tests/data/synthetic/toads_default_atlas_dispatch/transformed_shapes.nii.gz"
OUTPUT = ROOT / "tests/data/synthetic/toads_default_registration_handoff/atlas_shapes.nii.gz"


def main() -> None:
    src = nib.load(str(SOURCE))
    values = np.asarray(src.dataobj, dtype=np.float32)
    # The atlas template grid is deliberately larger than the subject grid.
    # Pad the folded priors by one zero voxel on each spatial face and shift
    # the RAS origin so the original support remains at the same world points.
    template = np.pad(values, ((1, 1), (1, 1), (1, 1), (0, 0)), mode="constant")
    affine = src.affine.copy()
    affine[:3, 3] -= affine[:3, :3] @ np.ones(3)
    image = nib.Nifti1Image(template, affine, src.header.copy())
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    nib.save(image, str(OUTPUT))
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
