"""Bundle the P18cyb field and P18cyc updated centroids for P18cyd."""

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_approximate_centroids import update_default_toads_approximate_centroids


ROOT = Path(__file__).parent
SOURCE = ROOT.parent / "toads_approximate_centroids"
REFERENCE = json.loads((SOURCE / "java_reference.json").read_text())


def load(name: str) -> np.ndarray:
    return np.asarray(nib.load(str(SOURCE / f"{name}.nii.gz")).dataobj).copy()


images = load("images").astype(np.float32)
field = load("field").astype(np.float32)
memberships = load("memberships").astype(np.float32)
best = load("best")
segmentation = load("segmentation")
centroids = np.loadtxt(SOURCE / "centroids.txt", dtype=np.float32)
exponents = np.loadtxt(SOURCE / "optimized_exponents.txt", dtype=np.float32)
updated = update_default_toads_approximate_centroids(
    images, memberships, best, centroids, field, exponents
).centroids
expected = np.asarray(REFERENCE["centroid_bits"], dtype=np.int32).view(np.float32).reshape(2, 11)
np.testing.assert_array_equal(updated.view(np.int32), expected.view(np.int32))

affine = nib.load(str(SOURCE / "images.nii.gz")).affine
for name, array in (
    ("images", images),
    ("field", field),
    ("memberships", memberships),
    ("best", best),
    ("segmentation", segmentation),
):
    nib.save(nib.Nifti1Image(array, affine), ROOT / f"{name}.nii.gz")
np.savetxt(ROOT / "updated_centroids.txt", updated, fmt="%.9g")
