"""Create the fixed-input folded RAS fixture used by P18cye."""

from pathlib import Path
import sys

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
sys.path.insert(0, str(REPO / "src"))
from pymacruise.toads_approximate_variances import compute_default_toads_approximate_variances

SOURCE = ROOT.parents[0] / "toads_approximate_centroids"
VARIANCE_SOURCE = ROOT.parents[0] / "toads_approximate_variances"
affine = nib.load(str(SOURCE / "images.nii.gz")).affine
images = np.asarray(nib.load(str(SOURCE / "images.nii.gz")).dataobj, dtype=np.float32)[..., :1]
field = np.asarray(nib.load(str(VARIANCE_SOURCE / "field.nii.gz")).dataobj, dtype=np.float32)[..., :1]
memberships = np.asarray(nib.load(str(VARIANCE_SOURCE / "memberships.nii.gz")).dataobj, dtype=np.float32)
best = np.asarray(nib.load(str(VARIANCE_SOURCE / "best.nii.gz")).dataobj, dtype=np.uint8)
centroids2 = np.loadtxt(VARIANCE_SOURCE / "updated_centroids.txt", dtype=np.float32)
centroids = centroids2[:1]
segmentation = np.asarray(nib.load(str(VARIANCE_SOURCE / "segmentation.nii.gz")).dataobj, dtype=np.uint8)
shape = images.shape[:3]
classes = 11

std = compute_default_toads_approximate_variances(
    images, memberships, best, segmentation, centroids, intensity_field=field
).standard_deviations
std[std <= 1e-7] = np.float32(0.1)

x, y, z = np.indices(shape, dtype=np.float32)
fold = np.sin((x + 0.35 * y) * np.float32(0.37)) + np.cos((z - 0.2 * y) * np.float32(0.43))
shapes = np.empty(shape + (classes,), dtype=np.float32)
direct = np.zeros(shape + (classes,), dtype=np.float32)
for k in range(classes):
    shapes[..., k] = np.clip(np.float32(0.42 + 0.19 * fold + 0.013 * k), 0.0, 1.0)
    direct[..., k] = (segmentation == k).astype(np.float32)

# Two exactly tied GM candidates at an interior voxel make strict lower-index
# bestIndex tie order observable; their six neighbors exclude both classes.
tie = (9, 9, 9)
centroids[0, 4] = centroids[0, 3]
std[0, 4] = std[0, 3]
shapes[tie + (4,)] = shapes[tie + (3,)]
direct[tie] = 0.0
direct[tie + (3,)] = 1.0
direct[tie + (4,)] = 1.0
for dx, dy, dz in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)):
    q = (tie[0]+dx, tie[1]+dy, tie[2]+dz)
    best[q] = np.asarray([1, 2, 5, 8], dtype=np.uint8)
    memberships[q] = np.float32(0.25)

# Dedicated source ZERO/INF candidate and a no-eligible-prior voxel.
zero = (8, 8, 8)
images[zero] = np.float32(0.0)
field[zero] = np.float32(1.7)
direct[zero] = 0.0
direct[zero + (0,)] = 1.0
shapes[zero] = 0.0
for dx, dy, dz in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)):
    q = (zero[0]+dx, zero[1]+dy, zero[2]+dz)
    best[q] = np.asarray([0, 0, 0, 0], dtype=np.uint8)
    memberships[q] = np.float32(0.0)
empty = (10, 10, 10)
direct[empty] = 0.0
shapes[empty] = 0.0

# A near-ZERO optimized residual gives the selected INF reciprocal branch.
near = (7, 7, 7)
images[near] = np.float32(centroids[0, 1] + 1.0e-6)
field[near] = 1.0
direct[near] = 0.0
direct[near + (1,)] = 1.0
shapes[near] = 0.0
for dx, dy, dz in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)):
    q = (near[0]+dx, near[1]+dy, near[2]+dz)
    best[q] = np.asarray([1, 1, 1, 1], dtype=np.uint8)
    memberships[q] = np.float32(0.0)

# Early interior voxel: all six neighbors still expose their old memberships
# when it is swept, allowing each face contribution to be perturbed separately.
sweep = (1, 1, 1)
direct[sweep] = 0.0
direct[sweep + (3,)] = 1.0
shapes[sweep] = 0.0
shapes[sweep + (3,)] = 0.4
for dx, dy, dz in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)):
    q = (sweep[0]+dx, sweep[1]+dy, sweep[2]+dz)
    best[q] = np.asarray([4, 5, 6, 7], dtype=np.uint8)
    memberships[q] = np.float32(0.25)

variance_weights = np.asarray([[1.0 + (k % 3) * 0.07 for k in range(classes)]], dtype=np.float32)
exponent_by_class = np.asarray([0.0] + [0.8, 0.8] + [1.0] * 5 + [1.0] * 3, dtype=np.float32)
modality = np.asarray([0.93], dtype=np.float32)
smooth_prior = np.asarray([0.61 + 0.037 * k for k in range(classes)], dtype=np.float32)
shape_prior = np.asarray([0.4 + 0.055 * k for k in range(classes)], dtype=np.float32)

def save(name: str, array: np.ndarray) -> None:
    nib.save(nib.Nifti1Image(np.asarray(array, dtype=np.float32), affine), str(ROOT / f"{name}.nii.gz"))

save("images", images)
save("field", field)
save("memberships", memberships)
save("best", best)
save("transformed_shapes", shapes)
save("direct_relation", direct)
for name, values in (("centroids", centroids), ("standard_deviations", std),
                     ("variance_weights", variance_weights), ("exponents", exponent_by_class),
                     ("modality", modality), ("smoothing_priors", smooth_prior),
                     ("shape_priors", shape_prior)):
    np.savetxt(ROOT / f"{name}.txt", values.reshape(-1)[None, :], fmt="%.9g")
(ROOT / "parameters.txt").write_text("0.12 0.73 2.0 0.028\n")
