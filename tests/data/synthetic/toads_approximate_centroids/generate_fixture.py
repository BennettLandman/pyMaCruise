"""Bundle a fixed centroid-update state using the Java-validated P18cyb fields."""

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_inhomogeneity_state import ToadsInhomogeneityFieldState
from pymacruise.toads_inhomogeneity_transition import update_toads_approximate_field_transition


ROOT = Path(__file__).parent
SOURCE = ROOT.parent / "toads_inhomogeneity_transition"


def load(name: str) -> np.ndarray:
    return np.asarray(nib.load(str(SOURCE / f"{name}.nii.gz")).dataobj).copy()


images = load("images").astype(np.float32)
memberships = load("memberships").astype(np.float32)
best = load("best")
segmentation = load("segmentation")
fields4 = load("attached_fields").astype(np.float32)
centroids = np.loadtxt(SOURCE / "centroids.txt", dtype=np.float32)
field_refs = [fields4[..., channel].copy() for channel in range(images.shape[3])]
state = ToadsInhomogeneityFieldState(channel_count=len(field_refs))
for channel, field in enumerate(field_refs):
    state.add_inhomogeneity_correction(field, channel)
attached = tuple(state.field_by_channel)
first = update_toads_approximate_field_transition(
    images,
    memberships,
    best,
    segmentation,
    centroids,
    attached,
    degree=1,
    polynomial_degree=3,
    previous_energy_ratio=0.4,
    low_energy_ratio=0.5,
)
second = update_toads_approximate_field_transition(
    images,
    memberships,
    best,
    segmentation,
    centroids,
    attached,
    degree=first.degree,
    polynomial_degree=3,
    previous_energy_ratio=0.4,
    low_energy_ratio=0.5,
)
assert first.degree == 2 and second.degree == 3

source_affine = nib.load(str(SOURCE / "images.nii.gz")).affine
for name, array in (
    ("images", images),
    ("memberships", memberships),
    ("best", best),
    ("segmentation", segmentation),
    ("field", np.stack(attached, axis=-1)),
):
    nib.save(nib.Nifti1Image(array, source_affine), ROOT / f"{name}.nii.gz")
np.savetxt(ROOT / "centroids.txt", centroids, fmt="%.9g")
np.savetxt(
    ROOT / "optimized_exponents.txt",
    np.asarray([[0.8, 1.0, 1.0], [1.0, 1.0, 1.0]], dtype=np.float32),
    fmt="%.9g",
)
