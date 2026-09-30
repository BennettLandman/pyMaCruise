"""Create the mixed-category LesionTOADS centroid reference phantom."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "lesion_mixed_centroids"


def make_inputs():
    nx, ny, nz, classes, channels = 6, 5, 4, 13, 2
    labels = np.empty((nx, ny, nz), dtype=np.int16)
    memberships = np.empty((nx, ny, nz, classes), dtype=np.float32)
    images = np.empty((nx, ny, nz, channels), dtype=np.float32)
    fields = np.empty_like(images)
    outlier = np.empty((nx, ny, nz), dtype=np.float32)
    lesion = np.empty_like(outlier)
    black = np.empty_like(outlier)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                labels[x, y, z] = (x + 2 * y + z) % 12  # class 12 is absent
                for k in range(classes):
                    term = np.float32(((x * 3 + y * 5 + z * 7 + k * 2) % 9) * np.float32(0.055))
                    memberships[x, y, z, k] = np.float32(np.float32(0.08) + term)
                for c in range(channels):
                    image_term = np.float32(((x * 11 + y * 7 + z * 3 + c * 5) % 23) * np.float32(0.12))
                    images[x, y, z, c] = np.float32(np.float32(0.25) + image_term)
                    field_term = np.float32(((x * 2 + y * 3 + z * 4 + c) % 5) * np.float32(0.11))
                    fields[x, y, z, c] = np.float32(np.float32(0.75) + field_term)
                outlier[x, y, z] = np.float32(((x + 3 * y + z) % 5) * np.float32(0.13))
                lesion_term = np.float32(((x + 2 * y + z) % 9) * np.float32(0.07))
                lesion[x, y, z] = np.float32(np.float32(0.04) + lesion_term)
                black_term = np.float32(((2 * x + y + 3 * z) % 9) * np.float32(0.075))
                black[x, y, z] = np.float32(np.float32(0.03) + black_term)
                if (x + y + z) % 6 == 0:
                    black[x, y, z] = lesion[x, y, z]
    # Find one voxel of each hard class, then force special auxiliary routes,
    # a tie, and a zero/below-ZERO case without changing the class coverage.
    positions = {
        k: tuple(int(v) for v in np.argwhere(labels == k)[0]) for k in range(12)
    }
    x, y, z = positions[7]
    memberships[x, y, z, 7], lesion[x, y, z], black[x, y, z] = 0.8, 0.8, 0.2
    x, y, z = positions[8]
    memberships[x, y, z, 8], lesion[x, y, z], black[x, y, z] = 0.3, 0.8, 0.8
    x, y, z = positions[9]
    memberships[x, y, z, 9], lesion[x, y, z], black[x, y, z] = 0.6, 0.6, 0.6
    x, y, z = positions[10]
    memberships[x, y, z, 10], lesion[x, y, z], black[x, y, z] = 0.3, 0.2, 0.9
    x, y, z = tuple(int(v) for v in np.argwhere(labels == 9)[1])
    memberships[x, y, z, 9] = np.float32(1e-6)
    lesion[x, y, z] = black[x, y, z] = np.float32(0.0)
    # Values in the fixture are consecutive 3-D frames: labels, memberships,
    # channels of images and fields, then outlier, lesion, and black-hole maps.
    return labels, memberships, images, fields, outlier, lesion, black


def main() -> None:
    labels, memberships, images, fields, outlier, lesion, black = make_inputs()
    frames = np.concatenate(
        [labels[..., None].astype(np.float32), memberships, images, fields,
         outlier[..., None], lesion[..., None], black[..., None]], axis=3
    )
    affine = np.eye(4)
    nii = nib.Nifti1Image(frames, affine)
    nii.set_qform(affine, code=1)
    nii.set_sform(affine, code=1)
    ROOT.mkdir(parents=True, exist_ok=True)
    nib.save(nii, str(ROOT / "mixed_centroid_inputs.nii.gz"))


if __name__ == "__main__":
    main()
