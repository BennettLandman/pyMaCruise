"""Generate the fixed, folded two-channel RAS input state for P39. """

from __future__ import annotations

import json
import hashlib
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1] / "tests/data/synthetic/lesion_initial_thinning_composition"
N, CLASSES = 17, 13


def fixture_state() -> tuple[np.ndarray, dict[str, object]]:
    x, y, z = np.indices((N, N, N), dtype=np.int32)
    fold = 8 + ((3 * x + 5 * y) % 5)
    ti = (x * 7 + y * 5 + z * 3 + ((x * y) % 3) * 2) % 16
    fi = (x * 3 + y * 7 + z * 5 + ((y * z) % 4) * 3) % 15
    folded = (z >= fold).astype(np.float32)
    t1 = np.float32(np.float32(0.4) + np.float32(0.025) * ti.astype(np.float32))
    t1 = np.float32(t1 + np.float32(0.018) * folded)
    flair = np.float32(np.float32(0.55) + np.float32(0.02) * fi.astype(np.float32))
    flair = np.float32(flair + np.float32(0.027) * folded)
    images = np.stack((t1, flair), axis=-1).astype(np.float32)

    k = np.arange(CLASSES, dtype=np.float32)
    shape = np.empty((N, N, N, CLASSES), dtype=np.float32)
    direct = np.empty_like(shape)
    for ix in range(N):
        for iy in range(N):
            for iz in range(N):
                shape[ix, iy, iz] = np.float32(
                    np.float32(0.4 + np.float32(0.02) * k)
                    + np.float32(0.001) * np.float32((ix + 2 * iy + 3 * iz) % 7)
                )
                direct[ix, iy, iz] = ((ix + 2 * iy + 3 * iz + np.arange(CLASSES)) % 2 == 0)

    segmentation = np.where(z < fold, 1, 12).astype(np.float32)
    border = ((x == 0) | (x == N - 1) | (y == 0) | (y == N - 1)
              | (z == 0) | (z == N - 1))
    skeleton = np.where(border, 7, -1).astype(np.float32)
    available = np.empty((N, N, N, CLASSES), dtype=np.float32)
    for cls in range(CLASSES):
        available[..., cls] = ((x + y + z + cls) % 2 == 0)
    ordering = np.zeros((N, N, N), dtype=np.float32)
    vent = np.where((x + y + z) % 4 == 0, -1, (x + 2 * y + z) % 3).astype(np.float32)
    gm = np.full((N, N, N), -1, dtype=np.float32)
    brainstem = np.where((x + z) % 5 == 0, -1, (2 * x + y + z) % 26).astype(np.float32)
    initial_memberships = np.zeros((N, N, N, CLASSES), dtype=np.float32)
    auxiliary = np.zeros((N, N, N), dtype=np.float32)
    packed = np.concatenate((
        images,
        shape,
        direct.astype(np.float32),
        segmentation[..., None], skeleton[..., None], available,
        ordering[..., None], vent[..., None], gm[..., None], brainstem[..., None],
        initial_memberships, auxiliary[..., None], auxiliary[..., None],
        auxiliary[..., None], auxiliary[..., None],
    ), axis=-1).astype(np.float32)
    names = (
        ["T1_MPRAGE", "FLAIR"]
        + [f"shape_{i}" for i in range(CLASSES)]
        + [f"direct_relation_{i}" for i in range(CLASSES)]
        + ["segmentation", "skeleton"]
        + [f"available_{i}" for i in range(CLASSES)]
        + ["ordering", "ventricle_distance", "gray_matter_distance", "brainstem_distance"]
        + [f"initial_membership_{i}" for i in range(CLASSES)]
        + ["initial_lesion", "initial_black_hole", "initial_outlier", "initial_lesion_relation"]
    )
    return packed, {
        "shape": [N, N, N, len(names)],
        "frame_names": names,
        "affine_ras_mm": np.eye(4).tolist(),
        "class_centroids": [[float(np.float32(0.2) + np.float32(0.04) * c) for c in range(CLASSES)],
                            [float(np.float32(0.3) + np.float32(0.03) * c) for c in range(CLASSES)]],
        "class_stddev": [[float(np.float32(0.3) + np.float32(0.01) * c) for c in range(CLASSES)],
                         [float(np.float32(0.4) + np.float32(0.005) * c) for c in range(CLASSES)]],
        "lesion_centroid": [0.5, 0.65],
        "lesion_stddev": [0.25, 0.3],
        "channel_weights": [1.0, 0.75],
        "shape_priors": [float(np.float32(0.5) + np.float32(0.05) * c) for c in range(CLASSES)],
        "smoothing_priors": [float(np.float32(0.8) + np.float32(0.025) * c) for c in range(CLASSES)],
        "outlier_energy": 0.1,
        "lesion_power_exponent": 1.2,
        "gm_flag_for_initial_membership": False,
        "object_types": ["lesion" if c == 10 else "object" for c in range(CLASSES)],
        "lesion_class": 10,
        "volume_ratio": 0.4,
        "thinning_limit": 1000.0,
        "use_lesion_max": False,
        "fold_definition": "segmentation is class 1 below z=8+((3*x+5*y) mod 5), class 12 above; T1 and FLAIR intensities include a fold-conditioned offset",
        "packed_sha256_float32le": hashlib.sha256(packed.astype("<f4", copy=False).tobytes()).hexdigest(),
        "image_sha256_float32le": hashlib.sha256(images.astype("<f4", copy=False).tobytes()).hexdigest(),
    }


def main() -> None:
    packed, spec = fixture_state()
    ROOT.mkdir(parents=True, exist_ok=True)
    image = nib.Nifti1Image(packed, np.eye(4))
    image.set_qform(np.eye(4), code=1)
    image.set_sform(np.eye(4), code=1)
    nib.save(image, ROOT / "initial_thinning_inputs_ras.nii.gz")
    (ROOT / "fixture.json").write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
