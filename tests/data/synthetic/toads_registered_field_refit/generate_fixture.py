"""Generate a compact RAS atlas/subject pair for registered field fitting."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).parent


def save(name: str, data: np.ndarray, affine: np.ndarray) -> None:
    image = nib.Nifti1Image(data, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    nib.save(image, str(ROOT / name))


def main() -> None:
    atlas_shape = (15, 13, 11)
    x, y, z = np.indices(atlas_shape, dtype=np.float32)
    topology = ((x.astype(np.int16) + 2 * y + 3 * z) % 11 + 1).astype(np.uint8)
    save("topology.nii.gz", topology, np.eye(4))
    structures = (
        ("Background", 1, "mask"), ("Csf1", 2, "optimized_csf"),
        ("Csf2", 3, "optimized_csf"), ("Gm1", 4, "optimized_gm"),
        ("Gm2", 5, "optimized_gm"), ("Gm3", 6, "optimized_gm"),
        ("Gm4", 7, "optimized_gm"), ("Gm5", 8, "optimized_gm"),
        ("Wm1", 9, "optimized_wm"), ("Wm2", 10, "optimized_wm"),
        ("Wm3", 11, "optimized_wm"),
    )
    lines = ["Structure Atlas File (edit at your own risks)", "", f"Structures:\t{len(structures)}"]
    lines.extend(f"{name}\t{label}\t{kind}" for name, label, kind in structures)
    lines.extend(("", "Topology Atlas:", "File:\ttopology.nii.gz", "Dimensions:\t15\t13\t11", "Resolutions:\t1\t1\t1", "", "Shape Atlas:\t1", "Dimensions:\t15\t13\t11", "Resolutions:\t1\t1\t1", "Orientations:\t1\t1\t6\t3"))
    for index, (name, _label, _kind) in enumerate(structures):
        shape_prior = (0.15 + 0.8 * np.exp(-((x - (2 + index % 7)) ** 2 + (y - (2 + (2 * index) % 7)) ** 2 + (z - (1 + (3 * index) % 6)) ** 2) / 24.0)).astype(np.float32)
        filename = f"shape-{index}.nii.gz"
        save(filename, shape_prior, np.eye(4))
        lines.extend((f"Structure: {name}", f"File:\t{filename}", "Dimensions:\t15\t13\t11", "Resolutions:\t1\t1\t1", "Center:\t0\t0\t0"))
    priors = "\t".join(f"{index / 10:.1f}" for index in range(11))
    lines.extend(("", "Intensity Atlas:\t0", f"T1_MPRAGE\t{priors}", ""))
    (ROOT / "atlas-11.txt").write_text("\n".join(lines), encoding="utf-8")

    subject_shape = (28, 27, 25)
    sx, sy, sz = np.indices(subject_shape, dtype=np.float32)
    mask = np.zeros(subject_shape, dtype=np.uint8)
    mask[5:23, 4:23, 3:22] = 1
    t1 = np.zeros(subject_shape, dtype=np.float32)
    inside = mask > 0
    t1[inside] = (20.0 + sx[inside] * 0.7 + sy[inside] * 0.31 + sz[inside] * 0.19 + ((sx[inside] + 2 * sy[inside] + sz[inside]) % 7) * 0.13).astype(np.float32)
    subject_affine = np.diag((0.8, 1.2, 1.5, 1.0))
    save("t1.nii.gz", t1, subject_affine)
    save("brain-mask.nii.gz", mask, subject_affine)
    (ROOT / "README.md").write_text(
        "# Registered TOADS field-fit integration fixture\n\n"
        "A deterministic default 11-class atlas and asymmetric subject T1/brain-mask pair, "
        "all stored as RAS `.nii.gz` images. The atlas is compact; the subject has an "
        "anisotropic affine and synthetic folded intensity structure. This fixture exercises "
        "atlas registration, crop/constructor state, and the degree-three correction fit; it "
        "is not a clinical image or whole-TOADS reference. Regenerate with `python generate_fixture.py`.\n"
    )


if __name__ == "__main__":
    main()
