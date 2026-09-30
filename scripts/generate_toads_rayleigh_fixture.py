"""Create true-RAS NIfTI scalar input pairs for Java Rayleigh-kernel tests."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


OUTPUT = Path(__file__).parents[1] / "tests/data/synthetic/toads_rayleigh_distance"


def _f32(bits: int) -> np.float32:
    return np.asarray([bits], dtype=np.uint32).view(np.float32)[0]


def main() -> None:
    tiny = _f32(1)
    largest = _f32(0x7F7FFFFF)
    nan = _f32(0x7FC00000)
    pos_inf = _f32(0x7F800000)
    neg_inf = _f32(0xFF800000)
    neg_zero = _f32(0x80000000)
    cases = [
        ("positive_zero", np.float32(0), np.float32(1)),
        ("negative_zero", neg_zero, np.float32(1)),
        ("negative_finite", np.float32(-1.25), np.float32(0.5)),
        ("ordinary_positive", np.float32(1.25), np.float32(0.5)),
        ("positive_val_zero_sigma", np.float32(1.25), np.float32(0)),
        ("positive_val_negative_zero_sigma", np.float32(1.25), neg_zero),
        ("positive_val_negative_sigma", np.float32(1.25), np.float32(-0.5)),
        ("subnormal_val", tiny, np.float32(1)),
        ("subnormal_sigma", np.float32(1), tiny),
        ("large_equal", largest, largest),
        ("large_val_small_sigma", largest, np.float32(0.5)),
        ("nan_val", nan, np.float32(1)),
        ("positive_infinity_val", pos_inf, np.float32(1)),
        ("negative_infinity_val", neg_inf, np.float32(1)),
        ("finite_val_infinite_sigma", np.float32(1), pos_inf),
        ("finite_val_nan_sigma", np.float32(1), nan),
        ("positive_infinity_pair", pos_inf, pos_inf),
    ]
    values = np.empty((len(cases), 1, 1, 2), dtype=np.float32)
    for index, (_, value, sigma) in enumerate(cases):
        values[index, 0, 0] = (value, sigma)
    affine = np.diag([0.8, 1.1, 1.4, 1.0]).astype(np.float64)
    affine[:3, 3] = (-6.4, -0.55, -0.7)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    image = nib.Nifti1Image(values, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    image.header.set_xyzt_units("mm")
    nib.save(image, OUTPUT / "rayleigh_inputs.nii.gz")
    (OUTPUT / "case_names.txt").write_text("\n".join(name for name, _, _ in cases) + "\n")
    print(f"wrote {len(cases)} value/sigma pairs to {OUTPUT}")


if __name__ == "__main__":
    main()
