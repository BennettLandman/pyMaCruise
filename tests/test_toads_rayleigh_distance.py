from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_rayleigh import toads_rayleigh_distance


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_rayleigh_distance"
SOURCE_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"


def _java_float_bits(values: np.ndarray) -> list[int]:
    canonical = np.asarray(values, dtype=np.float32).copy()
    canonical.view(np.uint32)[np.isnan(canonical)] = np.uint32(0x7FC00000)
    return canonical.view(np.int32).ravel().tolist()


def test_rayleigh_kernel_matches_pinned_java_ieee754_reference() -> None:
    image = nib.load(str(ROOT / "rayleigh_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), [0.8, 1.1, 1.4], rtol=0, atol=1e-6)
    inputs = np.asarray(image.dataobj, dtype=np.float32)
    assert inputs.shape == (17, 1, 1, 2)
    names = (ROOT / "case_names.txt").read_text().splitlines()
    assert len(names) == inputs.shape[0]

    reference = json.loads((ROOT / "java_reference.json").read_text())
    assert reference["source_sha256"] == SOURCE_SHA
    assert len(reference["cases"]) == len(names)
    values = inputs[:, 0, 0, 0]
    sigmas = inputs[:, 0, 0, 1]
    assert _java_float_bits(values) == [case["val_bits"] for case in reference["cases"]]
    assert _java_float_bits(sigmas) == [case["sig_bits"] for case in reference["cases"]]

    actual = toads_rayleigh_distance(values, sigmas)
    assert actual.dtype == np.float32
    assert _java_float_bits(actual) == [case["output_bits"] for case in reference["cases"]]

    by_name = dict(zip(names, actual.tolist(), strict=True))
    for name in ("positive_zero", "negative_zero", "negative_finite", "negative_infinity_val"):
        assert np.float32(by_name[name]).view(np.uint32) == np.uint32(0x3F800000)


def test_rayleigh_kernel_broadcasts_and_returns_float32_scalar() -> None:
    scalar = toads_rayleigh_distance(np.float32(1.25), np.float32(0.5))
    assert isinstance(scalar, np.float32)
    array = toads_rayleigh_distance(np.asarray([0.0, 1.25], dtype=np.float32), np.float32(0.5))
    assert array.dtype == np.float32
    np.testing.assert_array_equal(array[0], np.float32(1.0))
    np.testing.assert_array_equal(array[1], scalar)
