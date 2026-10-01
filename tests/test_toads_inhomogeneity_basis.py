"""Compare the selected TOADS polynomial-basis method excerpt to Java."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_inhomogeneity import toads_chebyshev_basis


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_inhomogeneity_basis"
JAVA = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "ToadsInhomogeneityBasisReferenceHarness.java"
HARNESS_SHA256 = "de3c8efa838db79e407846c8cee4d413fd0fc8ef8b77041371e9b5dcaaa53390"
SOURCE_SHA256 = "6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8"


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).tolist()


def test_basis_matches_source_formula_for_all_fixture_voxels_and_degrees() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert JAVA["source_sha256"] == SOURCE_SHA256
    assert JAVA["dimensions"] == [5, 4, 3]
    image = nib.load(str(ROOT / "grid.nii.gz"))
    assert image.shape == (5, 4, 3)
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), (.8, 1.1, 1.7), rtol=0, atol=1e-6)

    dimensions = image.shape
    for degree, reference in enumerate(JAVA["degrees"], start=1):
        assert reference["degree"] == degree
        actual = [
            _bits(toads_chebyshev_basis((x, y, z), dimensions, degree=degree))
            for x in range(dimensions[0])
            for y in range(dimensions[1])
            for z in range(dimensions[2])
        ]
        assert actual == reference["values"]


@pytest.mark.parametrize(
    ("point", "dimensions", "degree", "message"),
    [
        ((0, 0), (5, 4, 3), 1, "three integer indices"),
        ((0, 0, 0), (1, 4, 3), 1, "sizes >= 2"),
        ((5, 0, 0), (5, 4, 3), 1, "within"),
        ((0, 0, 0), (5, 4, 3), 0, "degree must"),
        ((0, 0, 0), (5, 4, 3), 4, "degree must"),
    ],
)
def test_basis_rejects_values_outside_the_supported_source_contract(point, dimensions, degree, message) -> None:
    with pytest.raises(ValueError, match=message):
        toads_chebyshev_basis(point, dimensions, degree=degree)
