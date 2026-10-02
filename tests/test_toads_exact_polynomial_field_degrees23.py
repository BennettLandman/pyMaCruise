from __future__ import annotations

import json
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_inhomogeneity_exact_fit import fit_toads_exact_chebyshev_field


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "toads_exact_polynomial_field_degrees23"
HARNESS = ROOT / "java_reference" / "ToadsExactPolynomialFieldReferenceHarness.java"
UPSTREAM = ROOT / "java_reference" / "upstream"


def _phantom():
    loaded = [nib.load(str(FIXTURE / f"{name}.nii.gz")) for name in ("image", "memberships", "segmentation", "initial_field")]
    for image in loaded:
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    arrays = [np.asarray(image.dataobj) for image in loaded]
    centers = np.asarray(json.loads((FIXTURE / "centroids.json").read_text()), dtype=np.float32)
    return arrays[0], arrays[1], arrays[2].astype(np.uint8), centers, arrays[3]


def _write_input(path: Path, labels, image, memberships, field, centers) -> None:
    path.write_bytes(b"".join((
        np.asarray(labels, dtype=np.int8).tobytes(order="C"),
        np.asarray(image, dtype=">f4").tobytes(order="C"),
        np.asarray(memberships, dtype=">f4").tobytes(order="C"),
        np.asarray(field, dtype=">f4").tobytes(order="C"),
        np.asarray(centers, dtype=">f4").tobytes(order="C"),
    )))


@pytest.mark.parametrize("degree", [2, 3])
def test_exact_polynomial_field_matches_selected_java_reference(degree: int, tmp_path: Path) -> None:
    image, memberships, segmentation, centers, initial_field = _phantom()
    input_path = tmp_path / f"degree{degree}.bin"
    _write_input(input_path, segmentation, image, memberships, initial_field, centers)
    build = tmp_path / "classes"
    build.mkdir()
    java_sources = sorted((UPSTREAM / "Jama").glob("*.java")) + sorted((UPSTREAM / "Jama" / "util").glob("*.java"))
    java_sources += [
        UPSTREAM / "gov/nih/mipav/model/structures/jama/JamaMatrix.java",
        UPSTREAM / "edu/jhmi/rad/medic/utilities/PowerTable.java",
        HARNESS,
    ]
    subprocess.run(["javac", "-d", str(build), *map(str, java_sources)], check=True, capture_output=True, text=True)
    completed = subprocess.run(
        ["java", "-cp", str(build), "ToadsExactPolynomialFieldReferenceHarness",
         str(input_path), *(str(n) for n in image.shape), "11", str(degree)],
        check=True, capture_output=True, text=True,
    )
    java = json.loads(completed.stdout)

    field = initial_field.copy()
    result = fit_toads_exact_chebyshev_field(
        image, memberships, segmentation, centers, field, degree=degree,
    )
    assert java["degree"] == degree
    assert int(result.correction_mask.sum()) == int(np.count_nonzero(segmentation))
    np.testing.assert_array_equal(result.normal_matrix.ravel(), np.asarray(java["normal"], dtype=np.float64))
    np.testing.assert_array_equal(result.right_hand_side, np.asarray(java["rhs"], dtype=np.float64))
    coefficients = np.asarray(java["coefficients"], dtype=np.float64)
    np.testing.assert_allclose(result.coefficients, coefficients, rtol=0, atol=3e-11)
    expected = np.asarray(java["field_bits"], dtype=np.int32).view(np.float32).reshape(image.shape)
    np.testing.assert_allclose(result.field, expected, rtol=0, atol=2e-6)
    np.testing.assert_array_equal(result.field[~result.correction_mask], 1.0)
    assert abs(float(np.mean(result.field[result.correction_mask], dtype=np.float64)) - 1.0) < 2e-6


@pytest.mark.parametrize("degree", [True, 0, 4, 1.5])
def test_exact_polynomial_field_rejects_unsupported_degree(degree) -> None:
    image, memberships, segmentation, centers, field = _phantom()
    with pytest.raises(ValueError, match="degree must be"):
        fit_toads_exact_chebyshev_field(
            image, memberships, segmentation, centers, field, degree=degree,
        )
