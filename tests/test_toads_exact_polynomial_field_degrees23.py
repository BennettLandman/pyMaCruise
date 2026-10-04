from __future__ import annotations

import json
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_inhomogeneity_exact_fit import fit_toads_exact_chebyshev_field
from pymacruise.toads_multichannel_exact_field_refit import toads_multichannel_exact_field_refit


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "toads_exact_polynomial_field_degrees23"
LESION_FIXTURE = ROOT / "data" / "synthetic" / "lesion_inhomogeneity_13class_refit"
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


@pytest.mark.parametrize("degree", [2, 3])
def test_multichannel_exact_polynomial_field_matches_java_in_source_order(degree: int, tmp_path: Path) -> None:
    image, memberships, segmentation, centers, field = _phantom()
    images = np.stack((image, np.asarray(image * np.float32(0.83) + np.float32(0.06), dtype=np.float32)), axis=-1)
    channel_centers = np.stack((centers, np.asarray(centers * np.float32(0.83) + np.float32(0.06), dtype=np.float32)))
    fields = np.stack((field, field), axis=-1)
    output = toads_multichannel_exact_field_refit(
        images, memberships, segmentation, channel_centers, fields, degree=degree,
    )
    assert output.channel_order == (0, 1)
    assert len(output.fits) == 2

    build = tmp_path / "classes"
    build.mkdir()
    java_sources = sorted((UPSTREAM / "Jama").glob("*.java")) + sorted((UPSTREAM / "Jama" / "util").glob("*.java"))
    java_sources += [UPSTREAM / "gov/nih/mipav/model/structures/jama/JamaMatrix.java", UPSTREAM / "edu/jhmi/rad/medic/utilities/PowerTable.java", HARNESS]
    subprocess.run(["javac", "-d", str(build), *map(str, java_sources)], check=True, capture_output=True, text=True)
    for channel in range(2):
        input_path = tmp_path / f"channel_{channel}_degree{degree}.bin"
        _write_input(input_path, segmentation, images[..., channel], memberships, fields[..., channel], channel_centers[channel])
        reference = json.loads(subprocess.run(
            ["java", "-cp", str(build), "ToadsExactPolynomialFieldReferenceHarness",
             str(input_path), *(str(n) for n in segmentation.shape), "11", str(degree)],
            check=True, capture_output=True, text=True,
        ).stdout)
        expected = np.asarray(reference["field_bits"], dtype=np.int32).view(np.float32).reshape(segmentation.shape)
        np.testing.assert_array_equal(output.fits[channel].normal_matrix.ravel(), np.asarray(reference["normal"]))
        np.testing.assert_array_equal(output.fits[channel].right_hand_side, np.asarray(reference["rhs"]))
        np.testing.assert_allclose(output.fields[..., channel], expected, rtol=0, atol=2e-6)
    np.testing.assert_array_equal(fields, np.stack((field, field), axis=-1))


def test_exact_polynomial_field_matches_java_for_13_class_lesion_profile(tmp_path: Path) -> None:
    loaded = [
        nib.load(str(LESION_FIXTURE / f"{name}.nii.gz"))
        for name in ("image", "memberships", "segmentation", "initial_field")
    ]
    for volume in loaded:
        assert nib.aff2axcodes(volume.affine) == ("R", "A", "S")
    image, lesion_memberships, lesion_segmentation, initial_field = (
        np.asarray(volume.dataobj) for volume in loaded
    )
    lesion_segmentation = lesion_segmentation.astype(np.uint8)
    lesion_centers = np.asarray(json.loads((LESION_FIXTURE / "centroids.json").read_text()), dtype=np.float32)
    descriptor_types = tuple(json.loads((LESION_FIXTURE / "class_types.json").read_text()))
    field = initial_field.copy()
    with pytest.raises(ValueError, match="class_types is required"):
        fit_toads_exact_chebyshev_field(
            image, lesion_memberships, lesion_segmentation, lesion_centers,
            initial_field.copy(), degree=2,
        )
    with pytest.raises(ValueError, match="one nonempty type name"):
        fit_toads_exact_chebyshev_field(
            image, lesion_memberships, lesion_segmentation, lesion_centers,
            initial_field.copy(), degree=2, class_types=("mask", "obj"),
        )
    result = fit_toads_exact_chebyshev_field(
        image, lesion_memberships, lesion_segmentation, lesion_centers, field,
        degree=2, class_types=descriptor_types,
    )

    input_path = tmp_path / "lesion_13class_degree2.bin"
    _write_input(input_path, lesion_segmentation, image, lesion_memberships, initial_field, lesion_centers)
    build = tmp_path / "classes"
    build.mkdir()
    java_sources = sorted((UPSTREAM / "Jama").glob("*.java")) + sorted((UPSTREAM / "Jama" / "util").glob("*.java"))
    java_sources += [UPSTREAM / "gov/nih/mipav/model/structures/jama/JamaMatrix.java", UPSTREAM / "edu/jhmi/rad/medic/utilities/PowerTable.java", HARNESS]
    subprocess.run(["javac", "-d", str(build), *map(str, java_sources)], check=True, capture_output=True, text=True)
    class_type_argument = ",".join(descriptor_types)
    reference = json.loads(subprocess.run(
        ["java", "-cp", str(build), "ToadsExactPolynomialFieldReferenceHarness",
         str(input_path), *(str(n) for n in image.shape), "13", "2", class_type_argument],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert reference["degree"] == 2
    assert int(result.correction_mask.sum()) == int(np.count_nonzero(lesion_segmentation))
    np.testing.assert_array_equal(result.normal_matrix.ravel(), np.asarray(reference["normal"]))
    np.testing.assert_array_equal(result.right_hand_side, np.asarray(reference["rhs"]))
    np.testing.assert_allclose(result.coefficients, np.asarray(reference["coefficients"]), rtol=0, atol=3e-11)
    expected = np.asarray(reference["field_bits"], dtype=np.int32).view(np.float32).reshape(image.shape)
    np.testing.assert_allclose(result.field, expected, rtol=0, atol=2e-6)
    np.testing.assert_array_equal(result.field[~result.correction_mask], 1.0)

    images = np.stack((image, np.asarray(image * np.float32(0.83) + np.float32(0.06), dtype=np.float32)), axis=-1)
    channel_centers = np.stack(
        (lesion_centers, np.asarray(lesion_centers * np.float32(0.83) + np.float32(0.06), dtype=np.float32))
    )
    channel_fields = np.stack((initial_field, initial_field), axis=-1)
    multichannel = toads_multichannel_exact_field_refit(
        images, lesion_memberships, lesion_segmentation, channel_centers, channel_fields,
        degree=2, class_types=descriptor_types,
    )
    assert multichannel.channel_order == (0, 1)
    for channel in range(2):
        if channel == 0:
            java_channel = reference
        else:
            channel_path = tmp_path / "lesion_13class_degree2_channel1.bin"
            _write_input(
                channel_path, lesion_segmentation, images[..., channel],
                lesion_memberships, initial_field, channel_centers[channel],
            )
            java_channel = json.loads(subprocess.run(
                ["java", "-cp", str(build), "ToadsExactPolynomialFieldReferenceHarness",
                 str(channel_path), *(str(n) for n in image.shape), "13", "2", class_type_argument],
                check=True, capture_output=True, text=True,
            ).stdout)
        fit = multichannel.fits[channel]
        np.testing.assert_array_equal(fit.normal_matrix.ravel(), np.asarray(java_channel["normal"]))
        np.testing.assert_array_equal(fit.right_hand_side, np.asarray(java_channel["rhs"]))
        expected_channel = np.asarray(java_channel["field_bits"], dtype=np.int32).view(np.float32).reshape(image.shape)
        np.testing.assert_allclose(multichannel.fields[..., channel], expected_channel, rtol=0, atol=2e-6)
    np.testing.assert_array_equal(channel_fields, np.stack((initial_field, initial_field), axis=-1))
