from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_alignment_startup import register_toads_from_atlas_assets
from pymacruise.toads_registered_field_refit import fit_registered_toads_correction_fields
from test_toads_registered_atlas_membership_integration import _synthetic_subject, _write_default_atlas


ROOT = Path(__file__).parent
HARNESS = ROOT / "java_reference" / "ToadsExactPolynomialFieldReferenceHarness.java"
UPSTREAM = ROOT / "java_reference" / "upstream"


def _field_input_bytes(labels, image, memberships, field, centers) -> bytes:
    return b"".join((
        np.asarray(labels, dtype=np.int8).tobytes(order="C"),
        np.asarray(image, dtype=">f4").tobytes(order="C"),
        np.asarray(memberships, dtype=">f4").tobytes(order="C"),
        np.asarray(field, dtype=">f4").tobytes(order="C"),
        np.asarray(centers, dtype=">f4").tobytes(order="C"),
    ))


def test_registered_startup_field_fit_matches_degree_three_java_reference(tmp_path: Path) -> None:
    descriptor = _write_default_atlas(tmp_path)
    t1, mask = _synthetic_subject()
    run = register_toads_from_atlas_assets(
        t1, mask, descriptor, levels=1, max_iterations=1,
        refinement_iterations=1, padding=2,
    )
    result = fit_registered_toads_correction_fields(run)
    images = np.moveaxis(run.template_refresh.recropped_images, 0, -1)
    memberships = run.segmentation_initialization.memberships
    segmentation = run.segmentation_constructor.segmentation
    centroids = run.segmentation_initialization.centroids

    build = tmp_path / "classes"
    build.mkdir()
    java_sources = sorted((UPSTREAM / "Jama").glob("*.java")) + sorted((UPSTREAM / "Jama" / "util").glob("*.java"))
    java_sources += [UPSTREAM / "gov/nih/mipav/model/structures/jama/JamaMatrix.java", UPSTREAM / "edu/jhmi/rad/medic/utilities/PowerTable.java", HARNESS]
    subprocess.run(["javac", "-d", str(build), *map(str, java_sources)], check=True, capture_output=True, text=True)
    input_path = tmp_path / "registered_exact_field.bin"
    input_path.write_bytes(_field_input_bytes(
        segmentation, images[..., 0], memberships, np.ones(segmentation.shape, dtype=np.float32), centroids[0],
    ))
    java = json.loads(subprocess.run(
        ["java", "-cp", str(build), "ToadsExactPolynomialFieldReferenceHarness",
         str(input_path), *(str(n) for n in segmentation.shape), "11", "3"],
        check=True, capture_output=True, text=True,
    ).stdout)

    assert result.degree == 3
    assert result.channel_order == (0,)
    assert result.stage_order == ("registered_startup.computeExact3DPolynomialField.degree3",)
    expected_field = np.asarray(java["field_bits"], dtype=np.int32).view(np.float32).reshape(segmentation.shape)
    np.testing.assert_array_equal(result.fits[0].normal_matrix.ravel(), np.asarray(java["normal"], dtype=np.float64))
    np.testing.assert_array_equal(result.fits[0].right_hand_side, np.asarray(java["rhs"], dtype=np.float64))
    np.testing.assert_allclose(result.fields[..., 0], expected_field, rtol=0, atol=2e-6)
    np.testing.assert_array_equal(result.fields[..., 0][segmentation == 0], 1.0)


@pytest.mark.parametrize("degree", [0, 4, True])
def test_registered_field_refit_rejects_unsupported_degree(degree) -> None:
    with pytest.raises(ValueError, match="degree must be"):
        fit_registered_toads_correction_fields(object(), degree=degree)
