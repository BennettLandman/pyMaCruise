from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_inhomogeneity_transition import (
    update_lesion_toads_exact_inhomogeneity_fields,
    update_lesion_toads_inhomogeneity_fields,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_inhomogeneity_transition"
HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsInhomogeneityTransitionReferenceHarness.java"
POLYNOMIAL_HARNESS = Path(__file__).parent / "java_reference" / "ToadsExactPolynomialFieldReferenceHarness.java"
POLYNOMIAL_UPSTREAM = Path(__file__).parent / "java_reference" / "upstream"
FIELD_FIXTURE = Path(__file__).parent / "data" / "synthetic" / "lesion_inhomogeneity_13class_refit"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text(encoding="utf-8"))
HARNESS_SHA256 = "9a4e519a5d98c397c141ed36bcd35788d5527a1ba46c30c65b469dff7f85c84c"
ALGORITHM_SHA256 = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
CORRECTION_SHA256 = "6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8"


def _bits(value: float | np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def test_degree_gate_and_channel_dispatch_match_selected_java_source(tmp_path: Path) -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert REFERENCE["algorithm_source_sha256"] == ALGORITHM_SHA256
    assert REFERENCE["correction_source_sha256"] == CORRECTION_SHA256
    assert REFERENCE["source"] == (
        "AlgorithmLesionToads.runAlgorithm degree gate + "
        "InhomogeneityCorrection2.computeExact3DPolynomialField"
    )

    java = shutil.which("java")
    javac = shutil.which("javac")
    if java is None or javac is None:
        pytest.skip("Java is needed to run the pinned source-excerpt harness")
    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run([javac, "-d", str(classes), str(HARNESS)], check=True, capture_output=True, text=True)
    actual_reference = json.loads(subprocess.run(
        [java, "-cp", str(classes), "LesionToadsInhomogeneityTransitionReferenceHarness"],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert actual_reference == REFERENCE

    # Synthetic two-channel fields resemble attached mutable correction maps.
    # The callback keeps refit implementation injectable: this test isolates
    # wrapper scheduling and channel order from the independent numerical fit.
    calls: list[tuple[int, int, bool]] = []

    def refit(channel: int, field: np.ndarray, degree: int):
        calls.append((channel, degree, field is original_refs[channel]))
        if degree:
            field[1:-1, 1:-1, 1:-1] += np.float32((channel + 1) * degree / 100.0)
        return {"channel": channel, "degree": degree}

    cases = {case["name"]: case for case in REFERENCE["cases"]}
    for name in (
        "below_threshold", "at_threshold", "above_threshold", "advance_one_degree",
        "hold_at_cap", "degree_zero_noop_dispatch", "nan_ratio_holds",
    ):
        case = cases[name]
        fields = [np.ones((7, 7, 7), dtype=np.float32) for _ in range(2)]
        original_refs = tuple(fields)
        calls.clear()
        result = update_lesion_toads_inhomogeneity_fields(
            fields,
            refit,
            degree=case["starting_degree"],
            polynomial_degree=3,
            previous_energy_ratio=np.asarray(case["energy_ratio_bits"], dtype=np.uint32).view(np.float32).item(),
            max_distance=np.asarray(case["max_distance_bits"], dtype=np.uint32).view(np.float32).item(),
        )
        assert _bits(result.low_energy_ratio) == case["low_ratio_bits"]
        assert result.degree == case["degree"]
        assert [(c, d) for c, d, _ in calls] == [tuple(map(int, item.split(":"))) for item in case["calls"].split(",")]
        assert all(same_ref for _, _, same_ref in calls)
        assert all(result.fields[c] is original_refs[c] for c in range(2))
        if name == "advance_one_degree":
            np.testing.assert_allclose(fields[0][1:-1, 1:-1, 1:-1], 1.02)
            np.testing.assert_allclose(fields[1][1:-1, 1:-1, 1:-1], 1.04)

    assert calls == [(0, 1, True), (1, 1, True)]


def test_invalid_degree_or_field_state_is_rejected_before_refit() -> None:
    field = np.ones((3, 3, 3), dtype=np.float32)
    calls: list[int] = []

    def refit(channel, attached, degree):
        calls.append(channel)
        return None

    with pytest.raises(ValueError, match="cannot exceed"):
        update_lesion_toads_inhomogeneity_fields(
            [field], refit, degree=2, polynomial_degree=1,
            previous_energy_ratio=0.0, max_distance=0.0001,
        )
    with pytest.raises(ValueError, match="float32"):
        update_lesion_toads_inhomogeneity_fields(
            [field.astype(np.float64)], refit, degree=0, polynomial_degree=1,
            previous_energy_ratio=0.0, max_distance=0.0001,
        )
    assert calls == []


def test_adaptive_transition_runs_13_class_exact_refits_in_channel_order(tmp_path: Path) -> None:
    volumes = [
        nib.load(str(FIELD_FIXTURE / f"{name}.nii.gz"))
        for name in ("image", "memberships", "segmentation", "initial_field")
    ]
    for volume in volumes:
        assert nib.aff2axcodes(volume.affine) == ("R", "A", "S")
    image, memberships, segmentation, initial_field = (
        np.asarray(volume.dataobj) for volume in volumes
    )
    segmentation = segmentation.astype(np.uint8)
    centers = np.asarray(json.loads((FIELD_FIXTURE / "centroids.json").read_text()), dtype=np.float32)
    class_types = tuple(json.loads((FIELD_FIXTURE / "class_types.json").read_text()))
    images = np.stack((image, np.asarray(image * np.float32(0.83) + np.float32(0.06), dtype=np.float32)), axis=-1)
    channel_centers = np.stack(
        (centers, np.asarray(centers * np.float32(0.83) + np.float32(0.06), dtype=np.float32))
    )
    fields = [initial_field.copy(), initial_field.copy()]

    result = update_lesion_toads_exact_inhomogeneity_fields(
        images, memberships, segmentation, channel_centers, fields, class_types,
        degree=1, polynomial_degree=2, previous_energy_ratio=0.005, max_distance=0.0001,
    )

    assert result.degree == 2
    assert result.low_energy_ratio == np.float32(0.01)
    assert result.stage_order[-2:] == (
        "computeCorrectionField.channel0.degree2",
        "computeCorrectionField.channel1.degree2",
    )
    assert result.fields[0] is fields[0]
    assert result.fields[1] is fields[1]
    assert len(result.fits) == 2

    build = tmp_path / "classes"
    build.mkdir()
    sources = sorted((POLYNOMIAL_UPSTREAM / "Jama").glob("*.java"))
    sources += sorted((POLYNOMIAL_UPSTREAM / "Jama" / "util").glob("*.java"))
    sources += [
        POLYNOMIAL_UPSTREAM / "gov/nih/mipav/model/structures/jama/JamaMatrix.java",
        POLYNOMIAL_UPSTREAM / "edu/jhmi/rad/medic/utilities/PowerTable.java",
        POLYNOMIAL_HARNESS,
    ]
    subprocess.run(["javac", "-d", str(build), *map(str, sources)], check=True, capture_output=True, text=True)
    type_argument = ",".join(class_types)
    for channel, fit in enumerate(result.fits):
        input_path = tmp_path / f"lesion_refit_channel{channel}.bin"
        input_path.write_bytes(b"".join((
            segmentation.astype(np.int8).tobytes(order="C"),
            images[..., channel].astype(">f4").tobytes(order="C"),
            memberships.astype(">f4").tobytes(order="C"),
            initial_field.astype(">f4").tobytes(order="C"),
            channel_centers[channel].astype(">f4").tobytes(order="C"),
        )))
        reference = json.loads(subprocess.run(
            ["java", "-cp", str(build), "ToadsExactPolynomialFieldReferenceHarness",
             str(input_path), *(str(n) for n in image.shape), "13", "2", type_argument],
            check=True, capture_output=True, text=True,
        ).stdout)
        np.testing.assert_array_equal(fit.normal_matrix.ravel(), np.asarray(reference["normal"]))
        np.testing.assert_array_equal(fit.right_hand_side, np.asarray(reference["rhs"]))
        expected = np.asarray(reference["field_bits"], dtype=np.int32).view(np.float32).reshape(image.shape)
        np.testing.assert_allclose(fields[channel], expected, rtol=0, atol=2e-6)
