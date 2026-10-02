from __future__ import annotations

import json
import subprocess
from dataclasses import replace
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_post_registration_next_thinning import toads_post_registration_next_thinning
from pymacruise.toads_post_thinning_exact_field import toads_post_thinning_exact_field_update
from pymacruise.toads_thinning import ToadsThinningPass
from test_toads_post_registration_next_thinning import _registered_branch


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "toads_exact_polynomial_field_degrees23"
HARNESS = ROOT / "java_reference" / "ToadsExactPolynomialFieldReferenceHarness.java"
UPSTREAM = ROOT / "java_reference" / "upstream"


def _load_inputs():
    arrays = [np.asarray(nib.load(str(FIXTURE / f"{name}.nii.gz")).dataobj) for name in (
        "image", "memberships", "segmentation", "initial_field",
    )]
    centers = np.asarray(json.loads((FIXTURE / "centroids.json").read_text()), dtype=np.float32)
    return arrays[0], arrays[1], arrays[2].astype(np.int16), arrays[3], centers


def _java_fit(tmp_path: Path, degree: int, labels, image, memberships, field, centers):
    input_path = tmp_path / f"degree{degree}.bin"
    input_path.write_bytes(b"".join((
        np.asarray(labels, dtype=np.int8).tobytes(order="C"),
        np.asarray(image, dtype=">f4").tobytes(order="C"),
        np.asarray(memberships, dtype=">f4").tobytes(order="C"),
        np.asarray(field, dtype=">f4").tobytes(order="C"),
        np.asarray(centers, dtype=">f4").tobytes(order="C"),
    )))
    build = tmp_path / f"classes{degree}"
    build.mkdir()
    sources = sorted((UPSTREAM / "Jama").glob("*.java")) + sorted((UPSTREAM / "Jama" / "util").glob("*.java"))
    sources += [UPSTREAM / "gov/nih/mipav/model/structures/jama/JamaMatrix.java", UPSTREAM / "edu/jhmi/rad/medic/utilities/PowerTable.java", HARNESS]
    subprocess.run(["javac", "-d", str(build), *map(str, sources)], check=True, capture_output=True, text=True)
    completed = subprocess.run(
        ["java", "-cp", str(build), "ToadsExactPolynomialFieldReferenceHarness",
         str(input_path), *(str(v) for v in labels.shape), "11", str(degree)],
        check=True, capture_output=True, text=True,
    )
    return json.loads(completed.stdout)


def _next_thinning_with_fixture(current_degree: int):
    base = _registered_branch()
    energy = base.handoff.membership_pass.objective
    base = _registered_branch(
        previous_energy=np.float32(energy * np.float32(1.005)),
        convergence_threshold=np.float32(0.0001),
    )
    assert np.float32(0.0001) < base.control.energy_ratio < np.float32(0.01)
    next_thinning = toads_post_registration_next_thinning(
        base, np.zeros_like(base.handoff.membership_pass.memberships, dtype=np.float32),
    )
    image, memberships, segmentation, _field, _centers = _load_inputs()
    membership_pass = replace(base.handoff.membership_pass, memberships=memberships)
    handoff = replace(base.handoff, membership_pass=membership_pass)
    branch = replace(base, handoff=handoff)
    thinning = replace(next_thinning.thinning, segmentation=segmentation)
    return replace(next_thinning, branch=branch, thinning=thinning)


@pytest.mark.parametrize("current_degree", [0, 1, 2])
def test_exact_field_ramp_updates_one_degree_per_eligible_iteration(current_degree: int, tmp_path: Path) -> None:
    image, memberships, segmentation, initial_field, centers = _load_inputs()
    next_thinning = _next_thinning_with_fixture(current_degree)
    field = initial_field.copy()
    result = toads_post_thinning_exact_field_update(
        next_thinning, image, centers, field, current_degree=current_degree,
    )
    selected_degree = current_degree + 1
    java = _java_fit(tmp_path, selected_degree, segmentation, image, memberships, initial_field, centers)
    expected = np.asarray(java["field_bits"], dtype=np.int32).view(np.float32).reshape(image.shape)
    assert result.previous_degree == current_degree
    assert result.degree == selected_degree
    assert result.fit is not None
    assert result.stage_order[-1] == f"next_iteration.computeExact3DPolynomialField.degree{selected_degree}"
    np.testing.assert_allclose(result.fit.field, expected, rtol=0, atol=2e-6)


def test_exact_field_ramp_degree_zero_is_source_noop() -> None:
    base = _registered_branch()
    branch = replace(base, control=replace(base.control, energy_ratio=np.float32(1.0)))
    next_thinning = toads_post_registration_next_thinning(
        branch, np.zeros_like(branch.handoff.membership_pass.memberships, dtype=np.float32),
    )
    image, _memberships, _segmentation, initial_field, centers = _load_inputs()
    field = initial_field.copy()
    before = field.copy()
    result = toads_post_thinning_exact_field_update(
        next_thinning, image, centers, field, current_degree=0,
    )
    assert result.degree == result.previous_degree == 0
    assert result.fit is None
    assert result.stage_order[-1] == "next_iteration.computeExact3DPolynomialField.degree0.noop"
    np.testing.assert_array_equal(field, before)
