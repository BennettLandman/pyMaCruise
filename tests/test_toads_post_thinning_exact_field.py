from __future__ import annotations

import json
from dataclasses import replace
import hashlib
import subprocess
from pathlib import Path

import numpy as np

from pymacruise.toads_post_registration_next_thinning import toads_post_registration_next_thinning
from pymacruise.toads_post_thinning_exact_field import toads_post_thinning_exact_field_degree1
from test_toads_post_registration_next_thinning import _registered_branch


ROOT = Path(__file__).parent
HARNESS_SCRIPT = ROOT.parent / "scripts/run_java_toads_post_thinning_exact_field_reference.sh"
INHOMOGENEITY_UPSTREAM_SOURCE_SHA256 = "6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8"
WRAPPER_UPSTREAM_SOURCE_SHA256 = "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
COMPILED_POWER_TABLE_COPY_SHA256 = "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"
VANDY_POWER_TABLE_UPSTREAM_SOURCE_SHA256 = "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678"
FORMULA_HARNESS_SHA256 = "dedffe89f33b4032feaf754edce9c7f5ce74c4a628212c2bc0ce80b007ffe26d"


def _input_bytes(labels, image, memberships, field, centroids) -> bytes:
    return b"".join(
        (
            np.asarray(labels, dtype=np.int8).tobytes(order="C"),
            np.asarray(image, dtype=">f4").tobytes(order="C"),
            np.asarray(memberships, dtype=">f4").tobytes(order="C"),
            np.asarray(field, dtype=">f4").tobytes(order="C"),
            np.asarray(centroids, dtype=">f4").tobytes(order="C"),
        )
    )


def _from_bits(values: list[int]) -> np.ndarray:
    return np.asarray(values, dtype=np.int32).view(np.float32)


def test_selected_post_thinning_exact_degree_one_fit_matches_java_formula_harness(tmp_path: Path) -> None:
    # previous/current are chosen to preserve the continuing branch while
    # producing a ratio in (0.0001, sqrt(0.0001)), which schedules degree 1.
    branch = _registered_branch()
    current_energy = branch.handoff.membership_pass.objective
    previous_energy = np.float32(current_energy * np.float32(1.005))
    branch = _registered_branch(
        previous_energy=previous_energy,
        convergence_threshold=np.float32(0.0001),
    )
    ratio = branch.control.energy_ratio
    assert np.float32(0.0001) < ratio < np.sqrt(np.float32(0.0001))
    # Expand the committed small-grid RAS state by nearest-neighbor replication
    # so stride-three exact fitting has a full-rank degree-one design. The
    # following Java harness independently evaluates this expanded fixed state.
    expand = lambda array: np.repeat(np.repeat(np.repeat(array, 3, axis=0), 3, axis=1), 3, axis=2)
    grown = replace(
        branch.growing,
        segmentation=expand(branch.growing.segmentation),
        skeleton=expand(branch.growing.skeleton),
        available=expand(branch.growing.available),
        ordering=expand(branch.growing.ordering),
        class_volumes=np.asarray(branch.growing.class_volumes * np.float32(27.0), dtype=np.float32),
    )
    membership_pass = replace(
        branch.handoff.membership_pass,
        class_data_energy=expand(branch.handoff.membership_pass.class_data_energy),
        memberships=expand(branch.handoff.membership_pass.memberships),
        denominator=expand(branch.handoff.membership_pass.denominator),
    )
    handoff = replace(branch.handoff, membership_pass=membership_pass)
    branch = replace(branch, handoff=handoff, growing=grown)
    next_thinning = toads_post_registration_next_thinning(
        branch,
        np.zeros_like(branch.handoff.membership_pass.memberships, dtype=np.float32),
    )

    nx, ny, nz = next_thinning.thinning.segmentation.shape
    x, y, z = np.meshgrid(
        np.linspace(0.0, 1.0, nx, dtype=np.float32),
        np.linspace(0.0, 1.0, ny, dtype=np.float32),
        np.linspace(0.0, 1.0, nz, dtype=np.float32),
        indexing="ij",
    )
    image = np.asarray(0.35 + 0.31 * x + 0.19 * y + 0.13 * z + 0.07 * x * y, dtype=np.float32)
    centroids = np.asarray([0.0, 0.20, 0.31, 0.39, 0.46, 0.52, 0.58, 0.63, 0.69, 0.75, 0.82], dtype=np.float32)
    initial_field = np.ones((nx, ny, nz), dtype=np.float32)
    java_input = tmp_path / "exact_degree1.bin"
    java_input.write_bytes(
        _input_bytes(
            next_thinning.thinning.segmentation,
            image,
            branch.handoff.membership_pass.memberships,
            initial_field,
            centroids,
        )
    )
    completed = subprocess.run(
        ["sh", str(HARNESS_SCRIPT), str(java_input), str(nx), str(ny), str(nz), "11"],
        check=True,
        capture_output=True,
        text=True,
    )
    java = json.loads(completed.stdout)

    field = initial_field.copy()
    result = toads_post_thinning_exact_field_degree1(
        next_thinning,
        image,
        centroids,
        field,
    )
    assert result.degree == java["degree"] == 1
    assert java["reference_kind"] == "formula_harness_not_source_class_execution"
    assert java["inhomogeneity_upstream_source_sha256"] == INHOMOGENEITY_UPSTREAM_SOURCE_SHA256
    assert java["wrapper_upstream_source_sha256"] == WRAPPER_UPSTREAM_SOURCE_SHA256
    assert java["compiled_power_table_copy_sha256"] == COMPILED_POWER_TABLE_COPY_SHA256
    assert java["vandy_power_table_upstream_source_sha256"] == VANDY_POWER_TABLE_UPSTREAM_SOURCE_SHA256
    harness_source = ROOT / "java_reference/ToadsPostThinningExactFieldReferenceHarness.java"
    assert hashlib.sha256(harness_source.read_bytes()).hexdigest() == FORMULA_HARNESS_SHA256
    np.testing.assert_allclose(result.fit.normal_matrix, java["normal"], rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(result.fit.right_hand_side, java["rhs"], rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(result.fit.coefficients, java["coefficients"], rtol=2e-6, atol=2e-8)
    java_field = _from_bits(java["field_bits"]).reshape((nx, ny, nz))
    np.testing.assert_allclose(result.fit.field, java_field, rtol=0, atol=2e-7)
    assert result.stage_order[-1] == "next_iteration.computeExact3DPolynomialField.degree1"
    assert np.all(result.fit.field[~result.fit.correction_mask] == 1.0)
    assert abs(float(np.mean(result.fit.field[result.fit.correction_mask])) - 1.0) < 2e-7
