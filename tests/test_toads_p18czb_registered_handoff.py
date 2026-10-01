from __future__ import annotations

import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.toads_control import toads_energy_ratio
from pymacruise.toads_default_membership_registration_control import (
    continue_default_membership_through_registration_control,
)
from pymacruise.toads_p18cza_carried_degree_iteration import (
    continue_p18cyz_thinning_through_carried_degree1_iteration,
)
from pymacruise.toads_p18cyx_next_loop_thinning import (
    continue_p18cyx_branch_through_next_thinning,
)
from test_toads_p18cyx_registered_cycle import _p18cyw_state
from test_toads_p18cyz_next_loop_thinning import _make_p18cyx
from test_toads_post_thinning_exact_iteration import _bits

ROOT = Path(__file__).parent
CASE = json.loads(
    (ROOT / "data/synthetic/toads_p18czb_registered_handoff/case.json").read_text()
)
RAS_FIXTURE = ROOT / "data" / "synthetic" / CASE["source_fixture"]
REGISTRATION_HARNESS = ROOT / "java_reference/run_toads_p18cza_carried_transform_registration_reference.sh"
BRANCH_HARNESS = ROOT / "java_reference/ToadsP18cytRegisteredIterationReferenceHarness.java"


def _atlas(image: np.ndarray) -> np.ndarray:
    nx, ny, nz = image.shape
    gx, gy, gz = np.meshgrid(
        np.linspace(-1, 1, nx, dtype=np.float32),
        np.linspace(-1, 1, ny, dtype=np.float32),
        np.linspace(-1, 1, nz, dtype=np.float32),
        indexing="ij",
    )
    atlas = np.empty((*image.shape, CASE["membership_class_count"]), dtype=np.float32)
    for k in range(atlas.shape[-1]):
        radius = (
            (gx - np.float32(-0.45 + 0.09 * k)) ** 2
            + 1.3 * (gy - np.float32(0.3 * np.sin(k))) ** 2
            + 0.8 * (gz - np.float32(0.25 * np.cos(k))) ** 2
        )
        atlas[..., k] = np.asarray(
            np.exp(-radius / np.float32(0.42)) * 0.75 + 0.05,
            dtype=np.float32,
        )
    return atlas


def test_p18cza_membership_uses_carried_transform_and_matches_java_registration_control(tmp_path: Path) -> None:
    previous, image = _p18cyw_state()
    cycle = _make_p18cyx(previous, image)
    prior_branch = cycle.registered_iteration.branch
    assert prior_branch.relations is not None
    assert prior_branch.growing is not None
    shapes = np.zeros_like(prior_branch.handoff.membership_pass.memberships, dtype=np.float32)
    p18cyz = continue_p18cyx_branch_through_next_thinning(cycle, shapes)
    prior_exact = previous.exact_iteration
    p18cza = continue_p18cyz_thinning_through_carried_degree1_iteration(
        p18cyz.next_thinning,
        image,
        np.asarray(prior_exact.centroids[0], dtype=np.float32).copy(),
        np.asarray(prior_exact.field_update.fit.field, dtype=np.float32).copy(),
        previous_degree=1,
        polynomial_degree=3,
        low_energy_ratio=0.01,
        category_exponents=np.asarray([0.8, 1.0, 1.0], dtype=np.float32),
        exponent_by_channel_class=np.asarray(
            [[1.0, 0.8, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]],
            dtype=np.float32,
        ),
        modality_weights=np.asarray([0.73], dtype=np.float32),
        class_variance_weights=np.asarray(
            [[0.81, 0.93, 1.17, 0.76, 0.88, 1.02, 1.14, 0.69, 0.96, 1.21, 0.73]],
            dtype=np.float32,
        ),
        transformed_shapes=shapes,
        closest_labels=prior_branch.relations.closest_labels,
        smoothing_priors=np.asarray(
            [0.65, 0.71, 0.84, 0.91, 1.02, 1.13, 0.76, 0.88, 1.05, 1.17, 0.82],
            dtype=np.float32,
        ),
        atlas_shape_priors=np.asarray(
            [0.10, 0.08, 0.09, 0.11, 0.12, 0.09, 0.10, 0.11, 0.07, 0.06, 0.07],
            dtype=np.float32,
        ),
        smoothing=0.043,
        prior_factor=0.64,
        prior_coefficient=0.16,
        prior_scale=0.05,
    )
    membership = p18cza.exact_iteration.membership_pass
    atlas = _atlas(image)
    prior_energy = np.float32(prior_branch.handoff.membership_pass.objective)
    carried_transform = np.asarray(cycle.registration_control.handoff.transform, dtype=np.float32)
    assert carried_transform.shape == (6,) and np.isfinite(carried_transform).all()

    source = nib.load(str(RAS_FIXTURE))
    assert nib.aff2axcodes(source.affine) == tuple(CASE["source_fixture_axis_codes"])
    affine = source.affine @ np.diag([1 / CASE["spatial_repeat_factor"]] * 3 + [1])
    assert nib.aff2axcodes(affine) == ("R", "A", "S")
    assert nib.Nifti1Image(membership.memberships, affine).shape == (
        *image.shape, CASE["membership_class_count"]
    )

    handoff = continue_default_membership_through_registration_control(
        membership,
        atlas,
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        initial_transform=tuple(float(value) for value in carried_transform),
        previous_energy=prior_energy,
        iteration=CASE["iteration"],
        maximum_iterations=CASE["maximum_iterations"],
        convergence_threshold=CASE["convergence_threshold"],
        do_growing=False,
        align_shapes=True,
        subsample=2,
        max_registration_iterations=5,
        registration_energy_tolerance=0.0,
        minimum_lambda=1e-3,
    )
    assert handoff.handoff.membership_pass is membership
    assert handoff.handoff.registration is not None
    assert handoff.handoff.transform.shape == (6,)
    assert handoff.control.energy_ratio == toads_energy_ratio(prior_energy, membership.objective)
    assert handoff.control.stop is False
    assert handoff.control.run_growing is True
    assert handoff.stage_order == (
        "computeExactMemberships",
        "energyRatio",
        "atlas.registerShapes() -> registerAllShapes()",
        "convergenceCheck",
        "iterationLimitCheck",
    )

    nx, ny, nz = image.shape
    mem_path = tmp_path / "p18czb_memberships.bin"
    atlas_path = tmp_path / "p18czb_atlas.bin"
    transform_path = tmp_path / "p18czb_incoming_transform.bin"
    mem_path.write_bytes(np.asarray(membership.memberships, dtype="<f4").tobytes(order="C"))
    atlas_path.write_bytes(np.asarray(atlas, dtype="<f4").tobytes(order="C"))
    transform_path.write_bytes(np.asarray(carried_transform, dtype="<f4").tobytes(order="C"))
    java_registration = json.loads(subprocess.run(
        [
            "sh", str(REGISTRATION_HARNESS), str(mem_path), str(atlas_path), str(transform_path),
            str(nx), str(ny), str(nz),
        ],
        check=True, capture_output=True, text=True,
    ).stdout.splitlines()[-1])
    assert java_registration["source"] == (
        "VandyCRUISE MultipleToadDeformableAtlas.updateShapeRegistration(rigid) + "
        "registerShapes() + private registerAllShapes()"
    )
    assert np.asarray(java_registration["incoming_transform_bits"], dtype=np.int32).view(np.float32).tolist() == carried_transform.tolist()
    python_registration = handoff.handoff.registration
    assert python_registration is not None
    java_transform = np.asarray(java_registration["transform_bits"], dtype=np.int32).view(np.float32)
    java_energy = np.asarray([java_registration["energy_bits"]], dtype=np.int32).view(np.float32)
    np.testing.assert_allclose(python_registration.transform, java_transform, rtol=3e-5, atol=3e-6)
    np.testing.assert_allclose([python_registration.energy], java_energy, rtol=3e-5, atol=3e-6)

    subprocess.run(["javac", "-d", str(tmp_path), str(BRANCH_HARNESS)], check=True, capture_output=True, text=True)
    ratio = np.float32(handoff.control.energy_ratio)
    branch_reference = json.loads(subprocess.run(
        [
            "java", "-cp", str(tmp_path), "ToadsP18cytRegisteredIterationReferenceHarness",
            str(_bits([prior_energy])[0]), str(_bits([membership.objective])[0]),
            str(_bits([CASE["convergence_threshold"]])[0]), str(CASE["iteration"]),
            str(CASE["maximum_iterations"]), "false", "true",
        ],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert branch_reference["ratio_bits"] == _bits([ratio])[0]
    assert branch_reference["stop"] is handoff.control.stop is False
    assert branch_reference["run_growing"] is handoff.control.run_growing is True
    assert branch_reference["source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert branch_reference["order"] == [
        "computeExactMemberships", "energyRatio", "atlas.registerShapes",
        "convergenceCheck", "iterationLimitCheck", "propagateGrowing",
        "computeRelations", "computeBoundaryDistances", "computeClosestLabels",
    ]
