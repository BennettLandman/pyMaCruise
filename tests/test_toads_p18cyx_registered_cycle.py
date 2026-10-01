from __future__ import annotations

import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.toads_p18cyx_registered_cycle import (
    continue_p18cyw_membership_through_registration_branch,
)
from test_toads_p18cyw_exact_iteration import _p18cyv_state
from pymacruise.toads_p18cyw_exact_iteration import (
    continue_p18cyv_thinning_through_exact_iteration_degree1,
)

ROOT = Path(__file__).parent
CASE = json.loads(
    (ROOT / "data/synthetic/toads_p18cyx_registered_cycle/case.json").read_text()
)
RAS_FIXTURE = ROOT / "data/synthetic/toads_post_registration_branch/membership_maps.nii.gz"
REGISTRATION_HARNESS = ROOT / "java_reference/generate_toads_p18cyp_actual_registration_reference.sh"
BRANCH_HARNESS = ROOT / "java_reference/ToadsP18cytRegisteredIterationReferenceHarness.java"


def _bits(values):
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).astype(int).tolist()


def _float32_from_bits(bits):
    return np.asarray([bits], dtype=np.int32).view(np.float32)[0]


def _p18cyw_state():
    p18cyt, p18cyv, image = _p18cyv_state()
    memberships = p18cyv.next_thinning.branch.handoff.membership_pass.memberships
    relations = p18cyv.registered_iteration.branch.relations
    assert relations is not None
    field = np.asarray(p18cyt.field_update.fit.field, dtype=np.float32).copy()
    initial_centroids = np.asarray(p18cyt.centroids[0], dtype=np.float32).copy()
    exponent_map = np.empty((1, 11), dtype=np.float32)
    exponent_map[0, 0] = 1.0
    exponent_map[0, 1:3] = 0.8
    exponent_map[0, 3:8] = 1.0
    exponent_map[0, 8:11] = 1.0
    exact = continue_p18cyv_thinning_through_exact_iteration_degree1(
        p18cyv,
        image,
        initial_centroids,
        field,
        category_exponents=np.asarray([0.8, 1.0, 1.0], dtype=np.float32),
        exponent_by_channel_class=exponent_map,
        modality_weights=np.asarray([0.73], dtype=np.float32),
        class_variance_weights=np.asarray(
            [[0.81, 0.93, 1.17, 0.76, 0.88, 1.02, 1.14, 0.69, 0.96, 1.21, 0.73]],
            dtype=np.float32,
        ),
        transformed_shapes=np.zeros_like(memberships, dtype=np.float32),
        closest_labels=relations.closest_labels,
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
    return exact, image


def test_p18cyw_membership_runs_second_loop_registration_and_continuing_branch(tmp_path: Path) -> None:
    exact, image = _p18cyw_state()
    source = nib.load(str(RAS_FIXTURE))
    assert nib.aff2axcodes(source.affine) == tuple(CASE["source_fixture_axis_codes"])
    affine = source.affine @ np.diag([1 / CASE["spatial_repeat_factor"]] * 3 + [1])
    assert nib.aff2axcodes(affine) == ("R", "A", "S")
    assert nib.Nifti1Image(image, affine).shape == image.shape

    iteration_result = exact.exact_iteration
    membership_pass = iteration_result.membership_pass
    labels = exact.p18cyv.next_thinning.thinning.segmentation
    assert membership_pass.memberships.shape == (*image.shape, CASE["membership_classes"])
    assert labels.shape == image.shape

    nx, ny, nz = image.shape
    gx, gy, gz = np.meshgrid(
        np.linspace(-1, 1, nx, dtype=np.float32),
        np.linspace(-1, 1, ny, dtype=np.float32),
        np.linspace(-1, 1, nz, dtype=np.float32),
        indexing="ij",
    )
    atlas = np.empty((*image.shape, 11), dtype=np.float32)
    for k in range(11):
        radius = (
            (gx - np.float32(-0.45 + 0.09 * k)) ** 2
            + 1.3 * (gy - np.float32(0.3 * np.sin(k))) ** 2
            + 0.8 * (gz - np.float32(0.25 * np.cos(k))) ** 2
        )
        atlas[..., k] = np.asarray(
            np.exp(-radius / np.float32(0.42)) * 0.75 + 0.05,
            dtype=np.float32,
        )

    prior_branch = exact.p18cyv.registered_iteration.branch
    prior_relations = prior_branch.relations
    assert prior_branch.growing is not None and prior_relations is not None
    previous_segmentation = exact.p18cyv.next_thinning.thinning
    previous_energy = np.float32(prior_branch.handoff.membership_pass.objective)
    parameters = dict(
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        iteration=CASE["iteration"],
        maximum_iterations=CASE["maximum_iterations"],
        convergence_threshold=CASE["convergence_threshold"],
        do_growing=False,
        segmentation=previous_segmentation.segmentation,
        skeleton=previous_segmentation.skeleton,
        ordering=prior_branch.growing.ordering,
        class_volumes=previous_segmentation.class_volumes,
        boundary_distances=prior_relations.boundary_distances,
        template_labels=np.asarray([0, 4, 7, 9, 11, 13, 17, 19, 21, 23, 25], dtype=np.int8),
        previous_classification=np.full(image.shape, -1, dtype=np.int8),
        relation_limit=1,
        initial_available=previous_segmentation.available,
        subsample=2,
        max_registration_iterations=5,
        registration_energy_tolerance=0.0,
        minimum_lambda=1e-3,
    )
    result = continue_p18cyw_membership_through_registration_branch(
        exact,
        atlas,
        **parameters,
    )
    assert result.registration_control.handoff.membership_pass is membership_pass
    assert result.registered_iteration.registration_control is result.registration_control
    assert result.registered_iteration.branch.control == result.registration_control.control
    assert result.registered_iteration.branch.control.stop is False
    assert result.registered_iteration.branch.growing is not None
    assert result.registered_iteration.branch.relations is not None
    assert result.registered_iteration.branch.classification is None
    assert result.stage_order[-8:] == (
        "energyRatio",
        "atlas.registerShapes() -> registerAllShapes()",
        "convergenceCheck",
        "iterationLimitCheck",
        "propagateGrowing",
        "computeRelations",
        "computeBoundaryDistances",
        "computeClosestLabels",
    )

    # Execute selected VandyCRUISE rigid registration on the exact P18cyw
    # membership/atlas arrays; this method harness uses the same five-iteration
    # cap, scale-one rigid route, offsets, and fixture's identity transform.
    membership_path = tmp_path / "p18cyw_memberships.bin"
    atlas_path = tmp_path / "p18cyw_atlas.bin"
    membership_path.write_bytes(np.asarray(membership_pass.memberships, dtype="<f4").tobytes(order="C"))
    atlas_path.write_bytes(np.asarray(atlas, dtype="<f4").tobytes(order="C"))
    java_registration = json.loads(subprocess.run(
        ["sh", str(REGISTRATION_HARNESS), str(membership_path), str(atlas_path), str(nx), str(ny), str(nz)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()[-1])
    assert java_registration["source"] == (
        "VandyCRUISE MultipleToadDeformableAtlas.updateShapeRegistration(rigid) + "
        "registerShapes() + private registerAllShapes()"
    )
    python_registration = result.registration_control.handoff.registration
    assert python_registration is not None
    java_transform = np.asarray(java_registration["transform_bits"], dtype=np.int32).view(np.float32)
    java_energy = np.asarray([java_registration["energy_bits"]], dtype=np.int32).view(np.float32)
    np.testing.assert_allclose(java_transform, python_registration.transform, rtol=3e-5, atol=3e-6)
    np.testing.assert_allclose(java_energy, [python_registration.energy], rtol=3e-5, atol=3e-6)

    # The selected Java wrapper excerpt independently checks the current
    # energy ratio and source branch order for this exact P18cyw objective.
    compile_result = subprocess.run(
        ["javac", "-d", str(tmp_path), str(BRANCH_HARNESS)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert not compile_result.stderr
    current_energy = np.float32(membership_pass.objective)
    threshold = np.float32(CASE["convergence_threshold"])
    branch_reference = json.loads(subprocess.run(
        [
            "java", "-cp", str(tmp_path), "ToadsP18cytRegisteredIterationReferenceHarness",
            str(int(np.asarray([previous_energy], dtype=np.float32).view(np.int32)[0])),
            str(int(np.asarray([current_energy], dtype=np.float32).view(np.int32)[0])),
            str(int(np.asarray([threshold], dtype=np.float32).view(np.int32)[0])),
            str(CASE["iteration"]),
            str(CASE["maximum_iterations"]),
            "false",
            "true",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout)
    assert branch_reference["ratio_bits"] == int(
        np.asarray([result.registration_control.control.energy_ratio], dtype=np.float32).view(np.int32)[0]
    )
    assert branch_reference["stop"] is result.registered_iteration.branch.control.stop is False
    assert branch_reference["run_growing"] is result.registered_iteration.branch.control.run_growing is True
    assert branch_reference["source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert branch_reference["order"] == [
        "computeExactMemberships", "energyRatio", "atlas.registerShapes",
        "convergenceCheck", "iterationLimitCheck", "propagateGrowing",
        "computeRelations", "computeBoundaryDistances", "computeClosestLabels",
    ]
