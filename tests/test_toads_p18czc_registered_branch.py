from __future__ import annotations

import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.toads_p18czc_registered_branch import (
    continue_p18cza_membership_through_registered_branch,
)
from pymacruise.toads_p18cza_carried_degree_iteration import (
    continue_p18cyz_thinning_through_carried_degree1_iteration,
)
from pymacruise.toads_p18cyx_next_loop_thinning import (
    continue_p18cyx_branch_through_next_thinning,
)
from test_toads_p18cyx_registered_cycle import _p18cyw_state
from test_toads_p18cyz_next_loop_thinning import _make_p18cyx
from test_toads_p18czb_registered_handoff import _atlas
from test_toads_post_thinning_exact_iteration import _bits

ROOT = Path(__file__).parent
CASE = json.loads(
    (ROOT / "data/synthetic/toads_p18czc_registered_branch/case.json").read_text()
)
RAS_FIXTURE = ROOT / "data" / "synthetic" / CASE["source_fixture"]
BRANCH_HARNESS = ROOT / "java_reference/ToadsP18cytRegisteredIterationReferenceHarness.java"
RELATION_HARNESS = ROOT / "java_reference/run_toads_p18czc_relation_refresh_reference.sh"


def _p18cza_state():
    previous, image = _p18cyw_state()
    cycle = _make_p18cyx(previous, image)
    prior_branch = cycle.registered_iteration.branch
    assert prior_branch.relations is not None and prior_branch.growing is not None
    shapes = np.zeros_like(prior_branch.handoff.membership_pass.memberships, dtype=np.float32)
    p18cyz = continue_p18cyx_branch_through_next_thinning(cycle, shapes)
    prior_exact = previous.exact_iteration
    exact = continue_p18cyz_thinning_through_carried_degree1_iteration(
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
    return previous, image, cycle, prior_branch, p18cyz, exact, shapes


def test_p18czb_continues_through_growth_and_java_referenced_relation_refresh(tmp_path: Path) -> None:
    previous, image, cycle, prior_branch, p18cyz, exact, shapes = _p18cza_state()
    source = nib.load(str(RAS_FIXTURE))
    assert nib.aff2axcodes(source.affine) == tuple(CASE["source_fixture_axis_codes"])
    affine = source.affine @ np.diag([1 / CASE["spatial_repeat_factor"]] * 3 + [1])
    assert nib.aff2axcodes(affine) == ("R", "A", "S")

    thinning = p18cyz.next_thinning.thinning
    boundary_state = prior_branch.relations.boundary_distances
    assert thinning.segmentation.shape == image.shape
    result = continue_p18cza_membership_through_registered_branch(
        exact,
        _atlas(image),
        previous_energy=np.float32(prior_branch.handoff.membership_pass.objective),
        initial_transform=tuple(float(v) for v in cycle.registration_control.handoff.transform),
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        iteration=CASE["iteration"],
        maximum_iterations=CASE["maximum_iterations"],
        convergence_threshold=CASE["convergence_threshold"],
        segmentation=thinning.segmentation,
        skeleton=thinning.skeleton,
        ordering=prior_branch.growing.ordering,
        class_volumes=thinning.class_volumes,
        boundary_distances=boundary_state,
        template_labels=np.asarray([0, 4, 7, 9, 11, 13, 17, 19, 21, 23, 25], dtype=np.int8),
        previous_classification=np.full(image.shape, -1, dtype=np.int8),
        relation_limit=CASE["relation_limit"],
        do_growing=False,
        initial_available=thinning.available,
        subsample=2,
        max_registration_iterations=5,
        registration_energy_tolerance=0.0,
        minimum_lambda=1e-3,
    )
    registration = result.registration_control
    branch = result.registered_iteration.branch
    assert registration.handoff.membership_pass is exact.exact_iteration.membership_pass
    assert registration.handoff.registration is not None
    assert branch.control == registration.control
    assert branch.control.stop is False and branch.control.run_growing is True
    assert branch.growing is not None and branch.relations is not None
    assert branch.classification is None
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

    subprocess.run(["javac", "-d", str(tmp_path), str(BRANCH_HARNESS)], check=True, capture_output=True, text=True)
    current_energy = np.float32(exact.exact_iteration.membership_pass.objective)
    previous_energy = np.float32(prior_branch.handoff.membership_pass.objective)
    control_reference = json.loads(subprocess.run(
        [
            "java", "-cp", str(tmp_path), "ToadsP18cytRegisteredIterationReferenceHarness",
            str(_bits([previous_energy])[0]), str(_bits([current_energy])[0]),
            str(_bits([CASE["convergence_threshold"]])[0]), str(CASE["iteration"]),
            str(CASE["maximum_iterations"]), "false", "true",
        ],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert control_reference["source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert control_reference["ratio_bits"] == _bits([registration.control.energy_ratio])[0]
    assert control_reference["stop"] is branch.control.stop is False
    assert control_reference["run_growing"] is branch.control.run_growing is True
    assert control_reference["order"] == [
        "computeExactMemberships", "energyRatio", "atlas.registerShapes",
        "convergenceCheck", "iterationLimitCheck", "propagateGrowing",
        "computeRelations", "computeBoundaryDistances", "computeClosestLabels",
    ]

    # Java recomputes the selected relation, boundary, and closest-label
    # methods on the exact post-growth labels, P18cza memberships, and carried
    # boundary-distance bytes produced/consumed by this branch.
    seg_path = tmp_path / "p18czc_segmentation.bin"
    mem_path = tmp_path / "p18czc_memberships.bin"
    boundary_path = tmp_path / "p18czc_boundary.bin"
    labels = np.asarray(branch.growing.segmentation, dtype=np.int8)
    memberships = np.asarray(registration.handoff.membership_pass.memberships, dtype=np.float32)
    seg_path.write_bytes(labels.tobytes(order="C"))
    mem_path.write_bytes(np.asarray(memberships, dtype="<f4").tobytes(order="C"))
    boundary_path.write_bytes(np.asarray(boundary_state, dtype=np.int8).tobytes(order="C"))
    nx, ny, nz = image.shape
    java_relation = json.loads(subprocess.run(
        [
            "sh", str(RELATION_HARNESS), str(seg_path), str(mem_path), str(boundary_path),
            str(nx), str(ny), str(nz), str(CASE["membership_class_count"]), str(CASE["relation_limit"]),
        ],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert java_relation["source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    relation = branch.relations
    np.testing.assert_array_equal(
        np.asarray(java_relation["relation_map"], dtype=np.int16).reshape(image.shape),
        relation.relation.relation_map,
    )
    np.testing.assert_array_equal(
        np.asarray(java_relation["relation_membership"], dtype=bool),
        relation.relation.relation_membership,
    )
    np.testing.assert_array_equal(
        np.asarray(java_relation["boundary_distances"], dtype=np.int8).reshape((*image.shape, 11)),
        relation.boundary_distances,
    )
    np.testing.assert_array_equal(
        np.asarray(java_relation["closest_labels"], dtype=np.int8).reshape((*image.shape, 4)),
        relation.closest_labels,
    )
