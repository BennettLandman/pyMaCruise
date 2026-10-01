from __future__ import annotations

import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.toads_p18czj_next_loop_thinning import (
    continue_p18czi_branch_through_next_loop_thinning,
)
from pymacruise.toads_p18czk_next_exact_iteration import (
    continue_p18czj_thinning_through_next_exact_iteration_degree1,
)
from pymacruise.toads_p18czl_terminal_branch import (
    continue_p18czk_membership_through_terminal_branch,
)
from test_toads_p18czb_registered_handoff import _atlas
from test_toads_p18cze_next_exact_iteration import _parameters
from test_toads_p18czj_next_loop_thinning import _p18czi_state
from test_toads_post_thinning_exact_iteration import _bits

ROOT = Path(__file__).parent
CASE = json.loads(
    (ROOT / "data/synthetic/toads_p18czl_terminal_branch/case.json").read_text()
)
RAS_FIXTURE = ROOT / "data" / "synthetic" / CASE["source_fixture"]
REGISTRATION_HARNESS = ROOT / "java_reference/run_toads_p18cza_carried_transform_registration_reference.sh"
BRANCH_HARNESS = ROOT / "java_reference/ToadsP18cytRegisteredIterationReferenceHarness.java"


def test_p18czk_membership_hits_java_referenced_iteration_limit(tmp_path: Path) -> None:
    image, p18czi = _p18czi_state()
    p18czj = continue_p18czi_branch_through_next_loop_thinning(
        p18czi,
        np.zeros_like(
            p18czi.registered_iteration.branch.handoff.membership_pass.memberships,
            dtype=np.float32,
        ),
    )
    memberships = p18czj.next_thinning.branch.handoff.membership_pass.memberships
    p18czk = continue_p18czj_thinning_through_next_exact_iteration_degree1(
        p18czj,
        image,
        polynomial_degree=3,
        low_energy_ratio=0.01,
        **_parameters(memberships.shape),
    )
    result = continue_p18czk_membership_through_terminal_branch(
        p18czk,
        _atlas(image),
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        iteration=CASE["iteration"],
        maximum_iterations=CASE["maximum_iterations"],
        convergence_threshold=CASE["convergence_threshold"],
        template_labels=np.asarray([0, 4, 7, 9, 11, 13, 17, 19, 21, 23, 25], dtype=np.int8),
        previous_classification=np.full(image.shape, -1, dtype=np.int8),
        relation_limit=1,
        subsample=2,
        max_registration_iterations=5,
        registration_energy_tolerance=0.0,
        minimum_lambda=1e-3,
    )

    source = nib.load(str(RAS_FIXTURE))
    assert nib.aff2axcodes(source.affine) == tuple(CASE["source_fixture_axis_codes"])
    affine = source.affine @ np.diag([1 / CASE["spatial_repeat_factor"]] * 3 + [1])
    assert nib.aff2axcodes(affine) == ("R", "A", "S")
    current = p18czk.exact_iteration.exact_iteration.membership_pass
    assert nib.Nifti1Image(current.memberships, affine).shape == (*image.shape, CASE["membership_class_count"])

    previous_branch = p18czi.registered_iteration.branch
    previous_energy = np.float32(previous_branch.handoff.membership_pass.objective)
    assert result.previous_energy == previous_energy
    assert result.registration_control.handoff.membership_pass is current
    registration = result.registration_control.handoff.registration
    assert registration is not None
    branch = result.registered_iteration.branch
    assert branch.control.stop is True
    assert branch.control.run_growing is False
    assert branch.growing is None and branch.relations is None
    assert branch.classification is not None
    assert result.stage_order[-6:] == (
        "computeExactMemberships",
        "energyRatio",
        "atlas.registerShapes() -> registerAllShapes()",
        "convergenceCheck",
        "iterationLimitCheck",
        "updateClassificationFromLabels",
    )

    nx, ny, nz = image.shape
    membership_path = tmp_path / "p18czl_memberships.bin"
    atlas_path = tmp_path / "p18czl_atlas.bin"
    transform_path = tmp_path / "p18czl_incoming_transform.bin"
    membership_path.write_bytes(np.asarray(current.memberships, dtype="<f4").tobytes(order="C"))
    atlas = _atlas(image)
    atlas_path.write_bytes(np.asarray(atlas, dtype="<f4").tobytes(order="C"))
    transform_path.write_bytes(np.asarray(result.initial_transform, dtype="<f4").tobytes(order="C"))
    java_registration = json.loads(subprocess.run(
        [
            "sh", str(REGISTRATION_HARNESS), str(membership_path), str(atlas_path),
            str(transform_path), str(nx), str(ny), str(nz),
        ],
        check=True, capture_output=True, text=True,
    ).stdout.splitlines()[-1])
    assert java_registration["source"] == (
        "VandyCRUISE MultipleToadDeformableAtlas.updateShapeRegistration(rigid) + "
        "registerShapes() + private registerAllShapes()"
    )
    python_transform = np.asarray(registration.transform, dtype=np.float32)
    java_transform = np.asarray(java_registration["transform_bits"], dtype=np.int32).view(np.float32)
    java_energy = np.asarray([java_registration["energy_bits"]], dtype=np.int32).view(np.float32)
    np.testing.assert_allclose(python_transform, java_transform, rtol=3e-5, atol=3e-6)
    np.testing.assert_allclose([registration.energy], java_energy, rtol=3e-5, atol=3e-6)

    subprocess.run(["javac", "-d", str(tmp_path), str(BRANCH_HARNESS)], check=True, capture_output=True, text=True)
    java_control = json.loads(subprocess.run(
        [
            "java", "-cp", str(tmp_path), "ToadsP18cytRegisteredIterationReferenceHarness",
            str(_bits([previous_energy])[0]), str(_bits([current.objective])[0]),
            str(_bits([CASE["convergence_threshold"]])[0]), str(CASE["iteration"]),
            str(CASE["maximum_iterations"]), "false", "true",
        ],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert java_control["source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert java_control["ratio_bits"] == _bits([result.registration_control.control.energy_ratio])[0]
    assert java_control["stop"] is branch.control.stop is True
    assert java_control["run_growing"] is branch.control.run_growing is False
    assert java_control["order"] == [
        "computeExactMemberships", "energyRatio", "atlas.registerShapes",
        "convergenceCheck", "iterationLimitCheck", "updateClassificationFromLabels",
    ]
