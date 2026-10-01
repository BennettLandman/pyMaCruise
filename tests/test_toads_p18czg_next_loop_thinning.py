from __future__ import annotations

import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.toads_p18czf_registered_branch import (
    continue_p18cze_membership_through_registered_branch,
)
from pymacruise.toads_p18cze_next_exact_iteration import (
    continue_p18czd_thinning_through_next_exact_iteration_degree1,
)
from pymacruise.toads_p18czg_next_loop_thinning import (
    continue_p18czf_branch_through_next_loop_thinning,
)
from test_toads_p18czb_registered_handoff import _atlas
from test_toads_p18czd_next_loop_thinning import _float_bits, _labels, _serialize
from test_toads_p18cze_next_exact_iteration import _parameters, _p18czd_state

ROOT = Path(__file__).parent
CASE = json.loads(
    (ROOT / "data/synthetic/toads_p18czg_next_loop_thinning/case.json").read_text()
)
RAS_FIXTURE = ROOT / "data" / "synthetic" / CASE["source_fixture"]
THINNING_SCRIPT = ROOT.parent / "scripts/run_java_toads_post_registration_next_thinning_reference.sh"
ORDER_HARNESS = ROOT / "java_reference/ToadsP18cyvNextIterationOrderReferenceHarness.java"


def test_p18czf_continuing_branch_feeds_java_referenced_next_thinning(tmp_path: Path) -> None:
    image, p18czc, p18czd = _p18czd_state()
    prior_branch = p18czc.registered_iteration.branch
    assert prior_branch.growing is not None and prior_branch.relations is not None
    memberships = p18czd.next_thinning.branch.handoff.membership_pass.memberships
    exact = continue_p18czd_thinning_through_next_exact_iteration_degree1(
        p18czd,
        image,
        polynomial_degree=3,
        low_energy_ratio=0.01,
        **_parameters(memberships.shape),
    )
    p18czf = continue_p18cze_membership_through_registered_branch(
        exact,
        _atlas(image),
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        iteration=CASE["iteration"],
        maximum_iterations=CASE["maximum_iterations"],
        convergence_threshold=CASE["convergence_threshold"],
        do_growing=False,
        template_labels=np.asarray([0, 4, 7, 9, 11, 13, 17, 19, 21, 23, 25], dtype=np.int8),
        previous_classification=np.full(image.shape, -1, dtype=np.int8),
        relation_limit=1,
        subsample=2,
        max_registration_iterations=5,
        registration_energy_tolerance=0.0,
        minimum_lambda=1e-3,
    )
    branch = p18czf.registered_iteration.branch
    assert branch.control.stop is False and branch.growing is not None and branch.relations is not None
    assert branch.classification is None

    source = nib.load(str(RAS_FIXTURE))
    assert nib.aff2axcodes(source.affine) == tuple(CASE["source_fixture_axis_codes"])
    affine = source.affine @ np.diag([1 / CASE["spatial_repeat_factor"]] * 3 + [1])
    assert nib.aff2axcodes(affine) == ("R", "A", "S")
    shapes = np.zeros_like(branch.handoff.membership_pass.memberships, dtype=np.float32)
    result = continue_p18czf_branch_through_next_loop_thinning(p18czf, shapes)
    thinning = result.next_thinning.thinning
    assert result.next_thinning.branch is branch
    assert result.stage_order[-9:] == (
        "energyRatio",
        "atlas.registerShapes() -> registerAllShapes()",
        "convergenceCheck",
        "iterationLimitCheck",
        "propagateGrowing",
        "computeRelations",
        "computeBoundaryDistances",
        "computeClosestLabels",
        "next_iteration.propagateThinning",
    )
    assert nib.aff2axcodes(nib.Nifti1Image(thinning.segmentation, affine).affine) == ("R", "A", "S")

    subprocess.run(
        ["javac", "-d", str(tmp_path), str(ORDER_HARNESS)],
        check=True, capture_output=True, text=True,
    )
    order = json.loads(subprocess.run(
        ["java", "-cp", str(tmp_path), "ToadsP18cyvNextIterationOrderReferenceHarness"],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert order["order"][-3:] == [
        "nextIteration.propagateThinning", "initThinningLabels", "propagateThinning",
    ]
    assert order["useApprox"] is False
    assert order["volumeRatio"] == CASE["volume_ratio"]
    assert order["limit"] == CASE["thinning_limit"]

    memberships = branch.handoff.membership_pass.memberships
    input_path = tmp_path / "p18czg_thinning_input.bin"
    input_path.write_bytes(_serialize(
        branch.growing.segmentation,
        memberships,
        shapes,
        branch.relations.relation.class_adjacency,
    ))
    java = json.loads(subprocess.run(
        [
            "sh", str(THINNING_SCRIPT), str(input_path),
            *(str(v) for v in image.shape), str(CASE["membership_class_count"]),
        ],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert java["source"] == (
        "selected VandyCRUISE MultipleToadSegmentation.initThinningLabels + "
        "propagateThinning; short constructor useApprox=false"
    )
    assert java["source_sha256"] == "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert java["shape"] == [*image.shape, CASE["membership_class_count"]]
    assert java["volume_ratio_bits"] == int(np.asarray([np.float32(CASE["volume_ratio"])]).view(np.int32)[0])
    assert java["limit_bits"] == int(np.asarray([np.float32(CASE["thinning_limit"])]).view(np.int32)[0])

    initial_volumes = np.zeros(CASE["membership_class_count"], dtype=np.float32)
    labels = branch.growing.segmentation
    for x in range(1, labels.shape[0] - 1):
        for y in range(1, labels.shape[1] - 1):
            for z in range(1, labels.shape[2] - 1):
                klass = int(labels[x, y, z])
                initial_volumes[klass] = np.float32(initial_volumes[klass] + memberships[x, y, z, klass])
    assert java["initial_volumes_bits"] == _float_bits(initial_volumes)
    assert java["minimum_volumes_bits"] == _float_bits(thinning.minimum_class_volumes)
    assert java["boundary_counts"] == thinning.boundary_counts.astype(int).tolist()
    assert java["processed"] == [
        [int(np.asarray([item.value], dtype=np.float32).view(np.int32)[0]), item.x, item.y, item.z, item.label]
        for item in thinning.processed_items
    ]
    assert java["segmentation"] == _labels(thinning.segmentation)
    assert java["skeleton"] == _labels(thinning.skeleton, empty=True)
    assert java["available"] == "".join("1" if value else "0" for value in thinning.available.ravel())
    assert java["final_volumes_bits"] == _float_bits(thinning.class_volumes)
