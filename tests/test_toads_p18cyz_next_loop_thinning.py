from __future__ import annotations

import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.toads_p18cyx_registered_cycle import continue_p18cyw_membership_through_registration_branch
from pymacruise.toads_p18cyx_next_loop_thinning import continue_p18cyx_branch_through_next_thinning
from test_toads_p18cyx_registered_cycle import _p18cyw_state

ROOT = Path(__file__).parent
CASE = json.loads(
    (ROOT / "data/synthetic/toads_p18cyz_next_loop_thinning/case.json").read_text()
)
RAS_FIXTURE = ROOT / "data/synthetic/toads_post_registration_branch/membership_maps.nii.gz"
THINNING_SCRIPT = ROOT.parent / "scripts/run_java_toads_post_registration_next_thinning_reference.sh"
ORDER_HARNESS = ROOT / "java_reference/ToadsP18cyvNextIterationOrderReferenceHarness.java"


def _serialize(segmentation, memberships, shapes, adjacency):
    return b"".join((
        np.asarray(segmentation, dtype=np.int8).tobytes(order="C"),
        np.asarray(memberships, dtype=">f4").tobytes(order="C"),
        np.asarray(shapes, dtype=">f4").tobytes(order="C"),
        np.asarray(adjacency, dtype=np.uint8).tobytes(order="C"),
    ))


def _labels(values, *, empty=False):
    return "".join(
        "." if empty and int(value) == -1 else
        (str(int(value)) if int(value) < 10 else chr(ord("a") + int(value) - 10))
        for value in values.ravel()
    )


def _float_bits(values):
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).astype(int).tolist()


def _make_p18cyx(exact, image):
    iteration_result = exact.exact_iteration
    membership_pass = iteration_result.membership_pass
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
            np.exp(-radius / np.float32(0.42)) * 0.75 + 0.05, dtype=np.float32
        )
    prior_branch = exact.p18cyv.registered_iteration.branch
    prior_relations = prior_branch.relations
    assert prior_branch.growing is not None and prior_relations is not None
    thinning = exact.p18cyv.next_thinning.thinning
    return continue_p18cyw_membership_through_registration_branch(
        exact,
        atlas,
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        iteration=4,
        maximum_iterations=8,
        convergence_threshold=0.0,
        do_growing=False,
        segmentation=thinning.segmentation,
        skeleton=thinning.skeleton,
        ordering=prior_branch.growing.ordering,
        class_volumes=thinning.class_volumes,
        boundary_distances=prior_relations.boundary_distances,
        template_labels=np.asarray([0, 4, 7, 9, 11, 13, 17, 19, 21, 23, 25], dtype=np.int8),
        previous_classification=np.full(image.shape, -1, dtype=np.int8),
        relation_limit=CASE["relation_limit"],
        initial_available=thinning.available,
        subsample=2,
        max_registration_iterations=5,
        registration_energy_tolerance=0.0,
        minimum_lambda=1e-3,
    )


def test_p18cyx_continuing_branch_feeds_java_referenced_next_thinning(tmp_path: Path) -> None:
    exact, image = _p18cyw_state()
    cycle = _make_p18cyx(exact, image)
    branch = cycle.registered_iteration.branch
    assert branch.control.stop is False and branch.growing is not None and branch.relations is not None

    source = nib.load(str(RAS_FIXTURE))
    assert nib.aff2axcodes(source.affine) == tuple(CASE["source_fixture_axis_codes"])
    affine = source.affine @ np.diag([1 / CASE["spatial_repeat_factor"]] * 3 + [1])
    assert nib.aff2axcodes(affine) == ("R", "A", "S")
    assert nib.Nifti1Image(branch.handoff.membership_pass.memberships, affine).shape == (
        *image.shape, CASE["membership_classes"],
    )

    shapes = np.zeros_like(branch.handoff.membership_pass.memberships, dtype=np.float32)
    result = continue_p18cyx_branch_through_next_thinning(cycle, shapes)
    assert result.registered_cycle is cycle
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
    thinning = result.next_thinning.thinning
    assert thinning.segmentation.shape == image.shape
    assert nib.aff2axcodes(nib.Nifti1Image(thinning.segmentation, affine).affine) == ("R", "A", "S")

    # The selected Java order harness identifies the source wrapper edge;
    # the numeric harness below executes the selected thinning reference.
    order_classes = subprocess.run(
        ["javac", "-d", str(tmp_path), str(ORDER_HARNESS)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert not order_classes.stderr
    order_reference = json.loads(subprocess.run(
        ["java", "-cp", str(tmp_path), "ToadsP18cyvNextIterationOrderReferenceHarness"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout)
    assert order_reference["order"][-3:] == [
        "nextIteration.propagateThinning", "initThinningLabels", "propagateThinning",
    ]
    assert order_reference["useApprox"] is False
    assert order_reference["volumeRatio"] == CASE["thinning_volume_ratio"]
    assert order_reference["limit"] == CASE["thinning_limit"]

    input_path = tmp_path / "p18cyz_thinning_input.bin"
    input_path.write_bytes(_serialize(
        branch.growing.segmentation,
        branch.handoff.membership_pass.memberships,
        shapes,
        branch.relations.relation.class_adjacency,
    ))
    java = json.loads(subprocess.run(
        [
            "sh", str(THINNING_SCRIPT), str(input_path),
            *(str(v) for v in image.shape), str(CASE["membership_classes"]),
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout)
    assert java["source"] == (
        "selected VandyCRUISE MultipleToadSegmentation.initThinningLabels + "
        "propagateThinning; short constructor useApprox=false"
    )
    assert java["source_sha256"] == "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert java["shape"] == [*image.shape, CASE["membership_classes"]]
    assert java["volume_ratio_bits"] == int(np.asarray([np.float32(0.5)]).view(np.int32)[0])
    assert java["limit_bits"] == int(np.asarray([np.float32(3.0)]).view(np.int32)[0])
    initial_volumes = np.zeros(CASE["membership_classes"], dtype=np.float32)
    labels = branch.growing.segmentation
    memberships = branch.handoff.membership_pass.memberships
    for x in range(1, labels.shape[0] - 1):
        for y in range(1, labels.shape[1] - 1):
            for z in range(1, labels.shape[2] - 1):
                k = int(labels[x, y, z])
                initial_volumes[k] = np.float32(initial_volumes[k] + memberships[x, y, z, k])
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
