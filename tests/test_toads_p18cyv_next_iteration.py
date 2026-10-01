from __future__ import annotations

import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.toads_default_membership_registration_control import (
    continue_default_membership_through_registration_control,
)
from pymacruise.toads_default_membership_registered_iteration import (
    continue_default_membership_registered_iteration,
)
from pymacruise.toads_p18cyv_next_iteration import (
    continue_p18cyv_iteration_through_next_thinning,
)
from test_toads_p18cyt_registered_iteration import _p18cyt_membership
from test_toads_post_thinning_exact_iteration import _continuation_state

ROOT = Path(__file__).parent
CASE = json.loads(
    (ROOT / "data/synthetic/toads_p18cyv_next_iteration/case.json").read_text()
)
RAS_FIXTURE = ROOT / "data/synthetic/toads_post_registration_branch/membership_maps.nii.gz"
NEXT_THINNING_SCRIPT = ROOT.parent / "scripts/run_java_toads_post_registration_next_thinning_reference.sh"
ORDER_HARNESS = ROOT / "java_reference/ToadsP18cyvNextIterationOrderReferenceHarness.java"


def _input_bytes(segmentation, memberships, shapes, adjacency):
    return b"".join((
        np.asarray(segmentation, dtype=np.int8).tobytes(order="C"),
        np.asarray(memberships, dtype=">f4").tobytes(order="C"),
        np.asarray(shapes, dtype=">f4").tobytes(order="C"),
        np.asarray(adjacency, dtype=np.uint8).tobytes(order="C"),
    ))


def _label_string(values, empty=False):
    return "".join(
        "." if empty and int(value) == -1 else
        (str(int(value)) if int(value) < 10 else chr(ord("a") + int(value) - 10))
        for value in values.ravel()
    )


def _float_bits(values):
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32).astype(int).tolist()


def test_actual_p18cyt_branch_flows_to_java_referenced_next_thinning(tmp_path: Path) -> None:
    next_thinning, image, centroids, field, closest = _continuation_state()
    p18cyt = _p18cyt_membership(next_thinning, image, centroids, field, closest)
    membership_pass = p18cyt.membership_pass
    assert image.ndim == 3 and membership_pass.memberships.shape[-1] == 11

    source = nib.load(str(RAS_FIXTURE))
    assert nib.aff2axcodes(source.affine) == tuple(CASE["source_fixture_axis_codes"])
    affine = source.affine @ np.diag([1 / CASE["spatial_repeat_factor"]] * 3 + [1])
    assert nib.aff2axcodes(affine) == ("R", "A", "S")
    assert nib.Nifti1Image(image, affine).shape == image.shape

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

    objective = np.float32(membership_pass.objective)
    previous_energy = np.float32(objective * np.float32(1.01))
    settings = dict(
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        initial_transform=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        previous_energy=previous_energy,
        iteration=3,
        maximum_iterations=5,
        convergence_threshold=1e-4,
        do_growing=False,
        align_shapes=False,
    )
    registration_control = continue_default_membership_through_registration_control(
        membership_pass, atlas, **settings,
    )

    source_branch = next_thinning.branch
    thinning = next_thinning.thinning
    old_relations = source_branch.relations
    assert old_relations is not None and source_branch.growing is not None
    boundary = np.repeat(
        np.repeat(
            np.repeat(
                old_relations.boundary_distances,
                CASE["spatial_repeat_factor"],
                axis=0,
            ),
            CASE["spatial_repeat_factor"],
            axis=1,
        ),
        CASE["spatial_repeat_factor"],
        axis=2,
    )
    registered_iteration = continue_default_membership_registered_iteration(
        registration_control,
        previous_energy=previous_energy,
        iteration=3,
        maximum_iterations=5,
        convergence_threshold=1e-4,
        segmentation=thinning.segmentation,
        skeleton=thinning.skeleton,
        ordering=source_branch.growing.ordering,
        class_volumes=thinning.class_volumes,
        boundary_distances=boundary,
        template_labels=np.asarray([0, 4, 7, 9, 11, 13, 17, 19, 21, 23, 25], dtype=np.int8),
        previous_classification=np.full(image.shape, -1, dtype=np.int8),
        relation_limit=CASE["relation_limit"],
        do_growing=False,
        initial_available=thinning.available,
    )
    assert registered_iteration.branch.control.stop is False
    assert registered_iteration.branch.growing is not None
    assert registered_iteration.branch.relations is not None

    shape_values = np.zeros_like(membership_pass.memberships, dtype=np.float32)
    result = continue_p18cyv_iteration_through_next_thinning(
        registered_iteration,
        shape_values,
    )
    next_pass = result.next_thinning.thinning
    assert result.next_thinning.branch is registered_iteration.branch
    assert result.stage_order[-5:] == (
        "propagateGrowing",
        "computeRelations",
        "computeBoundaryDistances",
        "computeClosestLabels",
        "next_iteration.propagateThinning",
    )
    assert next_pass.segmentation.shape == image.shape
    assert nib.aff2axcodes(nib.Nifti1Image(next_pass.segmentation, affine).affine) == ("R", "A", "S")

    # Java formula reference for the wrapper's next-iteration dispatch.
    compile_order = subprocess.run(
        ["javac", "-d", str(tmp_path), str(ORDER_HARNESS)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert not compile_order.stderr
    order_ref = json.loads(subprocess.run(
        ["java", "-cp", str(tmp_path), "ToadsP18cyvNextIterationOrderReferenceHarness"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout)
    assert order_ref["algorithm_source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert order_ref["segmentation_source_sha256"] == "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert order_ref["order"][-3:] == [
        "nextIteration.propagateThinning", "initThinningLabels", "propagateThinning",
    ]
    assert order_ref["useApprox"] is False
    assert order_ref["volumeRatio"] == CASE["thinning_volume_ratio"]
    assert order_ref["limit"] == CASE["thinning_limit"]

    # Run the existing selected Java initThinningLabels + propagateThinning
    # reference on precisely the Python branch's labels, memberships, shape,
    # and refreshed relation adjacency.
    branch = registered_iteration.branch
    assert branch.growing is not None and branch.relations is not None
    input_path = tmp_path / "p18cyv_next_thinning.bin"
    input_path.write_bytes(_input_bytes(
        branch.growing.segmentation,
        branch.handoff.membership_pass.memberships,
        shape_values,
        branch.relations.relation.class_adjacency,
    ))
    completed = subprocess.run(
        [
            "sh", str(NEXT_THINNING_SCRIPT), str(input_path),
            *(str(v) for v in image.shape), str(membership_pass.memberships.shape[-1]),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    java = json.loads(completed.stdout)
    assert java["source"] == (
        "selected VandyCRUISE MultipleToadSegmentation.initThinningLabels + "
        "propagateThinning; short constructor useApprox=false"
    )
    assert java["source_sha256"] == order_ref["segmentation_source_sha256"]
    assert java["shape"] == [*image.shape, 11]
    assert java["volume_ratio_bits"] == int(np.asarray([np.float32(0.5)]).view(np.int32)[0])
    assert java["limit_bits"] == int(np.asarray([np.float32(3.0)]).view(np.int32)[0])
    initial_volumes = np.zeros(11, dtype=np.float32)
    labels = branch.growing.segmentation
    memberships = branch.handoff.membership_pass.memberships
    for x in range(1, labels.shape[0] - 1):
        for y in range(1, labels.shape[1] - 1):
            for z in range(1, labels.shape[2] - 1):
                k = int(labels[x, y, z])
                initial_volumes[k] = np.float32(initial_volumes[k] + memberships[x, y, z, k])
    assert java["initial_volumes_bits"] == _float_bits(initial_volumes)
    assert java["minimum_volumes_bits"] == _float_bits(next_pass.minimum_class_volumes)
    assert java["boundary_counts"] == next_pass.boundary_counts.astype(int).tolist()
    assert java["processed"] == [
        [int(np.asarray([item.value], dtype=np.float32).view(np.int32)[0]), item.x, item.y, item.z, item.label]
        for item in next_pass.processed_items
    ]
    assert java["segmentation"] == _label_string(next_pass.segmentation)
    assert java["skeleton"] == _label_string(next_pass.skeleton, empty=True)
    assert java["available"] == "".join("1" if value else "0" for value in next_pass.available.ravel())
    assert java["final_volumes_bits"] == _float_bits(next_pass.class_volumes)
