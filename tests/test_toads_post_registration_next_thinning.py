from __future__ import annotations

import json
import struct
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_default_atlas_membership import ToadsDefaultAtlasMembershipPass
from pymacruise.toads_default_registration_handoff import ToadsDefaultRegistrationHandoff
from pymacruise.toads_post_registration import toads_post_registration_branch
from pymacruise.toads_post_registration_next_thinning import toads_post_registration_next_thinning
from pymacruise.toads_registration_sampling import ToadsRigidRegistrationResult


ROOT = Path(__file__).parent / "data/synthetic/toads_post_registration_branch"
JAVA_REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
STATE = json.loads((ROOT / "state.json").read_text())
HARNESS_SCRIPT = Path(__file__).parent.parent / "scripts/run_java_toads_post_registration_next_thinning_reference.sh"


def _load(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / name))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    return np.asarray(image.dataobj)


def _float32(bits: int) -> np.float32:
    return np.float32(struct.unpack("!f", struct.pack("!I", bits))[0])


def _registered_branch(
    previous_energy: np.float32 | None = None,
    convergence_threshold: np.float32 | None = None,
):
    membership_maps = _load("membership_maps.nii.gz").astype(np.float32)
    segmentation = _load("segmentation.nii.gz").astype(np.int16)
    skeleton = _load("skeleton.nii.gz").astype(np.int16)
    ordering = _load("ordering.nii.gz").astype(np.float32)
    class_volumes = _load("class_volumes.nii.gz").reshape(11).astype(np.float32)
    boundary_distances = _load("boundary_distances.nii.gz").astype(np.int8)
    previous_classification = _load("previous_classification.nii.gz").astype(np.int8)
    transform = np.asarray(STATE["registration_transform_bits"], dtype=np.uint32).view(np.float32)
    membership_pass = ToadsDefaultAtlasMembershipPass(
        class_data_energy=np.zeros_like(membership_maps),
        memberships=membership_maps.copy(),
        denominator=np.zeros(membership_maps.shape[:3], dtype=np.float32),
        objective=_float32(STATE["membership_objective_bits"]),
    )
    handoff = ToadsDefaultRegistrationHandoff(
        membership_pass=membership_pass,
        registration=ToadsRigidRegistrationResult(
            transform=transform.copy(), energy=np.float32(0.0), iterations=(),
            stop_reason="fixed synthetic registration output",
        ),
        transform=transform.copy(), rotation=np.eye(3, dtype=np.float32), shape_bounds=(),
    )
    case = JAVA_REFERENCE["cases"][0]
    assert case["name"] == "continue-grow-by-source" and case["stop"] is False
    branch = toads_post_registration_branch(
        handoff,
        previous_energy=_float32(case["previous_bits"]) if previous_energy is None else previous_energy,
        segmentation=segmentation,
        skeleton=skeleton,
        ordering=ordering,
        class_volumes=class_volumes,
        boundary_distances=boundary_distances,
        template_labels=np.asarray(STATE["template_labels"], dtype=np.int8),
        previous_classification=previous_classification,
        iteration=case["iteration"],
        maximum_iterations=case["maximum_iterations"],
        convergence_threshold=(
            _float32(case["threshold_bits"])
            if convergence_threshold is None
            else convergence_threshold
        ),
        relation_limit=1,
        do_growing=case["do_growing"],
    )
    return branch


def _label_string(values: np.ndarray, empty: bool = False) -> str:
    return "".join(
        "." if empty and int(value) == -1 else
        (str(int(value)) if int(value) < 10 else chr(ord("a") + int(value) - 10))
        for value in values.ravel()
    )


def _input_bytes(segmentation: np.ndarray, memberships: np.ndarray, shapes: np.ndarray,
                 adjacency: np.ndarray) -> bytes:
    return b"".join((
        np.asarray(segmentation, dtype=np.int8).tobytes(order="C"),
        np.asarray(memberships, dtype=">f4").tobytes(order="C"),
        np.asarray(shapes, dtype=">f4").tobytes(order="C"),
        np.asarray(adjacency, dtype=np.uint8).tobytes(order="C"),
    ))


def test_selected_registered_continue_branch_flows_into_java_referenced_next_thinning(tmp_path: Path) -> None:
    assert JAVA_REFERENCE["algorithm_source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    branch = _registered_branch()
    assert branch.growing is not None and branch.relations is not None
    assert branch.control.stop is False and branch.control.refresh_relations is True

    # The current atlas transform/samples and correction field are held fixed
    # at this bounded seam. This test advances only to the next thinning call.
    shapes = np.zeros_like(branch.handoff.membership_pass.memberships, dtype=np.float32)
    result = toads_post_registration_next_thinning(branch, shapes)
    assert result.branch is branch
    assert result.stage_order[-1] == "next_iteration.propagateThinning"

    labels = branch.growing.segmentation
    memberships = branch.handoff.membership_pass.memberships
    adjacency = branch.relations.relation.class_adjacency
    input_path = tmp_path / "next_thinning_input.bin"
    input_path.write_bytes(_input_bytes(labels, memberships, shapes, adjacency))
    completed = subprocess.run(
        ["sh", str(HARNESS_SCRIPT), str(input_path), *(str(v) for v in labels.shape), str(memberships.shape[-1])],
        check=True, capture_output=True, text=True,
    )
    java = json.loads(completed.stdout)
    thinning = result.thinning

    assert java["source"] == (
        "selected VandyCRUISE MultipleToadSegmentation.initThinningLabels + "
        "propagateThinning; short constructor useApprox=false"
    )
    assert java["source_sha256"] == "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert java["shape"] == [*labels.shape, memberships.shape[-1]]
    assert java["volume_ratio_bits"] == int(np.asarray([np.float32(0.5)]).view(np.int32)[0])
    assert java["limit_bits"] == int(np.asarray([np.float32(3.0)]).view(np.int32)[0])
    initial_volumes = np.zeros(memberships.shape[-1], dtype=np.float32)
    for x in range(1, labels.shape[0] - 1):
        for y in range(1, labels.shape[1] - 1):
            for z in range(1, labels.shape[2] - 1):
                k = int(labels[x, y, z])
                initial_volumes[k] += memberships[x, y, z, k]
    assert java["initial_volumes_bits"] == initial_volumes.view(np.int32).astype(int).tolist()
    assert java["minimum_volumes_bits"] == thinning.minimum_class_volumes.view(np.int32).astype(int).tolist()
    assert java["boundary_counts"] == thinning.boundary_counts.astype(int).tolist()
    assert java["processed"] == [
        [int(np.asarray(item.value, dtype=np.float32).view(np.int32)), item.x, item.y, item.z, item.label]
        for item in thinning.processed_items
    ]
    assert java["segmentation"] == _label_string(thinning.segmentation)
    assert java["skeleton"] == _label_string(thinning.skeleton, empty=True)
    assert java["available"] == "".join("1" if v else "0" for v in thinning.available.ravel())
    assert java["final_volumes_bits"] == thinning.class_volumes.view(np.int32).astype(int).tolist()
