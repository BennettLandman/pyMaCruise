from __future__ import annotations

import json
import struct
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_control import toads_classification_from_labels
from pymacruise.toads_default_atlas_membership import ToadsDefaultAtlasMembershipPass
from pymacruise.toads_default_registration_handoff import ToadsDefaultRegistrationHandoff
from pymacruise.toads_post_registration import toads_post_registration_branch
from pymacruise.toads_registration_sampling import ToadsRigidRegistrationResult
from pymacruise.toads_relations import toads_refresh_relation_state
from pymacruise.toads_thinning import toads_propagate_growing
from test_toads_default_registration_handoff import _arguments as _p18czh_arguments


ROOT = Path(__file__).parent / "data/synthetic/toads_post_registration_branch"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
STATE = json.loads((ROOT / "state.json").read_text())


def _load(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / name))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    return np.asarray(image.dataobj)


def _float32(bits: int) -> np.float32:
    return np.float32(struct.unpack("!f", struct.pack("!I", bits))[0])


def _arrays() -> dict[str, np.ndarray]:
    names = (
        "membership_maps.nii.gz", "segmentation.nii.gz", "skeleton.nii.gz",
        "ordering.nii.gz", "class_volumes.nii.gz", "boundary_distances.nii.gz",
        "previous_classification.nii.gz",
    )
    values = {name.removesuffix(".nii.gz"): _load(name) for name in names}
    values["membership_maps"] = values["membership_maps"].astype(np.float32)
    values["segmentation"] = values["segmentation"].astype(np.int16)
    values["skeleton"] = values["skeleton"].astype(np.int16)
    values["ordering"] = values["ordering"].astype(np.float32)
    values["class_volumes"] = values["class_volumes"].reshape(11).astype(np.float32)
    values["boundary_distances"] = values["boundary_distances"].astype(np.int8)
    values["previous_classification"] = values["previous_classification"].astype(np.int8)
    return values


def _handoff(values: dict[str, np.ndarray]) -> ToadsDefaultRegistrationHandoff:
    membership = ToadsDefaultAtlasMembershipPass(
        class_data_energy=np.zeros_like(values["membership_maps"], dtype=np.float32),
        memberships=values["membership_maps"].copy(),
        denominator=np.zeros(values["membership_maps"].shape[:3], dtype=np.float32),
        objective=_float32(STATE["membership_objective_bits"]),
    )
    transform = np.asarray(STATE["registration_transform_bits"], dtype=np.uint32).view(np.float32)
    registration = ToadsRigidRegistrationResult(
        transform=transform.copy(), energy=np.float32(0.0), iterations=(), stop_reason="fixed synthetic registration output"
    )
    return ToadsDefaultRegistrationHandoff(
        membership_pass=membership,
        registration=registration,
        transform=transform.copy(),
        rotation=np.eye(3, dtype=np.float32),
        shape_bounds=(),
    )


def _run(values: dict[str, np.ndarray], handoff: ToadsDefaultRegistrationHandoff, case: dict):
    template = np.asarray(STATE["template_labels"], dtype=np.int8)
    return toads_post_registration_branch(
        handoff,
        previous_energy=_float32(case["previous_bits"]),
        segmentation=values["segmentation"],
        skeleton=values["skeleton"],
        ordering=values["ordering"],
        class_volumes=values["class_volumes"],
        boundary_distances=values["boundary_distances"],
        template_labels=template,
        previous_classification=values["previous_classification"],
        iteration=case["iteration"],
        maximum_iterations=case["maximum_iterations"],
        convergence_threshold=_float32(case["threshold_bits"]),
        relation_limit=1,
        do_growing=case["do_growing"],
    )


def test_post_registration_control_and_branch_order_match_java_reference() -> None:
    assert REFERENCE["algorithm_source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert REFERENCE["post_registration_order"] == [
        "energyRatio", "registerShapesAlreadyComplete", "convergenceCheck",
        "iterationLimitCheck", "propagateGrowing", "classificationOrRelationRefresh",
    ]
    values = _arrays()
    handoff = _handoff(values)
    original = {name: array.copy() for name, array in values.items()}
    membership = handoff.membership_pass.memberships.copy()
    transform = handoff.transform.copy()

    for case in REFERENCE["cases"]:
        result = _run(values, handoff, case)
        ratio_bits = int(np.asarray([result.control.energy_ratio], dtype=np.float32).view(np.uint32)[0])
        assert ratio_bits == case["ratio_bits"]
        assert result.control.stop is case["stop"]
        assert result.control.run_growing is case["run_growing"]
        assert result.control.update_classification is case["update_classification"]
        assert result.control.refresh_relations is case["refresh_relations"]
        assert list(result.stage_order) == case["stage_order"]
        assert (result.growing is not None) is case["growing"]

        if result.growing is not None:
            expected_grow = toads_propagate_growing(
                values["segmentation"], values["skeleton"], membership,
                values["ordering"], values["class_volumes"],
            )
            np.testing.assert_array_equal(result.growing.segmentation, expected_grow.segmentation)
            np.testing.assert_array_equal(result.growing.skeleton, expected_grow.skeleton)
            np.testing.assert_array_equal(result.growing.available, expected_grow.available)
            np.testing.assert_array_equal(result.growing.ordering, expected_grow.ordering)
            np.testing.assert_array_equal(result.growing.class_volumes, expected_grow.class_volumes)
            assert len(result.growing.events) == len(expected_grow.events) == 1

        if result.control.stop:
            labels = result.growing.segmentation if result.growing else values["segmentation"]
            expected_classification, expected_changed = toads_classification_from_labels(
                labels, np.asarray(STATE["template_labels"], dtype=np.int8),
                values["previous_classification"],
            )
            np.testing.assert_array_equal(result.classification, expected_classification)
            assert result.changed_classification_voxels == expected_changed
            assert result.relations is None
        else:
            assert result.classification is None
            assert result.changed_classification_voxels is None
            expected_segmentation = result.growing.segmentation
            expected_relations = toads_refresh_relation_state(
                expected_segmentation,
                membership,
                values["boundary_distances"],
                1,
            )
            assert result.relations is not None
            np.testing.assert_array_equal(result.relations.relation.relation_map, expected_relations.relation.relation_map)
            np.testing.assert_array_equal(result.relations.relation.relation_membership, expected_relations.relation.relation_membership)
            np.testing.assert_array_equal(result.relations.boundary_distances, expected_relations.boundary_distances)
            np.testing.assert_array_equal(result.relations.closest_labels, expected_relations.closest_labels)

    for name, before in original.items():
        np.testing.assert_array_equal(values[name], before)
    np.testing.assert_array_equal(handoff.membership_pass.memberships, membership)
    np.testing.assert_array_equal(handoff.transform, transform)


def test_fixture_memberships_are_the_java_default_dispatch_output_and_true_ras() -> None:
    membership_reference = json.loads(
        (ROOT.parent / "toads_default_atlas_dispatch/java_reference.json").read_text()
    )
    case = membership_reference["cases"][STATE["membership_source_case"]]
    expected = np.asarray(case["membership_bits"], dtype=np.uint32)
    actual = _load("membership_maps.nii.gz").astype(np.float32).ravel()
    np.testing.assert_array_equal(actual.view(np.uint32), expected)
    assert int(np.asarray([case["objective_bits"]], dtype=np.uint32)[0]) == STATE["membership_objective_bits"]


def test_real_p18czh_handoff_flows_into_continuing_post_registration_branch() -> None:
    # Exercise the actual preceding composition, including its registered
    # transform, rather than only constructing a compatible fixed-state object.
    _, _, _, _, options, _ = _p18czh_arguments()
    from pymacruise.toads_default_registration_handoff import toads_default_atlas_registration_handoff

    handoff = toads_default_atlas_registration_handoff(**options)
    assert handoff.registration is not None
    assert not np.array_equal(handoff.transform, np.asarray(options["initial_transform"], dtype=np.float32))

    values = _arrays()
    # The integration intentionally supplies the actual P18czh membership map
    # and objective, while holding all post-membership segmentation state fixed.
    result = toads_post_registration_branch(
        handoff,
        previous_energy=handoff.membership_pass.objective,
        segmentation=values["segmentation"],
        skeleton=values["skeleton"],
        ordering=values["ordering"],
        class_volumes=values["class_volumes"],
        boundary_distances=values["boundary_distances"],
        template_labels=np.asarray(STATE["template_labels"], dtype=np.int8),
        previous_classification=values["previous_classification"],
        iteration=1,
        maximum_iterations=5,
        convergence_threshold=0.0,
        relation_limit=1,
        do_growing=False,
    )
    assert result.handoff is handoff
    assert result.control.stop is False
    assert result.control.run_growing is True
    assert result.control.refresh_relations is True
    assert result.growing is not None and result.relations is not None
    assert result.classification is None
    assert result.stage_order[-3:] == (
        "computeRelations", "computeBoundaryDistances", "computeClosestLabels"
    )
    np.testing.assert_array_equal(result.handoff.membership_pass.memberships, handoff.membership_pass.memberships)
    np.testing.assert_array_equal(result.handoff.transform, handoff.registration.transform)
