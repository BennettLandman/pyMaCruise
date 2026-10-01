from __future__ import annotations

import json
import struct
from pathlib import Path

import numpy as np

from pymacruise.toads_default_atlas_membership import ToadsDefaultAtlasMembershipPass
from pymacruise.toads_default_exact_iteration_transition import (
    toads_default_exact_iteration_transition,
)


ROOT = Path(__file__).parent / "data" / "synthetic"
REFERENCE = ROOT / "toads_default_exact_iteration_transition" / "java_reference.json"


def _float32_from_bits(value: int) -> np.float32:
    return np.float32(struct.unpack("!f", struct.pack("!I", value))[0])


def _inputs():
    size, classes = 9, 4
    segmentation = np.zeros((size, size, size), dtype=np.int16)
    skeleton = np.full_like(segmentation, -1)
    ordering = np.empty((size, size, size), dtype=np.float32)
    memberships = np.empty((size, size, size, classes), dtype=np.float32)
    for x in range(size):
        for y in range(size):
            for z in range(size):
                ordering[x, y, z] = np.float32(-((x * 13 + y * 7 + z * 3) % 19)) / np.float32(10.0)
                for k in range(classes):
                    memberships[x, y, z, k] = np.float32(
                        ((x * 5 + y * 7 + z * 11 + k * 3) % 29 + 1)
                    ) / np.float32(30.0)
                center = 2 + ((y + z) % 3)
                if 2 <= y <= 6 and 2 <= z <= 6 and abs(x - center) <= 1:
                    segmentation[x, y, z] = 1
                if 5 <= x <= 6 and 2 <= y <= 3 and 5 <= z <= 6:
                    segmentation[x, y, z] = 2
                if x == 1 and 2 <= y <= 6 and z == 7:
                    segmentation[x, y, z] = 3
    center = 2 + ((4 + 3) % 3)
    skeleton[center, 4, 3] = 1
    skeleton[5, 2, 5] = 2
    volumes = np.asarray([1, 2, 3, 4], dtype=np.float32)
    return segmentation, skeleton, memberships, ordering, volumes


def _membership_pass(memberships: np.ndarray, objective: np.float32) -> ToadsDefaultAtlasMembershipPass:
    spatial, classes = memberships.shape[:3], memberships.shape[-1]
    return ToadsDefaultAtlasMembershipPass(
        class_data_energy=np.zeros((*spatial, classes), dtype=np.float32),
        memberships=memberships,
        denominator=np.ones(spatial, dtype=np.float32),
        objective=objective,
    )


def test_exact_membership_output_composes_selected_continue_and_stop_branches() -> None:
    reference = json.loads(REFERENCE.read_text())
    growth_reference = json.loads((ROOT / "toads_growing_pass" / "java_reference.json").read_text())
    source_control = json.loads((ROOT / "toads_iteration_control" / "java_reference.json").read_text())
    assert reference["algorithm_source_sha256"] == source_control["source_sha256"]
    assert reference["segmentation_source_sha256"] == growth_reference["source_sha256"]
    assert reference["fixed_registration"] is True

    segmentation, skeleton, memberships, ordering, volumes = _inputs()
    growth_state = growth_reference["relation_state_after_growth"]
    shape = segmentation.shape
    classes = memberships.shape[-1]
    boundary_distances = np.asarray(
        growth_state["initial_boundary_distances"], dtype=np.int8
    ).reshape((*shape, classes))

    for case in reference["cases"]:
        result = toads_default_exact_iteration_transition(
            _membership_pass(memberships.copy(), _float32_from_bits(case["current_bits"])),
            _float32_from_bits(case["previous_bits"]),
            case["iteration"],
            case["maximum_iterations"],
            _float32_from_bits(case["threshold_bits"]),
            segmentation,
            skeleton,
            ordering,
            volumes,
            template_labels=(0, 11, 22, 33) if case["stop"] else None,
            previous_classification=np.full(shape, 99, dtype=np.int8) if case["stop"] else None,
            boundary_distances=boundary_distances if not case["stop"] else None,
            relation_limit=growth_state["limit"] if not case["stop"] else None,
            do_growing=case["do_growing"],
        )
        assert result.control.stop is case["stop"]
        assert result.control.run_growing is case["run_growing"]
        assert int(np.asarray(result.control.energy_ratio).view(np.uint32)) == case["ratio_bits"]

        if case["stop"]:
            assert case["source_order"] == (
                "membership,control,growing,updateClassificationFromLabels"
                if case["run_growing"]
                else "membership,control,updateClassificationFromLabels"
            )
            assert not result.control.refresh_relations
            if result.growing is not None:
                assert "".join(str(int(v)) for v in result.growing.segmentation.ravel()) == growth_reference["segmentation"]
                expected = np.asarray((0, 11, 22, 33), dtype=np.int8)[result.growing.segmentation]
            else:
                expected = np.asarray((0, 11, 22, 33), dtype=np.int8)[segmentation]
            np.testing.assert_array_equal(result.classification, expected)
        else:
            assert case["source_order"] == (
                "membership,control,growing,computeRelations,computeBoundaryDistances,computeClosestLabels"
            )
            assert result.control.refresh_relations and not result.control.update_classification
            assert result.growing.events
            assert "".join(str(int(v)) for v in result.growing.segmentation.ravel()) == growth_reference["segmentation"]
            state = result.relation_state
            np.testing.assert_array_equal(state.relation.relation_map.ravel(), growth_state["relation_map"])
            np.testing.assert_array_equal(
                state.relation.relation_membership.ravel().astype(np.int8),
                growth_state["relation_membership"],
            )
            np.testing.assert_array_equal(state.boundary_distances.ravel(), growth_state["boundary_distances"])
            np.testing.assert_array_equal(state.closest_labels.ravel(), growth_state["closest_labels"])
