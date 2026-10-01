from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes

from pymacruise.lesion_repeated_iteration import run_lesion_toads_repeated_iterations
from test_lesion_continuing_iteration import TYPES, _compose, _load, _options


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_repeated_iteration"
HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsRepeatedIterationReferenceHarness.java"
LOOP_STATE_FIXTURE = ROOT / "loop_states_ras.nii.gz"
WRAPPER_SHA = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"


def _bits(value: float) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def _java_init_thinning_class_volume(
    segmentation: np.ndarray,
    class_memberships: np.ndarray,
    lesion_membership: np.ndarray,
    black_hole_membership: np.ndarray,
    lesion_class: int,
) -> np.ndarray:
    """Mirror selected initThinningLabels float32 volume reset/scan order."""
    volume = np.zeros(class_memberships.shape[-1], dtype=np.float32)
    for x in range(1, segmentation.shape[0] - 1):
        for y in range(1, segmentation.shape[1] - 1):
            for z in range(1, segmentation.shape[2] - 1):
                current = int(segmentation[x, y, z])
                if current == lesion_class:
                    increment = np.float32(
                        np.float32(class_memberships[x, y, z, current] + lesion_membership[x, y, z])
                        + black_hole_membership[x, y, z]
                    )
                else:
                    increment = class_memberships[x, y, z, current]
                volume[current] = np.float32(volume[current] + increment)
    return volume


def test_selected_two_pass_loop_carries_memberships_relations_and_transform() -> None:
    reference = json.loads((ROOT / "java_reference.json").read_text())
    assert reference["source_sha256"] == WRAPPER_SHA
    assert reference["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()

    image = nib.load(str(Path(__file__).parent / "data/synthetic/lesion_continuing_iteration/iteration_inputs_ras.nii.gz"))
    assert aff2axcodes(image.affine) == ("R", "A", "S")
    frames, _ = _load()
    initial_bridge, initial_branch, _ = _compose(frames, use_field=False)
    starting_inputs = tuple(
        value.copy()
        for value in (
            initial_branch.segmentation,
            initial_branch.membership_result.membership.class_memberships,
            initial_branch.relations.relation_map,
            initial_branch.registration_handoff.transform,
            frames,
        )
    )
    initial_transform = np.asarray([0.03, -0.02, 0.01, 0.1, -0.1, 0.05], dtype=np.float32)
    result = run_lesion_toads_repeated_iterations(
        initial_branch,
        frames[..., :2],
        frames[..., 12:16],
        frames[..., 12:16],
        shape_priors=tuple(frames[..., 12 + index] for index in range(4)),
        initial_transform=initial_transform,
        image_resolutions_mm=(0.9, 1.2, 1.5),
        template_resolutions_mm=(0.9, 1.2, 1.5),
        bridge_options=_options(frames, use_field=False),
        volume_ratio=0.0,
        thinning_limit=1.0e6,
        classification=np.zeros(frames.shape[:3], dtype=np.int8),
        template_labels=np.asarray([10, 11, 12, 13], dtype=np.int8),
        object_types=TYPES,
        lesion_class=1,
        n_iter_max=2,
        max_distance=0.0,
        algorithm_mode="standard",
        relation_limit=1,
        init_align_iter=1,
    )

    assert len(result.passes) == 2
    assert result.stop_reason == "iteration_limit"
    assert [_bits(item.previous_energy) for item in result.passes] == reference["energy_bits"][:-1]
    assert [
        _bits(item.continuing.statistics_energy.membership.objective)
        for item in result.passes
    ] == reference["energy_bits"][1:]
    assert [
        _bits(item.registration_branch.branch.energy_ratio) for item in result.passes
    ] == reference["ratio_bits"]
    assert [
        item.registration_branch.branch.stop_reason or "continue" for item in result.passes
    ] == reference["stop_reasons"]
    assert [list(order) for order in result.stage_order] == reference["stage_order"]

    first, second = result.passes
    assert not first.registration_branch.branch.stop
    assert first.registration_branch.branch.relations is not None
    assert first.registration_branch.branch.growing is not None
    assert second.registration_branch.branch.stop
    assert second.registration_branch.branch.stop_reason == "iteration_limit"
    assert second.registration_branch.branch.growing is not None  # positive signed ceiling
    assert second.registration_branch.branch.relations is None
    assert second.registration_branch.branch.final_classification is not None
    assert first.registration_branch.branch.final_classification is None

    np.testing.assert_array_equal(first.initial_transform, initial_transform)
    np.testing.assert_array_equal(
        second.initial_transform, first.registration_branch.registration.transform
    )
    np.testing.assert_array_equal(
        second.starting_class_memberships,
        first.registration_branch.branch.membership_result.membership.class_memberships,
    )
    np.testing.assert_array_equal(
        second.starting_segmentation, first.registration_branch.branch.segmentation
    )
    np.testing.assert_array_equal(
        second.starting_relation_map,
        first.registration_branch.branch.relations.relation_map,
    )
    np.testing.assert_array_equal(first.starting_skeleton, initial_branch.growing.skeleton)
    np.testing.assert_array_equal(first.starting_available, initial_branch.growing.available)
    np.testing.assert_array_equal(first.starting_ordering, initial_branch.growing.ordering)
    np.testing.assert_array_equal(first.starting_class_volume, initial_branch.growing.class_volume)
    previous_growth = first.registration_branch.branch.growing
    assert previous_growth is not None
    np.testing.assert_array_equal(second.starting_skeleton, previous_growth.skeleton)
    np.testing.assert_array_equal(second.starting_available, previous_growth.available)
    np.testing.assert_array_equal(second.starting_ordering, previous_growth.ordering)
    np.testing.assert_array_equal(second.starting_class_volume, previous_growth.class_volume)
    for pass_state in result.passes:
        expected_volume = _java_init_thinning_class_volume(
            pass_state.starting_segmentation,
            pass_state.starting_class_memberships,
            pass_state.starting_lesion_membership,
            pass_state.starting_black_hole_membership,
            lesion_class=1,
        )
        np.testing.assert_array_equal(
            pass_state.continuing.thinning.initialization.class_volume,
            expected_volume,
        )
        assert not np.array_equal(
            pass_state.continuing.thinning.initialization.class_volume,
            pass_state.starting_class_volume,
        )

    loop_state_nii = nib.load(str(LOOP_STATE_FIXTURE))
    assert aff2axcodes(loop_state_nii.affine) == ("R", "A", "S")
    assert hashlib.sha256(LOOP_STATE_FIXTURE.read_bytes()).hexdigest() == reference[
        "loop_state_sha256"
    ]
    loop_states = np.asarray(loop_state_nii.dataobj, dtype=np.float32)
    assert loop_states.shape == (*frames.shape[:3], 14)
    for pass_index, pass_state in enumerate(result.passes):
        base = pass_index * 7
        np.testing.assert_array_equal(
            loop_states[..., base].astype(np.int16), pass_state.starting_segmentation
        )
        np.testing.assert_array_equal(
            loop_states[..., base + 1 : base + 5], pass_state.starting_class_memberships
        )
        np.testing.assert_array_equal(
            loop_states[..., base + 5], pass_state.starting_lesion_membership
        )
        np.testing.assert_array_equal(
            loop_states[..., base + 6], pass_state.starting_black_hole_membership
        )
        java_volume_bits = [
            int(value) & 0xFFFFFFFF
            for value in reference["thinning_initial_volume_bits"][pass_index]
        ]
        python_volume_bits = (
            pass_state.continuing.thinning.initialization.class_volume
            .view(np.uint32)
            .astype(np.uint64)
            .tolist()
        )
        assert python_volume_bits == java_volume_bits
    for value, before in zip(
        (
            initial_branch.segmentation,
            initial_branch.membership_result.membership.class_memberships,
            initial_branch.relations.relation_map,
            initial_branch.registration_handoff.transform,
            frames,
        ),
        starting_inputs,
    ):
        np.testing.assert_array_equal(value, before)

    # P39zzzzx requires repeatability for the same supplied initialized state,
    # including numerical state carried between passes, not just equal stop
    # decisions from the Java control-flow harness.
    repeated = run_lesion_toads_repeated_iterations(
        initial_branch,
        frames[..., :2],
        frames[..., 12:16],
        frames[..., 12:16],
        shape_priors=tuple(frames[..., 12 + index] for index in range(4)),
        initial_transform=initial_transform,
        image_resolutions_mm=(0.9, 1.2, 1.5),
        template_resolutions_mm=(0.9, 1.2, 1.5),
        bridge_options=_options(frames, use_field=False),
        volume_ratio=0.0,
        thinning_limit=1.0e6,
        classification=np.zeros(frames.shape[:3], dtype=np.int8),
        template_labels=np.asarray([10, 11, 12, 13], dtype=np.int8),
        object_types=TYPES,
        lesion_class=1,
        n_iter_max=2,
        max_distance=0.0,
        algorithm_mode="standard",
        relation_limit=1,
        init_align_iter=1,
    )
    assert repeated.stop_reason == result.stop_reason
    assert repeated.stage_order == result.stage_order
    assert len(repeated.passes) == len(result.passes)
    for actual_pass, repeated_pass in zip(result.passes, repeated.passes):
        assert _bits(actual_pass.previous_energy) == _bits(repeated_pass.previous_energy)
        assert _bits(
            actual_pass.continuing.statistics_energy.membership.objective
        ) == _bits(repeated_pass.continuing.statistics_energy.membership.objective)
        assert _bits(actual_pass.registration_branch.branch.energy_ratio) == _bits(
            repeated_pass.registration_branch.branch.energy_ratio
        )
        assert actual_pass.registration_branch.branch.stop_reason == (
            repeated_pass.registration_branch.branch.stop_reason
        )
        for actual_array, repeated_array in (
            (actual_pass.starting_segmentation, repeated_pass.starting_segmentation),
            (actual_pass.starting_class_memberships, repeated_pass.starting_class_memberships),
            (actual_pass.starting_relation_map, repeated_pass.starting_relation_map),
            (actual_pass.starting_skeleton, repeated_pass.starting_skeleton),
            (actual_pass.starting_available, repeated_pass.starting_available),
            (actual_pass.starting_ordering, repeated_pass.starting_ordering),
            (actual_pass.starting_class_volume, repeated_pass.starting_class_volume),
            (
                actual_pass.continuing.thinning.initialization.class_volume,
                repeated_pass.continuing.thinning.initialization.class_volume,
            ),
            (
                actual_pass.registration_branch.registration.transform,
                repeated_pass.registration_branch.registration.transform,
            ),
            (
                actual_pass.registration_branch.branch.segmentation,
                repeated_pass.registration_branch.branch.segmentation,
            ),
            (
                actual_pass.registration_branch.branch.membership_result.membership.class_memberships,
                repeated_pass.registration_branch.branch.membership_result.membership.class_memberships,
            ),
        ):
            np.testing.assert_array_equal(actual_array, repeated_array)

    negative = run_lesion_toads_repeated_iterations(
        initial_branch,
        frames[..., :2],
        frames[..., 12:16],
        frames[..., 12:16],
        shape_priors=tuple(frames[..., 12 + index] for index in range(4)),
        initial_transform=initial_transform,
        image_resolutions_mm=(0.9, 1.2, 1.5),
        template_resolutions_mm=(0.9, 1.2, 1.5),
        bridge_options=_options(frames, use_field=False),
        volume_ratio=0.0,
        thinning_limit=1.0e6,
        classification=np.zeros(frames.shape[:3], dtype=np.int8),
        template_labels=np.asarray([10, 11, 12, 13], dtype=np.int8),
        object_types=TYPES,
        lesion_class=1,
        n_iter_max=int(reference["negative_ceiling"]["n_iter_max"]),
        max_distance=0.0,
        algorithm_mode="standard",
        relation_limit=1,
        init_align_iter=1,
    )
    negative_final = negative.final_branch
    assert len(negative.passes) == 2
    assert negative.stop_reason == reference["negative_ceiling"]["stop_reason"]
    assert negative_final.stop is reference["negative_ceiling"]["stop"]
    assert negative_final.run_growing is reference["negative_ceiling"]["run_growing"]
    assert negative_final.growing is None
    assert negative_final.final_classification is not None
    assert negative_final.relations is None

    for case in reference["convergence_cases"]:
        convergence_threshold = float(
            np.asarray(case["max_distance_bits"], dtype=np.uint32).view(np.float32)
        )
        assert convergence_threshold > float(
            np.asarray(case["ratio_bits"], dtype=np.uint32).view(np.float32)
        )
        converged = run_lesion_toads_repeated_iterations(
            initial_branch,
            frames[..., :2],
            frames[..., 12:16],
            frames[..., 12:16],
            shape_priors=tuple(frames[..., 12 + index] for index in range(4)),
            initial_transform=initial_transform,
            image_resolutions_mm=(0.9, 1.2, 1.5),
            template_resolutions_mm=(0.9, 1.2, 1.5),
            bridge_options=_options(frames, use_field=False),
            volume_ratio=0.0,
            thinning_limit=1.0e6,
            classification=np.zeros(frames.shape[:3], dtype=np.int8),
            template_labels=np.asarray([10, 11, 12, 13], dtype=np.int8),
            object_types=TYPES,
            lesion_class=1,
            n_iter_max=case["n_iter_max"],
            max_distance=convergence_threshold,
            algorithm_mode="standard",
            relation_limit=1,
            init_align_iter=1,
        )
        assert len(converged.passes) == case["passes"] == 1
        convergence_branch = converged.final_branch
        assert converged.stop_reason == case["stop_reason"] == "convergence"
        assert convergence_branch.stop is case["stop"] is True
        assert convergence_branch.convergence_stop is True
        assert convergence_branch.iteration_limit_stop is False
        assert convergence_branch.run_growing is case["run_growing"]
        assert (convergence_branch.growing is not None) is case["run_growing"]
        assert convergence_branch.final_classification is not None
        assert convergence_branch.relations is None
        assert _bits(convergence_branch.energy_ratio) == case["ratio_bits"]
        assert [list(order) for order in converged.stage_order] == [case["stage_order"]]
