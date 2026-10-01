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
WRAPPER_SHA = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"


def _bits(value: float) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


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
