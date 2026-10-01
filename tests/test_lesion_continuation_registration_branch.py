from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_registration_branch_composition import finish_lesion_toads_iteration_after_membership
from test_lesion_continuing_iteration import TYPES, _compose, _load


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_continuation_registration_branch"
HARNESS = Path(__file__).parent / "java_reference" / "LesionContinuationRegistrationBranchReferenceHarness.java"
REGISTRATION_HARNESS = Path(__file__).parent / "java_reference" / "LesionContinuationRegistrationReferenceHarness.java"
WRAPPER_SHA = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
REGISTRATION_SHA = "718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e"
REGISTRATION_HARNESS_SHA = "c04529e906093a4eaa2c06662fb2be7be3efca0fdf9425b2b2d5ee0c5e6d45eb"


def test_p39v_memberships_feed_registration_and_selected_continuing_branch() -> None:
    reference = json.loads((ROOT / "java_reference.json").read_text())
    assert reference["source_sha256"] == WRAPPER_SHA
    assert reference["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    assert reference["registration_source_sha256"] == REGISTRATION_SHA
    assert reference["registration_wrapper_source_sha256"] == WRAPPER_SHA
    assert reference["registration_harness_sha256"] == REGISTRATION_HARNESS_SHA
    assert reference["registration_harness_sha256"] == hashlib.sha256(
        REGISTRATION_HARNESS.read_bytes()
    ).hexdigest()
    registration_reference = reference["registration_reference"]
    assert registration_reference["name"] == "p39v_membership"

    frames, _ = _load()
    previous_membership, previous_branch, continuing = _compose(frames, use_field=False)
    current = continuing.relation_membership.membership
    previous_energy = previous_membership.membership.objective
    assert int(np.asarray(previous_energy, dtype=np.float32).view(np.uint32)) == reference["previous_energy_bits"]
    assert int(np.asarray(current.objective, dtype=np.float32).view(np.uint32)) == reference["membership_energy_bits"]

    membership_before = current.class_memberships.copy()
    thinning_before = tuple(
        value.copy()
        for value in (
            continuing.thinning.run.segmentation,
            continuing.thinning.run.skeleton,
            continuing.thinning.run.ordering,
            continuing.thinning.run.class_volume,
        )
    )
    result = finish_lesion_toads_iteration_after_membership(
        continuing,
        previous_energy=float(previous_energy),
        shape_priors=tuple(frames[..., 12 + index] for index in range(4)),
        initial_transform=(0.03, -0.02, 0.01, 0.1, -0.1, 0.05),
        image_resolutions_mm=(0.9, 1.2, 1.5),
        template_resolutions_mm=(0.9, 1.2, 1.5),
        classification=np.zeros(frames.shape[:3], dtype=np.int8),
        template_labels=np.asarray([10, 11, 12, 13], dtype=np.int8),
        object_types=TYPES,
        lesion_class=1,
        iteration=1,
        n_iter_max=3,
        max_distance=0.01,
        algorithm_mode="standard",
        relation_limit=1,
    )

    assert result.registration.optimization is not None
    np.testing.assert_array_equal(result.registration.memberships, current.class_memberships)
    np.testing.assert_array_equal(current.class_memberships, membership_before)
    assert result.registration.optimization.iterations
    assert result.branch.registration_handoff is result.registration
    assert result.branch.membership_result is continuing.statistics_energy
    assert not result.branch.stop
    assert result.branch.run_growing
    assert result.branch.growing is not None
    assert result.branch.relations is not None
    assert int(result.branch.energy_ratio.view(np.uint32)) == reference["ratio_bits"]
    assert result.branch.stop is reference["stop"]
    assert result.branch.run_growing is reference["run_growing"]
    assert (result.branch.relations is not None) is reference["refresh_relations"]
    assert list(result.stage_order) == reference["stage_order"]
    assert result.registration.transform.shape == (6,)
    assert result.registration.rotation.shape == (3, 3)
    assert len(result.registration.shape_bounds) == len(TYPES)

    # The registration oracle consumes the actual P39v membership output and
    # the same folded shape priors passed at this composition boundary.
    np.testing.assert_array_equal(
        result.registration.transform,
        np.asarray(registration_reference["transform"], dtype=np.float32),
    )
    np.testing.assert_array_equal(
        result.registration.rotation,
        np.asarray(registration_reference["rotation"], dtype=np.float32),
    )
    for bounds, expected_bounds in zip(
        result.registration.shape_bounds, registration_reference["bounds"], strict=True
    ):
        assert bounds.lower_xyz == tuple(expected_bounds["lower"])
        assert bounds.upper_xyz == tuple(expected_bounds["upper"])
        assert bounds.support_count == expected_bounds["support_count"]
    optimization = result.registration.optimization
    assert optimization is not None
    assert optimization.stop_reason == registration_reference["stop_reason"]
    assert len(optimization.iterations) == len(registration_reference["iterations"])
    for step, expected_step in zip(
        optimization.iterations, registration_reference["iterations"], strict=True
    ):
        np.testing.assert_array_equal(
            step.trial_transforms,
            np.asarray(
                [trial["transform"] for trial in expected_step["trials"]],
                dtype=np.float32,
            ),
        )
        np.testing.assert_array_equal(
            step.trial_energies,
            np.asarray(
                [trial["energy"] for trial in expected_step["trials"]],
                dtype=np.float32,
            ),
        )
        np.testing.assert_array_equal(
            step.result.transform,
            np.asarray(expected_step["transform"], dtype=np.float32),
        )
        assert step.result.energy == np.float32(expected_step["energy"])
        assert step.result.lambda_value == np.float32(expected_step["lambda"])
        assert step.result.relative_change == np.float32(expected_step["relative_change"])
        assert step.result.accepted is expected_step["accepted"]
        assert step.result.rejected_trials == expected_step["rejected_trials"]

    # The P39v thinning state is passed through as the next wrapper state;
    # registration and growth must not mutate the supplied arrays.
    for value, before in zip(
        (
            continuing.thinning.run.segmentation,
            continuing.thinning.run.skeleton,
            continuing.thinning.run.ordering,
            continuing.thinning.run.class_volume,
        ),
        thinning_before,
    ):
        np.testing.assert_array_equal(value, before)
    assert result.branch.relations.relation_map.shape == frames.shape[:3]
