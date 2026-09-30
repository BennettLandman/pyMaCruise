from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes

from pymacruise.lesion_continuing_iteration import run_lesion_toads_continuing_iteration
from pymacruise.lesion_continuing_post_membership import (
    run_lesion_toads_continuing_post_membership_branch,
)
from pymacruise.lesion_post_registration import lesion_toads_post_registration_branch
from pymacruise.lesion_registration_handoff import register_lesion_membership_handoff
from pymacruise.lesion_statistics_membership_bridge import run_lesion_toads_statistics_membership_bridge


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_continuing_iteration"
REFERENCE = Path(__file__).parent / "data" / "synthetic" / "lesion_continuing_post_membership" / "java_reference.json"
HARNESS = Path(__file__).parent / "java_reference" / "LesionContinuingPostMembershipReferenceHarness.java"
WRAPPER_SHA = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
ATLAS_SHA = "718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e"
SEGMENTATION_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
TYPES = ("object", "lesion", "blackhole", "outmask")
INITIAL_CENTROIDS = np.asarray([[0.22, 0.34, 0.46, 0.58], [0.30, 0.42, 0.54, 0.66]], dtype=np.float32)
INITIAL_LESION_CENTROID = np.asarray([0.72, 0.83], dtype=np.float32)
INITIAL_BLACK_CENTROID = np.asarray([0.44, 0.55], dtype=np.float32)
SMOOTHING_PRIORS = np.asarray([0.30, 0.45, 0.60, 0.75], dtype=np.float32)
ATLAS_PRIORS = np.asarray([0.15, 0.25, 0.35, 0.45], dtype=np.float32)
CENTROID_OPTIONS = {
    "mode": "none",
    "modalities": ("T1_MPRAGE", "T2"),
    "lesion_power_exponents": np.asarray([1.0, 1.0], dtype=np.float32),
}
CHANNEL_OPTIONS = {"outlier_distance": 2.75}
STRUCTURES = ("Cerebrum-GM", "Lesion", "BlackHole", "Outmask")


def _inputs():
    image = nib.load(str(ROOT / "iteration_inputs_ras.nii.gz"))
    assert aff2axcodes(image.affine) == ("R", "A", "S")
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (15, 15, 15, 19)
    assert not np.array_equal(frames[..., 16], frames[::-1, ..., 16])
    bridge = run_lesion_toads_statistics_membership_bridge(
        frames[..., :2], frames[..., 5:9], frames[..., 4].astype(np.int16),
        INITIAL_CENTROIDS, frames[..., 9], frames[..., 10], frames[..., 11],
        INITIAL_LESION_CENTROID, INITIAL_BLACK_CENTROID, TYPES,
        fields=frames[..., 2:4], centroid_options=CENTROID_OPTIONS,
        channel_weight_options=CHANNEL_OPTIONS, structure_names=STRUCTURES,
        prior_scale=0.67, initial_black_hole_similarity=np.ones(4, dtype=np.float32),
        smoothing=0.23, class_smoothing_priors=SMOOTHING_PRIORS,
        shape_priors=ATLAS_PRIORS, class_shapes=frames[..., 12:16],
        relation_types=("direct",) * 4, prior_coefficient=0.31,
        outlier_distance=2.75, previous_energy=3.25, t1_channel=0,
    )
    shape_priors = tuple(frames[..., 12 + k] for k in range(4))
    first_registration = register_lesion_membership_handoff(
        bridge, shape_priors, np.zeros(6, dtype=np.float32),
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0), align_shapes=False,
    )
    first_branch = lesion_toads_post_registration_branch(
        bridge, first_registration,
        previous_energy=float(np.float32(bridge.membership.objective * np.float32(2.0))),
        segmentation=frames[..., 16].astype(np.int16), skeleton=frames[..., 17].astype(np.int16),
        ordering=frames[..., 18], class_volume=np.ones(4, dtype=np.float32),
        object_types=TYPES, classification=np.zeros(frames.shape[:3], dtype=np.int8),
        template_labels=np.asarray([10, 11, 12, 13], dtype=np.int8), lesion_class=1,
        iteration=1, n_iter_max=4, max_distance=0.0,
        algorithm_mode="standard", relation_limit=1,
    )
    assert not first_branch.stop and first_branch.growing is not None and first_branch.relations is not None
    options = {
        "object_types": TYPES, "fields": frames[..., 2:4],
        "centroid_options": CENTROID_OPTIONS, "channel_weight_options": CHANNEL_OPTIONS,
        "structure_names": STRUCTURES, "prior_scale": 0.67,
        "initial_black_hole_similarity": np.ones(4, dtype=np.float32),
        "smoothing": 0.23, "class_smoothing_priors": SMOOTHING_PRIORS,
        "shape_priors": ATLAS_PRIORS, "class_shapes": frames[..., 12:16],
        "prior_coefficient": 0.31, "outlier_distance": 2.75, "t1_channel": 0,
    }
    iteration_result = run_lesion_toads_continuing_iteration(
        first_branch, frames[..., :2], frames[..., 12:16], frames[..., 12:16],
        volume_ratio=0.0, thinning_limit=1.0e6, bridge_options=options,
    )
    return frames, shape_priors, bridge, first_branch, iteration_result


def _run(frames, shape_priors, previous_bridge, previous_branch, iteration_result, *, threshold, iteration=2, n_iter_max=4, align_shapes=False):
    return run_lesion_toads_continuing_post_membership_branch(
        previous_branch, iteration_result, shape_priors,
        previous_energy=float(previous_bridge.membership.objective),
        image_resolutions_mm=(1.0, 1.0, 1.0), template_resolutions_mm=(1.0, 1.0, 1.0),
        object_types=TYPES, classification=np.zeros(frames.shape[:3], dtype=np.int8),
        template_labels=np.asarray([10, 11, 12, 13], dtype=np.int8), lesion_class=1,
        iteration=iteration, n_iter_max=n_iter_max, max_distance=threshold,
        algorithm_mode="standard", relation_limit=1, align_shapes=align_shapes,
    )


def test_actual_p39v_output_registration_stop_continue_and_next_iteration_handoff() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["wrapper_sha256"] == WRAPPER_SHA
    assert reference["atlas_sha256"] == ATLAS_SHA
    assert reference["segmentation_sha256"] == SEGMENTATION_SHA
    assert reference["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    assert reference["order"] == [
        "propagateThinning", "computeCentroids", "computeVariances", "computeChannelWeights",
        "computeMemberships", "energyRatio", "registerShapes", "convergenceCheck",
        "iterationLimitCheck", "propagateGrowing", "classificationOrRelationRefresh",
    ]

    frames, shapes, prior_bridge, prior_branch, p39v = _inputs()
    original = [frames.copy(), prior_branch.segmentation.copy(), p39v.thinning.run.segmentation.copy(),
                p39v.statistics_energy.membership.class_memberships.copy()]

    # Active registration consumes the actual per-voxel relation-aware P39v map.
    continuing = _run(
        frames, shapes, prior_bridge, prior_branch, p39v,
        threshold=0.0, iteration=2, n_iter_max=4, align_shapes=True,
    )
    assert not continuing.stop and continuing.growing is not None and continuing.relations is not None
    assert continuing.registration_handoff.optimization is not None
    assert np.array_equal(continuing.registration_handoff.memberships, p39v.statistics_energy.membership.class_memberships)
    assert continuing.energy_ratio.view(np.uint32) == p39v.energy_ratio.view(np.uint32)
    class_one_dispatch = p39v.relation_membership.direct_dispatch[1:-1, 1:-1, 1:-1, 1]
    assert class_one_dispatch.any() and (~class_one_dispatch).any()
    assert continuing.stage_order == (
        *p39v.stage_order, "registerShapes", "convergenceCheck", "iterationLimitCheck",
        "propagateGrowing", "computeRelations",
    )
    assert continuing.growing.events
    assert any(event.topology_accepted is True for event in continuing.growing.events)
    assert any(event.topology_accepted is False for event in continuing.growing.events)
    assert np.array_equal(continuing.segmentation[[0, -1], :, :], p39v.thinning.run.segmentation[[0, -1], :, :])

    # Stop branch: convergence with negative max suppresses final growth; positive
    # max permits the source's final growth before classification transfer.
    stopped_no_growth = _run(
        frames, shapes, prior_bridge, prior_branch, p39v,
        threshold=10.0, iteration=2, n_iter_max=-4,
    )
    assert stopped_no_growth.stop and stopped_no_growth.convergence_stop
    assert stopped_no_growth.growing is None and stopped_no_growth.final_classification is not None
    assert stopped_no_growth.relations is None
    stopped_with_growth = _run(
        frames, shapes, prior_bridge, prior_branch, p39v,
        threshold=10.0, iteration=2, n_iter_max=4,
    )
    assert stopped_with_growth.stop and stopped_with_growth.run_growing
    assert stopped_with_growth.growing is not None and stopped_with_growth.final_classification is not None

    # The signed iteration ceiling stops in both cases, but only positive max
    # performs final growth. Continued output itself is a valid P39v predecessor.
    limit_off = _run(
        frames, shapes, prior_bridge, prior_branch, p39v,
        threshold=0.0, iteration=4, n_iter_max=-4,
    )
    limit_on = _run(
        frames, shapes, prior_bridge, prior_branch, p39v,
        threshold=0.0, iteration=4, n_iter_max=4,
    )
    assert limit_off.iteration_limit_stop and limit_off.stop and not limit_off.run_growing
    assert limit_on.iteration_limit_stop and limit_on.stop and limit_on.run_growing

    next_iteration = run_lesion_toads_continuing_iteration(
        continuing, frames[..., :2], frames[..., 12:16], frames[..., 12:16],
        volume_ratio=0.0, thinning_limit=1.0e6,
        bridge_options={
            "object_types": TYPES, "fields": frames[..., 2:4],
            "centroid_options": CENTROID_OPTIONS, "channel_weight_options": CHANNEL_OPTIONS,
            "structure_names": STRUCTURES, "prior_scale": 0.67,
            "initial_black_hole_similarity": np.ones(4, dtype=np.float32),
            "smoothing": 0.23, "class_smoothing_priors": SMOOTHING_PRIORS,
            "shape_priors": ATLAS_PRIORS, "class_shapes": frames[..., 12:16],
            "prior_coefficient": 0.31, "outlier_distance": 2.75, "t1_channel": 0,
        },
    )
    assert next_iteration.stage_order[0] == "propagateThinning"
    for actual, before in zip(
        (frames, prior_branch.segmentation, p39v.thinning.run.segmentation,
         p39v.statistics_energy.membership.class_memberships), original,
    ):
        np.testing.assert_array_equal(actual, before)


def test_java_excerpt_reference_controls_and_explicit_previous_energy() -> None:
    reference = json.loads(REFERENCE.read_text())
    cases = {case["name"]: case for case in reference["controls"]}
    assert cases["continue"]["stop"] is False and cases["continue"]["refresh_relations"] is True
    assert cases["threshold_stop_growth_off"]["run_growing"] is False
    assert cases["threshold_stop_growth_on"]["run_growing"] is True
    assert cases["iteration_stop_growth_off"]["final_outputs"] is True
    assert cases["iteration_stop_growth_on"]["run_growing"] is True

    frames, shapes, previous_bridge, previous_branch, p39v = _inputs()
    try:
        run_lesion_toads_continuing_post_membership_branch(
            previous_branch, p39v, shapes,
            previous_energy=float(previous_bridge.membership.objective) + 1.0,
            image_resolutions_mm=(1.0, 1.0, 1.0), template_resolutions_mm=(1.0, 1.0, 1.0),
            object_types=TYPES, classification=np.zeros(frames.shape[:3], dtype=np.int8),
            template_labels=np.asarray([10, 11, 12, 13], dtype=np.int8), lesion_class=1,
            iteration=2, n_iter_max=4, max_distance=0.0,
            algorithm_mode="standard", relation_limit=1, align_shapes=False,
        )
    except ValueError as exc:
        assert "previous_energy" in str(exc)
    else:
        raise AssertionError("mismatched previous energy was accepted")
