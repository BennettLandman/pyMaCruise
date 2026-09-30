from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes

from pymacruise.lesion_fixed_input_loop import run_lesion_toads_fixed_input_loop
from test_lesion_continuing_post_membership import _inputs


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_continuing_iteration"
REFERENCE = Path(__file__).parent / "data" / "synthetic" / "lesion_fixed_input_loop" / "java_reference.json"
HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsFixedInputLoopReferenceHarness.java"
WRAPPER_SHA = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
TYPES = ("object", "lesion", "blackhole", "outmask")
CENTROID_OPTIONS = {
    "mode": "none",
    "modalities": ("T1_MPRAGE", "T2"),
    "lesion_power_exponents": np.asarray([1.0, 1.0], dtype=np.float32),
}
CHANNEL_OPTIONS = {"outlier_distance": 2.75}
STRUCTURES = ("Cerebrum-GM", "Lesion", "BlackHole", "Outmask")
SMOOTHING_PRIORS = np.asarray([0.30, 0.45, 0.60, 0.75], dtype=np.float32)
SHAPE_PRIORS = np.asarray([0.15, 0.25, 0.35, 0.45], dtype=np.float32)


def _bridge_options(frames: np.ndarray) -> dict:
    return {
        "object_types": TYPES,
        "fields": frames[..., 2:4],
        "centroid_options": CENTROID_OPTIONS,
        "channel_weight_options": CHANNEL_OPTIONS,
        "structure_names": STRUCTURES,
        "prior_scale": 0.67,
        "initial_black_hole_similarity": np.ones(4, dtype=np.float32),
        "smoothing": 0.23,
        "class_smoothing_priors": SMOOTHING_PRIORS,
        "shape_priors": SHAPE_PRIORS,
        "class_shapes": frames[..., 12:16],
        "prior_coefficient": 0.31,
        "outlier_distance": 2.75,
        "t1_channel": 0,
    }


def _run(frames, shapes, branch, *, n_iter_max: int, max_distance: float):
    return run_lesion_toads_fixed_input_loop(
        branch,
        frames[..., :2],
        frames[..., 12:16],
        frames[..., 12:16],
        shapes,
        starting_iteration=1,
        n_iter_max=n_iter_max,
        max_distance=max_distance,
        relation_limit=1,
        volume_ratio=0.0,
        thinning_limit=1.0e6,
        bridge_options=_bridge_options(frames),
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        object_types=TYPES,
        classification=np.zeros(frames.shape[:3], dtype=np.int8),
        template_labels=np.asarray([10, 11, 12, 13], dtype=np.int8),
        lesion_class=1,
        align_shapes=False,
    )


def test_repeated_loop_matches_java_control_and_carries_state_until_stop() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["wrapper_sha256"] == WRAPPER_SHA
    assert reference["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    assert reference["iteration_order"] == [
        "propagateThinning", "computeCentroids", "computeVariances", "computeChannelWeights",
        "computeMemberships", "energyRatio", "registerShapes", "convergenceCheck",
        "iterationLimitCheck", "propagateGrowing", "classificationOrRelationRefresh",
    ]
    java_cases = {case["name"]: case for case in reference["cases"]}
    assert [item["iteration"] for item in java_cases["two_continues_then_converges"]["passes"]] == [2, 3, 4]
    assert [item["stop"] for item in java_cases["two_continues_then_converges"]["passes"]] == [False, False, True]
    assert java_cases["signed_ceiling_growth_on"]["passes"][-1]["run_growing"] is True
    assert java_cases["signed_ceiling_growth_off"]["passes"][-1]["run_growing"] is False

    fixture = nib.load(str(ROOT / "iteration_inputs_ras.nii.gz"))
    assert aff2axcodes(fixture.affine) == ("R", "A", "S")
    frames, shapes, initial_bridge, initial_branch, _ = _inputs()
    assert frames.shape[:3] == (15, 15, 15)
    initial_segmentation = initial_branch.segmentation.copy()
    initial_memberships = initial_branch.membership_result.membership.class_memberships.copy()

    growth_enabled = _run(frames, shapes, initial_branch, n_iter_max=4, max_distance=0.0)
    assert [entry.iteration for entry in growth_enabled.passes] == [2, 3, 4]
    assert [entry.post_membership_branch.stop for entry in growth_enabled.passes] == [False, False, True]
    assert growth_enabled.stop_reason == "iteration_limit"
    assert growth_enabled.final_branch.run_growing
    assert growth_enabled.final_branch.growing is not None
    assert growth_enabled.final_branch.final_classification is not None
    assert growth_enabled.final_branch.relations is None
    assert growth_enabled.passes[0].post_membership_branch.relations is not None
    assert growth_enabled.passes[1].post_membership_branch.relations is not None
    assert all(pass_result.continuing_iteration.stage_order[0] == "propagateThinning" for pass_result in growth_enabled.passes)
    assert all("registerShapes" in pass_result.post_membership_branch.stage_order for pass_result in growth_enabled.passes)

    growth_disabled = _run(frames, shapes, initial_branch, n_iter_max=-4, max_distance=0.0)
    assert [entry.iteration for entry in growth_disabled.passes] == [2, 3, 4]
    assert growth_disabled.stop_reason == "iteration_limit"
    assert growth_disabled.final_branch.stop and not growth_disabled.final_branch.run_growing
    assert growth_disabled.final_branch.growing is None
    assert growth_disabled.final_branch.final_classification is not None

    convergence = _run(frames, shapes, initial_branch, n_iter_max=4, max_distance=1.0e10)
    assert len(convergence.passes) == 1
    assert convergence.final_branch.convergence_stop
    assert convergence.stop_reason == "convergence"

    np.testing.assert_array_equal(initial_branch.segmentation, initial_segmentation)
    np.testing.assert_array_equal(
        initial_branch.membership_result.membership.class_memberships, initial_memberships
    )
    assert initial_bridge is not None  # fixture setup remains explicit and deterministic


def test_repeat_loop_rejects_inconsistent_initial_state_and_zero_ceiling() -> None:
    frames, shapes, _, branch, _ = _inputs()
    options = _bridge_options(frames)
    args = (
        branch, frames[..., :2], frames[..., 12:16], frames[..., 12:16], shapes,
    )
    kwargs = {
        "starting_iteration": 1, "n_iter_max": 4, "max_distance": 0.0,
        "relation_limit": 1, "volume_ratio": 0.0, "thinning_limit": 1.0e6,
        "bridge_options": options, "image_resolutions_mm": (1.0, 1.0, 1.0),
        "template_resolutions_mm": (1.0, 1.0, 1.0), "object_types": TYPES,
        "classification": np.zeros(frames.shape[:3], dtype=np.int8),
        "template_labels": np.asarray([10, 11, 12, 13], dtype=np.int8),
        "lesion_class": 1, "align_shapes": False,
    }
    try:
        run_lesion_toads_fixed_input_loop(*args, **{**kwargs, "n_iter_max": 0})
    except ValueError as exc:
        assert "n_iter_max" in str(exc)
    else:
        raise AssertionError("zero iteration ceiling was accepted")

