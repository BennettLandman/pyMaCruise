from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_growth import propagate_lesion_toads_thinning
from pymacruise.lesion_next_iteration_thinning import run_lesion_toads_next_iteration_thinning
from pymacruise.lesion_post_registration import LesionToadsPostRegistrationBranch


ROOT = Path(__file__).parent / "data/synthetic/lesion_toads_thinning_pass"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text(encoding="utf-8"))


def _signed_bits(value: np.float32 | float) -> int:
    bits = int(np.asarray(value, dtype=np.float32).view(np.uint32))
    return bits if bits < 0x80000000 else bits - 0x100000000


def _branch(frames: np.ndarray) -> LesionToadsPostRegistrationBranch:
    membership = SimpleNamespace(
        class_memberships=frames[..., 2:5],
        lesion_membership=frames[..., 5],
        black_hole_membership=frames[..., 6],
    )
    growth = SimpleNamespace(
        skeleton=frames[..., 1].astype(np.int16),
        available=frames[..., 10:13] > 0.5,
        ordering=frames[..., 13],
    )
    return LesionToadsPostRegistrationBranch(
        membership_result=SimpleNamespace(membership=membership),
        registration_handoff=SimpleNamespace(),
        energy_ratio=np.float32(0.5),
        convergence_stop=False,
        iteration_limit_stop=False,
        stop=False,
        stop_reason=None,
        run_growing=True,
        growing=growth,
        segmentation=frames[..., 0].astype(np.int16),
        final_classification=None,
        relations=SimpleNamespace(),
        stage_order=("propagateGrowing", "computeRelations"),
    )


def test_continuing_branch_handoff_runs_java_pinned_next_thinning_pass() -> None:
    image = nib.load(str(ROOT / "thinning_pass_inputs.nii.gz"))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    frames = np.asarray(image.dataobj, dtype=np.float32)
    prior = _branch(frames)
    actual = run_lesion_toads_next_iteration_thinning(
        prior,
        frames[..., 7:10],
        object_types=("object", "lesion", "object"),
        lesion_class=1,
        volume_ratio=0.4,
        limit=1000.0,
    )
    expected_direct = propagate_lesion_toads_thinning(
        prior.segmentation,
        prior.membership_result.membership.class_memberships,
        ("object", "lesion", "object"),
        prior.membership_result.membership.lesion_membership,
        prior.membership_result.membership.black_hole_membership,
        frames[..., 7:10],
        prior.growing.skeleton,
        prior.growing.available,
        prior.growing.ordering,
        lesion_class=1,
        volume_ratio=0.4,
        limit=1000.0,
    )
    assert [_signed_bits(value) for value in actual.initialization.class_volume] == REFERENCE["initial_volume_bits"]
    assert [_signed_bits(value) for value in actual.minimum_class_volume] == REFERENCE["minimum_bits"]
    assert [
        [_signed_bits(event.popped_item.value), event.popped_item.x, event.popped_item.y,
         event.popped_item.z, event.popped_item.label, event.processed]
        for event in actual.run.events
    ] == REFERENCE["events"]
    np.testing.assert_array_equal(actual.run.segmentation, expected_direct.run.segmentation)
    np.testing.assert_array_equal(actual.run.skeleton, expected_direct.run.skeleton)
    np.testing.assert_array_equal(actual.run.available, expected_direct.run.available)
    assert [_signed_bits(value) for value in actual.run.ordering.ravel()] == REFERENCE["ordering_bits"]
    assert [_signed_bits(value) for value in actual.run.class_volume] == REFERENCE["volume_bits"]


def test_continuing_branch_handoff_rejects_terminal_or_mismatched_state() -> None:
    frames = np.asarray(nib.load(str(ROOT / "thinning_pass_inputs.nii.gz")).dataobj, dtype=np.float32)
    prior = _branch(frames)
    with pytest.raises(ValueError, match="must continue"):
        run_lesion_toads_next_iteration_thinning(
            LesionToadsPostRegistrationBranch(**{**prior.__dict__, "stop": True}),
            frames[..., 7:10], object_types=("object", "lesion", "object"),
            lesion_class=1, volume_ratio=0.4, limit=1000.0,
        )
    with pytest.raises(ValueError, match="class count"):
        run_lesion_toads_next_iteration_thinning(
            prior, frames[..., 7:10], object_types=("object", "lesion"),
            lesion_class=1, volume_ratio=0.4, limit=1000.0,
        )
