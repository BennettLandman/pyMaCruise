from __future__ import annotations

from dataclasses import replace
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes

from pymacruise.lesion_post_registration import lesion_toads_post_registration_branch
from pymacruise.lesion_registration_handoff import register_lesion_membership_handoff
from test_lesion_statistics_membership_bridge import _run as _bridge


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_post_registration"
HARNESS = Path(__file__).parent / "java_reference" / "LesionPostRegistrationBranchReferenceHarness.java"
SOURCE_ROOT = (
    Path(__file__).parents[2]
    / "toads-cruise-scm-2026-09-26"
    / "VandyCRUISE"
    / "src"
    / "edu"
    / "jhmi"
    / "rad"
    / "medic"
)
SOURCE_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
WRAPPER_SHA = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
ATLAS_SHA = "718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e"
TYPES = ("object", "lesion", "blackhole", "outmask")


def _arrays() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, tuple[np.ndarray, ...]]:
    image = nib.load(str(ROOT / "branch_inputs_ras.nii.gz"))
    assert aff2axcodes(image.affine) == ("R", "A", "S")
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (6, 5, 4, 19)
    assert not np.array_equal(frames[..., 16], frames[::-1, ..., 16])
    assert np.array_equal(
        frames[..., :16],
        np.asarray(
            nib.load(
                str(Path(__file__).parent / "data/synthetic/lesion_statistics_membership_bridge/bridge_inputs_ras.nii.gz")
            ).dataobj,
            dtype=np.float32,
        ),
    )
    shapes = tuple(frames[..., 12 + index] for index in range(4))
    return (
        frames[..., 16].astype(np.int16),
        frames[..., 17].astype(np.int16),
        frames[..., 18].copy(),
        frames,
        shapes,
    )


def _f32(value: float) -> np.float32:
    return np.float32(value)


def _bits(value: float) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def _context(use_field: bool = False):
    bridge = _bridge(use_field)
    _, _, _, frames, shapes = _arrays()
    handoff = register_lesion_membership_handoff(
        bridge,
        shapes,
        (0.03, -0.02, 0.01, 0.1, -0.1, 0.05),
        image_resolutions_mm=(0.9, 1.2, 1.5),
        template_resolutions_mm=(0.9, 1.2, 1.5),
    )
    return bridge, handoff, frames


def _branch(bridge, handoff, segmentation, skeleton, ordering, *, previous, current, iteration, n_iter_max, threshold=0.1, mode="standard"):
    membership = replace(bridge.membership, objective=_f32(current))
    fixed_bridge = replace(bridge, membership=membership)
    return lesion_toads_post_registration_branch(
        fixed_bridge,
        handoff,
        previous_energy=previous,
        segmentation=segmentation,
        skeleton=skeleton,
        ordering=ordering,
        class_volume=np.ones(4, dtype=np.float32),
        object_types=TYPES,
        classification=np.zeros(segmentation.shape, dtype=np.int8),
        template_labels=np.asarray([10, 11, 12, 13], dtype=np.int8),
        lesion_class=1,
        iteration=iteration,
        n_iter_max=n_iter_max,
        max_distance=threshold,
        algorithm_mode=mode,
        relation_limit=1,
    )


def test_selected_java_branch_control_matrix_and_stage_order() -> None:
    ref = json.loads((ROOT / "java_reference.json").read_text())
    assert ref["source_sha256"] == WRAPPER_SHA
    assert ref["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    assert ref["post_registration_order"] == [
        "computeMemberships", "energyRatio", "registerShapes", "convergenceCheck",
        "iterationLimitCheck", "propagateGrowing", "classificationOrRelationRefresh",
    ]
    bridge, handoff, _ = _context()
    segmentation, skeleton, ordering, _, _ = _arrays()
    originals = [
        array.copy()
        for array in (
            segmentation,
            skeleton,
            ordering,
            handoff.memberships,
            bridge.membership.lesion_membership,
            bridge.membership.black_hole_membership,
            handoff.transform,
        )
    ]
    for case in ref["cases"]:
        previous, current = (3.25, 2.5) if "iteration_stop" in case["name"] or case["name"] in ("continue", "narrow_band_suppression") else (3.25, 3.20)
        iteration = 3 if "iteration_stop" in case["name"] else 1
        n_iter_max = -3 if case["name"].endswith("off") else 3
        mode = "narrow_band_competition" if case["name"] == "narrow_band_suppression" else "standard"
        actual = _branch(
            bridge, handoff, segmentation, skeleton, ordering,
            previous=previous, current=current, iteration=iteration,
            n_iter_max=n_iter_max, mode=mode,
        )
        assert _bits(actual.energy_ratio) == case["ratio_bits"], case["name"]
        assert actual.stop is case["stop"], case["name"]
        assert actual.run_growing is case["run_growing"], case["name"]
        assert (actual.final_classification is not None) is case["update_classification"]
        assert (actual.relations is not None) is case["refresh_relations"]
        expected = ["membership_energy", "atlas_alignment", "energy_ratio", "convergence_check", "iteration_limit_check"]
        if case["run_growing"]:
            expected.append("propagateGrowing")
        expected.extend(("updateClassificationFromLabels", "updateLesionClassifcation") if case["stop"] else ("computeRelations",))
        assert list(actual.stage_order) == expected
        if not actual.stop:
            assert actual.relations is not None
            assert not hasattr(actual.relations, "closest_labels")
            assert not hasattr(actual.relations, "boundary_distances")
    for value, before in zip(
        (
            segmentation,
            skeleton,
            ordering,
            handoff.memberships,
            bridge.membership.lesion_membership,
            bridge.membership.black_hole_membership,
            handoff.transform,
        ),
        originals,
    ):
        np.testing.assert_array_equal(value, before)


def test_folded_ras_growth_queue_has_accepted_rejected_and_empty_cases() -> None:
    bridge, handoff, frames = _context()
    segmentation = frames[..., 16].astype(np.int16)
    skeleton = frames[..., 17].astype(np.int16)
    ordering = frames[..., 18].astype(np.float32)
    before = tuple(value.copy() for value in (segmentation, skeleton, ordering))
    result = _branch(
        bridge, handoff, segmentation, skeleton, ordering,
        previous=3.25, current=2.5, iteration=1, n_iter_max=3,
    )
    assert result.growing is not None and result.relations is not None
    assert result.growing.events
    assert any(event.topology_accepted for event in result.growing.events)
    assert any(event.topology_accepted is False for event in result.growing.events)
    assert np.array_equal(result.growing.segmentation[[0, -1], :, :], segmentation[[0, -1], :, :])
    assert np.array_equal(result.growing.segmentation[:, [0, -1], :], segmentation[:, [0, -1], :])
    assert np.array_equal(result.growing.segmentation[:, :, [0, -1]], segmentation[:, :, [0, -1]])
    for value, original in zip((segmentation, skeleton, ordering), before):
        np.testing.assert_array_equal(value, original)

    empty = np.zeros_like(skeleton)
    no_queue = _branch(
        bridge, handoff, segmentation, empty, ordering,
        previous=3.25, current=2.5, iteration=1, n_iter_max=3,
    )
    assert no_queue.growing is not None
    assert no_queue.growing.events == ()
    assert no_queue.growing.remaining_queue == ()


def test_stopping_branch_transfers_labels_and_preserves_strict_subtype_tie() -> None:
    bridge, handoff, frames = _context()
    segmentation = np.zeros((6, 5, 4), dtype=np.int16)
    segmentation[2, 2, 1] = 1
    segmentation[2, 2, 2] = 1
    skeleton = np.zeros_like(segmentation)
    ordering = frames[..., 18].copy()
    membership = bridge.membership
    lesion = membership.lesion_membership.copy()
    black = membership.black_hole_membership.copy()
    lesion[2, 2, 1] = np.float32(0.9)
    black[2, 2, 1] = np.float32(0.9)  # strict source comparison gives lesion subtype on tie
    lesion[2, 2, 2] = np.float32(0.8)
    black[2, 2, 2] = np.float32(0.9)
    adjusted = replace(membership, objective=np.float32(3.20), lesion_membership=lesion, black_hole_membership=black)
    fixed_bridge = replace(bridge, membership=adjusted)
    result = lesion_toads_post_registration_branch(
        fixed_bridge, handoff, previous_energy=3.25, segmentation=segmentation,
        skeleton=skeleton, ordering=ordering, class_volume=np.ones(4, np.float32),
        object_types=TYPES,
        classification=np.zeros_like(segmentation, dtype=np.int8),
        template_labels=np.asarray([10, 11, 12, 13], dtype=np.int8), lesion_class=1,
        iteration=1, n_iter_max=-3, max_distance=0.1, algorithm_mode="standard",
        relation_limit=1,
    )
    assert result.stop and not result.run_growing
    assert result.stop_reason == "convergence"
    assert result.final_classification is not None
    assert result.final_classification.classification[2, 2, 1] == 11
    assert result.final_classification.lesion_classification[2, 2, 1] == 1
    assert result.final_classification.lesion_classification[2, 2, 2] == 2
    assert result.relations is None


def test_source_java_files_and_fixture_provenance_are_pinned() -> None:
    wrapper = SOURCE_ROOT / "algorithms" / "AlgorithmLesionToads.java"
    source = SOURCE_ROOT / "methods" / "LesionToadSegmentation.java"
    atlas = SOURCE_ROOT / "methods" / "LesionToadAtlas.java"
    assert hashlib.sha256(wrapper.read_bytes()).hexdigest() == WRAPPER_SHA
    assert hashlib.sha256(source.read_bytes()).hexdigest() == SOURCE_SHA
    assert hashlib.sha256(atlas.read_bytes()).hexdigest() == ATLAS_SHA
