from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes

from pymacruise.lesion_continuing_iteration import run_lesion_toads_continuing_iteration
from pymacruise.lesion_data_energy_maps import compute_lesion_toads_data_energy_maps
from pymacruise.lesion_growth import propagate_lesion_toads_thinning
from pymacruise.lesion_post_registration import lesion_toads_post_registration_branch
from pymacruise.lesion_registration_handoff import register_lesion_membership_handoff
from pymacruise.lesion_relation_aware_membership import run_lesion_toads_relation_aware_membership_pass
from pymacruise.lesion_relations import LesionRelationList
from pymacruise.lesion_statistics_membership_bridge import run_lesion_toads_statistics_membership_bridge


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_continuing_iteration"
HARNESS = Path(__file__).parent / "java_reference" / "LesionRelationAwareMembershipReferenceHarness.java"
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
WRAPPER_SHA = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
SEGMENTATION_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
ATLAS_SHA = "718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e"
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


def _load() -> tuple[np.ndarray, np.ndarray]:
    image = nib.load(str(ROOT / "iteration_inputs_ras.nii.gz"))
    assert aff2axcodes(image.affine) == ("R", "A", "S")
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (15, 15, 15, 19)
    assert not np.array_equal(frames[..., 16], frames[::-1, ..., 16])
    return frames, image.affine


def _bridge(frames: np.ndarray, use_field: bool):
    shapes = frames[..., 12:16]
    return run_lesion_toads_statistics_membership_bridge(
        frames[..., :2],
        frames[..., 5:9],
        frames[..., 4].astype(np.int16),
        INITIAL_CENTROIDS,
        frames[..., 9],
        frames[..., 10],
        frames[..., 11],
        INITIAL_LESION_CENTROID,
        INITIAL_BLACK_CENTROID,
        TYPES,
        fields=frames[..., 2:4] if use_field else None,
        centroid_options=CENTROID_OPTIONS,
        channel_weight_options=CHANNEL_OPTIONS,
        structure_names=STRUCTURES,
        prior_scale=0.67,
        initial_black_hole_similarity=np.ones(4, dtype=np.float32),
        smoothing=0.23,
        class_smoothing_priors=SMOOTHING_PRIORS,
        shape_priors=ATLAS_PRIORS,
        class_shapes=shapes,
        relation_types=("direct",) * 4,
        prior_coefficient=0.31,
        outlier_distance=2.75,
        previous_energy=3.25,
        t1_channel=0,
    )


def _options(frames: np.ndarray, use_field: bool) -> dict:
    return {
        "object_types": TYPES,
        "fields": frames[..., 2:4] if use_field else None,
        "centroid_options": CENTROID_OPTIONS,
        "channel_weight_options": CHANNEL_OPTIONS,
        "structure_names": STRUCTURES,
        "prior_scale": 0.67,
        "initial_black_hole_similarity": np.ones(4, dtype=np.float32),
        "smoothing": 0.23,
        "class_smoothing_priors": SMOOTHING_PRIORS,
        "shape_priors": ATLAS_PRIORS,
        "class_shapes": frames[..., 12:16],
        "relation_types": ("direct",) * 4,
        "prior_coefficient": 0.31,
        "outlier_distance": 2.75,
        "t1_channel": 0,
    }


def _compose(frames: np.ndarray, use_field: bool):
    previous_bridge = _bridge(frames, use_field)
    shapes = tuple(frames[..., 12 + k] for k in range(4))
    registration = register_lesion_membership_handoff(
        previous_bridge,
        shapes,
        (0, 0, 0, 0, 0, 0),
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        align_shapes=False,
    )
    energy = np.float32(previous_bridge.membership.objective)
    branch = lesion_toads_post_registration_branch(
        previous_bridge,
        registration,
        previous_energy=float(np.float32(energy * np.float32(2.0))),
        segmentation=frames[..., 16].astype(np.int16),
        skeleton=frames[..., 17].astype(np.int16),
        ordering=frames[..., 18],
        class_volume=np.ones(4, dtype=np.float32),
        object_types=TYPES,
        classification=np.zeros(frames.shape[:3], dtype=np.int8),
        template_labels=np.asarray([10, 11, 12, 13], dtype=np.int8),
        lesion_class=1,
        iteration=1,
        n_iter_max=3,
        max_distance=0.01,
        algorithm_mode="standard",
        relation_limit=1,
    )
    assert not branch.stop and branch.growing is not None and branch.relations is not None
    originals = [
        array.copy()
        for array in (
            frames,
            previous_bridge.membership.class_memberships,
            previous_bridge.membership.lesion_membership,
            previous_bridge.membership.black_hole_membership,
            registration.transform,
            branch.relations.relation_map,
            branch.relations.relation_membership,
            branch.growing.segmentation,
            branch.growing.skeleton,
            branch.growing.available,
            branch.growing.ordering,
            branch.growing.class_volume,
        )
    ]
    result = run_lesion_toads_continuing_iteration(
        branch,
        frames[..., :2],
        frames[..., 12:16],
        frames[..., 12:16],
        volume_ratio=0.0,
        thinning_limit=1.0e6,
        bridge_options=_options(frames, use_field),
    )
    for value, original in zip(
        (
            frames,
            previous_bridge.membership.class_memberships,
            previous_bridge.membership.lesion_membership,
            previous_bridge.membership.black_hole_membership,
            registration.transform,
            branch.relations.relation_map,
            branch.relations.relation_membership,
            branch.growing.segmentation,
            branch.growing.skeleton,
            branch.growing.available,
            branch.growing.ordering,
            branch.growing.class_volume,
        ),
        originals,
    ):
        np.testing.assert_array_equal(value, original)
    return previous_bridge, branch, result


def test_actual_continuing_handoff_thinning_statistics_energy_and_memberships() -> None:
    reference = json.loads((ROOT / "relation_dispatch_java_reference.json").read_text())
    assert reference["source_sha256"] == SEGMENTATION_SHA
    assert reference["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    assert reference["class1_dispatch"] == [False, True, True, False, False, True, True, False]
    assert reference["direct_bits"] == int(np.asarray(np.float32(0.5)).view(np.uint32))
    assert reference["prior_bits"] == int(np.asarray(np.float32(0.25)).view(np.uint32))
    frames, _ = _load()
    previous, branch, result = _compose(frames, use_field=False)

    assert result.stage_order == (
        "propagateThinning",
        "computeCentroids",
        "computeVariances",
        "computeChannelWeights",
        "computeMembershipsDataEnergy",
        "computeMembershipsRelationAwareSweep",
        "energyRatio",
    )
    assert result.statistics_energy.membership is result.relation_membership.membership
    assert int(result.statistics_energy.energy_ratio.view(np.uint32)) == int(result.energy_ratio.view(np.uint32))
    assert result.thinning.run.events
    assert any(event.topology_accepted is False for event in result.thinning.run.events)
    assert result.thinning.run.events[-1].popped_item.value < 0.0
    direct = result.relation_membership.direct_dispatch[1:-1, 1:-1, 1:-1, 1]
    assert direct.any() and (~direct).any()  # one class is direct and prior across the volume
    spatial = frames.shape[:3]
    assert result.relation_membership.membership.class_memberships.shape == (*spatial, 4)
    assert result.relation_membership.membership.lesion_membership.shape == spatial
    assert result.relation_membership.membership.black_hole_membership.shape == spatial
    assert np.array_equal(
        result.thinning.run.segmentation[[0, -1], :, :], branch.segmentation[[0, -1], :, :]
    )
    assert np.array_equal(
        result.relation_membership.membership.class_memberships[[0, -1], :, :],
        previous.membership.class_memberships[[0, -1], :, :],
    )
    empty_frontier = propagate_lesion_toads_thinning(
        np.zeros(frames.shape[:3], dtype=np.int16),
        previous.membership.class_memberships,
        TYPES,
        previous.membership.lesion_membership,
        previous.membership.black_hole_membership,
        frames[..., 12:16],
        np.full(frames.shape[:3], -1, dtype=np.int16),
        np.zeros((*frames.shape[:3], 4), dtype=np.bool_),
        frames[..., 18],
        lesion_class=1,
        volume_ratio=0.0,
        limit=1.0,
    )
    assert empty_frontier.initialization.inserted_items == ()
    assert empty_frontier.run.events == ()
    expected_energy = np.float32(
        np.float32(2.0)
        * np.float32(previous.membership.objective - result.relation_membership.membership.objective)
    )
    expected_energy = np.float32(expected_energy / np.float32(previous.membership.objective + result.relation_membership.membership.objective))
    expected_energy = np.float32(abs(expected_energy))
    assert int(result.energy_ratio.view(np.uint32)) == int(expected_energy.view(np.uint32))

    # P39u's raw weighted-residual map remains the input to the next update.
    stats = result.statistics_energy.statistics
    raw = compute_lesion_toads_data_energy_maps(
        frames[..., :2],
        stats.centroids.class_centroids[:, :3].T,
        stats.centroids.lesion_centroids,
        stats.centroids.black_hole_centroids,
        stats.channel_weights,
        ("object", "lesion", "blackhole"),
        use_field=False,
    )
    np.testing.assert_array_equal(result.statistics_energy.data_energies.class_energy, raw.class_energy)
    np.testing.assert_array_equal(result.statistics_energy.data_energies.lesion_energy, raw.lesion_energy)
    np.testing.assert_array_equal(result.statistics_energy.data_energies.black_hole_energy, raw.black_hole_energy)


def test_relation_dispatch_changes_membership_at_same_voxel_and_field_refresh_runs() -> None:
    frames, _ = _load()
    previous, branch, result = _compose(frames, use_field=True)
    relation = branch.relations
    assert relation is not None
    membership = previous.membership
    next_bridge = result.statistics_energy
    class_data = np.zeros_like(membership.class_memberships, dtype=np.float32)
    class_data[..., :3] = next_bridge.data_energies.class_energy
    class_data[..., 3] = next_bridge.outmask_class_energy

    def with_relation_value(value: bool):
        table = np.full_like(relation.relation_membership, value, dtype=np.bool_)
        forced = LesionRelationList(relation.relation_map.copy(), table, relation.diffusion_iterations)
        return run_lesion_toads_relation_aware_membership_pass(
            membership.class_memberships,
            membership.lesion_membership,
            membership.black_hole_membership,
            membership.outlier_membership,
            class_data,
            next_bridge.data_energies.lesion_energy,
            next_bridge.data_energies.black_hole_energy,
            next_bridge.outmask_data_energy,
            frames[..., 12:16],
            SMOOTHING_PRIORS,
            tuple(value.upper() for value in TYPES),
            forced,
            smoothing=0.23,
            shape_priors=ATLAS_PRIORS,
            class_similarities=next_bridge.similarity.class_similarity,
            lesion_similarities=next_bridge.similarity.lesion_similarity,
            black_hole_similarities=next_bridge.similarity.black_hole_similarity,
            prior_coefficient=0.31,
        )

    direct = with_relation_value(True).membership.class_memberships
    prior = with_relation_value(False).membership.class_memberships
    assert not np.array_equal(direct[1:-1, 1:-1, 1:-1], prior[1:-1, 1:-1, 1:-1])
    assert np.any(result.relation_membership.membership.lesion_membership != membership.lesion_membership)
    assert np.any(result.relation_membership.membership.black_hole_membership != membership.black_hole_membership)
    stats = next_bridge.statistics
    raw_field = compute_lesion_toads_data_energy_maps(
        frames[..., :2],
        stats.centroids.class_centroids[:, :3].T,
        stats.centroids.lesion_centroids,
        stats.centroids.black_hole_centroids,
        stats.channel_weights,
        ("object", "lesion", "blackhole"),
        use_field=True,
        intensity_fields=frames[..., 2:4],
    )
    raw_unfielded = compute_lesion_toads_data_energy_maps(
        frames[..., :2],
        stats.centroids.class_centroids[:, :3].T,
        stats.centroids.lesion_centroids,
        stats.centroids.black_hole_centroids,
        stats.channel_weights,
        ("object", "lesion", "blackhole"),
        use_field=False,
    )
    np.testing.assert_array_equal(next_bridge.data_energies.class_energy, raw_field.class_energy)
    assert not np.array_equal(raw_field.class_energy, raw_unfielded.class_energy)


def test_single_voxel_relation_aware_membership_matches_java_bits() -> None:
    reference = json.loads((ROOT / "relation_dispatch_java_reference.json").read_text())
    expected = reference["single_voxel"]
    shape = (3, 3, 3)
    classes = 4
    initial = np.full((*shape, classes), np.float32(0.2), dtype=np.float32)
    lesion = np.full(shape, np.float32(0.1), dtype=np.float32)
    black = np.full(shape, np.float32(0.1), dtype=np.float32)
    outlier = np.full(shape, np.float32(0.25), dtype=np.float32)
    class_data = np.zeros((*shape, classes), dtype=np.float32)
    class_data[1, 1, 1] = (2.0, 3.0, 4.0, 5.0)
    lesion_data = np.zeros(shape, dtype=np.float32)
    black_data = np.zeros(shape, dtype=np.float32)
    outmask_data = np.zeros(shape, dtype=np.float32)
    lesion_data[1, 1, 1] = 6.0
    black_data[1, 1, 1] = 7.0
    outmask_data[1, 1, 1] = 8.0
    shape_values = np.ones((*shape, classes), dtype=np.float32)
    shape_values[1, 1, 1] = (0.30, 0.25, 0.20, 0.45)
    relation_map = np.zeros(shape, dtype=np.int16)
    relation_map[1, 1, 1] = 1
    table = np.asarray([[True, False, True, False], [True, True, False, False]])
    relation = LesionRelationList(relation_map, table)
    result = run_lesion_toads_relation_aware_membership_pass(
        initial,
        lesion,
        black,
        outlier,
        class_data,
        lesion_data,
        black_data,
        outmask_data,
        shape_values,
        np.zeros(classes, dtype=np.float32),
        ("OBJECT", "LESION", "BLACKHOLE", "OUTMASK"),
        relation,
        smoothing=0.0,
        shape_priors=np.zeros(classes, dtype=np.float32),
        class_similarities=np.zeros((classes, classes), dtype=np.float32),
        lesion_similarities=np.zeros(classes, dtype=np.float32),
        black_hole_similarities=np.zeros(classes, dtype=np.float32),
        prior_coefficient=0.0,
    )
    membership = result.membership
    assert result.direct_dispatch[1, 1, 1].tolist() == [True, True, False, False]
    bits = lambda values: np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()
    assert bits(membership.class_memberships[1, 1, 1]) == expected["class_bits"]
    assert bits([membership.lesion_membership[1, 1, 1]])[0] == expected["lesion_bits"]
    assert bits([membership.black_hole_membership[1, 1, 1]])[0] == expected["black_hole_bits"]
    assert bits([membership.outlier_membership[1, 1, 1]])[0] == expected["outlier_bits"]
    assert bits([membership.denominator[1, 1, 1]])[0] == expected["denominator_bits"]
    assert bits([membership.objective])[0] == expected["objective_bits"]


def test_selected_java_sources_are_pinned() -> None:
    wrapper = SOURCE_ROOT / "algorithms" / "AlgorithmLesionToads.java"
    segmentation = SOURCE_ROOT / "methods" / "LesionToadSegmentation.java"
    atlas = SOURCE_ROOT / "methods" / "LesionToadAtlas.java"
    assert hashlib.sha256(wrapper.read_bytes()).hexdigest() == WRAPPER_SHA
    assert hashlib.sha256(segmentation.read_bytes()).hexdigest() == SEGMENTATION_SHA
    assert hashlib.sha256(atlas.read_bytes()).hexdigest() == ATLAS_SHA
