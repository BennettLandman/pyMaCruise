from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes

from pymacruise.lesion_statistics_membership_bridge import (
    run_lesion_toads_statistics_membership_bridge,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_statistics_membership_bridge"
HARNESS = Path(__file__).parent / "java_reference" / "LesionStatisticsMembershipBridgeReferenceHarness.java"
SOURCE_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
WRAPPER_SHA = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
TYPES = ("OBJECT", "LESION", "BLACKHOLE", "OUTMASK")
RELATIONS = ("direct", "prior", "direct", "prior")


def _bits(values: np.ndarray) -> list:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def _inputs() -> tuple[np.ndarray, ...]:
    image = nib.load(str(ROOT / "bridge_inputs_ras.nii.gz"))
    assert aff2axcodes(image.affine) == ("R", "A", "S")
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (6, 5, 4, 16)
    return (
        frames[..., :2],
        frames[..., 5:9],
        frames[..., 4].astype(np.int16),
        frames[..., 9],
        frames[..., 10],
        frames[..., 11],
        frames[..., 2:4],
        frames[..., 12:16],
    )


def _run(use_field: bool):
    images, memberships, labels, lesion, black, outlier, fields, shapes = _inputs()
    return run_lesion_toads_statistics_membership_bridge(
        images,
        memberships,
        labels,
        np.asarray([[0.22, 0.34, 0.46, 0.58], [0.30, 0.42, 0.54, 0.66]], dtype=np.float32),
        lesion,
        black,
        outlier,
        np.asarray([0.72, 0.83], dtype=np.float32),
        np.asarray([0.44, 0.55], dtype=np.float32),
        TYPES,
        fields=fields if use_field else None,
        centroid_options={
            "mode": "none",
            "modalities": ("T1_MPRAGE", "T2"),
            "lesion_power_exponents": np.asarray([1.0, 1.0], dtype=np.float32),
        },
        channel_weight_options={"outlier_distance": 2.75},
        structure_names=("Cerebrum-GM", "Lesion", "BlackHole", "Outmask"),
        prior_scale=0.67,
        initial_black_hole_similarity=np.ones(4, dtype=np.float32),
        smoothing=0.23,
        class_smoothing_priors=np.asarray([0.30, 0.45, 0.60, 0.75], dtype=np.float32),
        shape_priors=np.asarray([0.15, 0.25, 0.35, 0.45], dtype=np.float32),
        class_shapes=shapes,
        relation_types=RELATIONS,
        prior_coefficient=0.31,
        outlier_distance=2.75,
        previous_energy=3.25,
        t1_channel=0,
    )


def test_stats_to_membership_bridge_has_expected_stage_shapes_and_order() -> None:
    images, memberships, labels, lesion, black, outlier, fields, _ = _inputs()
    originals = tuple(a.copy() for a in (images, memberships, labels, lesion, black, outlier, fields))
    for use_field in (False, True):
        result = _run(use_field)
        assert result.statistics.centroids.class_centroids.shape == (2, 4)
        assert result.statistics.variances.class_standard_deviation.shape == (2, 4)
        assert result.statistics.channel_weights.shape == (2,)
        assert result.data_energies.class_energy.shape == (6, 5, 4, 3)
        assert result.outmask_class_energy.shape == (6, 5, 4)
        assert result.membership.class_memberships.shape == (6, 5, 4, 4)
        assert np.isfinite(result.energy_ratio)
        assert set(np.unique(labels)) == {0, 1, 2, 3}
        assert np.all(result.membership.denominator[1:-1, 1:-1, 1:-1] > 0)
    for value, before in zip((images, memberships, labels, lesion, black, outlier, fields), originals):
        np.testing.assert_array_equal(value, before)


def test_stats_to_membership_bridge_matches_java_intermediates_and_final_maps() -> None:
    reference = json.loads((ROOT / "java_reference.json").read_text())
    assert reference["source_sha256"] == SOURCE_SHA
    assert reference["wrapper_sha256"] == WRAPPER_SHA
    assert reference["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    for mode, use_field in (("no_field", False), ("with_field", True)):
        actual = _run(use_field)
        expected = reference[mode]
        values = {
            "centroids_bits": actual.statistics.centroids.class_centroids,
            "lesion_centroid_bits": actual.statistics.centroids.lesion_centroids,
            "black_hole_centroid_bits": actual.statistics.centroids.black_hole_centroids,
            "class_std_bits": actual.statistics.variances.class_standard_deviation,
            "lesion_std_bits": actual.statistics.variances.lesion_standard_deviation,
            "black_hole_std_bits": actual.statistics.variances.black_hole_standard_deviation,
            "weights_bits": actual.statistics.channel_weights,
            "class_similarity_bits": actual.similarity.class_similarity,
            "lesion_similarity_bits": actual.similarity.lesion_similarity,
            "raw_class_energy_bits": actual.data_energies.class_energy,
            "raw_lesion_energy_bits": actual.data_energies.lesion_energy,
            "raw_black_hole_energy_bits": actual.data_energies.black_hole_energy,
            "outmask_energy_bits": actual.outmask_class_energy,
            "class_memberships_bits": actual.membership.class_memberships,
            "lesion_membership_bits": actual.membership.lesion_membership,
            "black_hole_membership_bits": actual.membership.black_hole_membership,
            "outlier_membership_bits": actual.membership.outlier_membership,
            "denominator_bits": actual.membership.denominator,
            "objective_energy_bits": actual.membership.class_energy,
            "lesion_objective_energy_bits": actual.membership.lesion_energy,
            "black_hole_objective_energy_bits": actual.membership.black_hole_energy,
        }
        for key, values_array in values.items():
            actual_bits = np.asarray(values_array, dtype=np.float32).view(np.uint32).reshape(-1).astype(np.uint64).tolist()
            expected_bits = expected[key]
            if key == "raw_class_energy_bits":
                expected_bits = np.asarray(expected_bits, dtype=np.uint32).reshape(6, 5, 4, 4)[..., :3].reshape(-1).astype(np.uint64).tolist()
            assert actual_bits == expected_bits, f"{mode} {key}"
        assert int(np.asarray(actual.membership.objective).view(np.uint32)) == (expected["objective_bits"] & 0xFFFFFFFF)
        assert int(np.asarray(actual.energy_ratio).view(np.uint32)) == (expected["energy_ratio_bits"] & 0xFFFFFFFF)


def test_stats_to_membership_bridge_source_hash_pins() -> None:
    reference = json.loads((ROOT / "java_reference.json").read_text())
    assert reference["source_sha256"] == SOURCE_SHA
    assert reference["wrapper_sha256"] == WRAPPER_SHA
    assert reference["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
