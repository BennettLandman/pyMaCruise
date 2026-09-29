from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_data_terms import compute_lesion_toads_unfielded_data_energies
from pymacruise.lesion_membership import (
    convert_lesion_toads_class_membership,
    modulate_lesion_toads_blackhole_memberships_by_distances,
    modulate_lesion_toads_membership_by_distances,
    normalize_lesion_toads_voxel_memberships,
    update_lesion_toads_prior_energies,
    update_lesion_toads_smoothing_energies,
)


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_fielded_data_energy_membership/java_reference.json"


def _from_bits(value: int) -> np.float32:
    return np.asarray(value & 0xFFFFFFFF, dtype=np.uint32).view(np.float32)[()]


def _bits(value: float | np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def _shape_factor(shape: np.float32) -> np.float32:
    twice = np.float32(np.float32(2.0) * shape)
    return np.float32(1.0) if np.float32(1.0) < twice else twice


def test_fielded_data_energy_flows_through_membership_pipeline() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert len(reference["cases"]) == 2
    for case in reference["cases"]:
        target = case["target_type"]
        types = ["OBJECT", target]
        relations = ["direct", case["relation_type"]]
        image = np.asarray([_from_bits(v) for v in case["image_bits"]])
        intensity_field = np.asarray([_from_bits(v) for v in case["field_bits"]])
        centroids = np.asarray([[_from_bits(v) for v in row] for row in case["centroid_bits"]])
        lesion_centroid = np.asarray([_from_bits(v) for v in case["lesion_centroid_bits"]])
        black_hole_centroid = np.asarray([_from_bits(v) for v in case["black_hole_centroid_bits"]])
        weights = np.asarray([_from_bits(v) for v in case["channel_weight_bits"]])
        data = []
        for index, kind in enumerate(types):
            object_type = "object" if kind == "OBJECT" else kind.lower()
            data.append(
                compute_lesion_toads_unfielded_data_energies(
                    image,
                    centroids[index],
                    lesion_centroid,
                    black_hole_centroid,
                    weights,
                    object_type=object_type,
                    use_field=True,
                    intensity_field=intensity_field,
                )
            )
        assert [[_bits(v) for v in (e.class_energy, e.lesion_energy, e.black_hole_energy)] for e in data] == case["class_data_bits"], case

        shapes = np.asarray([_from_bits(v) for v in case["shape_values_bits"]])
        shape_priors = np.asarray([_from_bits(v) for v in case["shape_prior_bits"]])
        similarities = np.asarray([[_from_bits(v) for v in row] for row in case["similarity_bits"]])
        lesion_sim = np.asarray([_from_bits(v) for v in case["lesion_similarity_bits"]])
        black_hole_sim = np.asarray([_from_bits(v) for v in case["black_hole_similarity_bits"]])
        neighbors = np.asarray([[_from_bits(v) for v in row] for row in case["neighbor_class_bits"]])
        nlesion = np.asarray([_from_bits(v) for v in case["neighbor_lesion_bits"]])
        nblack = np.asarray([_from_bits(v) for v in case["neighbor_black_hole_bits"]])

        memberships: list[np.float32] = []
        energies: list[np.float32] = []
        lesion = np.float32(0.0)
        black_hole = np.float32(0.0)
        relation_value = np.float32(0.0)
        lesion_energy = np.float32(0.0)
        black_hole_energy = np.float32(0.0)
        for index, kind in enumerate(types):
            is_target = index == 1
            smooth = update_lesion_toads_smoothing_energies(
                neighbors,
                nlesion,
                nblack,
                class_index=index,
                class_type=kind,
                smoothing_prior=float(_from_bits(case["prior_bits"])),
                smoothing=float(_from_bits(case["smoothing_bits"])),
                class_membership_data_term=float(data[index].class_energy),
                lesion_membership_data_term=(float(data[index].lesion_energy) if is_target else 0.0),
                black_hole_membership_data_term=(float(data[index].black_hole_energy) if is_target else 0.0),
                outmask_membership_data_term=0.0,
                class_energy=float(data[index].class_energy),
                lesion_energy=(float(data[index].lesion_energy) if is_target else 0.0),
                black_hole_energy=(float(data[index].black_hole_energy) if is_target else 0.0),
            )
            prior = update_lesion_toads_prior_energies(
                shape_priors,
                similarities[index],
                lesion_sim,
                black_hole_sim,
                shapes,
                class_index=index,
                class_type=kind,
                prior_coefficient=float(_from_bits(case["prior_coefficient_bits"])),
                class_membership_data_term=float(smooth.class_membership_data_term),
                lesion_membership_data_term=float(smooth.lesion_membership_data_term),
                black_hole_membership_data_term=float(smooth.black_hole_membership_data_term),
                outmask_membership_data_term=0.0,
                class_energy=float(smooth.class_energy),
                lesion_energy=float(smooth.lesion_energy),
                black_hole_energy=float(smooth.black_hole_energy),
                has_lesion=True,
                lesion_class_index=1,
                lesion_class_type=target,
            )
            memberships.append(
                convert_lesion_toads_class_membership(
                    float(prior.class_membership_data_term),
                    float(shapes[index]),
                    relation_type=relations[index],
                )
            )
            energies.append(prior.class_energy)
            if not is_target:
                continue
            lesion = convert_lesion_toads_class_membership(
                float(prior.lesion_membership_data_term),
                float(shapes[index]),
                relation_type=relations[index],
            )
            relation_value = np.float32(1.0) if relations[index] == "direct" else _shape_factor(shapes[index])
            lesion_energy, black_hole_energy = prior.lesion_energy, prior.black_hole_energy
            distance_args = {
                "ventricle_distance": case["ventricle_distance"],
                "gray_matter_distance": case["gray_matter_distance"],
                "brainstem_distance": case["brainstem_distance"],
                "max_ventricle_distance": case["max_ventricle_distance"],
                "max_gray_matter_distance": case["max_gray_matter_distance"],
                "max_brainstem_distance": case["max_brainstem_distance"],
                "gm_flag": True,
                "spread": float(_from_bits(case["spread_bits"])),
            }
            if target == "LESION":
                modulated = modulate_lesion_toads_membership_by_distances(
                    float(lesion), float(relation_value), **distance_args
                )
                lesion, relation_value = modulated.membership, modulated.relation
            else:
                black_hole = convert_lesion_toads_class_membership(
                    float(prior.black_hole_membership_data_term),
                    float(shapes[index]),
                    relation_type=relations[index],
                )
                modulated = modulate_lesion_toads_blackhole_memberships_by_distances(
                    float(lesion), float(black_hole), float(relation_value), **distance_args
                )
                lesion, black_hole, relation_value = (
                    modulated.lesion_membership,
                    modulated.black_hole_membership,
                    modulated.relation,
                )

        normalized = normalize_lesion_toads_voxel_memberships(
            np.asarray(memberships),
            np.asarray(energies),
            types,
            lesion_membership=float(lesion),
            black_hole_membership=float(black_hole),
            outlier_membership=0.0,
            lesion_energy=float(lesion_energy),
            black_hole_energy=float(black_hole_energy),
        )
        assert [_bits(v) for v in normalized.class_memberships] == case["class_memberships_bits"], case
        assert [_bits(v) for v in energies] == case["class_energies_bits"], case
        assert _bits(normalized.lesion_membership) == case["lesion_membership_bits"], case
        assert _bits(normalized.black_hole_membership) == case["black_hole_membership_bits"], case
        assert _bits(relation_value) == case["lesion_relation_bits"], case
        assert _bits(normalized.denominator) == case["denominator_bits"], case
        assert _bits(normalized.energy) == case["local_energy_bits"], case
