from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_membership import (
    convert_lesion_toads_class_membership,
    modulate_lesion_toads_blackhole_memberships_by_distances,
    modulate_lesion_toads_membership_by_distances,
    normalize_lesion_toads_voxel_memberships,
    update_lesion_toads_prior_energies,
    update_lesion_toads_smoothing_energies,
)


REFERENCE = (
    Path(__file__).parent
    / "data/synthetic/lesion_toads_smoothing_prior_membership/java_reference.json"
)


def _from_bits(value: int) -> np.float32:
    return np.asarray(value & 0xFFFFFFFF, dtype=np.uint32).view(np.float32)[()]


def _bits(value: float | np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def _shape_factor(shape: np.float32) -> np.float32:
    twice = np.float32(np.float32(2.0) * shape)
    return np.float32(1.0) if np.float32(1.0) < twice else twice


def test_smoothing_prior_membership_composition_matches_selected_java() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert len(reference["cases"]) == 2

    for case in reference["cases"]:
        target = case["target_type"]
        types = ["OBJECT", target]
        relations = ["direct", case["relation_type"]]
        shapes = [_from_bits(v) for v in case["shape_values_bits"]]
        priors = np.asarray([_from_bits(v) for v in case["shape_priors_bits"]])
        similarities = np.asarray(
            [[_from_bits(v) for v in row] for row in case["similarities_bits"]]
        )
        lesion_sim = np.asarray([_from_bits(v) for v in case["lesion_similarities_bits"]])
        black_hole_sim = np.asarray([_from_bits(v) for v in case["black_hole_similarities_bits"]])
        neighbors = np.asarray([[_from_bits(v) for v in row] for row in case["neighbor_class_bits"]])
        neighbor_lesion = np.asarray([_from_bits(v) for v in case["neighbor_lesion_bits"]])
        neighbor_black_hole = np.asarray([_from_bits(v) for v in case["neighbor_black_hole_bits"]])
        class_terms = [_from_bits(v) for v in case["class_terms_bits"]]
        class_energy_inputs = [_from_bits(v) for v in case["class_energies_bits"]]
        class_memberships: list[np.float32] = []
        class_energies: list[np.float32] = []
        lesion = np.float32(0.0)
        black_hole = np.float32(0.0)
        lesion_relation = np.float32(0.0)
        lesion_energy = np.float32(0.0)
        black_hole_energy = np.float32(0.0)

        for index, kind in enumerate(types):
            is_target = index == 1
            smoothed = update_lesion_toads_smoothing_energies(
                neighbors,
                neighbor_lesion,
                neighbor_black_hole,
                class_index=index,
                class_type=kind,
                smoothing_prior=float(_from_bits(case["prior_bits"])),
                smoothing=float(_from_bits(case["smoothing_bits"])),
                class_membership_data_term=float(class_terms[index]),
                lesion_membership_data_term=(float(_from_bits(case["lesion_term_bits"])) if is_target else 0.0),
                black_hole_membership_data_term=(float(_from_bits(case["black_hole_term_bits"])) if is_target else 0.0),
                outmask_membership_data_term=0.0,
                class_energy=float(class_energy_inputs[index]),
                lesion_energy=(float(_from_bits(case["lesion_energy_input_bits"])) if is_target else 0.0),
                black_hole_energy=(float(_from_bits(case["black_hole_energy_input_bits"])) if is_target else 0.0),
            )
            prior_updated = update_lesion_toads_prior_energies(
                priors,
                similarities[index],
                lesion_sim,
                black_hole_sim,
                np.asarray(shapes),
                class_index=index,
                class_type=kind,
                prior_coefficient=float(_from_bits(case["prior_coefficient_bits"])),
                class_membership_data_term=float(smoothed.class_membership_data_term),
                lesion_membership_data_term=float(smoothed.lesion_membership_data_term),
                black_hole_membership_data_term=float(smoothed.black_hole_membership_data_term),
                outmask_membership_data_term=0.0,
                class_energy=float(smoothed.class_energy),
                lesion_energy=float(smoothed.lesion_energy),
                black_hole_energy=float(smoothed.black_hole_energy),
                has_lesion=case["has_lesion"],
                lesion_class_index=1,
                lesion_class_type=target,
            )
            class_memberships.append(
                convert_lesion_toads_class_membership(
                    float(prior_updated.class_membership_data_term),
                    float(shapes[index]),
                    relation_type=relations[index],
                )
            )
            class_energies.append(prior_updated.class_energy)
            if not is_target:
                continue

            lesion = convert_lesion_toads_class_membership(
                float(prior_updated.lesion_membership_data_term),
                float(shapes[index]),
                relation_type=relations[index],
            )
            lesion_relation = (
                np.float32(1.0) if relations[index] == "direct" else _shape_factor(shapes[index])
            )
            lesion_energy = prior_updated.lesion_energy
            black_hole_energy = prior_updated.black_hole_energy
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
                    float(lesion), float(lesion_relation), **distance_args
                )
                lesion, lesion_relation = modulated.membership, modulated.relation
            else:
                black_hole = convert_lesion_toads_class_membership(
                    float(prior_updated.black_hole_membership_data_term),
                    float(shapes[index]),
                    relation_type=relations[index],
                )
                modulated = modulate_lesion_toads_blackhole_memberships_by_distances(
                    float(lesion), float(black_hole), float(lesion_relation), **distance_args
                )
                lesion = modulated.lesion_membership
                black_hole = modulated.black_hole_membership
                lesion_relation = modulated.relation

        normalized = normalize_lesion_toads_voxel_memberships(
            np.asarray(class_memberships),
            np.asarray(class_energies),
            types,
            lesion_membership=float(lesion),
            black_hole_membership=float(black_hole),
            outlier_membership=0.0,
            lesion_energy=float(lesion_energy),
            black_hole_energy=float(black_hole_energy),
        )
        assert [_bits(v) for v in normalized.class_memberships] == case["class_memberships_bits"], case
        assert [_bits(v) for v in class_energies] == case["class_energy_bits"], case
        assert _bits(normalized.lesion_membership) == case["lesion_membership_bits"], case
        assert _bits(normalized.black_hole_membership) == case["black_hole_membership_bits"], case
        assert _bits(lesion_relation) == case["lesion_relation_bits"], case
        assert _bits(normalized.denominator) == case["denominator_bits"], case
        assert _bits(normalized.energy) == case["energy_bits"], case
