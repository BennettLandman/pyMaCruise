from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_membership import (
    convert_lesion_toads_class_membership,
    modulate_lesion_toads_blackhole_memberships_by_distances,
    modulate_lesion_toads_membership_by_distances,
    normalize_lesion_toads_voxel_memberships,
    update_lesion_toads_smoothing_energies,
)


REFERENCE = (
    Path(__file__).parent
    / "data/synthetic/lesion_toads_smoothing_membership_integration/java_reference.json"
)


def _from_bits(value: int) -> np.float32:
    return np.asarray(value & 0xFFFFFFFF, dtype=np.uint32).view(np.float32)[()]


def _bits(value: float | np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def _shape_factor(shape: np.float32) -> np.float32:
    twice = np.float32(np.float32(2.0) * shape)
    return np.float32(1.0) if np.float32(1.0) < twice else twice


def test_smoothing_to_membership_composition_matches_selected_java() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert len(reference["cases"]) == 2

    for case in reference["cases"]:
        target_type = case["target_type"]
        types = ["OBJECT", target_type]
        relations = ["direct", case["relation_type"]]
        shapes = [_from_bits(v) for v in case["shape_bits"]]
        neighbor_classes = np.asarray(
            [[_from_bits(v) for v in face] for face in case["neighbor_class_bits"]]
        )
        neighbor_lesion = np.asarray([_from_bits(v) for v in case["neighbor_lesion_bits"]])
        neighbor_black_hole = np.asarray([_from_bits(v) for v in case["neighbor_black_hole_bits"]])
        class_data_terms = [_from_bits(v) for v in case["class_term_bits"]]
        class_energy_inputs = [_from_bits(v) for v in case["class_energy_input_bits"]]
        memberships: list[np.float32] = []
        class_energies: list[np.float32] = []
        lesion_membership = np.float32(0.0)
        black_hole_membership = np.float32(0.0)
        lesion_relation = np.float32(0.0)
        lesion_energy = np.float32(0.0)
        black_hole_energy = np.float32(0.0)

        for index, kind in enumerate(types):
            is_target = index == 1
            smoothed = update_lesion_toads_smoothing_energies(
                neighbor_classes,
                neighbor_lesion,
                neighbor_black_hole,
                class_index=index,
                class_type=kind,
                smoothing_prior=float(_from_bits(case["smoothing_prior_bits"])),
                smoothing=float(_from_bits(case["smoothing_bits"])),
                class_membership_data_term=float(class_data_terms[index]),
                lesion_membership_data_term=(float(_from_bits(case["lesion_term_bits"])) if is_target else 0.0),
                black_hole_membership_data_term=(float(_from_bits(case["black_hole_term_bits"])) if is_target else 0.0),
                outmask_membership_data_term=0.0,
                class_energy=float(class_energy_inputs[index]),
                lesion_energy=(float(_from_bits(case["lesion_energy_input_bits"])) if is_target else 0.0),
                black_hole_energy=(float(_from_bits(case["black_hole_energy_input_bits"])) if is_target else 0.0),
            )
            memberships.append(
                convert_lesion_toads_class_membership(
                    float(smoothed.class_membership_data_term),
                    float(shapes[index]),
                    relation_type=relations[index],
                )
            )
            class_energies.append(smoothed.class_energy)
            if not is_target:
                continue

            lesion_membership = convert_lesion_toads_class_membership(
                float(smoothed.lesion_membership_data_term),
                float(shapes[index]),
                relation_type=relations[index],
            )
            lesion_relation = (
                np.float32(1.0)
                if relations[index] == "direct"
                else _shape_factor(shapes[index])
            )
            lesion_energy = smoothed.lesion_energy
            black_hole_energy = smoothed.black_hole_energy
            distance_args = {
                "ventricle_distance": case["ventricle_distance"],
                "gray_matter_distance": case["gray_matter_distance"],
                "brainstem_distance": case["brainstem_distance"],
                "max_ventricle_distance": case["max_ventricle_distance"],
                "max_gray_matter_distance": case["max_gray_matter_distance"],
                "max_brainstem_distance": case["max_brainstem_distance"],
                "gm_flag": case["gm_flag"],
                "spread": float(_from_bits(case["spread_bits"])),
            }
            if target_type == "LESION":
                result = modulate_lesion_toads_membership_by_distances(
                    float(lesion_membership), float(lesion_relation), **distance_args
                )
                lesion_membership, lesion_relation = result.membership, result.relation
            else:
                black_hole_membership = convert_lesion_toads_class_membership(
                    float(smoothed.black_hole_membership_data_term),
                    float(shapes[index]),
                    relation_type=relations[index],
                )
                result = modulate_lesion_toads_blackhole_memberships_by_distances(
                    float(lesion_membership),
                    float(black_hole_membership),
                    float(lesion_relation),
                    **distance_args,
                )
                lesion_membership = result.lesion_membership
                black_hole_membership = result.black_hole_membership
                lesion_relation = result.relation

        normalized = normalize_lesion_toads_voxel_memberships(
            np.asarray(memberships),
            np.asarray(class_energies),
            types,
            lesion_membership=float(lesion_membership),
            black_hole_membership=float(black_hole_membership),
            outlier_membership=0.0,
            lesion_energy=float(lesion_energy),
            black_hole_energy=float(black_hole_energy),
        )
        assert [_bits(v) for v in normalized.class_memberships] == case["class_membership_bits"], case
        assert [_bits(v) for v in class_energies] == case["class_energy_bits"], case
        assert _bits(normalized.lesion_membership) == case["lesion_membership_bits"], case
        assert _bits(normalized.black_hole_membership) == case["black_hole_membership_bits"], case
        assert _bits(lesion_relation) == case["lesion_relation_bits"], case
        assert _bits(normalized.denominator) == case["denominator_bits"], case
        assert _bits(normalized.energy) == case["energy_bits"], case
