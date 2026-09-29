from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_membership import (
    convert_lesion_toads_class_membership,
    modulate_lesion_toads_blackhole_memberships_by_distances,
    modulate_lesion_toads_membership_by_distances,
    normalize_lesion_toads_voxel_memberships,
)


REFERENCE = (
    Path(__file__).parent
    / "data/synthetic/lesion_toads_membership_integration/java_reference.json"
)


def _from_bits(value: int) -> np.float32:
    return np.asarray(value & 0xFFFFFFFF, dtype=np.uint32).view(np.float32)[()]


def _bits(value: float | np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def _shape_factor(shape: np.float32) -> np.float32:
    twice = np.float32(np.float32(2.0) * shape)
    return np.float32(1.0) if np.float32(1.0) < twice else twice


def test_selected_membership_kernels_compose_like_java_voxel_formula() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert len(reference["cases"]) == 5

    for case in reference["cases"]:
        types = case["types"]
        relations = case["relations"]
        shapes = [_from_bits(v) for v in case["shape_bits"]]
        class_energies = [_from_bits(v) for v in case["class_energy_bits"]]
        lesion_energies = [_from_bits(v) for v in case["lesion_energy_bits"]]
        black_hole_energies = [_from_bits(v) for v in case["black_hole_energy_bits"]]
        class_memberships = []
        lesion_membership = np.float32(0.0)
        black_hole_membership = np.float32(0.0)
        lesion_relation = np.float32(0.0)
        lesion_energy = np.float32(0.0)
        black_hole_energy = np.float32(0.0)

        for index, kind in enumerate(types):
            relation = relations[index]
            shape = shapes[index]
            class_memberships.append(
                convert_lesion_toads_class_membership(
                    float(class_energies[index]),
                    float(shape),
                    relation_type=relation,
                )
            )
            if kind in {"LESION", "OPTIMIZED_LESION"}:
                lesion_membership = convert_lesion_toads_class_membership(
                    float(lesion_energies[index]),
                    float(shape),
                    relation_type=relation,
                )
                lesion_relation = (
                    np.float32(1.0) if relation == "direct" else _shape_factor(shape)
                )
                lesion_energy = lesion_energies[index]
                result = modulate_lesion_toads_membership_by_distances(
                    float(lesion_membership),
                    float(lesion_relation),
                    ventricle_distance=case["ventricle_distance"],
                    gray_matter_distance=case["gray_matter_distance"],
                    brainstem_distance=case["brainstem_distance"],
                    max_ventricle_distance=case["max_ventricle_distance"],
                    max_gray_matter_distance=case["max_gray_matter_distance"],
                    max_brainstem_distance=case["max_brainstem_distance"],
                    gm_flag=case["gm_flag"],
                    spread=float(_from_bits(case["spread_bits"])),
                )
                lesion_membership, lesion_relation = result.membership, result.relation
            elif kind == "BLACKHOLE":
                lesion_membership = convert_lesion_toads_class_membership(
                    float(lesion_energies[index]), float(shape), relation_type=relation
                )
                black_hole_membership = convert_lesion_toads_class_membership(
                    float(black_hole_energies[index]), float(shape), relation_type=relation
                )
                lesion_relation = (
                    np.float32(1.0) if relation == "direct" else _shape_factor(shape)
                )
                lesion_energy = lesion_energies[index]
                black_hole_energy = black_hole_energies[index]
                result = modulate_lesion_toads_blackhole_memberships_by_distances(
                    float(lesion_membership),
                    float(black_hole_membership),
                    float(lesion_relation),
                    ventricle_distance=case["ventricle_distance"],
                    gray_matter_distance=case["gray_matter_distance"],
                    brainstem_distance=case["brainstem_distance"],
                    max_ventricle_distance=case["max_ventricle_distance"],
                    max_gray_matter_distance=case["max_gray_matter_distance"],
                    max_brainstem_distance=case["max_brainstem_distance"],
                    gm_flag=case["gm_flag"],
                    spread=float(_from_bits(case["spread_bits"])),
                )
                lesion_membership = result.lesion_membership
                black_hole_membership = result.black_hole_membership
                lesion_relation = result.relation

        normalized = normalize_lesion_toads_voxel_memberships(
            np.asarray(class_memberships),
            np.asarray(class_energies),
            types,
            lesion_membership=float(lesion_membership),
            black_hole_membership=float(black_hole_membership),
            outlier_membership=0.0,
            lesion_energy=float(lesion_energy),
            black_hole_energy=float(black_hole_energy),
        )
        assert [_bits(v) for v in normalized.class_memberships] == case["class_membership_bits"], case
        assert _bits(normalized.lesion_membership) == case["lesion_membership_bits"], case
        assert _bits(normalized.black_hole_membership) == case["black_hole_membership_bits"], case
        assert _bits(lesion_relation) == case["lesion_relation_bits"], case
        assert _bits(normalized.denominator) == case["denominator_bits"], case
        assert _bits(normalized.energy) == case["energy_bits"], case
