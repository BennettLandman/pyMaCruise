from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_membership import (
    modulate_lesion_toads_blackhole_memberships_by_distances,
    modulate_lesion_toads_membership_by_distances,
    normalize_lesion_toads_voxel_memberships,
    update_lesion_toads_outmask_membership,
)


REFERENCE = (
    Path(__file__).parent
    / "data/synthetic/lesion_toads_auxiliary_membership_integration/java_reference.json"
)


def _from_bits(value: int) -> np.float32:
    return np.asarray(value & 0xFFFFFFFF, dtype=np.uint32).view(np.float32)[()]


def _bits(value: float | np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def test_lesion_blackhole_and_outmask_auxiliaries_compose_like_java() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert len(reference["cases"]) == 2
    for case in reference["cases"]:
        types = ["OBJECT", case["auxiliary_type"], "OUTMASK"]
        memberships = np.asarray([_from_bits(v) for v in case["class_membership_bits"]])
        energies = np.asarray([_from_bits(v) for v in case["class_energy_bits"]])
        lesion = _from_bits(case["lesion_bits"])
        black_hole = _from_bits(case["black_hole_bits"])
        relation = _from_bits(case["lesion_relation_bits"])
        outlier = _from_bits(case["outlier_bits"])
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
        if case["auxiliary_type"] == "LESION":
            result = modulate_lesion_toads_membership_by_distances(
                float(lesion), float(relation), **distance_args
            )
            lesion, relation = result.membership, result.relation
        else:
            result = modulate_lesion_toads_blackhole_memberships_by_distances(
                float(lesion),
                float(black_hole),
                float(relation),
                **distance_args,
            )
            lesion = result.lesion_membership
            black_hole = result.black_hole_membership
            relation = result.relation

        outmask = update_lesion_toads_outmask_membership(
            float(memberships[2]),
            float(_from_bits(case["outmask_value_bits"])),
            float(outlier),
            relation_type=case["outmask_relation"],
            shape_prior=float(_from_bits(case["outmask_shape_bits"])),
        )
        memberships[2] = outmask.class_membership
        outlier = outmask.outlier_membership
        normalized = normalize_lesion_toads_voxel_memberships(
            memberships,
            energies,
            types,
            lesion_membership=float(lesion),
            black_hole_membership=float(black_hole),
            outlier_membership=float(outlier),
            lesion_energy=float(_from_bits(case["lesion_energy_bits"])),
            black_hole_energy=float(_from_bits(case["black_hole_energy_bits"])),
        )
        assert [_bits(v) for v in normalized.class_memberships] == case["class_memberships_output_bits"], case
        assert _bits(normalized.lesion_membership) == case["lesion_output_bits"], case
        assert _bits(normalized.black_hole_membership) == case["black_hole_output_bits"], case
        assert _bits(relation) == case["lesion_relation_output_bits"], case
        assert _bits(normalized.outlier_membership) == case["outlier_output_bits"], case
        assert _bits(normalized.denominator) == case["denominator_bits"], case
        assert _bits(normalized.energy) == case["local_energy_bits"], case
