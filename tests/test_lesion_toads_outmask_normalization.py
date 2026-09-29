from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_membership import (
    normalize_lesion_toads_voxel_memberships,
    update_lesion_toads_outmask_membership,
)


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_outmask_normalization/java_reference.json"


def _from_bits(value: int) -> np.float32:
    return np.asarray(value & 0xFFFFFFFF, dtype=np.uint32).view(np.float32)[()]


def _bits(value: float | np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def test_outmask_mixture_composes_with_shared_normalization() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert len(reference["cases"]) == 3
    for case in reference["cases"]:
        memberships = [_from_bits(v) for v in case["class_values_bits"]]
        outlier = _from_bits(case["initial_outlier_bits"])
        outmask_values = [_from_bits(v) for v in case["outmask_values_bits"]]
        shapes = [_from_bits(v) for v in case["shape_bits"]]
        for index, kind in enumerate(case["types"]):
            if kind == "OUTMASK":
                result = update_lesion_toads_outmask_membership(
                    float(memberships[index]),
                    float(outmask_values[index]),
                    float(outlier),
                    relation_type=case["relations"][index],
                    shape_prior=float(shapes[index]),
                )
                memberships[index] = result.class_membership
                outlier = result.outlier_membership
        normalized = normalize_lesion_toads_voxel_memberships(
            np.asarray(memberships),
            np.asarray([_from_bits(v) for v in case["energy_bits"]]),
            case["types"],
            lesion_membership=0.0,
            black_hole_membership=0.0,
            outlier_membership=float(outlier),
            lesion_energy=0.0,
            black_hole_energy=0.0,
        )
        assert [_bits(v) for v in normalized.class_memberships] == case["class_memberships_bits"], case
        assert _bits(normalized.outlier_membership) == case["outlier_membership_bits"], case
        assert _bits(normalized.denominator) == case["denominator_bits"], case
        assert _bits(normalized.energy) == case["local_energy_bits"], case
