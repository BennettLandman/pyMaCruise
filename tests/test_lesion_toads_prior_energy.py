from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_membership import update_lesion_toads_prior_energies


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_prior_energy/java_reference.json"


def _from_bits(value: int) -> np.float32:
    return np.asarray(value & 0xFFFFFFFF, dtype=np.uint32).view(np.float32)[()]


def _bits(value: float | np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def test_prior_energy_terms_match_selected_java_float32() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert len(reference["cases"]) == 6

    for case in reference["cases"]:
        result = update_lesion_toads_prior_energies(
            np.asarray([_from_bits(v) for v in case["shape_priors_bits"]]),
            np.asarray([_from_bits(v) for v in case["class_similarities_bits"]]),
            np.asarray([_from_bits(v) for v in case["lesion_similarities_bits"]]),
            np.asarray([_from_bits(v) for v in case["black_hole_similarities_bits"]]),
            np.asarray([_from_bits(v) for v in case["shape_values_bits"]]),
            class_index=case["class_index"],
            class_type=case["class_type"],
            prior_coefficient=float(_from_bits(case["prior_coefficient_bits"])),
            class_membership_data_term=float(_from_bits(case["initial_class_term_bits"])),
            lesion_membership_data_term=float(_from_bits(case["initial_lesion_term_bits"])),
            black_hole_membership_data_term=float(_from_bits(case["initial_black_hole_term_bits"])),
            outmask_membership_data_term=float(_from_bits(case["initial_outmask_term_bits"])),
            class_energy=float(_from_bits(case["initial_class_energy_bits"])),
            lesion_energy=float(_from_bits(case["initial_lesion_energy_bits"])),
            black_hole_energy=float(_from_bits(case["initial_black_hole_energy_bits"])),
            has_lesion=case["has_lesion"],
            lesion_class_index=case["lesion_class_index"],
            lesion_class_type=case["lesion_class_type"],
        )
        assert _bits(result.class_prior) == case["class_prior_bits"], case
        assert _bits(result.lesion_prior) == case["lesion_prior_bits"], case
        assert _bits(result.black_hole_prior) == case["black_hole_prior_bits"], case
        assert _bits(result.class_membership_data_term) == case["class_term_bits"], case
        assert _bits(result.lesion_membership_data_term) == case["lesion_term_bits"], case
        assert _bits(result.black_hole_membership_data_term) == case["black_hole_term_bits"], case
        assert _bits(result.outmask_membership_data_term) == case["outmask_term_bits"], case
        assert _bits(result.class_energy) == case["class_energy_bits"], case
        assert _bits(result.lesion_energy) == case["lesion_energy_bits"], case
        assert _bits(result.black_hole_energy) == case["black_hole_energy_bits"], case
