from __future__ import annotations

import json
import hashlib
from pathlib import Path

import numpy as np

from pymacruise.toads_math import toads_powertable_lookup


REFERENCE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_powertable"
    / "java_reference.json"
)
POWER_TABLE_SHA256 = "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"
JAVA_SOURCE = (
    Path(__file__).parent
    / "java_reference"
    / "upstream"
    / "edu"
    / "jhmi"
    / "rad"
    / "medic"
    / "utilities"
    / "PowerTable.java"
)


def _normalize_source_inverse_energy(power: np.ndarray) -> np.ndarray:
    energy = np.asarray(power, dtype=np.float32)
    inverse = np.full(3, np.float32(1e10), dtype=np.float32)
    np.divide(1.0, energy, out=inverse, where=energy > np.float32(1e-10))
    denominator = np.float32(0.0)
    for value in inverse:
        denominator = np.float32(denominator + value)
    return np.asarray(inverse / denominator, dtype=np.float32)


def test_numpy_power_substitution_error_against_selected_java_table() -> None:
    fixture = json.loads(REFERENCE.read_text())
    assert hashlib.sha256(JAVA_SOURCE.read_bytes()).hexdigest() == POWER_TABLE_SHA256
    assert fixture["source_sha256"] == POWER_TABLE_SHA256
    assert len(fixture["cases"]) == 1042
    assert fixture["input_range"] == [0, 4]
    assert fixture["step"] == float(np.float32(1e-6))
    assert fixture["seed"] == 1804

    exponents = np.asarray(fixture["exponents"], dtype=np.float32)
    direct_power_error = []
    resolved_power_error = []
    direct_membership_error = []
    compatible_membership_error = []
    direct_winner_mismatches = 0
    compatible_winner_mismatches = 0
    compatible_worst = None
    max_relative_power_error = 0.0
    for case in fixture["cases"]:
        p = np.asarray(case["p"], dtype=np.float32)
        java_power = np.asarray(case["power"], dtype=np.float64)
        direct_power = np.power(p.astype(np.float64), exponents.astype(np.float64))
        compatible_power = np.asarray(
            [
                toads_powertable_lookup(
                    p[k], float(exponents[k]), cached_exponent=float(exponents[k])
                ).item()
                for k in range(3)
            ],
            dtype=np.float64,
        )
        direct_power_error.extend(np.abs(direct_power - java_power))
        nonzero = java_power != 0.0
        max_relative_power_error = max(
            max_relative_power_error,
            float(
                np.max(
                    np.abs((direct_power[nonzero] - java_power[nonzero]) / java_power[nonzero])
                )
            ),
        )
        resolved = p >= np.float32(0.01)
        if np.any(resolved):
            resolved_power_error.extend(np.abs(direct_power[resolved] - java_power[resolved]))

        java_membership = np.asarray(case["membership"], dtype=np.float32)
        direct_membership = _normalize_source_inverse_energy(direct_power)
        compatible_membership = _normalize_source_inverse_energy(compatible_power)
        direct_membership_error.extend(np.abs(java_membership - direct_membership))
        compatible_membership_error.extend(np.abs(java_membership - compatible_membership))
        if compatible_worst is None or max(np.abs(java_membership - compatible_membership)) > compatible_worst[0]:
            compatible_worst = (
                float(max(np.abs(java_membership - compatible_membership))),
                case,
                compatible_power.tolist(),
                compatible_membership.tolist(),
            )
        direct_winner_mismatches += int(
            np.argmax(java_membership) != np.argmax(direct_membership)
        )
        compatible_winner_mismatches += int(
            np.argmax(java_membership) != np.argmax(compatible_membership)
        )

    assert np.isfinite(direct_power_error).all()
    assert np.isfinite(direct_membership_error).all()
    assert np.isfinite(compatible_membership_error).all()
    # Unquantized NumPy power is close as a scalar, but lookup rounding near
    # zero materially changes normalized class memberships and three winners.
    assert max(direct_power_error) < 1e-5
    assert max(resolved_power_error) < 1e-6
    assert np.isclose(max_relative_power_error, 0.4256508272473058, rtol=0, atol=1e-12)
    assert np.isclose(max(direct_membership_error), 0.484713018, rtol=0, atol=1e-7)
    assert direct_winner_mismatches == 3
    # NumPy power with the source's table-index rounding matches Java after
    # the source's float32 inverse-energy and normalization steps.
    assert max(compatible_membership_error) < 1e-6, f"worst quantized case: {compatible_worst}"
    assert compatible_winner_mismatches == 0


def test_powertable_lookup_uses_direct_power_for_cache_misses_and_out_of_range_inputs() -> None:
    q = float(np.float32(0.8))
    cache_miss = np.asarray([0.25, 0.3, 4.25], dtype=np.float64)
    np.testing.assert_array_equal(
        toads_powertable_lookup(cache_miss, q + 0.01, cached_exponent=q),
        np.power(cache_miss, q + 0.01),
    )
    in_range = np.asarray([0.2, 1.1], dtype=np.float64)
    out_of_range = toads_powertable_lookup(in_range, q, cached_exponent=q, minimum=0.1, maximum=1.0)
    assert out_of_range[1] == np.power(in_range[1], q)
