from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np

from pymacruise.lesion_energy_ratio import lesion_toads_energy_ratio


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_energy_ratio"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "LesionEnergyRatioReferenceHarness.java"
HARNESS_SHA = "24abbf20c5ec3f4ab7a22e23898b623cd7e340de07e768ffb00d2fa0f9d3c106"
WRAPPER_SHA = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
NUMERICS_SHA = "13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb"


def _bits(value: float | np.float32) -> str:
    return f"{int(np.asarray(value, dtype=np.float32).view(np.uint32)):08x}"


def _float_from_bits(value: str) -> np.float32:
    return np.asarray(int(value, 16), dtype=np.uint32).view(np.float32)[()]


def test_energy_ratio_matches_selected_java_raw_bits(tmp_path: Path) -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA
    assert REFERENCE["wrapper_sha256"] == WRAPPER_SHA
    assert REFERENCE["numerics_sha256"] == NUMERICS_SHA
    assert REFERENCE["source_call_order"] == [
        "Eprev=energy",
        "energy=computeMemberships()",
        "Eratio=Numerics.abs(2.0f*(Eprev-energy)/(Eprev+energy))",
    ]

    synthetic_root = Path(__file__).parent / "data" / "synthetic"
    prior_reference = json.loads(
        (synthetic_root / REFERENCE["prior_objective_reference"]).read_text()
    )
    current_reference = json.loads(
        (synthetic_root / REFERENCE["current_objective_reference"]).read_text()
    )
    main_case = next(case for case in REFERENCE["cases"] if case["case"] == "p39ak_to_p39am_nan")
    assert main_case["previous_energy_bits"] == f"{prior_reference['objective_bits']:08x}"
    assert main_case["current_energy_bits"] == f"{current_reference['objective_bits']:08x}"

    javac = shutil.which("javac")
    java = shutil.which("java")
    assert javac is not None and java is not None, "a JDK is required for the pinned Java harness"
    subprocess.run([javac, "-d", str(tmp_path), str(HARNESS)], check=True, capture_output=True, text=True)
    completed = subprocess.run(
        [java, "-cp", str(tmp_path), "LesionEnergyRatioReferenceHarness"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert json.loads(completed.stdout) == REFERENCE

    for case in REFERENCE["cases"]:
        result = lesion_toads_energy_ratio(
            _float_from_bits(case["previous_energy_bits"]),
            _float_from_bits(case["current_energy_bits"]),
        )
        for field in (
            "previous_energy",
            "current_energy",
            "difference",
            "numerator",
            "denominator",
            "quotient",
            "energy_ratio",
        ):
            assert _bits(getattr(result, field)) == case[f"{field}_bits"]

    # NumPy's abs preserves a positive quiet-NaN sign here; Java's selected
    # Numerics.abs(float) falls through to unary negation and flips that sign.
    java_nan_ratio = next(case for case in REFERENCE["cases"] if case["case"] == "p39ak_to_p39am_nan")
    assert java_nan_ratio["energy_ratio_bits"] == "ffc00000"
    assert _bits(np.abs(_float_from_bits("7fc00000"))) == "7fc00000"
    assert next(case for case in REFERENCE["cases"] if case["case"] == "equal")["energy_ratio_bits"] == "80000000"
