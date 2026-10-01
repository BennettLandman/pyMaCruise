from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_initial_loop_entry import seed_lesion_toads_initial_loop
from pymacruise.lesion_standard13_volume import run_lesion_toads_standard13_volume
from test_lesion_standard13_volume import _state


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_initial_loop_entry"
HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsInitialLoopEntryReferenceHarness.java"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
WRAPPER_SHA = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
HARNESS_SHA = "feeb6ab95df14d6f1839fbf42f657940b6db4ba0d0607467cb501895c7843f13"


def _bits(value: float | np.float32) -> int:
    return int(np.asarray(value, dtype=np.float32).view(np.uint32))


def test_initial_standard13_memberships_seed_selected_main_loop_gate() -> None:
    state = _state()
    initial = run_lesion_toads_standard13_volume(**state)
    volume_reference = json.loads(
        (
            Path(__file__).parent
            / "data"
            / "synthetic"
            / "lesion_standard13_volume"
            / "java_reference.json"
        ).read_text()
    )
    assert _bits(initial.objective) == volume_reference["objective_bits"]
    assert REFERENCE["wrapper_sha256"] == WRAPPER_SHA
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA
    assert REFERENCE["source_call_order"] == [
        "computeMemberships",
        "setGMflag",
        "Eprev=energy",
        "Eratio=1.0f",
        "stop=false",
        "iter=0",
        "if(iterations==0) stop=true",
        "while(!stop && !threadStopped && isWorking())",
        "iter++",
        "propagateThinning",
    ]

    maps_before = tuple(
        value.copy()
        for value in (
            initial.class_memberships,
            initial.lesion_membership,
            initial.black_hole_membership,
            initial.outlier_membership,
            initial.lesion_relation,
        )
    )
    for case in REFERENCE["cases"]:
        result = seed_lesion_toads_initial_loop(
            initial,
            iterations=case["iterations"],
            thread_stopped=case["thread_stopped"],
            is_working=case["is_working"],
        )
        assert result.initial_membership is initial
        assert _bits(result.previous_energy) == case["previous_energy_bits"]
        assert _bits(result.energy_ratio) == case["energy_ratio_bits"]
        assert result.gm_flag is case["gm_flag_after_initial_memberships"]
        assert result.iteration == case["iteration"]
        assert result.stop is case["stop"]
        assert result.enter_loop is case["enter_loop"]
        assert result.next_stage == case["next_stage"]
        assert list(result.stage_order[:5]) == [
            "computeMemberships_initial",
            "setGMflag",
            "Eprev=initial_energy",
            "Eratio=1.0f",
            "initialize_stop_and_iteration",
        ]
        assert ("propagateThinning_pending" in result.stage_order) is result.enter_loop

    for actual, before in zip(
        (
            initial.class_memberships,
            initial.lesion_membership,
            initial.black_hole_membership,
            initial.outlier_membership,
            initial.lesion_relation,
        ),
        maps_before,
    ):
        np.testing.assert_array_equal(actual, before)
