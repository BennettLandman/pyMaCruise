from __future__ import annotations

import json
import struct
from pathlib import Path

import numpy as np

from pymacruise.toads_control import toads_iteration_control


REFERENCE = Path(__file__).parent / "data/synthetic/toads_iteration_control/java_reference.json"


def _float32_from_bits(value: int) -> np.float32:
    return np.float32(struct.unpack("!f", struct.pack("!I", value))[0])


def test_iteration_control_matches_selected_java_post_membership_flow() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    for case in reference["cases"]:
        result = toads_iteration_control(
            _float32_from_bits(case["previous_bits"]),
            _float32_from_bits(case["current_bits"]),
            case["iteration"],
            case["maximum_iterations"],
            _float32_from_bits(case["threshold_bits"]),
            do_growing=case["do_growing"],
        )
        assert int(np.asarray(result.energy_ratio).view(np.uint32)) == case["ratio_bits"]
        assert result.stop is case["stop"]
        assert result.run_growing is case["run_growing"]
        assert result.update_classification is case["update_classification"]
        assert result.refresh_relations is case["refresh_relations"]
