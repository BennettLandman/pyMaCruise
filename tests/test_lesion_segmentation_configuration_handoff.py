from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_segmentation_configuration_handoff import (
    LesionToadsSegmentationConfiguration,
    configure_lesion_toads_segmentation,
)


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "lesion_segmentation_configuration_handoff"
HARNESS = ROOT / "java_reference" / "LesionToadsSegmentationConfigurationHandoffReferenceHarness.java"
REFERENCE = FIXTURE / "java_reference.json"
SOURCE_HASHES = {
    "segmentation": "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99",
    "wrapper": "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb",
    "preprocess": "72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08",
}


def _bits(values: np.ndarray) -> list[int]:
    return values.view(np.uint32).astype(np.uint64).tolist()


def test_configuration_handoff_matches_java_setters_and_order() -> None:
    ref = json.loads(REFERENCE.read_text())
    fixture = json.loads((FIXTURE / "input.json").read_text())
    assert ref["source_sha256"] == SOURCE_HASHES
    assert ref["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    assert ref["source_methods"] == [
        "setLesionWeight/reSetLesionWeight",
        "setIntensityMax",
        "setIntensityScale",
        "setGMflag",
    ]
    assert ref["call_order"] == [
        "lesion_weight_branch",
        "setIntensityMax",
        "setIntensityScale",
        "setGMflag",
    ]
    intensity_max = np.asarray(fixture["intensity_max_bits"], dtype=np.uint32).view(np.float32)
    intensity_scale = np.asarray(fixture["intensity_scale_bits"], dtype=np.uint32).view(np.float32)
    assert _bits(intensity_max) == ref["max_input_bits"]
    assert _bits(intensity_scale) == ref["scale_input_bits"]
    for lesion_weight in (False, True):
        scenario = ref["true_case" if lesion_weight else "false_case"]
        state = LesionToadsSegmentationConfiguration()
        configure_lesion_toads_segmentation(
            state,
            lesion_weight=lesion_weight,
            intensity_max=intensity_max,
            intensity_scale=intensity_scale,
        )
        assert state.use_lesion_weight is scenario["use_lesion_weight"] is lesion_weight
        assert state.gm_flag is scenario["gm_flag"] is True
        assert scenario["intensity_max_same_reference"] is True
        assert scenario["intensity_scale_same_reference"] is True
        assert state.intensity_max is intensity_max
        assert state.intensity_scale is intensity_scale
        assert _bits(state.intensity_max) == scenario["intensity_max_bits"]
        assert _bits(state.intensity_scale) == scenario["intensity_scale_bits"]
