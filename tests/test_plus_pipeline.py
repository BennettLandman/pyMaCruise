from __future__ import annotations

import json

import pytest

from pymacruise.plus_pipeline import PlusPipelineOptions


def test_plus_pipeline_options_match_selected_matlab_entrypoint() -> None:
    options = PlusPipelineOptions()
    assert options.as_dict() == {
        "lesion_toads_output": "dura removal inputs",
        "lesion_toads_correct": False,
        "wm_iso_level": 0.5,
        "gvf_lambda": 0.2,
        "update_enhance": 0,
        "generate_intermediate": True,
        "use_updated_csf_in_ace": 3,
        "boundary_erosion_method": "3D",
        "boundary_overlook_labels": [120, 121, 202, 203, 156, 157],
        "refine_topology_by_maseg": True,
    }


def test_plus_pipeline_options_serialize_as_json_provenance() -> None:
    options = PlusPipelineOptions()
    assert json.loads(json.dumps(options.as_dict())) == options.as_dict()


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("lesion_toads_output", "hard segmentation", "dura removal inputs"),
        ("wm_iso_level", 1.0, "strictly between"),
        ("gvf_lambda", -0.1, "nonnegative"),
        ("update_enhance", 2, "0 or 1"),
        ("use_updated_csf_in_ace", 4, "one of"),
        ("boundary_erosion_method", "2D", "unsupported"),
    ],
)
def test_plus_pipeline_options_reject_invalid_values(
    field: str, value: object, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        PlusPipelineOptions(**{field: value})
