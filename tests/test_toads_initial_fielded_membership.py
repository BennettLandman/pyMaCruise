from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.toads_initial_fielded_membership import (
    toads_initial_default_atlas_membership_after_field_attachment,
)
from test_toads_default_atlas_dispatch import _bits, _fixture


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_initial_fielded_membership"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "ToadsInitialFieldedMembershipReferenceHarness.java"
HARNESS_SHA256 = "731a9e94b8d8c3aaa1e166d5385b8fe132cd55a9fceb92a4e1e003893ba94c26"


def test_initial_11_class_membership_follows_field_attachment_order() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    data, parameters, cases, dispatch_reference = _fixture()
    assert REFERENCE["source_sha256"] == {
        "MedicAlgorithmToads08.java": "d00a7ad609ab0c6381334248b53f4daadc070161af67ef853e2a552dd008d115",
        "AlgorithmMultipleToads.java": "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d",
        "MultipleToadSegmentation.java": "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884",
    }
    selected = next(
        case for case in cases
        if (case["field"], case["psf"], case["smoothing"], case["prior"])
        == (True, True, True, True)
    )
    assert REFERENCE["membership_reference"]["fixture"] == "toads_default_atlas_dispatch"
    assert dispatch_reference["toads_source_sha256"] == REFERENCE["source_sha256"]["MultipleToadSegmentation.java"]

    fields = [data["field"][..., channel] for channel in range(data["images"].shape[-1])]
    originals = [field.copy() for field in fields]
    result = toads_initial_default_atlas_membership_after_field_attachment(
        image_values=data["images"],
        initial_memberships=data["initial_memberships"],
        transformed_shapes=data["transformed_shapes"],
        closest_labels=data["closest_labels"],
        smoothing=.043,
        prior_coefficient=.16,
        **parameters,
        channel_fields=fields,
        psf_bias=data["psf_bias"],
    )

    assert result.field_state.use_field is True
    for channel, field in enumerate(fields):
        assert result.field_state.field_by_channel[channel] is field
        np.testing.assert_array_equal(field, originals[channel])
    assert tuple(stage.removeprefix("fixed:") for stage in result.stage_order) == tuple(
        REFERENCE["stage_order"]
    )
    assert _bits(result.membership_pass.class_data_energy) == selected["class_energy_bits"]
    assert _bits(result.membership_pass.memberships) == selected["membership_bits"]
    assert _bits(result.membership_pass.denominator) == selected["denominator_bits"]
    assert _bits(result.membership_pass.objective) == [selected["objective_bits"]]


def test_initial_fielded_membership_requires_a_field_for_each_channel() -> None:
    data, parameters, _, _ = _fixture()
    with np.testing.assert_raises_regex(ValueError, "one field per image channel"):
        toads_initial_default_atlas_membership_after_field_attachment(
            image_values=data["images"],
            initial_memberships=data["initial_memberships"],
            transformed_shapes=data["transformed_shapes"],
            closest_labels=data["closest_labels"],
            smoothing=.043,
            prior_coefficient=.16,
            channel_fields=[data["field"][..., 0]],
            psf_bias=data["psf_bias"],
            **parameters,
        )
