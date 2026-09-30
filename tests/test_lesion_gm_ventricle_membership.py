from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_gm_ventricle_membership import lesion_toads_gm_ventricle_memberships


FIXTURE = Path(__file__).parent / "data" / "synthetic" / "lesion_gm_ventricle_membership" / "java_reference.json"
HARNESS = Path(__file__).parent / "java_reference" / "LesionGmVentricleMembershipReferenceHarness.java"
HARNESS_SHA256 = "074407cd41bc58963bc0bd642496bbb847d5c476453e37fa1c293277646328a1"
FIXTURE_SHA256 = "1957985e281c5f9f7b9aa03d32c15e786ff9d8ea5da9ffede2d0a9a08c46fcf4"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"


def _floats_from_bits(values: list[int]) -> np.ndarray:
    return np.asarray(values, dtype=np.uint32).view(np.float32)


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def test_selected_gm_ventricle_membership_branch_matches_java_excerpt() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == FIXTURE_SHA256
    reference = json.loads(FIXTURE.read_text())
    assert reference["source_method"] == (
        "LesionToadSegmentation.computeMemberships ordinary GM class branch and normalization"
    )
    assert reference["source_sha256"] == SOURCE_SHA256

    for case in reference["cases"]:
        actual = lesion_toads_gm_ventricle_memberships(
            case["class_names"],
            _floats_from_bits(case["num_bits"]),
            _floats_from_bits(case["shape_bits"]),
            relation_type=case["relation"],
            has_lesion=case["has_lesion"],
            ventricle_distance=case["ventricle_distance"],
            max_ventricle_distance=case["max_ventricle_distance"],
        )
        assert _bits(actual) == case["normalized_membership_bits"], case["name"]
        total = np.float32(0.0)
        for value in actual:
            total = np.float32(total + value)
        assert int(np.asarray(total).view(np.uint32)) == case["normalized_sum_bits"], case["name"]

    cases = {case["name"]: case for case in reference["cases"]}
    assert cases["distance_zero"]["normalized_membership_bits"][0] == 0
    assert cases["zero_maximum_falls_through"]["normalized_membership_bits"] == (
        cases["no_lesion_falls_through"]["normalized_membership_bits"]
    )


def test_selected_gm_ventricle_membership_validates_inputs() -> None:
    valid = dict(
        relation_type="DIRECT", has_lesion=True, ventricle_distance=2, max_ventricle_distance=8
    )
    with np.testing.assert_raises_regex(ValueError, "one value per class"):
        lesion_toads_gm_ventricle_memberships(["GM"], np.asarray([0.5, 0.2]), np.asarray([0.5]), **valid)
    with np.testing.assert_raises_regex(ValueError, "DIRECT or PRIOR"):
        lesion_toads_gm_ventricle_memberships(
            ["GM"], np.asarray([0.5]), np.asarray([0.5]), **(valid | {"relation_type": "OTHER"})
        )
