from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_atlas_priors import load_lesion_atlas_priors
from pymacruise.lesion_descriptor_membership import (
    map_lesion_descriptor_membership_types,
    run_lesion_toads_descriptor_membership_pass,
)


ATLAS = Path(__file__).parents[1] / "data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"
FIXTURE = Path(__file__).parent / "data/synthetic/lesion_descriptor_membership"
HARNESS = Path(__file__).parent / "java_reference/LesionDescriptorMembershipReferenceHarness.java"
SOURCE = Path(__file__).parents[1] / ".." / "toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java"
HARNESS_SHA = "bb2102a6e99edd0f94e9a15971b06e2e78bf827a62a8bb47a90566a352031d0f"
SOURCE_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"


def test_bundled_descriptor_drives_thirteen_class_membership_normalization() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA
    reference = json.loads((FIXTURE / "java_reference.json").read_text())
    assert reference["source_sha256"] == SOURCE_SHA
    descriptor = load_lesion_atlas_priors(str(ATLAS), ("T1", "FLAIR"))
    topologies = tuple(structure.topology for structure in descriptor.structures)
    assert list(topologies) == reference["class_types"]
    class_types = map_lesion_descriptor_membership_types(topologies)
    assert class_types.count("OBJECT") == 9
    assert class_types.count("MASK") == 3  # source MASK plus two ordinary optimized-CSF updates
    assert class_types.count("LESION") == 1

    shape = (3, 3, 3)
    classes = len(class_types)
    memberships = np.full((*shape, classes), 0.05, dtype=np.float32)
    lesion = np.full(shape, 0.05, dtype=np.float32)
    zero_map = np.zeros(shape, dtype=np.float32)
    class_data = np.ones((*shape, classes), dtype=np.float32)
    lesion_data = np.zeros(shape, dtype=np.float32)
    lesion_data[1, 1, 1] = 2.0
    shapes = np.ones((*shape, classes), dtype=np.float32)
    vector = np.zeros(classes, dtype=np.float32)
    matrix = np.zeros((classes, classes), dtype=np.float32)

    result = run_lesion_toads_descriptor_membership_pass(
        memberships,
        lesion,
        zero_map,
        zero_map,
        class_data,
        lesion_data,
        zero_map,
        zero_map,
        shapes,
        np.ones(classes, dtype=np.float32),
        descriptor_types=topologies,
        relation_types=("direct",) * classes,
        smoothing=0.0,
        shape_priors=vector,
        class_similarities=matrix,
        lesion_similarities=vector,
        black_hole_similarities=vector,
        prior_coefficient=0.0,
    )

    # Each class's raw term is 1, so each inverse term is 1. The LESION
    # auxiliary raw term is 2, contributing 0.5 to the shared denominator.
    np.testing.assert_allclose(result.class_memberships[1, 1, 1], 1.0 / 13.5, atol=1e-7)
    assert result.lesion_membership[1, 1, 1] == pytest.approx(0.5 / 13.5, abs=1e-7)
    assert result.denominator[1, 1, 1] == pytest.approx(13.5, abs=1e-6)
    assert result.class_memberships[1, 1, 1].view(np.uint32).tolist() == reference["class_membership_bits"]
    assert result.lesion_membership[1, 1, 1].view(np.uint32) == reference["lesion_membership_bits"]
    assert result.denominator[1, 1, 1].view(np.uint32) == reference["denominator_bits"]
    assert np.isfinite(result.objective)
    np.testing.assert_array_equal(result.class_memberships[0], memberships[0])


@pytest.mark.parametrize("bad", [("obj", "unknown", "lesion"), ("obj", "mask")])
def test_descriptor_membership_mapping_rejects_unknown_or_missing_lesion(bad: tuple[str, ...]) -> None:
    with pytest.raises(ValueError):
        map_lesion_descriptor_membership_types(bad)
