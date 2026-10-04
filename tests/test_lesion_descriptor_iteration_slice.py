from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_atlas_priors import load_lesion_atlas_priors
from pymacruise.lesion_descriptor_iteration_slice import run_lesion_toads_descriptor_iteration_slice


ROOT = Path(__file__).parent
ATLAS = ROOT.parent / "data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"
FIXTURE = ROOT / "data/synthetic/lesion_descriptor_iteration_slice"
HARNESS = ROOT / "java_reference/LesionDescriptorIterationSliceReferenceHarness.java"
SOURCE = ROOT.parent / ".." / "toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java"
HARNESS_SHA = "ebae0608757c6af9113612019e56524ae1bf9cf0e3c824f3d3a78395ac87722a"
SOURCE_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"


def test_descriptor_slice_matches_java_active_data_to_membership_transition() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA
    reference = json.loads((FIXTURE / "java_reference.json").read_text())
    descriptor = load_lesion_atlas_priors(str(ATLAS), ("T1", "FLAIR"))
    topologies = tuple(item.topology for item in descriptor.structures)
    assert list(topologies) == reference["class_types"]
    classes = len(topologies)
    shape = (3, 3, 3)
    images = np.zeros((*shape, 2), dtype=np.float32)
    images[1, 1, 1] = (0.0, 2.0)
    centroids = np.arange(classes * 2, dtype=np.float32).reshape(classes, 2) / np.float32(8.0)
    memberships = np.full((*shape, classes), 0.05, dtype=np.float32)
    lesion_membership = np.full(shape, 0.05, dtype=np.float32)
    zero_map = np.zeros(shape, dtype=np.float32)
    class_shapes = np.ones((*shape, classes), dtype=np.float32)
    class_vector = np.zeros(classes, dtype=np.float32)

    result = run_lesion_toads_descriptor_iteration_slice(
        images,
        memberships,
        topologies,
        centroids,
        np.asarray((0.25, 0.75), dtype=np.float32),
        np.asarray((0.5, 1.25), dtype=np.float32),
        np.asarray((0.5, 1.5), dtype=np.float32),
        class_shapes,
        class_vector,
        class_vector,
        np.zeros((classes, classes), dtype=np.float32),
        class_vector,
        class_vector,
        ("direct",) * classes,
        outlier_distance=3.25,
        smoothing=0.0,
        prior_coefficient=0.0,
        initial_lesion_membership=lesion_membership,
        initial_black_hole_membership=zero_map,
        initial_outlier_membership=zero_map,
    )

    assert result.stage_order == (
        "computeActiveClassDataEnergies",
        "updateDescriptorMemberships",
        "normalizeMemberships",
    )
    assert result.data_energies.class_energy[1, 1, 1].view(np.uint32).tolist() == reference["class_energy_bits"]
    assert result.membership.class_memberships[1, 1, 1].view(np.uint32).tolist() == reference["class_membership_bits"]
    assert result.data_energies.lesion_energy[1, 1, 1].view(np.uint32) == reference["lesion_energy_bits"]
    assert result.membership.lesion_membership[1, 1, 1].view(np.uint32) == reference["lesion_membership_bits"]
    assert result.membership.denominator[1, 1, 1].view(np.uint32) == reference["denominator_bits"]
    np.testing.assert_array_equal(result.membership.class_memberships[0], memberships[0])
