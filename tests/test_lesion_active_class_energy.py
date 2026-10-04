from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_active_class_energy import compute_lesion_toads_active_class_energies
from pymacruise.lesion_atlas_priors import load_lesion_atlas_priors
from pymacruise.lesion_data_terms import compute_lesion_toads_unfielded_data_energies


ATLAS = Path(__file__).parents[1] / "data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"
FIXTURE = Path(__file__).parent / "data/synthetic/lesion_active_class_energy"
HARNESS = Path(__file__).parent / "java_reference/LesionActiveClassEnergyReferenceHarness.java"
SOURCE = Path(__file__).parents[1] / ".." / "toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java"
HARNESS_SHA = "95eac136b52ef9c6bc00670455c6563c14a9fd48010a60d99a4ac6577d3900c9"
SOURCE_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"


def test_bundled_thirteen_class_active_dispatch_uses_source_fallback_terms() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA
    reference = json.loads((FIXTURE / "java_reference.json").read_text())
    assert reference["source_sha256"] == SOURCE_SHA
    descriptor = load_lesion_atlas_priors(str(ATLAS), ("T1", "FLAIR"))
    class_types = tuple(structure.topology for structure in descriptor.structures)
    assert list(class_types) == reference["class_types"]
    assert len(class_types) == 13
    assert class_types.count("obj") == 9
    assert class_types.count("mask") == 1
    assert class_types.count("optimized_csf") == 2
    assert class_types.count("lesion") == 1

    images = np.zeros((5, 4, 3, 2), dtype=np.float32)
    images[2, 2, 1] = (0.0, 2.0)
    centroids = np.arange(26, dtype=np.float32).reshape(13, 2) / np.float32(8)
    weights = np.asarray((0.5, 1.5), dtype=np.float32)
    lesion_centroid = np.asarray((0.25, 0.75), dtype=np.float32)
    black_hole_centroid = np.asarray((0.5, 1.25), dtype=np.float32)
    outlier = np.float32(3.25)

    result = compute_lesion_toads_active_class_energies(
        images,
        centroids,
        lesion_centroid,
        black_hole_centroid,
        weights,
        class_types,
        outlier_distance=float(outlier),
    )
    assert result.class_energy[2, 2, 1].view(np.uint32).tolist() == reference["class_energy_bits"]
    assert int(result.lesion_energy[2, 2, 1].view(np.uint32)) == reference["lesion_energy_bits"]
    assert int(result.black_hole_energy[2, 2, 1].view(np.uint32)) == reference["black_hole_energy_bits"]

    # Java computeMemberships currently falls through from optimized CSF to
    # `num = outlier`; it does not run the old commented PowerTable block.
    optimized = [i for i, kind in enumerate(class_types) if kind == "optimized_csf"]
    np.testing.assert_array_equal(result.class_energy[2, 2, 1, optimized], outlier)
    mask_index = class_types.index("mask")
    np.testing.assert_array_equal(result.class_energy[2, 2, 1, mask_index], weights[1])

    object_index = class_types.index("obj")
    expected_object = np.float32(0.0)
    for channel in range(2):
        difference = np.float32(images[2, 2, 1, channel] - centroids[object_index, channel])
        expected_object = np.float32(
            expected_object + np.float32(np.float32(weights[channel] * difference) * difference)
        )
    assert result.class_energy[2, 2, 1, object_index] == expected_object

    lesion_index = class_types.index("lesion")
    expected_lesion = np.float32(0.0)
    expected_lesion_data = np.float32(0.0)
    for channel in range(2):
        difference = np.float32(images[2, 2, 1, channel] - centroids[lesion_index, channel])
        expected_lesion = np.float32(
            expected_lesion + np.float32(np.float32(weights[channel] * difference) * difference)
        )
        difference = np.float32(images[2, 2, 1, channel] - lesion_centroid[channel])
        expected_lesion_data = np.float32(
            expected_lesion_data + np.float32(np.float32(weights[channel] * difference) * difference)
        )
    assert result.class_energy[2, 2, 1, lesion_index] == expected_lesion
    assert result.lesion_energy[2, 2, 1] == expected_lesion_data
    assert not np.any(result.black_hole_energy)

    boundary = np.ones(images.shape[:3], dtype=bool)
    boundary[1:-1, 1:-1, 1:-1] = False
    np.testing.assert_array_equal(result.class_energy[boundary], 0.0)
    np.testing.assert_array_equal(result.lesion_energy[boundary], 0.0)


def test_active_class_energy_rejects_unknown_source_type_and_missing_field() -> None:
    args = (
        np.ones((3, 3, 3, 1), dtype=np.float32),
        np.ones((1, 1), dtype=np.float32),
        np.ones(1, dtype=np.float32),
        np.ones(1, dtype=np.float32),
        np.ones(1, dtype=np.float32),
    )
    with pytest.raises(ValueError, match="not handled"):
        compute_lesion_toads_active_class_energies(
            *args, ("new_source_type",), outlier_distance=1.0
        )
    with pytest.raises(ValueError, match="required"):
        compute_lesion_toads_active_class_energies(
            *args, ("object",), outlier_distance=1.0, use_field=True
        )


def test_active_object_branch_multiplies_by_the_supplied_intensity_field() -> None:
    images = np.full((3, 3, 3, 2), 2.0, dtype=np.float32)
    fields = np.full_like(images, 0.5)
    result = compute_lesion_toads_active_class_energies(
        images,
        np.zeros((1, 2), dtype=np.float32),
        np.zeros(2, dtype=np.float32),
        np.zeros(2, dtype=np.float32),
        np.ones(2, dtype=np.float32),
        ("obj",),
        outlier_distance=7.0,
        use_field=True,
        intensity_fields=fields,
    )
    # field * image = 1 in both channels, so the two source residual terms sum to 2.
    assert result.class_energy[1, 1, 1, 0] == np.float32(2.0)


def test_active_lesion_branch_uses_the_selected_flair_power_table_term() -> None:
    image = np.asarray((0.3, 1.8), dtype=np.float32)
    class_centroid = np.asarray((0.5, 0.9), dtype=np.float32)
    lesion_centroid = np.asarray((0.2, 0.4), dtype=np.float32)
    black_hole_centroid = np.zeros(2, dtype=np.float32)
    weights = np.asarray((0.5, 1.5), dtype=np.float32)
    exponent = np.asarray((1.0, 0.5), dtype=np.float32)
    expected = compute_lesion_toads_unfielded_data_energies(
        image,
        class_centroid,
        lesion_centroid,
        black_hole_centroid,
        weights,
        object_type="lesion",
        use_lesion_weight=True,
        flair_channel=1,
        lesion_power_exponent=float(exponent[1]),
        cached_lesion_power_exponent=float(exponent[1]),
    )
    images = np.zeros((3, 3, 3, 2), dtype=np.float32)
    images[1, 1, 1] = image
    result = compute_lesion_toads_active_class_energies(
        images,
        class_centroid[None, :],
        lesion_centroid,
        black_hole_centroid,
        weights,
        ("lesion",),
        outlier_distance=2.0,
        use_lesion_weight=True,
        flair_channel=1,
        lesion_power_exponents=exponent,
    )
    assert result.class_energy[1, 1, 1, 0] == expected.class_energy
    assert result.lesion_energy[1, 1, 1] == expected.lesion_energy


def test_active_lesion_weighting_requires_a_valid_flair_exponent() -> None:
    images = np.ones((3, 3, 3, 2), dtype=np.float32)
    centroids = np.zeros((1, 2), dtype=np.float32)
    auxiliary = np.zeros(2, dtype=np.float32)
    weights = np.ones(2, dtype=np.float32)
    with pytest.raises(ValueError, match="exponents are required"):
        compute_lesion_toads_active_class_energies(
            images, centroids, auxiliary, auxiliary, weights, ("lesion",),
            outlier_distance=1.0, use_lesion_weight=True, flair_channel=1,
        )
    with pytest.raises(ValueError, match="valid channel index"):
        compute_lesion_toads_active_class_energies(
            images, centroids, auxiliary, auxiliary, weights, ("lesion",),
            outlier_distance=1.0, use_lesion_weight=True, flair_channel=2,
            lesion_power_exponents=np.ones(2, dtype=np.float32),
        )
