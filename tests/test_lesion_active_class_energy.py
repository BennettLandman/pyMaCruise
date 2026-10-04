from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_active_class_energy import compute_lesion_toads_active_class_energies
from pymacruise.lesion_atlas_priors import load_lesion_atlas_priors


ATLAS = Path(__file__).parents[1] / "data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"


def test_bundled_thirteen_class_active_dispatch_uses_source_fallback_terms() -> None:
    descriptor = load_lesion_atlas_priors(str(ATLAS), ("T1", "FLAIR"))
    class_types = tuple(structure.topology for structure in descriptor.structures)
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
