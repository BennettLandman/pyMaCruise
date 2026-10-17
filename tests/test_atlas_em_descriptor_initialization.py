from __future__ import annotations

from pathlib import Path

import numpy as np

from pymacruise.atlas_em_preprocessing import (
    AtlasEMNormalizedImages,
    initialize_atlas_em_centroids,
    initialize_atlas_em_memberships,
    initialize_atlas_em_memberships_from_descriptor,
)


def test_descriptor_intensity_profiles_initialize_memberships_in_source_class_order(
    tmp_path: Path,
) -> None:
    descriptor = tmp_path / "atlas.txt"
    descriptor.write_text(
        "\n".join(
            (
                "Structure Atlas File (edit at your own risks)", "", "Structures:\t2",
                "Background\t1\tmask", "Cortex\t2\tclosed", "", "Topology Atlas:",
                "File:\ttopology.nii.gz", "Dimensions:\t5\t6\t7", "Resolutions:\t1\t1\t1",
                "", "Shape Atlas:\t4", "Dimensions:\t5\t6\t7", "Resolutions:\t1\t1\t1",
                "Orientations:\t1\t1\t6\t3", "Structure:\tBackground", "File:\tbackground.nii.gz",
                "Structure:\tCortex", "File:\tcortex.nii.gz", "", "Intensity Atlas:\t10",
                "T1_MPRAGE\t0.2\t0.8",
            )
        ) + "\n",
        encoding="utf-8",
    )
    values = np.linspace(0.0, 1.5, 4 * 5 * 6, dtype=np.float32).reshape(4, 5, 6, 1)
    normalized = AtlasEMNormalizedImages(
        values,
        np.array([10.0], dtype=np.float32),
        np.array([20.0], dtype=np.float32),
    )

    result = initialize_atlas_em_memberships_from_descriptor(
        normalized, descriptor, ["MPRAGE"]
    )

    expected_centroids = initialize_atlas_em_centroids(
        normalized.robust_minima,
        normalized.robust_maxima,
        np.array([[0.2, 0.8]], dtype=np.float32),
    )
    expected_memberships = initialize_atlas_em_memberships(values, expected_centroids)
    assert result.class_names == ("Background", "Cortex")
    np.testing.assert_array_equal(result.centroids, expected_centroids)
    np.testing.assert_array_equal(result.memberships, expected_memberships)
    assert result.memberships.shape == (4, 5, 6, 2)
