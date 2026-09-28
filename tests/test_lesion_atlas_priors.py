from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_atlas_priors import (
    LesionAtlasDescriptorError,
    load_lesion_atlas_priors,
)


ROOT = Path(__file__).parents[1]
ATLAS_ROOT = ROOT / "data" / "atlases" / "lesionToads-atlas-2012"
SYNTHETIC = Path(__file__).parent / "data" / "synthetic" / "lesion_atlas_priors"


def test_minimal_lesion_atlas_preserves_structure_order_and_modality_order() -> None:
    priors = load_lesion_atlas_priors(
        SYNTHETIC / "minimal-lesion-atlas.txt", ["FLAIR", "MPRAGE"]
    )

    assert [structure.name for structure in priors.structures] == [
        "Background",
        "Lesion-WM",
        "Cortex",
    ]
    assert priors.modalities == ("FLAIR", "MPRAGE")
    np.testing.assert_array_equal(
        priors.intensity_priors,
        np.asarray([[0.1, 1.2, 0.3], [0.0, 0.4, 0.8]], dtype=np.float32),
    )
    np.testing.assert_array_equal(
        priors.lesion_priors, np.asarray([1.2, 0.4], dtype=np.float32)
    )
    np.testing.assert_array_equal(
        priors.black_hole_priors, np.asarray([0.6, 0.1], dtype=np.float32)
    )
    assert not priors.intensity_priors.flags.writeable


def test_bundled_2012_descriptor_matches_selected_t1_flair_priors() -> None:
    priors = load_lesion_atlas_priors(
        ATLAS_ROOT / "cruise-atlas-12obj-lesiontoads2012.txt",
        ["T1_MPRAGE", "FLAIR"],
    )

    assert [structure.name for structure in priors.structures] == [
        "Background",
        "Sulcal-CSF",
        "Ventricles",
        "Cerebellum-GM",
        "Cerebrum-GM",
        "Caudate",
        "Thalamus",
        "Putamen",
        "Brainstem",
        "Cerebellum-WM",
        "OuterVentricular-WM",
        "InterVentricular-WM",
        "Frontal-WM",
    ]
    np.testing.assert_array_equal(
        priors.intensity_priors[0],
        np.asarray(
            [0.0, 0.1, 0.1, 0.7, 0.7, 0.6, 0.6, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0],
            dtype=np.float32,
        ),
    )
    np.testing.assert_array_equal(
        priors.intensity_priors[1],
        np.asarray(
            [0.0, 0.5, 0.5, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7, 0.7, 0.7, 1.2, 0.75],
            dtype=np.float32,
        ),
    )
    np.testing.assert_array_equal(
        priors.lesion_priors, np.asarray([0.4, 1.2], dtype=np.float32)
    )
    np.testing.assert_array_equal(
        priors.black_hole_priors, np.asarray([0.1, 0.6], dtype=np.float32)
    )


def test_lesion_atlas_parser_rejects_incomplete_selected_modality(tmp_path) -> None:
    source = (SYNTHETIC / "minimal-lesion-atlas.txt").read_text()
    descriptor = tmp_path / "incomplete.txt"
    descriptor.write_text(source.replace("FLAIR\t0.1\t1.2\t0.3\n", ""))

    with pytest.raises(LesionAtlasDescriptorError, match="No intensity row.*FLAIR"):
        load_lesion_atlas_priors(descriptor, ["T1_MPRAGE", "FLAIR"])


def test_lesion_atlas_parser_rejects_row_with_wrong_structure_count(tmp_path) -> None:
    source = (SYNTHETIC / "minimal-lesion-atlas.txt").read_text()
    descriptor = tmp_path / "wrong-count.txt"
    descriptor.write_text(source.replace("T1_MPRAGE\t0.0\t0.4\t0.8", "T1_MPRAGE\t0.0\t0.4"))

    with pytest.raises(LesionAtlasDescriptorError, match="expected 3"):
        load_lesion_atlas_priors(descriptor, ["T1_MPRAGE"])
