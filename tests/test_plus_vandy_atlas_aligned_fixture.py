"""Regression checks for the atlas-aligned plus synthetic tissue labels."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


FIXTURE_DIR = Path(__file__).parent / "data/synthetic/plus_vandy_atlas_aligned"
SPEC = importlib.util.spec_from_file_location(
    "plus_vandy_atlas_aligned_generator", FIXTURE_DIR / "generate_fixture.py"
)
assert SPEC is not None and SPEC.loader is not None
GENERATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GENERATOR)


def test_atlas_topology_maps_cerebral_and_non_cortical_tissues_separately() -> None:
    topology = np.array([1, 5, 6, 14, 15, 16, 17, 18, 22, 24, 25, 26, 27, 80])

    actual = GENERATOR._maseg_from_topology(topology)

    np.testing.assert_array_equal(
        actual,
        np.array([0, 5, 5, 75, 100, 75, 75, 75, 75, 75, 44, 44, 44, 0], dtype=np.uint16),
    )


def test_atlas_aligned_fixture_contains_required_tissue_codes() -> None:
    image = GENERATOR.nib.load(str(FIXTURE_DIR / "maseg.nii.gz"))
    labels = np.asanyarray(image.dataobj)

    assert {0, 5, 44, 75, 100}.issubset(set(np.unique(labels).tolist()))
    assert GENERATOR.nib.aff2axcodes(image.affine) == ("R", "A", "S")
