from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_atlas import ToadsAtlasError, load_toads_intensity_variance_factors


ROOT = Path(__file__).parent.parent
FIXTURE = Path(__file__).parent / "data" / "synthetic" / "toads_intensity_variance_factors"
DESCRIPTOR = FIXTURE / "mini-atlas.txt"
JAVA_REFERENCE = FIXTURE / "java_reference.json"


def test_variance_factors_match_java_reference_and_keep_channel_order() -> None:
    expected = json.loads(JAVA_REFERENCE.read_text())
    actual = load_toads_intensity_variance_factors(
        DESCRIPTOR, expected["modalities"]
    )

    assert actual.shape == (len(expected["modalities"]), expected["class_count"])
    assert actual.dtype == np.float32
    np.testing.assert_array_equal(actual, np.asarray(expected["factors"], dtype=np.float32))


def test_variance_factor_parser_does_not_read_intensity_means() -> None:
    actual = load_toads_intensity_variance_factors(DESCRIPTOR, ["T1_MPRAGE"])
    np.testing.assert_array_equal(actual[0], np.asarray([0.25, 1.5, 2.75], dtype=np.float32))
    assert not np.array_equal(actual[0], np.asarray([11, 12, 13], dtype=np.float32))


def test_variance_factors_default_to_ones_without_section_or_for_known_unlisted_modality(
    tmp_path: Path,
) -> None:
    no_section = tmp_path / "no-variance-section.txt"
    no_section.write_text("Structures: 2\nCSF 1 closed\nGM 2 closed\nIntensity Atlas: 0\n")
    actual = load_toads_intensity_variance_factors(no_section, ["T1", "PD"])
    np.testing.assert_array_equal(actual, np.ones((2, 2), dtype=np.float32))


def test_variance_parser_rejects_malformed_class_count(tmp_path: Path) -> None:
    malformed = tmp_path / "malformed-atlas.txt"
    text = DESCRIPTOR.read_text().replace("T1_MPRAGE 0.25 1.5 2.75", "T1_MPRAGE 0.25 1.5")
    malformed.write_text(text)
    with pytest.raises(ToadsAtlasError, match="expected 3"):
        load_toads_intensity_variance_factors(malformed, ["MPRAGE"])
