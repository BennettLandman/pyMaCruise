from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_atlas import ToadsAtlasError, load_toads_intensity_priors


ROOT = Path(__file__).parent.parent
REFERENCE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_atlas_priors"
    / "java_reference.json"
)
BASE_DESCRIPTOR = ROOT / "data" / "atlases" / "cruise-atlas-10obj-toads2008.txt"


def test_selected_base_t1_priors_match_java_descriptor_harness() -> None:
    expected = json.loads(REFERENCE.read_text())
    actual = load_toads_intensity_priors(BASE_DESCRIPTOR, ["T1_MPRAGE"])

    assert actual.shape == (1, expected["class_count"])
    assert expected["modalities"] == ["T1_MPRAGE"]
    np.testing.assert_array_equal(actual, np.asarray(expected["priors"], dtype=np.float32))


def test_source_modality_aliases_and_requested_channel_order() -> None:
    result = load_toads_intensity_priors(
        BASE_DESCRIPTOR, ["MPRAGE", "T1", "FLAIR"]
    )
    assert result.shape == (3, 11)
    np.testing.assert_array_equal(
        result[0], load_toads_intensity_priors(BASE_DESCRIPTOR, ["T1_MPRAGE"])[0]
    )
    np.testing.assert_array_equal(
        result[1], load_toads_intensity_priors(BASE_DESCRIPTOR, ["T1_SPGR"])[0]
    )
    np.testing.assert_array_equal(
        result[2], np.asarray([0.0, 0.5, 0.5, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8], dtype=np.float32)
    )


def test_prior_parser_rejects_missing_or_malformed_profiles(tmp_path: Path) -> None:
    descriptor = tmp_path / "mini-atlas.txt"
    descriptor.write_text(
        "Structures: 2\n"
        "Intensity Atlas: 0\n"
        "T1_MPRAGE 0.1 0.9\n"
        "End of atlas\n"
    )
    with pytest.raises(ToadsAtlasError, match="no explicit intensity profile"):
        load_toads_intensity_priors(descriptor, ["FLAIR"])

    descriptor.write_text(
        "Structures: 2\n"
        "Intensity Atlas: 0\n"
        "T1_MPRAGE 0.1\n"
        "End of atlas\n"
    )
    with pytest.raises(ToadsAtlasError, match="expected 2"):
        load_toads_intensity_priors(descriptor, ["T1_MPRAGE"])


def test_prior_parser_rejects_unknown_modality() -> None:
    with pytest.raises(ToadsAtlasError, match="Unsupported TOADS image modality"):
        load_toads_intensity_priors(BASE_DESCRIPTOR, ["T1_TURBO_UNKNOWN"])
