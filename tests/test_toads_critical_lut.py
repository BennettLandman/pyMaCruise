from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_critical_lut import (
    LUT_SIZE,
    critical186_is_regular,
    key_from_neighborhood,
    load_critical186_lut,
)


JAVA_REFERENCE = (
    Path(__file__).parent
    / "data"
    / "synthetic"
    / "toads_critical_lut"
    / "java_reference.json"
)


def test_runtime_lut_is_source_sized_read_only_and_has_expected_entries() -> None:
    table = load_critical186_lut()
    assert table.shape == (LUT_SIZE,)
    assert table.dtype == np.uint8
    assert not table.flags.writeable
    assert table[0] == 0
    assert table[1] == 1  # Java bit 0: (0, 0, 0)
    assert table[1 << 12] == 1  # Java bit 12: (1, 1, 0)


def test_key_encoding_matches_java_x_major_y_z_order_and_omits_center() -> None:
    neighborhood = np.zeros((3, 3, 3), dtype=bool)
    neighborhood[0, 0, 0] = True
    neighborhood[1, 1, 0] = True
    neighborhood[1, 1, 1] = True  # Center never contributes a bit.
    assert key_from_neighborhood(neighborhood) == (1 | (1 << 12))

    neighborhood[:] = False
    neighborhood[2, 2, 2] = True
    assert key_from_neighborhood(neighborhood) == (1 << 25)


@pytest.mark.parametrize(
    "neighborhood",
    [np.zeros((3, 3), dtype=bool), np.full((3, 3, 3), 2, dtype=np.uint8)],
)
def test_key_rejects_wrong_shape_or_nonbinary_values(neighborhood: np.ndarray) -> None:
    with pytest.raises(ValueError):
        key_from_neighborhood(neighborhood)


def test_loader_rejects_wrong_payload_length(tmp_path) -> None:
    malformed = tmp_path / "short.raw.gz"
    malformed.write_bytes(gzip.compress(b"short payload"))
    with pytest.raises(ValueError, match="expected 67108864"):
        load_critical186_lut(malformed)


def test_loader_rejects_invalid_gzip(tmp_path) -> None:
    malformed = tmp_path / "not-gzip.raw.gz"
    malformed.write_bytes(b"not a gzip stream")
    with pytest.raises(ValueError, match="Invalid compressed"):
        load_critical186_lut(malformed)


def test_runtime_lut_matches_selected_java_key_and_regular_lookup() -> None:
    reference = json.loads(JAVA_REFERENCE.read_text())
    assert reference["source"] == "VandyCRUISE CriticalPointLUT"
    assert len(reference["cases"]) == 68
    for case in reference["cases"]:
        neighborhood = np.fromiter(
            (value == "1" for value in case["neighborhood"]),
            dtype=np.bool_,
            count=27,
        ).reshape((3, 3, 3))
        assert key_from_neighborhood(neighborhood) == case["key"]
        assert critical186_is_regular(neighborhood) is case["regular"]


def test_critical186_regular_predicate_ignores_center_and_checks_lut_shape() -> None:
    neighborhood = np.zeros((3, 3, 3), dtype=bool)
    neighborhood[0, 0, 0] = True
    expected = critical186_is_regular(neighborhood)
    neighborhood[1, 1, 1] = True
    assert critical186_is_regular(neighborhood) is expected
