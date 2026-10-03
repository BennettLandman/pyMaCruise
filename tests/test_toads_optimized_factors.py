from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_atlas import ToadsAtlasError, load_toads_optimized_factors


ROOT = Path(__file__).parent.parent
DESCRIPTOR = ROOT / "tests/data/synthetic/toads_intensity_variance_factors/mini-atlas.txt"


def test_optimized_factors_read_source_modality_rows_and_keep_default_ones(
    tmp_path: Path,
) -> None:
    descriptor = tmp_path / "optimized-factors.txt"
    descriptor.write_text(
        DESCRIPTOR.read_text()
        + "\nOptimized factors:\nT1_MPRAGE 0.8 1.0 1.0\nPD 1.2 1.3 1.4\n"
    )

    actual = load_toads_optimized_factors(descriptor, ["MPRAGE", "T2", "PD"])

    expected = np.asarray(
        [[0.8, 1.0, 1.0], [1.0, 1.0, 1.0], [1.2, 1.3, 1.4]],
        dtype=np.float32,
    )
    np.testing.assert_array_equal(actual, expected)


def test_optimized_factor_parser_rejects_malformed_profile(tmp_path: Path) -> None:
    descriptor = tmp_path / "malformed-optimized-factors.txt"
    descriptor.write_text(
        DESCRIPTOR.read_text() + "\nOptimized factors:\nT1_MPRAGE 0.8 1.0\n"
    )

    with pytest.raises(ToadsAtlasError, match="expected 3"):
        load_toads_optimized_factors(descriptor, ["MPRAGE"])
