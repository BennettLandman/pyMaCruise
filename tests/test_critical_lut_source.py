from __future__ import annotations

import gzip
import hashlib
from importlib.resources import files

from pymacruise.topology import topology_numbers_618


LUT_RESOURCE = files("pymacruise").joinpath("data", "critical186LUT.raw.gz")
LUT_SHA256 = "cbf1e339d78a0a2c673c939ee2c8a30ee3a35bfd9603d07ed2d48864bb0e88be"
LUT_SIZE = 1 << 26


def _key_pattern(key: int) -> list[list[list[bool]]]:
    """Decode the Java keyFromPattern order (x-major, y, z; center omitted)."""
    pattern = [[[False for _ in range(3)] for _ in range(3)] for _ in range(3)]
    bit = 0
    for x in range(3):
        for y in range(3):
            for z in range(3):
                if (x, y, z) == (1, 1, 1):
                    continue
                pattern[x][y][z] = bool((key >> bit) & 1)
                bit += 1
    return pattern


def test_selected_toads_lut_asset_matches_source_and_java_key_order() -> None:
    compressed = LUT_RESOURCE.read_bytes()
    assert hashlib.sha256(compressed).hexdigest() == LUT_SHA256
    lut = gzip.decompress(compressed)
    assert len(lut) == LUT_SIZE

    # Java's key bits 0 and 12 address (0,0,0) and (1,1,0), respectively.
    assert lut[0] == 0
    assert lut[1 << 12] == 1
    pattern = _key_pattern(1 << 12)
    assert pattern[1][1][0]
    assert sum(sum(sum(row) for row in plane) for plane in pattern) == 1


def test_toads_lut_is_not_the_p25_6_18_simple_point_predicate() -> None:
    # A single corner is ignored by TOADS' 18-neighbor connectivity, while
    # the P25 6/18 helper rejects the same key. Preserve this discrepancy as
    # an explicit guard until the exhaustive 18/6 comparison is completed.
    key = 1
    assert gzip.decompress(LUT_RESOURCE.read_bytes())[key] == 1
    assert topology_numbers_618(_key_pattern(key)) == (0, 1)
