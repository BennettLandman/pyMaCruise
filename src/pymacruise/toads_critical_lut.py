"""Runtime access to TOADS' source 18/6 critical-point lookup table.

The table is distinct from the P25 6/18 topology predicate in ``topology``.
Its 2**26 entries index the 26 neighbors around a center voxel, packed with
the Java ``CriticalPointLUT.keyFromPattern(boolean[][][])`` convention:
x-major, then y, then z, with the center omitted.

The LUT is the local digital-topology primitive behind the homeomorphic
multi-object constraint discussed by Bazin & Pham (2007, 2008). The copied
selected-source table and lookup define its behavior; a lookup by itself is
not the full multi-label topology test. See ``docs/papers/toads.md``.
"""

from __future__ import annotations

from functools import lru_cache
import gzip
from importlib.resources import files
from pathlib import Path
from typing import TypeAlias

import numpy as np

Neighborhood: TypeAlias = np.ndarray
LUT_SIZE = 1 << 26


def key_from_neighborhood(neighborhood: Neighborhood) -> int:
    """Return the Java 26-bit key for a binary 3x3x3 neighborhood.

    Bits are assigned in the order ``(x, y, z)`` produced by three nested
    loops over ``range(3)``; ``(1, 1, 1)`` is skipped. The center value is
    ignored, as it is by Java's boolean-array key method.

    The input must contain booleans or numeric 0/1 values. This helper only
    packs the neighborhood; it does not decide whether the point is simple.
    """
    values = np.asarray(neighborhood)
    if values.shape != (3, 3, 3):
        raise ValueError("A TOADS LUT neighborhood must have shape (3, 3, 3)")
    if values.dtype != np.bool_:
        try:
            valid = np.isfinite(values).all() and np.isin(values, (0, 1)).all()
        except TypeError as exc:
            raise ValueError("Neighborhood values must be boolean or 0/1") from exc
        if not valid:
            raise ValueError("Neighborhood values must be boolean or 0/1")

    key = 0
    bit = 0
    for x in range(3):
        for y in range(3):
            for z in range(3):
                if (x, y, z) == (1, 1, 1):
                    continue
                if bool(values[x, y, z]):
                    key |= 1 << bit
                bit += 1
    return key


def critical186_is_regular(
    neighborhood: Neighborhood,
) -> bool:
    """Query the selected VandyCRUISE 18/6 LUT for one binary neighborhood.

    This is only the source `CriticalPointLUT.get(keyFromPattern(...))`
    predicate. Multi-label topology preservation also checks the changed
    object and its neighboring label relations; this helper does not perform
    that higher-level segmentation check.
    """
    key = key_from_neighborhood(neighborhood)
    table = load_critical186_lut()
    return bool(table[key] == 1)


def load_critical186_lut(path: str | Path | None = None) -> np.ndarray:
    """Load the bundled source LUT, or a specified compressed LUT file.

    Returns a read-only ``uint8`` array indexed directly by
    :func:`key_from_neighborhood`. A table value is ``1`` for a regular
    neighborhood and ``0`` otherwise, matching Java's ``list[n] == 1`` test.
    A custom path is provided for validation and source-table comparisons.
    """
    if path is None:
        return _load_bundled_lut()
    return _read_lut(Path(path))


@lru_cache(maxsize=1)
def _load_bundled_lut() -> np.ndarray:
    resource = files("pymacruise").joinpath("data", "critical186LUT.raw.gz")
    try:
        compressed = resource.read_bytes()
    except OSError as exc:
        raise RuntimeError("Bundled TOADS critical186 LUT is unavailable") from exc
    return _decode_lut(compressed, str(resource))


def _read_lut(path: Path) -> np.ndarray:
    try:
        compressed = path.read_bytes()
    except OSError as exc:
        raise ValueError(f"Cannot read TOADS critical186 LUT: {path}") from exc
    return _decode_lut(compressed, str(path))


def _decode_lut(compressed: bytes, source: str) -> np.ndarray:
    try:
        payload = gzip.decompress(compressed)
    except (OSError, EOFError) as exc:
        raise ValueError(f"Invalid compressed TOADS critical186 LUT: {source}") from exc
    if len(payload) != LUT_SIZE:
        raise ValueError(
            f"TOADS critical186 LUT has {len(payload)} entries; expected {LUT_SIZE}: {source}"
        )
    table = np.frombuffer(payload, dtype=np.uint8)
    if np.any((table != 0) & (table != 1)):
        raise ValueError(f"TOADS critical186 LUT contains values other than 0 and 1: {source}")
    return table
