"""JSON-safe encoding for IEEE floating-point values recorded as evidence."""

from __future__ import annotations

import math


def json_number(value: float) -> float | str:
    """Keep finite values numeric and encode NaN or infinity explicitly.

    Strict JSON has no NaN or Infinity literals. The returned strings can be
    parsed with ``float(value)`` when a checkpoint is hydrated.
    """
    number = float(value)
    if math.isnan(number):
        return "NaN"
    if math.isinf(number):
        return "Infinity" if number > 0 else "-Infinity"
    return number
