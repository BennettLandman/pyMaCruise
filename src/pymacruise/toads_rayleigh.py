"""Source-faithful scalar Rayleigh distance kernel from VandyCRUISE TOADS."""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike, NDArray

# Match Java's expression evaluation: log(double pi / float 2), rounded to
# float, minus (float pi / float 2), with each float operation rounded.
_LOGPI2MPI2 = np.float32(
    np.float32(math.log(math.pi / 2.0))
    - np.float32(np.float32(math.pi) / np.float32(2.0))
)


def toads_rayleigh_distance(val: ArrayLike, sig: ArrayLike) -> NDArray[np.float32] | np.float32:
    """Evaluate the selected Java ``rayleighDistance(float,float)`` kernel.

    The result preserves Java float32 operation order, including the special
    ``val <= 0`` return of exactly ``1``. Inputs broadcast as NumPy arrays;
    scalar inputs return a NumPy float32 scalar. NaNs, infinities, signed zero,
    and IEEE overflow/underflow follow Java's float/double conversion points.
    """

    values, sigmas = np.broadcast_arrays(
        np.asarray(val, dtype=np.float32), np.asarray(sig, dtype=np.float32)
    )
    result = np.empty(values.shape, dtype=np.float32)
    nonpositive = values <= np.float32(0.0)
    result[nonpositive] = np.float32(1.0)
    active = ~nonpositive

    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        v = values[active]
        s = sigmas[active]
        # Java parses val*val/sig*sig left-to-right; force float32 rounding at
        # each operation rather than allowing a wider intermediate.
        square = np.multiply(v, v, dtype=np.float32)
        scaled = np.divide(square, s, dtype=np.float32)
        gaussian_term = np.multiply(scaled, s, dtype=np.float32)
        ratio = np.divide(v, s, dtype=np.float32)
        log_term = np.asarray(
            [np.float32(math.log(float(item))) if item > 0 else np.float32(np.nan) for item in ratio],
            dtype=np.float32,
        )
        twice_log = np.multiply(np.float32(2.0), log_term, dtype=np.float32)
        after_log = np.subtract(gaussian_term, twice_log, dtype=np.float32)
        result[active] = np.add(after_log, _LOGPI2MPI2, dtype=np.float32)

    if result.ndim == 0:
        return np.float32(result[()])
    return result
