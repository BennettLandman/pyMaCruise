from __future__ import annotations

import numpy as np
import pytest

from pymacruise.toads_registered_loop import (
    toads_effective_prior_coefficient,
    toads_effective_smoothing,
)


def test_selected_wrapper_coefficients_follow_java_float32_division_order() -> None:
    smoothing = np.float32(toads_effective_smoothing(0.2, 11))
    prior = np.float32(toads_effective_prior_coefficient(2.0, 11))

    java_smoothing = np.float32(np.float32(np.float32(0.2) / np.float32(6.0)) / np.float32(11.0))
    java_prior = np.float32(np.float32(2.0) / np.float32(11.0))
    python_double_smoothing = np.float32(0.2 / 6.0 / 11.0)

    assert smoothing.view(np.uint32) == java_smoothing.view(np.uint32)
    assert prior.view(np.uint32) == java_prior.view(np.uint32)
    assert smoothing.view(np.uint32) != python_double_smoothing.view(np.uint32)


@pytest.mark.parametrize("class_count", [0, -1, 1.5, True])
def test_effective_coefficients_reject_invalid_class_counts(class_count: object) -> None:
    with pytest.raises(ValueError, match="class_count"):
        toads_effective_smoothing(0.2, class_count)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="class_count"):
        toads_effective_prior_coefficient(2.0, class_count)  # type: ignore[arg-type]


@pytest.mark.parametrize("value", [-0.1, float("nan"), float("inf")])
def test_effective_coefficients_reject_invalid_wrapper_values(value: float) -> None:
    with pytest.raises(ValueError, match="raw_smoothing"):
        toads_effective_smoothing(value, 11)
    with pytest.raises(ValueError, match="raw_coefficient"):
        toads_effective_prior_coefficient(value, 11)
