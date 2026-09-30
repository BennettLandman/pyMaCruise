"""Selected Java check for the TOADS MASK binary mismatch energy."""

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_math import toads_mask_intensity_energy


REFERENCE = Path(__file__).parent / "data/synthetic/toads_mask_intensity_energy/java_reference.json"


def _values(bits: list[int], shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32).reshape(shape)


def _bits(values: np.ndarray) -> list[int]:
    return values.astype(np.float32, copy=False).ravel().view(np.int32).tolist()


def test_mask_intensity_energy_matches_selected_java_bits() -> None:
    reference = json.loads(REFERENCE.read_text())
    shape = tuple(reference["shape"])
    classes = reference["classes"]
    channels = shape[-1]
    images = _values(reference["images_bits"], shape)
    modality_weights = _values(reference["modality_weights_bits"], (channels,))
    class_variances = _values(reference["class_variance_bits"], (channels, classes))
    image_variances = _values(reference["image_variances_bits"], (channels,))
    copies = [value.copy() for value in (images, modality_weights, class_variances, image_variances)]

    result = toads_mask_intensity_energy(
        images, modality_weights, class_variances, image_variances
    )

    assert result.shape == shape[:3] + (classes,)
    assert result.dtype == np.float32
    assert _bits(result) == reference["energy_bits"]
    for actual, copied in zip((images, modality_weights, class_variances, image_variances), copies, strict=True):
        np.testing.assert_array_equal(actual, copied)


def test_mask_intensity_energy_rejects_malformed_weights() -> None:
    images = np.ones((2, 2, 2, 2), dtype=np.float32)
    with pytest.raises(ValueError, match="class_variance_weights"):
        toads_mask_intensity_energy(images, np.ones(2), np.ones((3, 2)), np.ones(2))
    invalid = np.ones((2, 2, 2, 2), dtype=np.float32)
    invalid[0, 0, 0, 0] = np.inf
    with pytest.raises(ValueError, match="image_values"):
        toads_mask_intensity_energy(invalid, np.ones(2), np.ones((2, 1)), np.ones(2))
