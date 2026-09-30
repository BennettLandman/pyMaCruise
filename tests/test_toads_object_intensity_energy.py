"""Selected Java comparison for the TOADS ordinary OBJECT data term."""

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_math import toads_object_intensity_energy


REFERENCE = Path(__file__).parent / "data/synthetic/toads_object_intensity_energy/java_reference.json"


def _values(bits: list[int], shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32).reshape(shape)


def _bits(values: np.ndarray) -> list[int]:
    return values.astype(np.float32, copy=False).ravel().view(np.int32).tolist()


def test_object_intensity_energy_matches_selected_java_bits() -> None:
    reference = json.loads(REFERENCE.read_text())
    for case in reference["cases"]:
        shape = tuple(case["shape"])
        classes = case["classes"]
        images = _values(case["images"], shape)
        centroids = _values(case["centroids"], (shape[-1], classes))
        modality_weights = _values(case["modality_weights"], (shape[-1],))
        class_weights = _values(case["class_variance_weights"], (shape[-1], classes))
        image_variances = _values(case["image_variances"], (shape[-1],))
        field = None if case["field"] is None else _values(case["field"], shape)
        bias = None if case["psf_bias"] is None else _values(case["psf_bias"], shape)
        inputs = [images.copy(), centroids.copy(), modality_weights.copy(), class_weights.copy(), image_variances.copy()]
        if field is not None:
            inputs.append(field.copy())
        if bias is not None:
            inputs.append(bias.copy())

        result = toads_object_intensity_energy(
            images,
            centroids,
            modality_weights,
            class_weights,
            image_variances,
            intensity_field=field,
            psf_bias=bias,
        )

        assert result.shape == shape[:3] + (classes,)
        assert result.dtype == np.float32
        assert _bits(result) == case["energy_bits"]
        for original, copied in zip(
            [images, centroids, modality_weights, class_weights, image_variances]
            + ([] if field is None else [field])
            + ([] if bias is None else [bias]),
            inputs,
            strict=True,
        ):
            np.testing.assert_array_equal(original, copied)


def test_object_intensity_energy_rejects_invalid_shapes_and_nonfinite_values() -> None:
    images = np.ones((2, 2, 2, 2), dtype=np.float32)
    centroids = np.ones((2, 3), dtype=np.float32)
    modality = np.ones(2, dtype=np.float32)
    class_weights = np.ones((2, 3), dtype=np.float32)
    image_variances = np.ones(2, dtype=np.float32)
    with pytest.raises(ValueError, match="centroids"):
        toads_object_intensity_energy(images, np.ones((3, 2)), modality, class_weights, image_variances)
    invalid_field = np.ones_like(images)
    invalid_field[0, 0, 0, 0] = np.inf
    with pytest.raises(ValueError, match="intensity_field"):
        toads_object_intensity_energy(
            images, centroids, modality, class_weights, image_variances, intensity_field=invalid_field
        )
    invalid_bias = np.ones_like(images)
    invalid_bias[0, 0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="psf_bias"):
        toads_object_intensity_energy(
            images, centroids, modality, class_weights, image_variances, psf_bias=invalid_bias
        )
