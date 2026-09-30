"""Selected Java comparison for the TOADS OUTLIER data terms."""

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_math import toads_outlier_intensity_energy


REFERENCE = Path(__file__).parent / "data/synthetic/toads_outlier_intensity_energy/java_reference.json"


def _bits(values: np.ndarray) -> list[int]:
    return values.astype(np.float32, copy=False).ravel().view(np.int32).tolist()


def _values(bits: list[int], shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32).reshape(shape)


def test_outlier_intensity_energies_match_selected_java_bits() -> None:
    reference = json.loads(REFERENCE.read_text())
    shape = tuple(reference["shape"])
    classes = reference["classes"]
    for case in reference["cases"]:
        images = _values(case["images_bits"], shape)
        centroids = _values(case["centroids_bits"], (shape[-1], classes))
        modality = _values(case["modality_bits"], (shape[-1],))
        class_weights = _values(case["class_weight_bits"], (shape[-1], classes))
        image_variances = _values(case["image_variance_bits"], (shape[-1],))
        selected_field = None if case["field_bits"] is None else _values(case["field_bits"], shape)
        selected_bias = None if case["bias_bits"] is None else _values(case["bias_bits"], shape)
        residual, outlier = toads_outlier_intensity_energy(
            images,
            centroids,
            modality,
            class_weights,
            image_variances,
            case["outlier"],
            intensity_field=selected_field,
            psf_bias=selected_bias,
        )
        assert _bits(residual) == case["residual_bits"]
        assert _bits(outlier) == case["outlier_bits"]


def test_outlier_intensity_energy_rejects_nonfinite_scalar() -> None:
    with pytest.raises(ValueError, match="outlier_value"):
        toads_outlier_intensity_energy(
            np.ones((1, 1, 1, 1), dtype=np.float32),
            np.zeros((1, 1), dtype=np.float32),
            np.ones(1, dtype=np.float32),
            np.ones((1, 1), dtype=np.float32),
            np.ones(1, dtype=np.float32),
            float("nan"),
        )
