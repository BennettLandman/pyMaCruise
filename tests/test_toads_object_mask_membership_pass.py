"""Selected Java comparison for composed OBJECT/MASK membership passes."""

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_math import toads_object_mask_membership_pass, toads_object_membership_pass


REFERENCE = Path(__file__).parent / "data/synthetic/toads_object_mask_membership_pass/java_reference.json"


def _values(bits: list[int], shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32).reshape(shape)


def _bits(values: np.ndarray) -> list[int]:
    return values.astype(np.float32, copy=False).ravel().view(np.int32).tolist()


def test_object_mask_membership_pass_matches_selected_java() -> None:
    reference = json.loads(REFERENCE.read_text())
    spatial_shape = tuple(reference["shape"][:3])
    channels = reference["shape"][3]
    classes = reference["classes"]
    energy_shape = (*spatial_shape, classes)
    images = _values(reference["images_bits"], (*spatial_shape, channels))
    centroids = _values(reference["centroids_bits"], (channels, classes))
    modality_weights = _values(reference["modality_weights_bits"], (channels,))
    class_variances = _values(reference["class_variance_bits"], (channels, classes))
    image_variances = _values(reference["image_variances_bits"], (channels,))
    initial = _values(reference["initial_memberships_bits"], energy_shape)
    shapes = _values(reference["transformed_shapes_bits"], energy_shape)
    closest = np.asarray(reference["closest_labels"], dtype=np.int8).reshape((*spatial_shape, 4))
    similarity = _values(reference["similarity_bits"], (classes, classes))
    smoothing_priors = _values(reference["smoothing_prior_bits"], (classes,))
    atlas_priors = _values(reference["atlas_prior_bits"], (classes,))

    for case in reference["cases"]:
        field = None if case["intensity_field"] is None else _values(case["intensity_field"], (*spatial_shape, channels))
        bias = None if case["psf_bias"] is None else _values(case["psf_bias"], (*spatial_shape, channels))
        initial_before = initial.copy()
        result, objective = toads_object_mask_membership_pass(
            images,
            initial,
            centroids,
            modality_weights,
            class_variances,
            image_variances,
            shapes,
            closest,
            smoothing_priors,
            similarity,
            atlas_priors,
            reference["smoothing"],
            reference["prior_factor"],
            reference["prior_coefficient"],
            class_types=reference["class_types"],
            intensity_field=field,
            psf_bias=bias,
        )
        assert _bits(result) == case["updated_bits"]
        assert _bits(np.asarray([objective])) == [case["objective_bits"]]
        assert result.dtype == np.float32
        np.testing.assert_array_equal(initial, initial_before)
        np.testing.assert_array_equal(result[0], initial[0])
        np.testing.assert_array_equal(result[:, -1], initial[:, -1])
        np.testing.assert_array_equal(result[:, :, 0], initial[:, :, 0])


def test_object_membership_pass_rejects_mismatched_model_grid() -> None:
    shape = (3, 3, 3)
    images = np.ones((*shape, 1), dtype=np.float32)
    memberships = np.ones((*shape, 2), dtype=np.float32) / 2
    with pytest.raises(ValueError, match="centroids"):
        toads_object_membership_pass(
            images,
            memberships,
            np.ones((2, 2), dtype=np.float32),
            np.ones(1, dtype=np.float32),
            np.ones((1, 2), dtype=np.float32),
            np.ones(1, dtype=np.float32),
            np.ones_like(memberships),
            np.zeros((*shape, 4), dtype=np.int8),
            np.ones(2, dtype=np.float32),
            np.zeros((2, 2), dtype=np.float32),
            np.ones(2, dtype=np.float32),
            0.1,
            0.5,
            0.2,
        )
