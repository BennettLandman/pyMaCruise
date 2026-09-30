"""Selected Java comparison for composed OBJECT/MASK/OUTLIER membership."""

import json
from pathlib import Path

import numpy as np

from pymacruise.toads_math import toads_object_mask_outlier_membership_pass


REFERENCE = Path(__file__).parent / "data/synthetic/toads_object_mask_outlier_membership_pass/java_reference.json"


def _values(bits: list[int], shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32).reshape(shape)


def _bits(values: np.ndarray) -> list[int]:
    return values.astype(np.float32, copy=False).ravel().view(np.int32).tolist()


def test_object_mask_outlier_membership_pass_matches_selected_java() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    spatial = tuple(reference["shape"])
    classes = reference["classes"]
    class_shape = (*spatial, classes)
    images = _values(reference["images_bits"], (*spatial, 1))
    field = _values(reference["field_bits"], (*spatial, 1))
    bias = _values(reference["bias_bits"], (*spatial, 1))
    initial = _values(reference["initial_memberships_bits"], class_shape)
    shapes = _values(reference["shape_bits"], class_shape)
    labels = np.asarray(reference["closest_labels"], dtype=np.int8).reshape((*spatial, 4))
    initial_outlier = _values(reference["initial_outlier_bits"], spatial)
    initial_before = initial.copy()
    outlier_before = initial_outlier.copy()
    updated, outlier, objective = toads_object_mask_outlier_membership_pass(
        images,
        initial,
        initial_outlier,
        _values(reference["centroids_bits"], (1, classes)),
        _values([reference["modality_weight_bits"]], (1,)),
        _values(reference["variance_bits"], (1, classes)),
        _values([reference["image_variance_bits"]], (1,)),
        reference["outlier_value"],
        shapes,
        labels,
        _values(reference["smoothing_priors_bits"], (classes,)),
        _values(reference["similarity_bits"], (classes, classes)),
        _values(reference["atlas_priors_bits"], (classes,)),
        reference["smoothing"],
        reference["prior_factor"],
        reference["prior_coefficient"],
        class_types=reference["class_types"],
        intensity_field=field,
        psf_bias=bias,
    )
    assert _bits(updated) == reference["updated_bits"]
    assert _bits(outlier) == reference["updated_outlier_bits"]
    assert _bits(np.asarray([objective])) == [reference["objective_bits"]]
    np.testing.assert_array_equal(initial, initial_before)
    np.testing.assert_array_equal(initial_outlier, outlier_before)

