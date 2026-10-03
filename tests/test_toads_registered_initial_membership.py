from __future__ import annotations

from types import SimpleNamespace
from pathlib import Path

import numpy as np
import pytest

from pymacruise.toads_default_atlas_membership import (
    DEFAULT_ATLAS_CLASS_TYPES,
    toads_default_atlas_membership_pass,
)
from pymacruise.toads_initial_statistics_membership import (
    compute_toads_initial_statistics_membership,
)
from pymacruise.toads_registered_initial_membership import (
    compute_registered_toads_initial_membership,
)
from pymacruise.toads_atlas import load_toads_intensity_variance_factors
from test_toads_default_atlas_dispatch import _fixture


def _alignment_run():
    data, parameters, _cases, _reference = _fixture()
    data = dict(data)
    parameters = dict(parameters)
    data["images"] = data["images"][..., :1].copy()
    parameters["centroids"] = parameters["centroids"][:1].copy()
    parameters["exponent_by_channel_class"] = parameters["exponent_by_channel_class"][:1].copy()
    parameters["modality_weights"] = parameters["modality_weights"][:1].copy()
    shape = data["images"].shape[:3]
    x, y, z = np.indices(shape)
    segmentation = ((x + 2 * y + 3 * z) % 11).astype(np.int8)
    channel_first_images = np.moveaxis(data["images"], -1, 0)
    result = SimpleNamespace(
        atlas=SimpleNamespace(
            descriptor=Path(__file__).parent.parent / "data/atlases/cruise-atlas-10obj-toads2008.txt",
            structures=tuple(
                SimpleNamespace(topology=kind.lower())
                for kind in DEFAULT_ATLAS_CLASS_TYPES
            )
        ),
        template_refresh=SimpleNamespace(recropped_images=channel_first_images),
        segmentation_initialization=SimpleNamespace(
            memberships=data["initial_memberships"],
            centroids=parameters["centroids"],
        ),
        segmentation_constructor=SimpleNamespace(segmentation=segmentation),
        relation_state=SimpleNamespace(closest_labels=data["closest_labels"]),
        transformed_shape_priors=data["transformed_shapes"],
        modality="T1_MPRAGE",
    )
    return result, data, parameters, segmentation


def test_registered_startup_arrays_feed_initial_statistics_and_membership() -> None:
    run, data, parameters, segmentation = _alignment_run()
    inputs = dict(
        exponent_by_channel_class=parameters["exponent_by_channel_class"],
        modality_weights=parameters["modality_weights"],
        smoothing_priors=parameters["smoothing_priors"],
        atlas_shape_priors=parameters["atlas_shape_priors"],
        smoothing=0.043,
        prior_factor=parameters["prior_factor"],
        prior_coefficient=0.16,
    )

    result = compute_registered_toads_initial_membership(run, **inputs)

    expected = compute_toads_initial_statistics_membership(
        image_values=data["images"],
        initial_memberships=data["initial_memberships"],
        centroids=parameters["centroids"],
        segmentation=segmentation,
        class_types=DEFAULT_ATLAS_CLASS_TYPES,
        class_variance_weights=load_toads_intensity_variance_factors(
            run.atlas.descriptor, [run.modality]
        ),
        exponent_by_channel_class=parameters["exponent_by_channel_class"],
        modality_weights=parameters["modality_weights"],
        transformed_shapes=data["transformed_shapes"],
        closest_labels=data["closest_labels"],
        smoothing_priors=parameters["smoothing_priors"],
        atlas_shape_priors=parameters["atlas_shape_priors"],
        smoothing=0.043,
        prior_factor=parameters["prior_factor"],
        prior_coefficient=0.16,
    )
    np.testing.assert_array_equal(
        result.class_variances.standard_deviations,
        expected.class_variances.standard_deviations,
    )
    np.testing.assert_array_equal(
        result.image_variance.image_variance_weights,
        expected.image_variance.image_variance_weights,
    )
    np.testing.assert_array_equal(
        result.centroid_similarity.similarity,
        expected.centroid_similarity.similarity,
    )
    np.testing.assert_array_equal(
        result.membership_pass.memberships, expected.membership_pass.memberships
    )
    np.testing.assert_array_equal(
        result.membership_pass.denominator, expected.membership_pass.denominator
    )
    np.testing.assert_array_equal(result.membership_pass.objective, expected.membership_pass.objective)


def test_registered_startup_membership_rejects_nondefault_descriptor_order() -> None:
    run, _data, parameters, _segmentation = _alignment_run()
    run.atlas.structures = run.atlas.structures[:-1]
    with pytest.raises(ValueError, match="default 11-class atlas order"):
        compute_registered_toads_initial_membership(
            run,
            exponent_by_channel_class=parameters["exponent_by_channel_class"],
            modality_weights=parameters["modality_weights"],
            smoothing_priors=parameters["smoothing_priors"],
            atlas_shape_priors=parameters["atlas_shape_priors"],
            smoothing=0.043,
            prior_factor=parameters["prior_factor"],
            prior_coefficient=0.16,
        )


def test_registered_startup_membership_rejects_missing_run_state() -> None:
    with pytest.raises(TypeError, match="alignment_run must contain"):
        compute_registered_toads_initial_membership(
            object(),
            exponent_by_channel_class=np.ones((1, 11), dtype=np.float32),
            modality_weights=np.ones(1, dtype=np.float32),
            smoothing_priors=np.ones(11, dtype=np.float32),
            atlas_shape_priors=np.ones(11, dtype=np.float32),
            smoothing=0.01,
            prior_factor=0.2,
            prior_coefficient=0.1,
        )
