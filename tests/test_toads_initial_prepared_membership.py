from __future__ import annotations

import numpy as np

from pymacruise.toads_default_atlas_membership import toads_default_atlas_membership_pass
from pymacruise.toads_initial_prepared_membership import (
    toads_prepared_initial_default_atlas_membership,
)
from pymacruise.toads_initial_relation_state import toads_initialize_relation_state
from test_toads_default_atlas_dispatch import _bits, _fixture


def test_startup_composition_attaches_then_prepares_then_updates_memberships() -> None:
    data, parameters, _cases, _reference = _fixture()
    segmentation = data["folded_ribbon"].astype(np.int8)
    channel_fields = [data["field"][..., channel] for channel in range(2)]
    initial_fields = [field.copy() for field in channel_fields]

    result = toads_prepared_initial_default_atlas_membership(
        image_values=data["images"],
        segmentation=segmentation,
        initial_memberships=data["initial_memberships"],
        transformed_shapes=data["transformed_shapes"],
        smoothing=0.043,
        prior_coefficient=0.16,
        **parameters,
        channel_fields=channel_fields,
        limit=2,
        psf_bias=data["psf_bias"],
    )

    assert result.stage_order == (
        "addInhomogeneityCorrection[channel=0]",
        "addInhomogeneityCorrection[channel=1]",
        "computeRelations",
        "computeBoundaryDistances",
        "computeClosestLabels",
        "fixed:computeVariances",
        "fixed:computeImageVariance+centroidSimilarity+transformedShapes",
        "computeExactMemberships",
    )
    assert result.field_state.use_field is True
    for channel, field in enumerate(channel_fields):
        assert result.field_state.field_by_channel[channel] is field
        np.testing.assert_array_equal(field, initial_fields[channel])

    expected_relations = toads_initialize_relation_state(
        segmentation, data["initial_memberships"], limit=2
    )
    np.testing.assert_array_equal(
        result.relation_state.relation.relation_map,
        expected_relations.relation.relation_map,
    )
    np.testing.assert_array_equal(
        result.relation_state.boundary_distances,
        expected_relations.boundary_distances,
    )
    np.testing.assert_array_equal(
        result.relation_state.closest_labels,
        expected_relations.closest_labels,
    )

    direct_membership = toads_default_atlas_membership_pass(
        image_values=data["images"],
        initial_memberships=data["initial_memberships"],
        transformed_shapes=data["transformed_shapes"],
        closest_labels=result.relation_state.closest_labels,
        smoothing=0.043,
        prior_coefficient=0.16,
        intensity_field=data["field"],
        psf_bias=data["psf_bias"],
        **parameters,
    )
    assert _bits(result.membership_pass.class_data_energy) == _bits(
        direct_membership.class_data_energy
    )
    assert _bits(result.membership_pass.memberships) == _bits(direct_membership.memberships)
    assert _bits(result.membership_pass.denominator) == _bits(direct_membership.denominator)
    assert _bits(result.membership_pass.objective) == _bits(direct_membership.objective)
