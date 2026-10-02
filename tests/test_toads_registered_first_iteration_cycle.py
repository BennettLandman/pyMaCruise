from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_alignment_startup import register_toads_from_atlas_assets
from pymacruise.toads_registered_first_exact_iteration import (
    continue_registered_unity_startup_through_first_exact_iteration,
)
from pymacruise.toads_registered_first_iteration_cycle import (
    continue_registered_first_exact_iteration_through_next_thinning,
)
from pymacruise.toads_registered_unity_field_startup import (
    compute_registered_toads_unity_field_startup,
)
from pymacruise.toads_shape_volume import transform_toads_shape_priors
from pymacruise.toads_control import toads_energy_ratio
from pymacruise.toads_relations import toads_refresh_relation_state
from pymacruise.toads_thinning import toads_propagate_growing, toads_propagate_thinning


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "toads_registered_first_iteration_cycle"
ATLAS_DESCRIPTOR = ROOT / "data" / "synthetic" / "toads_registered_field_refit" / "atlas-11.txt"


def test_registered_first_exact_loop_reaches_rigid_branch_and_next_thinning() -> None:
    run = register_toads_from_atlas_assets(
        nib.load(str(FIXTURE / "t1.nii.gz")),
        nib.load(str(FIXTURE / "brain-mask.nii.gz")),
        ATLAS_DESCRIPTOR,
        levels=1,
        max_iterations=1,
        refinement_iterations=1,
        padding=2,
    )
    settings = dict(
        exponent_by_channel_class=np.full((1, 11), 2.0, dtype=np.float32),
        modality_weights=np.ones(1, dtype=np.float32),
        smoothing_priors=np.ones(11, dtype=np.float32),
        atlas_shape_priors=np.ones(11, dtype=np.float32),
        smoothing=np.float32(0.2 / 6.0 / 11.0),
        prior_factor=np.float32(0.2),
        prior_coefficient=np.float32(2.0 / 11.0),
        prior_scale=np.float32(0.05),
    )
    startup = compute_registered_toads_unity_field_startup(run, **settings)
    exact = continue_registered_unity_startup_through_first_exact_iteration(
        run,
        startup,
        category_exponents=np.asarray([[1.0, 1.0, 1.0]], dtype=np.float32),
        class_variance_weights=run.initial_variances.standard_deviations,
        **settings,
    )

    result = continue_registered_first_exact_iteration_through_next_thinning(
        run,
        exact,
        maximum_iterations=2,
        convergence_threshold=0.0,
        relation_limit=1,
        do_growing=False,
        subsample=2,
        max_registration_iterations=1,
        registration_energy_tolerance=0.0,
        minimum_lambda=1e-3,
    )

    assert result.exact_iteration is exact
    assert result.registration_control.handoff.membership_pass is exact.membership_pass
    assert result.registration_control.handoff.registration is not None
    expected_ratio = toads_energy_ratio(
        startup.statistics_membership.membership_pass.objective,
        exact.membership_pass.objective,
    )
    assert result.registration_control.control.energy_ratio == expected_ratio
    assert result.registered_iteration.branch.control.energy_ratio == expected_ratio
    assert not result.registered_iteration.branch.control.stop
    assert result.registered_iteration.branch.growing is not None
    assert result.registered_iteration.branch.relations is not None
    assert result.registered_iteration.branch.classification is None
    assert result.next_thinning.segmentation.shape == exact.membership_pass.memberships.shape[:3]
    assert result.next_transformed_shape_priors.shape == exact.membership_pass.memberships.shape
    assert result.stage_order[-11:] == (
        "iteration1.computeExactMemberships",
        "energyRatio",
        "atlas.registerShapes() -> registerAllShapes()",
        "convergenceCheck",
        "iterationLimitCheck",
        "propagateGrowing",
        "computeRelations",
        "computeBoundaryDistances",
        "computeClosestLabels",
        "iteration2.initThinningLabels",
        "iteration2.propagateThinning",
    )

    handoff = result.registration_control.handoff
    expected_shapes = transform_toads_shape_priors(
        np.stack(
            tuple(run.atlas_data.shapes[item.name] for item in run.atlas.structures),
            axis=-1,
        ),
        image_dimensions_xyz=exact.membership_pass.memberships.shape[:3],
        image_resolutions_mm=run.alignment.setup.signed_image_resolutions_mm,
        template_resolutions_mm=run.alignment.setup.template_resolutions_mm,
        transform=tuple(float(value) for value in handoff.transform),
        class_bounds=handoff.shape_bounds,
    )
    np.testing.assert_array_equal(result.next_transformed_shape_priors, expected_shapes)
    direct_growing = toads_propagate_growing(
        startup.thinning.segmentation,
        startup.thinning.skeleton,
        exact.membership_pass.memberships,
        np.zeros(startup.thinning.segmentation.shape, dtype=np.float32),
        startup.thinning.class_volumes,
        initial_available=startup.thinning.available,
    )
    np.testing.assert_array_equal(
        result.registered_iteration.branch.growing.segmentation,
        direct_growing.segmentation,
    )
    np.testing.assert_array_equal(
        result.registered_iteration.branch.growing.ordering,
        direct_growing.ordering,
    )
    direct_relations = toads_refresh_relation_state(
        direct_growing.segmentation,
        exact.membership_pass.memberships,
        run.relation_state.boundary_distances,
        1,
    )
    np.testing.assert_array_equal(
        result.registered_iteration.branch.relations.relation.relation_map,
        direct_relations.relation.relation_map,
    )
    np.testing.assert_array_equal(
        result.registered_iteration.branch.relations.closest_labels,
        direct_relations.closest_labels,
    )
    expected_thinning = toads_propagate_thinning(
        direct_growing.segmentation,
        exact.membership_pass.memberships,
        direct_relations.relation.class_adjacency,
        expected_shapes,
        volume_ratio=0.5,
        limit=3.0,
    )
    np.testing.assert_array_equal(result.next_thinning.segmentation, expected_thinning.segmentation)
    np.testing.assert_array_equal(result.next_thinning.skeleton, expected_thinning.skeleton)
    assert result.next_thinning.processed_items == expected_thinning.processed_items
