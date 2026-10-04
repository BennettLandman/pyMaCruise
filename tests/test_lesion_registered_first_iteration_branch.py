from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_next_iteration_branch import continue_after_next_iteration_registration
from pymacruise.lesion_next_iteration_registration import register_after_next_gm_membership_update
from pymacruise.lesion_registered_first_iteration import (
    run_lesion_toads_registered_first_iteration,
)
from pymacruise.lesion_registered_first_iteration_branch import (
    run_lesion_toads_registered_first_iteration_branch,
)
from pymacruise.lesion_registered_initial_membership import (
    LesionRegisteredInitialMembership,
    run_lesion_toads_registered_initial_membership,
)
from pymacruise.lesion_standard13_voxel import STANDARD13_TYPES


FIXTURE = Path(__file__).parent / "data" / "synthetic" / "lesion_registered_initial_membership"


def _initial_state() -> tuple[
    LesionRegisteredInitialMembership,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    images = np.asanyarray(nib.load(FIXTURE / "images_ras.nii.gz").dataobj).astype(np.float32)
    labels = np.asanyarray(nib.load(FIXTURE / "constructor_labels_ras.nii.gz").dataobj).astype(np.int16)
    shapes = np.asanyarray(nib.load(FIXTURE / "transformed_shapes_ras.nii.gz").dataobj).astype(np.float32)
    centroids = np.stack(
        (np.linspace(0.2, 0.7, 13), np.linspace(0.25, 0.8, 13))
    ).astype(np.float32)
    lesion_centroids = np.asarray((0.57, 0.63), dtype=np.float32)
    black_centroids = np.asarray((0.34, 0.4), dtype=np.float32)
    initial = run_lesion_toads_registered_initial_membership(
        images,
        labels,
        centroids,
        lesion_centroids,
        black_centroids,
        shapes,
        labels.copy(),
        np.zeros((*labels.shape, 13), dtype=np.bool_),
        np.zeros(labels.shape, dtype=np.float32),
        iterations=1,
        first_thinning_limit=1000.0,
    )
    return initial, images, labels, shapes, centroids, np.stack((lesion_centroids, black_centroids))


def test_registered_first_iteration_composes_registration_and_terminal_branch() -> None:
    initial, images, labels, shapes, centroids, auxiliary = _initial_state()
    lesion_centroids, black_centroids = auxiliary
    shape_priors = np.ones(13, dtype=np.float32)
    smoothing_priors = np.ones(13, dtype=np.float32)
    transition = run_lesion_toads_registered_first_iteration(
        initial,
        images,
        centroids,
        lesion_centroids,
        black_centroids,
        shapes,
        shape_priors,
        smoothing_priors,
        centroid_options={
            "mode": "none",
            "modalities": ("T1_MPRAGE", "FLAIR"),
            "lesion_power_exponents": (0.4, 1.2),
            "tissue_power_exponents": np.ones((2, 13), dtype=np.float32),
        },
        channel_weight_options={"outlier_distance": 0.1, "norm_type": "Identity"},
        outlier_energy=0.1,
        lesion_power_exponent=1.2,
    )
    transform = np.zeros(6, dtype=np.float32)
    branch_inputs = {
        "image_resolutions_mm": (1.0, 1.0, 1.0),
        "template_resolutions_mm": (1.0, 1.0, 1.0),
        "object_types": STANDARD13_TYPES,
        "classification": labels,
        "template_labels": np.arange(13, dtype=np.int16),
        "lesion_class": 10,
        "iteration": 1,
        "n_iter_max": 99,
        "max_distance": 1.0e10,
        "algorithm_mode": "diffeomorphism",
        "relation_limit": 3,
        "align_shapes": True,
        "init_align_iter": 1,
    }
    result = run_lesion_toads_registered_first_iteration_branch(
        initial,
        transition,
        tuple(shapes[..., class_index] for class_index in range(13)),
        transform,
        **branch_inputs,
    )

    direct_registration = register_after_next_gm_membership_update(
        transition.membership,
        tuple(shapes[..., class_index] for class_index in range(13)),
        transform,
        image_resolutions_mm=branch_inputs["image_resolutions_mm"],
        template_resolutions_mm=branch_inputs["template_resolutions_mm"],
        align_shapes=True,
        init_align_iter=1,
    )
    thinning = initial.first_thinning.thinning.run
    direct_branch = continue_after_next_iteration_registration(
        direct_registration,
        segmentation=thinning.segmentation,
        skeleton=thinning.skeleton,
        ordering=thinning.ordering,
        class_volume=thinning.class_volume,
        object_types=STANDARD13_TYPES,
        classification=labels,
        template_labels=np.arange(13, dtype=np.int16),
        lesion_class=10,
        iteration=1,
        n_iter_max=99,
        max_distance=1.0e10,
        algorithm_mode="diffeomorphism",
        relation_limit=3,
    )

    np.testing.assert_array_equal(
        result.registration.registration.transform,
        direct_registration.registration.transform,
    )
    np.testing.assert_array_equal(
        result.branch.branch.final_classification.classification,
        direct_branch.branch.final_classification.classification,
    )
    assert result.branch.branch.stop is True
    assert result.branch.branch.final_classification.classification.shape == labels.shape
    assert result.stage_order.index("registerShapes") < result.stage_order.index(
        "convergence_check"
    )
    assert result.stage_order[-3:] == (
        "propagateGrowing",
        "updateClassificationFromLabels",
        "updateLesionClassifcation",
    )
    np.testing.assert_array_equal(transform, np.zeros(6, dtype=np.float32))
