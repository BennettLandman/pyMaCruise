from __future__ import annotations

import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.toads_default_membership_registration_control import (
    continue_default_membership_through_registration_control,
)
from pymacruise.toads_default_membership_registered_iteration import (
    continue_default_membership_registered_iteration,
)
from pymacruise.toads_post_thinning_exact_iteration import (
    toads_post_thinning_exact_iteration_degree1,
)
from test_toads_post_thinning_exact_iteration import _continuation_state

ROOT = Path(__file__).parent
RAS_FIXTURE = ROOT / "data/synthetic/toads_post_registration_branch/membership_maps.nii.gz"
JAVA_REFERENCE = ROOT / "java_reference/ToadsP18cytRegisteredIterationReferenceHarness.java"


def _p18cyt_membership(next_thinning, image, centroids, field, closest):
    """Run the actual one-channel degree-one post-thinning membership pass."""
    shape = np.zeros((*image.shape, 11), dtype=np.float32)
    exponents = np.empty((1, 11), dtype=np.float32)
    exponents[0, 0] = 1.0
    exponents[0, 1:3] = 0.8
    exponents[0, 3:] = 1.0
    return toads_post_thinning_exact_iteration_degree1(
        next_thinning,
        image,
        centroids,
        field,
        category_exponents=np.asarray([0.8, 1.0, 1.0], dtype=np.float32),
        exponent_by_channel_class=exponents,
        modality_weights=np.asarray([0.73], dtype=np.float32),
        class_variance_weights=np.asarray(
            [[0.81, 0.93, 1.17, 0.76, 0.88, 1.02, 1.14, 0.69, 0.96, 1.21, 0.73]],
            dtype=np.float32,
        ),
        transformed_shapes=shape,
        closest_labels=closest,
        smoothing_priors=np.asarray(
            [0.65, 0.71, 0.84, 0.91, 1.02, 1.13, 0.76, 0.88, 1.05, 1.17, 0.82],
            dtype=np.float32,
        ),
        atlas_shape_priors=np.asarray(
            [0.10, 0.08, 0.09, 0.11, 0.12, 0.09, 0.10, 0.11, 0.07, 0.06, 0.07],
            dtype=np.float32,
        ),
        smoothing=0.043,
        prior_factor=0.64,
        prior_coefficient=0.16,
        prior_scale=0.05,
    )


def test_p18cyt_registration_state_composes_with_selected_branch_and_ras_case(tmp_path: Path) -> None:
    next_thinning, image, centroids, field, closest = _continuation_state()
    p18cyt = _p18cyt_membership(next_thinning, image, centroids, field, closest)
    membership = p18cyt.membership_pass
    assert image.ndim == 3
    assert membership.memberships.shape == (*image.shape, 11)

    # The selected source phantom is committed as RAS NIfTI. The synthetic
    # continuation expands each source voxel 3x on all axes; its derived
    # affine retains the same R/A/S directions and corresponding spacing.
    source = nib.load(str(RAS_FIXTURE))
    assert nib.aff2axcodes(source.affine) == ("R", "A", "S")
    expanded_affine = source.affine @ np.diag([1 / 3, 1 / 3, 1 / 3, 1])
    assert nib.aff2axcodes(expanded_affine) == ("R", "A", "S")
    assert nib.Nifti1Image(membership.memberships, expanded_affine).shape == membership.memberships.shape

    nx, ny, nz = image.shape
    gx, gy, gz = np.meshgrid(
        np.linspace(-1, 1, nx, dtype=np.float32),
        np.linspace(-1, 1, ny, dtype=np.float32),
        np.linspace(-1, 1, nz, dtype=np.float32),
        indexing="ij",
    )
    atlas = np.empty((*image.shape, 11), dtype=np.float32)
    for k in range(11):
        radius = (
            (gx - np.float32(-0.45 + 0.09 * k)) ** 2
            + 1.3 * (gy - np.float32(0.3 * np.sin(k))) ** 2
            + 0.8 * (gz - np.float32(0.25 * np.cos(k))) ** 2
        )
        atlas[..., k] = np.asarray(
            np.exp(-radius / np.float32(0.42)) * 0.75 + 0.05, dtype=np.float32
        )

    objective = np.float32(membership.objective)
    previous_energy = np.float32(objective * np.float32(1.01))
    settings = dict(
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        initial_transform=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        previous_energy=previous_energy,
        iteration=3,
        maximum_iterations=5,
        convergence_threshold=1e-4,
        do_growing=False,
        align_shapes=False,
    )
    registration_control = continue_default_membership_through_registration_control(
        membership, atlas, **settings,
    )

    # Inputs to the wrapper branch come from the preceding synthetic thinning
    # state. P18cyt supplies this pass's membership maps.
    branch_source = next_thinning.branch
    thinning = next_thinning.thinning
    relation_state = branch_source.relations
    assert relation_state is not None
    # The synthetic state builder expands the image/growing fields 3x in
    # spatial dimensions; carry its prior boundary-distance state through the
    # same nearest-voxel replication before the next selected refresh.
    boundary_distances = np.repeat(
        np.repeat(np.repeat(relation_state.boundary_distances, 3, axis=0), 3, axis=1),
        3,
        axis=2,
    )
    ordering = branch_source.growing.ordering.copy()
    template_labels = np.asarray([0, 4, 7, 9, 11, 13, 17, 19, 21, 23, 25], dtype=np.int8)
    previous_classification = np.full(image.shape, -1, dtype=np.int8)
    arguments = dict(
        previous_energy=previous_energy,
        iteration=3,
        maximum_iterations=5,
        convergence_threshold=1e-4,
        segmentation=thinning.segmentation,
        skeleton=thinning.skeleton,
        ordering=ordering,
        class_volumes=thinning.class_volumes,
        boundary_distances=boundary_distances,
        template_labels=template_labels,
        previous_classification=previous_classification,
        relation_limit=1,
        do_growing=False,
        initial_available=thinning.available,
    )
    result = continue_default_membership_registered_iteration(registration_control, **arguments)
    assert result.registration_control is registration_control
    assert result.branch.handoff.membership_pass is membership
    assert result.branch.control == registration_control.control
    assert result.branch.growing is not None and result.branch.relations is not None
    assert result.branch.classification is None
    assert result.stage_order == (
        "computeExactMemberships",
        "energyRatio",
        "convergenceCheck",
        "iterationLimitCheck",
        "propagateGrowing",
        "computeRelations",
        "computeBoundaryDistances",
        "computeClosestLabels",
    )

    compile_result = subprocess.run(
        ["javac", "-d", str(tmp_path), str(JAVA_REFERENCE)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert not compile_result.stderr
    bits = lambda value: str(int(np.asarray([value], dtype=np.float32).view(np.int32)[0]))
    java_result = subprocess.run(
        [
            "java", "-cp", str(tmp_path), "ToadsP18cytRegisteredIterationReferenceHarness",
            bits(previous_energy), bits(objective), bits(1e-4), "3", "5", "false", "false",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    reference = json.loads(java_result.stdout)
    assert reference["source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert reference["source"] == "AlgorithmMultipleToads selected loop after computeMemberships"
    assert reference["ratio_bits"] == int(registration_control.control.energy_ratio.view(np.int32))
    assert reference["stop"] is result.branch.control.stop is False
    assert reference["run_growing"] is result.branch.control.run_growing is True
    assert reference["order"] == [
        "computeExactMemberships", "energyRatio", "convergenceCheck", "iterationLimitCheck",
        "propagateGrowing", "computeRelations", "computeBoundaryDistances", "computeClosestLabels",
    ]

    growing = result.branch.growing
    relations = result.branch.relations
    assert growing is not None and relations is not None
    grown_image = nib.Nifti1Image(growing.segmentation, expanded_affine)
    assert nib.aff2axcodes(grown_image.affine) == ("R", "A", "S")
    assert relations.relation.relation_map.shape == image.shape
    assert relations.boundary_distances.shape == (*image.shape, 11)
    assert relations.closest_labels.shape == (*image.shape, 4)
