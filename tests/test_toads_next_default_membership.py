from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_default_atlas_membership import toads_default_atlas_membership_pass
from pymacruise.toads_default_registration_handoff import toads_default_atlas_registration_handoff
from pymacruise.toads_next_default_membership import toads_next_default_membership_pass
from pymacruise.toads_post_registration import toads_post_registration_branch
from test_toads_default_atlas_dispatch import _fixture as default_fixture
from test_toads_default_registration_handoff import _arguments as registration_arguments
from test_toads_post_registration_branch import STATE as CZI_STATE, _arrays as czi_arrays


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_next_default_membership"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
SUBJECT_SHAPE = (25, 25, 25)


def _bits(value: np.ndarray | np.float32) -> list[int]:
    return np.asarray(value, dtype=np.float32).ravel().view(np.int32).tolist()


def _from_bits(bits: list[int], shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(bits, dtype=np.int32).view(np.float32).reshape(shape)


def _repeat_to_shape(values: np.ndarray) -> np.ndarray:
    indices = tuple(np.arange(size) % values.shape[axis] for axis, size in enumerate(SUBJECT_SHAPE))
    return values[np.ix_(*indices)]


def _case():
    data, parameters, _, _, options, subject_affine = registration_arguments()
    source_data, _, _, _ = default_fixture()
    data = {name: _repeat_to_shape(value) for name, value in source_data.items()}
    options = dict(options)
    for option_name, data_name in (
        ("image_values", "images"), ("initial_memberships", "initial_memberships"),
        ("transformed_shapes", "transformed_shapes"), ("closest_labels", "closest_labels"),
        ("intensity_field", "field"), ("psf_bias", "psf_bias"),
    ):
        options[option_name] = data[data_name]
    template = nib.load(str(ROOT / "template_atlas_shapes.nii.gz"))
    assert nib.aff2axcodes(template.affine) == ("R", "A", "S")
    options["atlas_shapes"] = np.asarray(template.dataobj, dtype=np.float32)

    handoff = toads_default_atlas_registration_handoff(**options)
    state = {
        name: _repeat_to_shape(value)
        if value.ndim >= 3 and value.shape[:3] == (9, 8, 7)
        else value
        for name, value in czi_arrays().items()
    }
    state["class_volumes"] = np.sum(
        np.where(
            np.eye(11, dtype=bool)[state["segmentation"]],
            handoff.membership_pass.memberships,
            np.float32(0.0),
        )[1:-1, 1:-1, 1:-1],
        axis=(0, 1, 2),
        dtype=np.float32,
    )
    branch = toads_post_registration_branch(
        handoff,
        previous_energy=handoff.membership_pass.objective,
        segmentation=state["segmentation"],
        skeleton=state["skeleton"],
        ordering=state["ordering"],
        class_volumes=state["class_volumes"],
        boundary_distances=state["boundary_distances"],
        template_labels=np.asarray(CZI_STATE["template_labels"], dtype=np.int8),
        previous_classification=state["previous_classification"],
        iteration=1,
        maximum_iterations=5,
        convergence_threshold=0.0,
        relation_limit=1,
        do_growing=False,
    )
    kwargs = {
        "image_values": data["images"],
        "centroids": parameters["centroids"],
        "exponent_by_channel_class": parameters["exponent_by_channel_class"],
        "modality_weights": parameters["modality_weights"],
        "class_variance_weights": parameters["class_variance_weights"],
        "image_variances": parameters["image_variances"],
        "smoothing_priors": parameters["smoothing_priors"],
        "centroid_similarity": parameters["centroid_similarity"],
        "atlas_shape_priors": parameters["atlas_shape_priors"],
        "atlas_shapes": options["atlas_shapes"],
        "image_resolutions_mm": options["image_resolutions_mm"],
        "template_resolutions_mm": options["template_resolutions_mm"],
        "smoothing": options["smoothing"],
        "prior_factor": options["prior_factor"],
        "prior_coefficient": options["prior_coefficient"],
        "volume_ratio": 0.2,
        "thinning_limit": 1.0,
        "intensity_field": data["field"],
        "psf_bias": data["psf_bias"],
    }
    return data, options, state, handoff, branch, kwargs


def _java_popped(result) -> list[list[int]]:
    return [
        [int(np.asarray([item.value], dtype=np.float32).view(np.int32)[0]), item.x, item.y, item.z, item.label]
        for item in result.thinning.processed_items
    ]


def test_actual_registered_continuation_thins_then_matches_java_membership() -> None:
    assert REFERENCE["algorithm_source_sha256"] == "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert REFERENCE["wrapper_source_sha256"] == "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"
    assert REFERENCE["powertable_source_sha256"] == "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678"
    assert REFERENCE["shape_xyz_channels"] == [25, 25, 25, 11]
    assert REFERENCE["template_xyz_channels"] == [11, 10, 9, 11]
    assert REFERENCE["stage_order"] == [
        "continuing_post_registration_result", "propagateThinning", "computeExactMemberships"
    ]

    image = nib.load(str(ROOT / "subject_images.nii.gz"))
    ribbon = nib.load(str(ROOT / "folded_ribbon.nii.gz"))
    template = nib.load(str(ROOT / "template_atlas_shapes.nii.gz"))
    assert image.shape == (25, 25, 25, 2)
    assert nib.aff2axcodes(image.affine) == nib.aff2axcodes(ribbon.affine) == nib.aff2axcodes(template.affine) == ("R", "A", "S")
    assert template.shape[:3] == (11, 10, 9) and template.shape[:3] != image.shape[:3]
    ribbon_values = np.asarray(ribbon.dataobj)
    assert ribbon_values.sum() > 20
    assert not np.array_equal(ribbon_values, ribbon_values[::-1, :, :])

    data, options, state, handoff, branch, kwargs = _case()
    assert branch.control.stop is False and branch.growing is not None and branch.relations is not None
    assert np.asarray(handoff.transform, dtype=np.float32).view(np.uint32).astype(int).tolist() == REFERENCE["registered_transform_bits"]
    assert not np.array_equal(handoff.transform, options["initial_transform"])
    assert branch.growing.segmentation.shape == SUBJECT_SHAPE

    # Snapshot every caller-owned array used by this composition.
    before = {name: value.copy() for name, value in data.items()}
    options_before = {name: value.copy() for name, value in options.items() if isinstance(value, np.ndarray)}
    segmentation_before = branch.growing.segmentation.copy()
    memberships_before = handoff.membership_pass.memberships.copy()
    relations_before = branch.relations.relation.class_adjacency.copy()
    closest_before = branch.relations.closest_labels.copy()

    result = toads_next_default_membership_pass(branch, **kwargs)
    assert result.stage_order == tuple(REFERENCE["stage_order"])
    assert len(result.thinning.processed_items) == REFERENCE["queue_pop_count"]

    # Java's real 100-voxel minimum permits transitions for this 25^3 phantom;
    # compare source heap trace and the resulting segmentation/skeleton state.
    java_queue = REFERENCE["queue"]
    assert len(result.thinning.processed_items) == len(java_queue["popped"]) > 1000
    assert _java_popped(result) == java_queue["popped"]
    assert _bits(result.thinning.class_volumes) == java_queue["volumes_bits"]
    assert _bits(result.thinning.minimum_class_volumes) == java_queue["minimum_bits"]
    assert result.thinning.boundary_counts.tolist() == java_queue["boundary_counts"]
    assert result.thinning.segmentation.ravel().tolist() == java_queue["segmentation"]
    # Java and Python both carry the post-growing skeleton into this pass.
    # The source initializer rewrites the interior and retains the border.
    java_skeleton = np.asarray(java_queue["skeleton"], dtype=np.int16).reshape(
        result.thinning.skeleton.shape
    )
    np.testing.assert_array_equal(result.thinning.skeleton, java_skeleton)
    assert result.thinning.available.ravel().tolist() == java_queue["available"]
    assert np.count_nonzero(result.thinning.segmentation != branch.growing.segmentation) > 0
    assert np.count_nonzero(result.thinning.segmentation != branch.growing.segmentation) == REFERENCE["segmentation_label_update_count"] == 978
    assert np.count_nonzero(result.thinning.skeleton == -1) > np.count_nonzero(
        branch.growing.skeleton == -1
    )

    for actual, key in (
        (result.membership_pass.class_data_energy, "class_energy_bits"),
        (result.membership_pass.memberships, "membership_bits"),
        (result.membership_pass.denominator, "denominator_bits"),
    ):
        assert _bits(actual) == REFERENCE[key]
    assert _bits(result.membership_pass.objective) == [REFERENCE["objective_bits"]]
    assert result.transformed_shapes.shape == SUBJECT_SHAPE + (11,)
    np.testing.assert_array_equal(result.transformed_shapes, nib.load(str(ROOT / "current_transform_shapes.nii.gz")).get_fdata(dtype=np.float32))

    # Membership borders remain the original maps and denominator borders are zero.
    for axis in range(3):
        for side in (0, -1):
            selector = [slice(None)] * 4
            selector[axis] = side
            np.testing.assert_array_equal(result.membership_pass.memberships[tuple(selector)], memberships_before[tuple(selector)])
            den_selector = tuple(selector[:3])
            assert not np.any(result.membership_pass.denominator[den_selector])
            seg_selector = [slice(None)] * 3
            seg_selector[axis] = side
            np.testing.assert_array_equal(result.thinning.segmentation[tuple(seg_selector)], branch.growing.segmentation[tuple(seg_selector)])

    # The refreshed closest-label map is consumed materially by the next sweep.
    stale = toads_default_atlas_membership_pass(
        image_values=kwargs["image_values"],
        initial_memberships=memberships_before,
        centroids=kwargs["centroids"],
        exponent_by_channel_class=kwargs["exponent_by_channel_class"],
        modality_weights=kwargs["modality_weights"],
        class_variance_weights=kwargs["class_variance_weights"],
        image_variances=kwargs["image_variances"],
        transformed_shapes=result.transformed_shapes,
        closest_labels=options["closest_labels"],
        smoothing_priors=kwargs["smoothing_priors"],
        centroid_similarity=kwargs["centroid_similarity"],
        atlas_shape_priors=kwargs["atlas_shape_priors"],
        smoothing=kwargs["smoothing"],
        prior_factor=kwargs["prior_factor"],
        prior_coefficient=kwargs["prior_coefficient"],
        intensity_field=kwargs["intensity_field"],
        psf_bias=kwargs["psf_bias"],
    )
    assert not np.array_equal(stale.memberships, result.membership_pass.memberships)
    assert stale.objective != result.membership_pass.objective

    for name, original in before.items():
        np.testing.assert_array_equal(data[name], original)
    for name, original in options_before.items():
        np.testing.assert_array_equal(options[name], original)
    np.testing.assert_array_equal(branch.growing.segmentation, segmentation_before)
    np.testing.assert_array_equal(handoff.membership_pass.memberships, memberships_before)
    np.testing.assert_array_equal(branch.relations.relation.class_adjacency, relations_before)
    np.testing.assert_array_equal(branch.relations.closest_labels, closest_before)
