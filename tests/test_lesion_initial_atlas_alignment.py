import json
import hashlib
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_initial_atlas_alignment import (
    initialize_lesion_atlas_alignment,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_initial_atlas_alignment"
HARNESS = Path(__file__).parent / "java_reference" / "LesionInitialAtlasAlignmentReferenceHarness.java"
ALGORITHM_SHA = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
ATLAS_SHA = "718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e"
TRANSFORM_SHA = "22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac"
HARNESS_SHA = "7c2a5213d29476c74e522b95698587bff5dc7b1a0a68a858ad87ac915166d307"


def _read_fixture():
    metadata = json.loads((ROOT / "inputs.json").read_text())
    image_obj = nib.load(str(ROOT / metadata["input_nifti"]))
    assert nib.aff2axcodes(image_obj.affine) == ("R", "A", "S")
    channels_last = np.asarray(image_obj.dataobj, dtype=np.float32)
    images = np.transpose(channels_last, (3, 0, 1, 2))
    shapes = []
    for filename in metadata["shape_niftis"]:
        shape_obj = nib.load(str(ROOT / filename))
        assert nib.aff2axcodes(shape_obj.affine) == ("R", "A", "S")
        shapes.append(np.asarray(shape_obj.dataobj, dtype=np.float32))
    centroids = np.asarray(metadata["centroids_channels_classes"], dtype=np.float32)
    transform = np.asarray(metadata["preprocessing_transform"], dtype=np.float32)
    bounds = tuple(tuple(pair) for pair in metadata["crop_bounds_source_coordinates"])
    return metadata, images, centroids, shapes, transform, bounds


def test_initial_memberships_seed_selected_rigid_atlas_pyramid_in_source_order():
    metadata, images, centroids, shapes, transform, bounds = _read_fixture()
    reference = json.loads((ROOT / "java_reference.json").read_text())
    input_copies = (images.copy(), centroids.copy(), transform.copy(), [s.copy() for s in shapes])

    result = initialize_lesion_atlas_alignment(
        images,
        centroids,
        shapes,
        transform,
        crop_bounds=bounds,
        pad_border=metadata["pad_border"],
        image_resolutions_mm=tuple(metadata["image_resolutions_mm"]),
        template_resolutions_mm=tuple(metadata["template_resolutions_mm"]),
        init_align_iterations=metadata["init_align_iterations"],
        pyramid_levels=metadata["pyramid_levels"],
    )

    assert reference["algorithm_source_sha256"] == ALGORITHM_SHA
    assert reference["atlas_source_sha256"] == ATLAS_SHA
    assert reference["transform_source_sha256"] == TRANSFORM_SHA
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA
    assert result.call_order == tuple(reference["call_order"])
    assert result.initialization_iterations == reference["init_align_iterations"] == 50
    assert result.pyramid_levels == reference["levels"] == 4
    assert result.registration_scale_factor == reference["subsample"] == 2
    assert result.coupling_factor.view(np.int32) == reference["coupling_factor_bits"]
    assert tuple(level.scale for level in result.registration.levels) == tuple(
        reference["pyramid_scales"]
    ) == (8, 4, 2)
    assert list(result.preprocessing_transform.view(np.int32)) == reference[
        "preprocessing_transform_bits"
    ]
    assert list(result.registration.levels[0].initial_transform.view(np.int32)) == reference[
        "post_update_transform_bits"
    ]
    # P39o's selected pyramid is reused directly: the preprocessing seed is
    # the first scale's input, and each finer scale starts from the preceding
    # optimizer result without resetting the transform.
    carried_transform = result.preprocessing_transform
    for level in result.registration.levels:
        np.testing.assert_array_equal(level.initial_transform, carried_transform)
        assert level.optimization.stop_reason
        for step in level.optimization.iterations:
            assert step.trial_transforms.shape[1:] == (6,)
            assert step.trial_energies.shape == (step.trial_transforms.shape[0],)
        carried_transform = level.optimization.transform
    np.testing.assert_array_equal(result.registration.transform, carried_transform)
    expected_memberships = np.asarray(reference["membership_bits"], dtype=np.int32)
    expected_memberships = expected_memberships.view(np.float32).reshape(
        images.shape[1:] + (centroids.shape[1],)
    )
    np.testing.assert_array_equal(
        result.memberships.memberships.view(np.int32), expected_memberships.view(np.int32)
    )
    for before, after in zip(input_copies[:3], (images, centroids, transform), strict=True):
        np.testing.assert_array_equal(after, before)
    for before, after in zip(input_copies[3], shapes, strict=True):
        np.testing.assert_array_equal(after, before)
