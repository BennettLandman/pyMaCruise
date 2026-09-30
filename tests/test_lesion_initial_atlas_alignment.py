import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_initial_atlas_alignment import (
    initialize_lesion_atlas_alignment,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_initial_atlas_alignment"


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
