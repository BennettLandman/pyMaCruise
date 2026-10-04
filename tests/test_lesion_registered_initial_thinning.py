from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes

from pymacruise.lesion_atlas_priors import (
    load_lesion_atlas_priors,
    load_lesion_atlas_shape_priors,
)
from pymacruise.lesion_initial_alignment_bounds import compute_lesion_initial_alignment_bounds
from pymacruise.lesion_initial_atlas_alignment import (
    initialize_lesion_atlas_alignment_from_descriptor,
)
from pymacruise.lesion_initial_segmentation_state import run_loaded_template_first_thinning
from pymacruise.lesion_registered_initial_thinning import run_registered_lesion_initial_thinning
from pymacruise.lesion_registered_shapes import transform_registered_lesion_shape_priors
from pymacruise.lesion_transformed_template import transform_lesion_toads_template


ATLAS = (
    Path(__file__).parents[1]
    / "data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"
)


def test_registered_loaded_atlas_composes_transform_constructor_and_first_thinning() -> None:
    priors = load_lesion_atlas_priors(str(ATLAS), ("T1", "FLAIR"))
    shapes = load_lesion_atlas_shape_priors(str(ATLAS))
    center = tuple(size // 2 for size in next(iter(shapes.shape_maps.values())).shape)
    slices = tuple(slice(value - 5, value + 6) for value in center)
    atlas_patch = np.stack(
        [shapes.shape_maps[item.name][slices] for item in priors.structures], axis=-1
    ).astype(np.float32)
    dominant = np.argmax(atlas_patch, axis=-1)
    padded_images = priors.intensity_priors[:, dominant].copy()
    source_memberships = padded_images.copy()
    source_shape = atlas_patch.copy()

    alignment = initialize_lesion_atlas_alignment_from_descriptor(
        str(ATLAS),
        padded_images,
        np.zeros(6, dtype=np.float32),
        crop_bounds=((0, 9), (0, 9), (0, 9)),
        pad_border=1,
        image_resolutions_mm=(0.9375,) * 3,
        template_resolutions_mm=(0.9375,) * 3,
        init_align_iterations=0,
        pyramid_levels=2,
    )
    result = run_registered_lesion_initial_thinning(
        alignment,
        ATLAS,
        image_resolutions_mm=(0.9375,) * 3,
        template_resolutions_mm=(0.9375,) * 3,
        volume_ratio=0.0,
        limit=10.0,
    )

    topology = nib.load(str(result.topology.path))
    assert aff2axcodes(topology.affine) == ("R", "A", "S")
    assert result.transformed_template.shape == result.transformed_shapes.shape[:3] == (11, 11, 11)
    assert result.transformed_shapes.shape[-1] == 13
    assert result.stage_order == (
        "preprocess.initialMemberships",
        'atlas.initShapeRegistration(memberships,"rigid",initAlignIter,levels,1.0f)',
        "atlas.setTransform(preprocess.getTransform())",
        'atlas.updateShapeRegistration(memberships,"rigid",initAlignIter,levels,1.0f)',
        "atlas.registerShapesPyramid",
        "atlas.computeTransformedShapeBoundingBox",
        "atlas.computeTransformedTemplate.afterInitialRegistration",
        "atlas.getFastTransformedShape",
        "MultipleToadSegmentation.constructor:initializeClassState",
        "setGMflag",
        "initThinningLabels",
        "propagateThinning",
    )
    assert result.excluded_wrapper_stages == (
        "preprocess.uncropImages",
        "preprocess.findCroppingBoundaries",
        "preprocess.cropImages",
        "preprocess.estimateInitialTransform",
    )

    ordered_shapes = tuple(shapes.shape_maps[item.name] for item in priors.structures)
    direct_bounds = compute_lesion_initial_alignment_bounds(
        alignment.alignment,
        ordered_shapes,
        image_resolutions_mm=(0.9375,) * 3,
        template_resolutions_mm=(0.9375,) * 3,
    )
    direct_template = transform_lesion_toads_template(
        np.asarray(result.topology.image.dataobj),
        direct_bounds.transform,
        (11, 11, 11),
        image_resolution_mm=(0.9375,) * 3,
        template_resolution_mm=(0.9375,) * 3,
        background_label=priors.structures[0].label,
    )
    direct_shapes = transform_registered_lesion_shape_priors(
        direct_bounds,
        ordered_shapes,
        image_dimensions_xyz=(11, 11, 11),
        image_resolutions_mm=(0.9375,) * 3,
        template_resolutions_mm=(0.9375,) * 3,
    )
    lesion_class = next(
        index for index, item in enumerate(priors.structures) if item.topology == "lesion"
    )
    direct = run_loaded_template_first_thinning(
        alignment,
        direct_template,
        direct_shapes,
        np.zeros((11, 11, 11), dtype=np.float32),
        object_types=tuple(item.topology for item in priors.structures),
        lesion_class=lesion_class,
        volume_ratio=0.0,
        limit=10.0,
    )

    np.testing.assert_array_equal(result.registration_bounds.transform, direct_bounds.transform)
    np.testing.assert_array_equal(result.transformed_template, direct_template)
    np.testing.assert_array_equal(result.transformed_shapes, direct_shapes)
    np.testing.assert_array_equal(
        result.constructor_state.segmentation, direct.state.segmentation
    )
    np.testing.assert_array_equal(
        result.constructor_state.lesion_membership,
        np.zeros((11, 11, 11), dtype=np.float32),
    )
    np.testing.assert_array_equal(
        result.first_thinning.thinning.run.segmentation,
        direct.first_thinning.thinning.run.segmentation,
    )
    np.testing.assert_array_equal(
        result.first_thinning.thinning.run.skeleton,
        direct.first_thinning.thinning.run.skeleton,
    )
    np.testing.assert_array_equal(
        result.first_thinning.thinning.run.class_volume,
        direct.first_thinning.thinning.run.class_volume,
    )
    np.testing.assert_array_equal(padded_images, source_memberships)
    np.testing.assert_array_equal(atlas_patch, source_shape)
