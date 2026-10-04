from __future__ import annotations

from dataclasses import replace
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_atlas_priors import (
    LesionAtlasShapePriors,
    LesionAtlasTopology,
    load_lesion_atlas_priors,
)
from pymacruise.lesion_initial_atlas_alignment import (
    LesionLoadedAtlasInitialAlignment,
    initialize_lesion_atlas_alignment,
)
from pymacruise.lesion_initial_alignment_bounds import compute_lesion_initial_alignment_bounds
from pymacruise.lesion_initial_segmentation_state import initialize_lesion_segmentation_state
from pymacruise.lesion_preprocessing_crop import crop_lesion_toads_preprocessing_inputs
from pymacruise.lesion_preprocessing_uncrop import uncrop_lesion_toads_preprocessing_inputs
from pymacruise.lesion_registered_recrop import run_registered_lesion_recrop
from pymacruise.lesion_registered_shapes import transform_registered_lesion_shape_priors
from pymacruise.lesion_transformed_template import transform_lesion_toads_template


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_preprocessing_crop_rigid"
ATLAS = (
    Path(__file__).parents[1]
    / "data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"
)


def test_registered_recrop_uses_preprocessing_transform_and_new_atlas_grid() -> None:
    metadata = json.loads((ROOT / "inputs.json").read_text())
    image = nib.load(str(ROOT / metadata["image_nifti"]))
    source_topology_image = nib.load(str(ROOT / metadata["topology_nifti"]))
    channels = np.asarray(image.dataobj, dtype=np.float32).transpose(3, 0, 1, 2)
    topology = np.asarray(source_topology_image.dataobj).copy()
    # Remap the crop fixture's two non-background labels to the descriptor's
    # source labels so constructor state has no undeclared class values.
    topology[topology == 2] = 5
    topology[topology == 3] = 6
    topology_image = nib.Nifti1Image(topology, source_topology_image.affine)
    assert nib.aff2axcodes(image.affine) == nib.aff2axcodes(topology_image.affine) == ("R", "A", "S")
    channels_before = channels.copy()
    topology_before = topology.copy()

    resolutions = tuple(metadata["image_resolution_mm"])
    template_resolutions = tuple(metadata["template_resolution_mm"])
    signs = tuple(metadata["orientation_signs"])
    padding = 2
    initial_crop = crop_lesion_toads_preprocessing_inputs(
        channels,
        topology,
        np.asarray(metadata["initial_transform"], dtype=np.float32),
        image_resolution_mm=resolutions,
        template_resolution_mm=template_resolutions,
        orientation_signs=signs,
        padding=padding,
    )

    priors = load_lesion_atlas_priors(str(ATLAS), ("T1_MPRAGE", "FLAIR"))
    # The atlas shape arrays are a small asymmetric synthetic patch so this
    # test exercises the entire post-registration grid without materializing
    # the 172x216x164 production atlas maps.
    shape_maps = {}
    for class_index, structure in enumerate(priors.structures):
        shape = np.zeros(topology.shape, dtype=np.float32)
        x0 = 1 + class_index % 3
        y0 = 1 + (class_index // 3) % 3
        z0 = 1 + (class_index // 9)
        shape[x0 : x0 + 5, y0 : y0 + 5, z0 : z0 + 5] = np.float32(
            0.4 + class_index / 50.0
        )
        shape_maps[structure.name] = shape
    shape_priors = LesionAtlasShapePriors(
        structures=priors.structures,
        shape_maps=shape_maps,
        files={name: ROOT / "synthetic-shape.nii.gz" for name in shape_maps},
    )
    shape_maps_before = {name: values.copy() for name, values in shape_maps.items()}
    inclusive = initial_crop.crop.bounds
    exclusive_membership_bounds = tuple(
        (inclusive[2 * axis], inclusive[2 * axis + 1] + 1) for axis in range(3)
    )
    alignment_result = initialize_lesion_atlas_alignment(
        initial_crop.crop.images,
        priors.intensity_priors,
        tuple(shape_maps[item.name] for item in priors.structures),
        initial_crop.crop.transform,
        crop_bounds=exclusive_membership_bounds,
        pad_border=padding,
        image_resolutions_mm=resolutions,
        template_resolutions_mm=template_resolutions,
        init_align_iterations=0,
        pyramid_levels=1,
    )
    # A fixed atlas-optimizer result is intentionally distinct from the
    # preprocessor transform used for the recrop location.
    atlas_transform = alignment_result.registration.transform.copy()
    atlas_transform[3:] += np.asarray((2.0, -3.0, 1.5), dtype=np.float32)
    alignment_result = replace(
        alignment_result,
        registration=replace(alignment_result.registration, transform=atlas_transform),
    )
    alignment = LesionLoadedAtlasInitialAlignment(priors, shape_priors, alignment_result)
    atlas_transform_before = atlas_transform.copy()
    topology_result = LesionAtlasTopology(ROOT / metadata["topology_nifti"], topology_image)
    reference_path = ROOT / "post_registration_recrop_java_reference.json"
    reference = json.loads(reference_path.read_text())
    assert hashlib.sha256(reference_path.read_bytes()).hexdigest() == (
        "85c670e8bed7aca159148b52869721d2df89642b0d98e8af1ee9aaae85454e39"
    )
    harness = Path(__file__).parent / "java_reference" / "LesionPostRegistrationRecropReferenceHarness.java"
    assert hashlib.sha256(harness.read_bytes()).hexdigest() == (
        "57b6cce2d236de401fe3f92291f13656f9e6617f5816dab89756786f8f6594f3"
    )

    original_dimensions = tuple(int(value) for value in channels.shape[1:])
    expected_uncropped, expected_uncropped_transform = uncrop_lesion_toads_preprocessing_inputs(
        initial_crop.crop.images,
        original_dimensions,
        initial_crop.crop.bounds,
        padding,
        initial_crop.crop.transform,
        image_resolution_mm=resolutions,
        orientation_signs=signs,
    )
    expected_recrop = crop_lesion_toads_preprocessing_inputs(
        expected_uncropped,
        topology,
        expected_uncropped_transform,
        image_resolution_mm=resolutions,
        template_resolution_mm=template_resolutions,
        orientation_signs=signs,
        padding=padding,
    )
    incorrectly_atlas_seeded_recrop = crop_lesion_toads_preprocessing_inputs(
        expected_uncropped,
        topology,
        atlas_transform,
        image_resolution_mm=resolutions,
        template_resolution_mm=template_resolutions,
        orientation_signs=signs,
        padding=padding,
    )

    result = run_registered_lesion_recrop(
        alignment,
        initial_crop,
        topology_result,
        original_dimensions,
        image_resolutions_mm=resolutions,
        template_resolutions_mm=template_resolutions,
        orientation_signs=signs,
    )

    assert list(initial_crop.crop.bounds) == reference["first_bounds"]
    assert list(initial_crop.crop.images.shape[1:]) == reference["first_dimensions"]
    assert list(initial_crop.crop.transform.view(np.int32)) == reference["first_transform_bits"]
    first_crop_bytes = np.asarray(initial_crop.crop.images, dtype="<f4", order="C").tobytes()
    assert hashlib.sha256(first_crop_bytes).hexdigest() == reference["first_crop_sha256"]
    assert list(result.uncropped_preprocessing_transform.view(np.int32)) == reference[
        "restored_transform_bits"
    ]
    assert list(result.preprocessing_crop.crop.bounds) == reference["second_bounds"]
    assert list(result.preprocessing_crop.crop.images.shape[1:]) == reference["second_dimensions"]
    assert list(result.preprocessing_crop.crop.transform.view(np.int32)) == reference[
        "second_transform_bits"
    ]
    second_crop_bytes = np.asarray(
        result.preprocessing_crop.crop.images, dtype="<f4", order="C"
    ).tobytes()
    assert hashlib.sha256(second_crop_bytes).hexdigest() == reference["second_crop_sha256"]
    assert result.preprocessing_crop.crop.bounds == expected_recrop.crop.bounds
    assert result.preprocessing_crop.crop.bounds != incorrectly_atlas_seeded_recrop.crop.bounds
    assert result.preprocessing_crop.crop.images.shape == expected_recrop.crop.images.shape
    np.testing.assert_array_equal(
        result.preprocessing_crop.crop.images.view(np.int32),
        expected_recrop.crop.images.view(np.int32),
    )
    np.testing.assert_array_equal(
        result.preprocessing_crop.crop.transform.view(np.int32),
        expected_recrop.crop.transform.view(np.int32),
    )
    np.testing.assert_array_equal(
        result.uncropped_images.view(np.int32), expected_uncropped.view(np.int32)
    )
    np.testing.assert_array_equal(
        result.uncropped_preprocessing_transform.view(np.int32),
        expected_uncropped_transform.view(np.int32),
    )
    recrop_dimensions = result.preprocessing_crop.crop.images.shape[1:]
    assert result.transformed_template.shape == recrop_dimensions
    assert result.transformed_shapes.shape == (*recrop_dimensions, 13)
    assert result.constructor_state.segmentation.shape == recrop_dimensions
    direct_bounds = compute_lesion_initial_alignment_bounds(
        alignment_result,
        tuple(shape_maps[item.name] for item in priors.structures),
        image_dimensions_xyz=recrop_dimensions,
        image_resolutions_mm=resolutions,
        template_resolutions_mm=template_resolutions,
    )
    assert direct_bounds.shape_bounds == result.registration_bounds.shape_bounds
    direct_template = transform_lesion_toads_template(
        topology,
        atlas_transform,
        recrop_dimensions,
        image_resolution_mm=resolutions,
        template_resolution_mm=template_resolutions,
        background_label=priors.structures[0].label,
    )
    direct_shapes = transform_registered_lesion_shape_priors(
        direct_bounds,
        tuple(shape_maps[item.name] for item in priors.structures),
        image_dimensions_xyz=recrop_dimensions,
        image_resolutions_mm=resolutions,
        template_resolutions_mm=template_resolutions,
    )
    np.testing.assert_array_equal(result.transformed_template, direct_template)
    np.testing.assert_array_equal(result.transformed_shapes, direct_shapes)
    direct_constructor = initialize_lesion_segmentation_state(
        direct_template, tuple(item.label for item in priors.structures)
    )
    np.testing.assert_array_equal(
        result.constructor_state.segmentation, direct_constructor.segmentation
    )
    np.testing.assert_array_equal(
        result.registration_bounds.transform.view(np.int32), atlas_transform.view(np.int32)
    )
    assert not np.array_equal(
        result.registration_bounds.transform.view(np.int32),
        result.preprocessing_crop.crop.transform.view(np.int32),
    )
    assert result.stage_order[-5:] == (
        "preprocess.cropImages:preprocessingTransform",
        "atlas.setImageInfo:recroppedDimensions",
        "atlas.computeTransformedTemplate:optimizedAtlasTransform",
        "preprocess.updateTransformedTemplate",
        "MultipleToadSegmentation.constructor:recroppedGrid",
    )
    assert result.excluded_wrapper_stages == (
        "preprocess.normalizeImages",
        "preprocess.estimateInitialTransform",
        "preprocess.updateTransformedTemplate centroid refresh",
        "LesionToadSegmentation.computeMemberships",
        "initThinningLabels",
        "propagateThinning",
    )
    np.testing.assert_array_equal(channels, channels_before)
    np.testing.assert_array_equal(topology, topology_before)
    np.testing.assert_array_equal(alignment_result.registration.transform, atlas_transform_before)
    for name, values in shape_maps.items():
        np.testing.assert_array_equal(values, shape_maps_before[name])

    left_affine = topology_image.affine.copy()
    left_affine[0, :3] *= -1.0
    left_affine[0, 3] += (topology.shape[0] - 1) * resolutions[0]
    left_topology = nib.Nifti1Image(np.flip(topology, axis=0), left_affine)
    assert nib.aff2axcodes(left_topology.affine) == ("L", "A", "S")
    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        left_result = run_registered_lesion_recrop(
            alignment,
            initial_crop,
            LesionAtlasTopology(ROOT / metadata["topology_nifti"], left_topology),
            original_dimensions,
            image_resolutions_mm=resolutions,
            template_resolutions_mm=template_resolutions,
            orientation_signs=signs,
        )
    assert nib.aff2axcodes(left_result.topology.image.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(left_result.topology.image.dataobj, topology)
    assert left_result.preprocessing_crop.crop.bounds == result.preprocessing_crop.crop.bounds
