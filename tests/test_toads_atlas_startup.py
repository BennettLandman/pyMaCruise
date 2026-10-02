from __future__ import annotations

import inspect

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_atlas import load_toads_lesion_priors
from pymacruise.toads_alignment_startup import register_toads_from_atlas_assets


def _write_atlas(root):
    shape = (15, 13, 11)
    x, y, z = np.indices(shape, dtype=np.float32)
    topology = np.ones(shape, dtype=np.uint8)
    region = (x - 7) ** 2 + 1.2 * (y - 6) ** 2 + 0.8 * (z - 5) ** 2 < 20
    topology[region] = 15
    background = np.ones(shape, dtype=np.float32)
    gm = np.exp(-((x - 7) ** 2 + (y - 6) ** 2 + (z - 5) ** 2) / 18).astype(np.float32)
    affine = np.eye(4)
    for name, values in (
        ("topology", topology),
        ("background", background),
        ("gm", gm),
    ):
        nib.save(nib.Nifti1Image(values, affine), root / f"{name}.nii.gz")

    descriptor = root / "atlas.txt"
    descriptor.write_text(
        "\n".join(
            (
                "Structure Atlas File (edit at your own risks)",
                "",
                "Structures:\t2",
                "Background\t1\tmask",
                "CerebralGM\t15\toptimized_gm",
                "",
                "Topology Atlas:",
                "File:\ttopology.nii.gz",
                "Dimensions:\t15\t13\t11",
                "Resolutions:\t1\t1\t1",
                "",
                "Shape Atlas:\t1",
                "Dimensions:\t15\t13\t11",
                "Resolutions:\t1\t1\t1",
                "Orientations:\t1\t1\t6\t3",
                "Structure: Background",
                "File:\tbackground.nii.gz",
                "Dimensions:\t15\t13\t11",
                "Resolutions:\t1\t1\t1",
                "Center:\t0\t0\t0",
                "Structure: CerebralGM",
                "File:\tgm.nii.gz",
                "Dimensions:\t15\t13\t11",
                "Resolutions:\t1\t1\t1",
                "Center:\t0\t0\t0",
                "",
                "Intensity Atlas:\t0",
                "T1_MPRAGE\t0.0\t0.8",
            )
        )
        + "\n",
        encoding="utf-8",
    )
    return descriptor


def _synthetic_subject():
    shape = (28, 27, 25)
    x, y, z = np.indices(shape, dtype=np.float32)
    mask = np.zeros(shape, dtype=np.uint8)
    mask[5:23, 4:23, 3:22] = 1
    values = np.zeros(shape, dtype=np.float32)
    inside = mask > 0
    values[inside] = (
        20.0
        + x[inside] * 0.7
        + y[inside] * 0.31
        + z[inside] * 0.19
        + ((x[inside] + 2 * y[inside] + z[inside]) % 7) * 0.13
    )
    affine = np.diag((0.8, 1.2, 1.5, 1.0))
    return nib.Nifti1Image(values, affine), nib.Nifti1Image(mask, affine)


def test_atlas_assets_compose_into_ras_startup_and_registration(tmp_path) -> None:
    parameters = inspect.signature(register_toads_from_atlas_assets).parameters
    assert parameters["levels"].default == 4
    assert parameters["max_iterations"].default == 50
    assert parameters["modality"].default == "T1_MPRAGE"

    descriptor = _write_atlas(tmp_path)
    t1, mask = _synthetic_subject()

    result = register_toads_from_atlas_assets(
        t1,
        mask,
        descriptor,
        levels=2,
        max_iterations=1,
        refinement_iterations=1,
    )

    assert result.atlas.descriptor == descriptor.resolve()
    assert result.modality == "T1_MPRAGE"
    assert result.atlas_data.source_axis_codes == {
        "topology.nii.gz": ("R", "A", "S"),
        "background.nii.gz": ("R", "A", "S"),
        "gm.nii.gz": ("R", "A", "S"),
    }
    assert result.alignment.startup.initial_memberships.shape[-1] == 2
    assert result.alignment.setup.signed_image_resolutions_mm == (-0.8, -1.2, 1.5)
    assert result.alignment.setup.template_resolutions_mm == (1.0, 1.0, 1.0)
    assert result.alignment.registration.levels[0].scale == 2
    assert np.isfinite(result.alignment.registration.transform).all()
    assert result.template_refresh.call_order == (
        "uncropImages",
        "findCroppingBoundaries",
        "cropImages",
        "setImageInfo",
        "computeTransformedTemplate",
        "updateTransformedTemplate",
    )
    assert result.template_refresh.recropped_images.shape[1:] == result.template_refresh.transformed_topology.shape
    assert np.isfinite(result.template_refresh.recropped_images).all()
    assert result.post_refresh_registration.call_order == (
        "segmentation.initCentroids",
        "preprocess.initialMemberships",
        "atlas.updateShapeRegistration",
        "atlas.registerShapes()",
    )
    initialized = result.segmentation_initialization
    assert initialized.call_order == (
        "segmentation.initCentroids",
        "segmentation.setIntensityMax",
        "segmentation.setIntensityScale",
        "segmentation.setLesionCentroid",
        "preprocess.initialMemberships",
    )
    assert initialized.approximate_class_count == 2
    assert initialized.use_approximate_memberships is False
    np.testing.assert_array_equal(initialized.centroids, initialized.prior_centroids)
    np.testing.assert_array_equal(initialized.memberships, result.post_refresh_registration.memberships)
    np.testing.assert_allclose(
        initialized.intensity_max,
        [result.alignment.startup.intensity.robust_maximum],
    )
    np.testing.assert_allclose(
        initialized.intensity_scale,
        [
            result.alignment.startup.intensity.robust_maximum
            - result.alignment.startup.intensity.robust_minimum
        ],
    )
    np.testing.assert_allclose(
        initialized.lesion_centroid,
        [
            result.alignment.startup.intensity.robust_minimum
            / (
                result.alignment.startup.intensity.robust_maximum
                - result.alignment.startup.intensity.robust_minimum
            )
        ],
    )
    constructor = result.segmentation_constructor
    assert constructor.call_order == (
        "segmentation.constructor.mapTemplateLabels",
        "segmentation.constructor.initializeSegmentationAndSkeleton",
        "segmentation.constructor.initializeRelationAndBoundaryState",
    )
    np.testing.assert_array_equal(constructor.classification_labels, result.template_refresh.transformed_topology)
    assert set(np.unique(constructor.segmentation)).issubset({0, 1})
    np.testing.assert_array_equal(constructor.segmentation, constructor.skeleton)
    np.testing.assert_array_equal(constructor.segmentation, constructor.relation_map)
    assert constructor.available.shape == (*constructor.segmentation.shape, 2)
    assert not np.any(constructor.available)
    assert not np.any(constructor.boundary_distances)
    assert int(constructor.class_volumes.sum()) == constructor.segmentation.size
    assert result.relation_state.relation.class_count == 2
    assert result.relation_state.boundary_distances.shape == (*constructor.segmentation.shape, 2)
    assert result.relation_state.closest_labels.shape == (*constructor.segmentation.shape, 4)
    assert result.post_refresh_registration.memberships.shape == (
        *result.template_refresh.recropped_images.shape[1:], 2
    )
    assert result.post_refresh_registration.registered_classes.tolist() == [False, True]
    assert np.isfinite(result.post_refresh_registration.final_transform).all()


def test_atlas_lesion_priors_follow_source_missing_and_present_rules(tmp_path) -> None:
    descriptor = _write_atlas(tmp_path)
    np.testing.assert_array_equal(load_toads_lesion_priors(descriptor, ["MPRAGE"]), [0.0])
    descriptor.write_text(
        descriptor.read_text(encoding="utf-8") + "\nLesions:\nT1_MPRAGE\t0.35\tFLAIR\t0.8\n",
        encoding="utf-8",
    )
    np.testing.assert_array_equal(
        load_toads_lesion_priors(descriptor, ["MPRAGE", "FLAIR", "T2"]),
        np.asarray([0.35, 0.8, 0.0], dtype=np.float32),
    )
    descriptor.write_text(
        descriptor.read_text(encoding="utf-8").replace("T1_MPRAGE\t0.35", "T1_MPRAGE\tNaN"),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Non-finite lesion prior"):
        load_toads_lesion_priors(descriptor, ["MPRAGE"])


def test_segmentation_constructor_rejects_unknown_template_labels() -> None:
    from pymacruise.toads_segmentation_constructor import (
        initialize_toads_segmentation_constructor_state,
    )

    labels = np.zeros((3, 4, 5), dtype=np.uint8)
    labels[1, 1, 1] = 9
    with pytest.raises(ValueError, match="absent from the atlas descriptor"):
        initialize_toads_segmentation_constructor_state(labels, [1, 2])
