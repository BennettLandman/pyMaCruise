from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.atlas_em_preprocessing import (
    AtlasEMNormalizedImages,
    crop_atlas_em_images,
)
from pymacruise.atlas_em_registration import (
    ATLAS_EM_PLUGIN_DEFAULT_INITIAL_REGISTRATION_ITERATIONS,
    ATLAS_EM_PLUGIN_DEFAULT_DEMONS_SCALE,
    ATLAS_EM_PLUGIN_DEFAULT_DEMONS_SMOOTHING,
    ATLAS_EM_PLUGIN_DEFAULT_LEVELS,
    ATLAS_EM_PLUGIN_DEFAULT_POST_RECROP_REGISTRATION_ITERATIONS,
    initialize_register_recrop_register_atlas_em_from_descriptor,
    register_atlas_em_rigid_single_scale,
    register_atlas_em_descriptor_pyramid,
)
from pymacruise.atlas_em_execution import (
    finish_atlas_em_registered_execution,
    restore_atlas_em_native_nifti_products,
    save_atlas_em_native_nifti_products,
)
from pymacruise.toads_registration_sampling import register_toads_atlas_rigid
from pymacruise.toads_registration_sampling import register_toads_atlas_rigid


def test_atlas_em_selected_plugin_registration_defaults_are_recorded() -> None:
    assert ATLAS_EM_PLUGIN_DEFAULT_LEVELS == 4
    assert ATLAS_EM_PLUGIN_DEFAULT_INITIAL_REGISTRATION_ITERATIONS == 50
    assert ATLAS_EM_PLUGIN_DEFAULT_POST_RECROP_REGISTRATION_ITERATIONS == 1
    assert ATLAS_EM_PLUGIN_DEFAULT_DEMONS_SMOOTHING == 1.0
    assert ATLAS_EM_PLUGIN_DEFAULT_DEMONS_SCALE == 2.0


def test_descriptor_loaded_atlas_reaches_a_scheduled_registration_update(tmp_path: Path) -> None:
    atlas_shape = (5, 6, 7)
    x, y, z = np.indices(atlas_shape, dtype=np.float32)
    prior = np.exp(-((x - 2.0) ** 2 / 2.5 + (y - 2.5) ** 2 / 4.0 + (z - 3.0) ** 2 / 5.0))
    topology = np.where(prior > 0.4, 2, 1).astype(np.uint8)
    affine = np.diag((0.8, 0.9, 1.2, 1.0))
    nib.save(nib.Nifti1Image(topology, affine), tmp_path / "topology.nii.gz")
    nib.save(nib.Nifti1Image(np.ones(atlas_shape, dtype=np.float32), affine), tmp_path / "background.nii.gz")
    nib.save(nib.Nifti1Image(prior.astype(np.float32), affine), tmp_path / "cortex.nii.gz")
    descriptor = tmp_path / "atlas.txt"
    descriptor.write_text(
        "\n".join(
            (
                "Structure Atlas File (edit at your own risks)",
                "",
                "Structures:\t2",
                "Background\t1\tmask",
                "Cortex\t2\tclosed",
                "",
                "Topology Atlas:",
                "File:\ttopology.nii.gz",
                "Dimensions:\t5\t6\t7",
                "Resolutions:\t0.8\t0.9\t1.2",
                "",
                "Shape Atlas:\t4",
                "Dimensions:\t5\t6\t7",
                "Resolutions:\t0.8\t0.9\t1.2",
                "Orientations:\t1\t1\t6\t3",
                "Structure:\tBackground",
                "File:\tbackground.nii.gz",
                "Structure:\tCortex",
                "File:\tcortex.nii.gz",
                "Registered Shapes (0/1: true/false)",
                "0\t1",
                "0\t1",
            )
        ) + "\n",
        encoding="utf-8",
    )

    subject_shape = (8, 8, 8)
    sx, sy, sz = np.indices(subject_shape, dtype=np.float32)
    cortical_membership = np.exp(-((sx - 4.5) ** 2 / 4.0 + (sy - 4.0) ** 2 / 6.0 + (sz - 3.5) ** 2 / 8.0))
    memberships = np.stack((1.0 - cortical_membership, cortical_membership), axis=-1).astype(np.float32)

    result = register_atlas_em_descriptor_pyramid(
        memberships,
        str(descriptor),
        levels=2,
        image_resolutions_mm=(0.8, 0.9, 1.2),
        initial_transform=(0.0,) * 6,
        registration_iterations=1,
    )

    assert result.atlas.class_names == ("Background", "Cortex")
    assert result.atlas.registered_classes.tolist() == [False, True]
    assert result.registration.transform.shape == (6,)
    assert len(result.registration.levels) == 1
    assert result.registration.levels[0].scale == 2
    assert len(result.registration.levels[0].result.iterations) == 1


def test_descriptor_registration_completes_native_grid_recrop_handoff(tmp_path: Path) -> None:
    shape = (12, 12, 12)
    x, y, z = np.indices(shape, dtype=np.float32)
    foreground = (x - 5.5) ** 2 + (y - 5.5) ** 2 + (z - 5.5) ** 2 < 17
    image = np.zeros((*shape, 1), dtype=np.float32)
    image[..., 0][foreground] = 0.5 + 0.02 * x[foreground]
    normalized = AtlasEMNormalizedImages(
        image, np.array([0.0], np.float32), np.array([1.0], np.float32)
    )
    crop = crop_atlas_em_images(
        image,
        (2, 2, 2),
        (9, 9, 9),
        border=2,
        signed_resolutions_mm=(1.0, 1.0, 1.0),
        transform=(0.0,) * 6,
    )
    background = np.ones(shape, dtype=np.float32)
    cortex = np.exp(-((x - 5.0) ** 2 + (y - 6.0) ** 2 + (z - 5.0) ** 2) / 18.0).astype(np.float32)
    affine = np.eye(4)
    nib.save(nib.Nifti1Image(np.where(cortex > 0.4, 2, 1).astype(np.uint8), affine), tmp_path / "topology.nii.gz")
    nib.save(nib.Nifti1Image(background, affine), tmp_path / "background.nii.gz")
    nib.save(nib.Nifti1Image(cortex, affine), tmp_path / "cortex.nii.gz")
    descriptor = tmp_path / "atlas.txt"
    descriptor.write_text(
        "\n".join(
            (
                "Structure Atlas File (edit at your own risks)", "", "Structures:\t2",
                "Background\t1\tmask", "Cortex\t2\tclosed", "", "Topology Atlas:",
                "File:\ttopology.nii.gz", "Dimensions:\t12\t12\t12", "Resolutions:\t1\t1\t1",
                "", "Shape Atlas:\t4", "Dimensions:\t12\t12\t12", "Resolutions:\t1\t1\t1",
                "Orientations:\t1\t1\t6\t3", "Structure:\tBackground", "File:\tbackground.nii.gz",
                "Structure:\tCortex", "File:\tcortex.nii.gz", "Registered Shapes:", "0\t1", "0\t1",
                "", "Intensity Atlas:\t10", "T1_MPRAGE\t0.2\t0.8",
            )
        ) + "\n",
        encoding="utf-8",
    )
    result = initialize_register_recrop_register_atlas_em_from_descriptor(
        crop,
        normalized,
        str(descriptor),
        ("T1_MPRAGE",),
        background_ratio=0.1,
        levels=2,
        initial_registration_iterations=1,
        post_recrop_registration_iterations=1,
    )

    assert len(result.initial_registration.levels) == 1
    assert result.initial_registration.levels[0].scale == 2
    assert len(result.initial_registration.levels[0].result.iterations) == 1
    assert result.initialization.class_names == ("Background", "Cortex")
    assert result.initialization.memberships.shape == (*crop.images.shape[:3], 2)
    assert np.isfinite(result.initialization.memberships).all()
    assert result.post_recrop_initialization.memberships.shape == (*result.recropped.images.shape[:3], 2)
    assert len(result.post_recrop_registration.iterations) == 1
    assert np.isfinite(result.post_recrop_registration.transform).all()
    finished = finish_atlas_em_registered_execution(
        result, str(descriptor), max_iterations=1
    )
    assert finished.mask.shape == result.recropped.images.shape[:3]
    assert finished.priors.shape == result.post_recrop_initialization.memberships.shape
    assert finished.iterations.iterations == 1
    assert finished.segmentation.shape == result.recropped.images.shape[:3]
    assert finished.dura_inputs.shape == (*result.recropped.images.shape[:3], 3)
    assert np.isfinite(finished.dura_inputs).all()
    reference_t1 = nib.Nifti1Image(image[..., 0], affine)
    native_products = restore_atlas_em_native_nifti_products(
        reference_t1, result, finished
    )
    for product in (
        native_products.memberships,
        native_products.segmentation,
        native_products.internal_mask,
        native_products.dura_inputs,
    ):
        assert product.shape[:3] == shape
        assert nib.aff2axcodes(product.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(product.affine, reference_t1.affine)
    paths = save_atlas_em_native_nifti_products(native_products, tmp_path / "outputs")
    assert set(paths) == {"memberships", "segmentation", "internal_mask", "dura_inputs"}
    assert all(path.name.endswith(".nii.gz") and path.is_file() for path in paths.values())
    assert result.recropped.original_dimensions_xyz == shape
    assert result.recropped.images.shape[3] == 1
    assert np.isfinite(result.recropped.images).all()


def test_atlas_em_post_recrop_single_scale_update_runs_source_faithful() -> None:
    shape = (7, 7, 7)
    x, y, z = np.indices(shape, dtype=np.float32)
    prior = np.exp(-((x - 3.2) ** 2 / 6.0 + (y - 2.8) ** 2 / 8.0 + (z - 3.5) ** 2 / 5.0))
    prior *= (0.75 + 0.025 * x + 0.015 * y).astype(np.float32)
    atlas = np.stack((np.ones(shape, dtype=np.float32), prior.astype(np.float32)), axis=-1)
    memberships = np.stack((0.25 + 0.01 * z, 0.3 + 0.02 * x + 0.01 * y), axis=-1).astype(np.float32)

    result = register_atlas_em_rigid_single_scale(
        memberships,
        atlas,
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        initial_transform=(0.01, -0.015, 0.02, 0.0, 0.0, 0.0),
        registration_iterations=1,
        registered_classes=np.array([False, True]),
    )

    assert len(result.iterations) == 1
    assert result.iterations[0].rejected_trials <= 11
    assert np.isfinite(result.energy)
    assert np.isfinite(result.transform).all()
    assert result.stop_reason in {"iteration_limit", "energy_tolerance"}

    expected = register_toads_atlas_rigid(
        memberships,
        atlas,
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        initial_transform=(0.01, -0.015, 0.02, 0.0, 0.0, 0.0),
        class_bounds_low_xyz=np.zeros((2, 3), dtype=np.int64),
        class_bounds_high_xyz=np.full((2, 3), 7, dtype=np.int64),
        max_iterations=1,
        registered_classes=np.array([False, True]),
        offset=0,
        subsample=2,
        energy_tolerance=1e-6,
        minimum_lambda=1e-3,
        source_faithful=True,
        allow_nonfinite_trials=True,
        scale=1.0,
        precomputed_affine=True,
    )
    np.testing.assert_array_equal(result.transform, expected.transform)
    np.testing.assert_equal(result.energy, expected.energy)
    assert len(result.iterations) == len(expected.iterations)
    for actual, source in zip(result.iterations, expected.iterations, strict=True):
        np.testing.assert_array_equal(actual.transform, source.transform)
        np.testing.assert_equal(actual.energy, source.energy)
        np.testing.assert_equal(actual.lambda_value, source.lambda_value)
        np.testing.assert_equal(actual.relative_change, source.relative_change)
        assert actual.accepted == source.accepted
        assert actual.rejected_trials == source.rejected_trials

    expected = register_toads_atlas_rigid(
        memberships,
        atlas,
        image_resolutions_mm=(1.0, 1.0, 1.0),
        template_resolutions_mm=(1.0, 1.0, 1.0),
        initial_transform=(0.01, -0.015, 0.02, 0.0, 0.0, 0.0),
        class_bounds_low_xyz=np.zeros((2, 3), dtype=np.int64),
        class_bounds_high_xyz=np.full((2, 3), 7, dtype=np.int64),
        max_iterations=1,
        registered_classes=np.array([False, True]),
        offset=0,
        subsample=2,
        energy_tolerance=1e-6,
        minimum_lambda=1e-3,
        source_faithful=True,
        allow_nonfinite_trials=True,
        scale=1.0,
        precomputed_affine=True,
    )
    np.testing.assert_array_equal(result.transform, expected.transform)
    np.testing.assert_equal(result.energy, expected.energy)
    assert len(result.iterations) == len(expected.iterations)
    for actual, source in zip(result.iterations, expected.iterations, strict=True):
        np.testing.assert_array_equal(actual.transform, source.transform)
        np.testing.assert_equal(actual.energy, source.energy)
        np.testing.assert_equal(actual.lambda_value, source.lambda_value)
        np.testing.assert_equal(actual.relative_change, source.relative_change)
        assert actual.accepted == source.accepted
        assert actual.rejected_trials == source.rejected_trials
