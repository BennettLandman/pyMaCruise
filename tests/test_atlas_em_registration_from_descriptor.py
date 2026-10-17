from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.atlas_em_preprocessing import (
    AtlasEMNormalizedImages,
    crop_atlas_em_images,
)
from pymacruise.atlas_em_registration import (
    register_and_recrop_atlas_em_from_descriptor,
    register_atlas_em_descriptor_pyramid,
)


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
            )
        ) + "\n",
        encoding="utf-8",
    )
    cortical_membership = np.exp(-((x - 5.0) ** 2 + (y - 6.0) ** 2 + (z - 5.0) ** 2) / 18.0)
    memberships = np.stack(
        (1.0 - cortical_membership, cortical_membership), axis=-1
    ).astype(np.float32)

    result = register_and_recrop_atlas_em_from_descriptor(
        crop,
        memberships,
        normalized,
        str(descriptor),
        background_ratio=0.1,
        levels=2,
        registration_iterations=1,
    )

    assert len(result.registration.levels) == 1
    assert result.registration.levels[0].scale == 2
    assert len(result.registration.levels[0].result.iterations) == 1
    assert result.recropped.original_dimensions_xyz == shape
    assert result.recropped.images.shape[3] == 1
    assert np.isfinite(result.recropped.images).all()
