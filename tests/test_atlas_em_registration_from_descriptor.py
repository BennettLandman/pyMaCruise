from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.atlas_em_registration import register_atlas_em_descriptor_pyramid


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
