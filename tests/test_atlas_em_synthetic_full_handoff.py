from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.atlas_em_execution import (
    finish_atlas_em_registered_execution,
    restore_atlas_em_native_nifti_products,
)
from pymacruise.atlas_em_preprocessing import (
    AtlasEMNormalizedImages,
    crop_atlas_em_images,
)
from pymacruise.atlas_em_registration import (
    initialize_register_recrop_register_atlas_em_from_descriptor,
)
from pymacruise.plus_atlas_em_dura import run_dura_from_atlas_em_products


def test_synthetic_descriptor_to_atlasem_products_to_dura(tmp_path: Path) -> None:
    shape = (12, 12, 12)
    x, y, z = np.indices(shape, dtype=np.float32)
    radius = np.sqrt((x - 5.5) ** 2 + (y - 5.5) ** 2 + (z - 5.5) ** 2)
    foreground = radius < 4.6
    t1_values = np.where(foreground, 0.35 + 0.025 * x + 0.015 * y, 0.0).astype(np.float32)
    images = t1_values[..., None]
    normalized = AtlasEMNormalizedImages(
        images, np.array([0.0], dtype=np.float32), np.array([1.0], dtype=np.float32)
    )
    crop = crop_atlas_em_images(
        images,
        (2, 2, 2),
        (9, 9, 9),
        border=2,
        signed_resolutions_mm=(1.0, 1.0, 1.0),
        transform=(0.0,) * 6,
    )

    gm = np.exp(-((radius - 3.5) ** 2) / 2.5).astype(np.float32)
    wm = np.exp(-(radius**2) / 18.0).astype(np.float32)
    background = np.ones(shape, dtype=np.float32)
    affine = np.eye(4)
    nib.save(nib.Nifti1Image(np.where(radius < 2.4, 3, np.where(radius < 4.6, 2, 1)).astype(np.uint8), affine), tmp_path / "topology.nii.gz")
    nib.save(nib.Nifti1Image(background, affine), tmp_path / "background.nii.gz")
    nib.save(nib.Nifti1Image(gm, affine), tmp_path / "gm.nii.gz")
    nib.save(nib.Nifti1Image(wm, affine), tmp_path / "wm.nii.gz")
    descriptor = tmp_path / "atlas.txt"
    descriptor.write_text(
        "\n".join(
            (
                "Structure Atlas File (edit at your own risks)", "", "Structures:\t3",
                "Background\t1\tmask", "Cerebrum-GM\t2\tclosed", "Cerebrum-WM\t3\tclosed",
                "", "Topology Atlas:", "File:\ttopology.nii.gz", "Dimensions:\t12\t12\t12",
                "Resolutions:\t1\t1\t1", "", "Shape Atlas:\t4", "Dimensions:\t12\t12\t12",
                "Resolutions:\t1\t1\t1", "Orientations:\t1\t1\t6\t3",
                "Structure:\tBackground", "File:\tbackground.nii.gz", "Structure:\tCerebrum-GM",
                "File:\tgm.nii.gz", "Structure:\tCerebrum-WM", "File:\twm.nii.gz",
                "Registered Shapes:", "0\t1\t1", "0\t1\t1", "", "Intensity Atlas:\t10",
                "T1_MPRAGE\t0.2\t0.5\t0.8",
            )
        ) + "\n",
        encoding="utf-8",
    )

    registered = initialize_register_recrop_register_atlas_em_from_descriptor(
        crop,
        normalized,
        str(descriptor),
        ("T1_MPRAGE",),
        background_ratio=0.1,
        levels=2,
        initial_registration_iterations=1,
        post_recrop_registration_iterations=1,
    )
    execution = finish_atlas_em_registered_execution(
        registered, str(descriptor), max_iterations=1
    )
    reference_t1 = nib.Nifti1Image(t1_values, affine)
    products = restore_atlas_em_native_nifti_products(reference_t1, registered, execution)
    dura = run_dura_from_atlas_em_products(
        reference_t1, products.dura_inputs, max_iterations=1
    )

    assert execution.iterations.iterations == 1
    assert products.dura_inputs.shape == (*shape, 3)
    assert np.any(np.asarray(products.dura_inputs.dataobj)[..., 1] > 0)
    assert dura.dura.evolution.completed_outer_passes == 1
    assert nib.aff2axcodes(dura.dura.stripped_t1.affine) == ("R", "A", "S")
    assert np.isfinite(np.asarray(dura.dura.stripped_t1.dataobj)).all()
