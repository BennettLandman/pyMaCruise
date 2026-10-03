from __future__ import annotations

import nibabel as nib
import numpy as np

from pymacruise.toads_alignment_startup import register_toads_from_atlas_assets
from pymacruise.toads_atlas import load_toads_intensity_variance_factors
from pymacruise.toads_default_atlas_membership import toads_default_atlas_membership_pass
from pymacruise.toads_registered_initial_membership import (
    compute_registered_toads_initial_membership,
)


_CLASS_TYPES = (
    "mask",
    "optimized_csf",
    "optimized_csf",
    "optimized_gm",
    "optimized_gm",
    "optimized_gm",
    "optimized_gm",
    "optimized_gm",
    "optimized_wm",
    "optimized_wm",
    "optimized_wm",
)


def _write_default_atlas(root):
    shape = (15, 13, 11)
    x, y, z = np.indices(shape, dtype=np.float32)
    topology = ((x.astype(np.int16) + 2 * y + 3 * z) % 11 + 1).astype(np.uint8)
    nib.save(nib.Nifti1Image(topology, np.eye(4)), root / "topology.nii.gz")
    structures = [
        ("Background", 1, "mask"),
        ("Csf1", 2, "optimized_csf"),
        ("Csf2", 3, "optimized_csf"),
        ("Gm1", 4, "optimized_gm"),
        ("Gm2", 5, "optimized_gm"),
        ("Gm3", 6, "optimized_gm"),
        ("Gm4", 7, "optimized_gm"),
        ("Gm5", 8, "optimized_gm"),
        ("Wm1", 9, "optimized_wm"),
        ("Wm2", 10, "optimized_wm"),
        ("Wm3", 11, "optimized_wm"),
    ]
    lines = [
        "Structure Atlas File (edit at your own risks)",
        "",
        f"Structures:\t{len(structures)}",
        *(f"{name}\t{label}\t{kind}" for name, label, kind in structures),
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
    ]
    for index, (name, _label, _kind) in enumerate(structures):
        shape_prior = (
            0.15
            + 0.8
            * np.exp(
                -(
                    (x - (2 + index % 7)) ** 2
                    + (y - (2 + (2 * index) % 7)) ** 2
                    + (z - (1 + (3 * index) % 6)) ** 2
                )
                / 24.0
            )
        ).astype(np.float32)
        path = f"shape-{index}.nii.gz"
        nib.save(nib.Nifti1Image(shape_prior, np.eye(4)), root / path)
        lines.extend(
            (
                f"Structure: {name}",
                f"File:\t{path}",
                "Dimensions:\t15\t13\t11",
                "Resolutions:\t1\t1\t1",
                "Center:\t0\t0\t0",
            )
        )
    priors = "\t".join(f"{i / 10:.1f}" for i in range(11))
    lines.extend(("", "Intensity Atlas:\t0", f"T1_MPRAGE\t{priors}", ""))
    descriptor = root / "atlas-11.txt"
    descriptor.write_text("\n".join(lines), encoding="utf-8")
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


def test_registered_startup_adapter_uses_actual_default_atlas_run(tmp_path) -> None:
    descriptor = _write_default_atlas(tmp_path)
    t1, mask = _synthetic_subject()
    run = register_toads_from_atlas_assets(
        t1,
        mask,
        descriptor,
        levels=1,
        max_iterations=1,
        refinement_iterations=1,
        padding=2,
    )

    shape = run.segmentation_initialization.memberships.shape
    exponents = np.full((1, 11), 2.0, dtype=np.float32)
    weights = np.ones(1, dtype=np.float32)
    smoothing_priors = np.ones(11, dtype=np.float32)
    atlas_shape_priors = np.ones(11, dtype=np.float32)
    settings = dict(
        exponent_by_channel_class=exponents,
        modality_weights=weights,
        smoothing_priors=smoothing_priors,
        atlas_shape_priors=atlas_shape_priors,
        smoothing=np.float32(0.2 / 6.0 / 11.0),
        prior_factor=0.2,
        prior_coefficient=np.float32(2.0 / 11.0),
    )

    composed = compute_registered_toads_initial_membership(run, **settings)
    image_values = np.moveaxis(run.template_refresh.recropped_images, 0, -1)
    direct = toads_default_atlas_membership_pass(
        image_values=image_values,
        initial_memberships=run.segmentation_initialization.memberships,
        centroids=run.segmentation_initialization.centroids,
        class_variance_weights=load_toads_intensity_variance_factors(
            descriptor, [run.modality]
        ),
        image_variances=run.initial_image_variance.image_variance_weights,
        transformed_shapes=run.transformed_shape_priors,
        closest_labels=run.relation_state.closest_labels,
        centroid_similarity=run.initial_centroid_similarity.similarity,
        class_types=_CLASS_TYPES,
        **{key: value for key, value in settings.items() if key != "prior_scale"},
    )

    assert run.transformed_shape_priors.shape == shape
    assert composed.class_variances.standard_deviations.shape == (1, 11)
    assert composed.image_variance.image_variance_weights.shape == (1,)
    np.testing.assert_array_equal(
        composed.membership_pass.memberships, direct.memberships
    )
    np.testing.assert_array_equal(
        composed.membership_pass.denominator, direct.denominator
    )
    np.testing.assert_array_equal(composed.membership_pass.objective, direct.objective)
