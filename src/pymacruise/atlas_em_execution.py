"""Compose the selected AtlasEM registration result with its fixed-prior EM tail."""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

import nibabel as nib
import numpy as np
from numpy.typing import NDArray

from pymacruise.atlas_em_iteration import (
    AtlasEMIterationResult,
    run_atlas_em_fixed_prior_iterations,
)
from pymacruise.atlas_em_preprocessing import (
    build_atlas_em_internal_mask,
    transform_atlas_em_rigid_priors,
)
from pymacruise.atlas_em_registration import AtlasEMTwoStageRegistrationResult
from pymacruise.atlas_em_preprocessing import uncrop_atlas_em_images
from pymacruise.atlas_em_statistics import (
    compute_atlas_em_atlas_prior_memberships,
    compute_atlas_em_diagonal_covariances,
)
from pymacruise.atlas_em_outputs import (
    generate_atlas_em_dura_inputs,
    generate_atlas_em_segmentation,
)
from pymacruise.atlas_em_descriptor import parse_atlas_em_descriptor
from pymacruise.image_io import normalize_image_to_ras, save_nifti_ras


@dataclass(frozen=True)
class AtlasEMRegisteredExecutionResult:
    """Transformed priors, EM results, tissue indices and Dura inputs."""

    mask: NDArray[np.bool_]
    priors: NDArray[np.float32]
    initial_memberships: NDArray[np.float32]
    initial_metrics: NDArray[np.float32]
    iterations: AtlasEMIterationResult
    segmentation: NDArray[np.int8]
    dura_inputs: NDArray[np.float32]
    class_names: tuple[str, ...]


@dataclass(frozen=True)
class AtlasEMNiftiProducts:
    """Native-grid RAS NIfTI products restored from the AtlasEM crop."""

    memberships: nib.Nifti1Image
    segmentation: nib.Nifti1Image
    internal_mask: nib.Nifti1Image
    dura_inputs: nib.Nifti1Image


def finish_atlas_em_registered_execution(
    registration: AtlasEMTwoStageRegistrationResult,
    descriptor_path: str,
    *,
    max_iterations: int = 99,
    max_membership_change: float = 0.01,
    smoothing: float = 10.0,
) -> AtlasEMRegisteredExecutionResult:
    """Run AtlasEM's initialization membership update, fixed-prior loop, outputs.

    The selected wrapper computes diagonal covariance estimates, transforms
    and normalizes the registered shape priors, performs an initial Atlas-mode
    membership update, then iterates centroid/covariance/posterior updates
    until the source maximum-change threshold or cap. This function consumes
    the preceding two-stage descriptor registration output. It does not run
    the pre-LesionTOADS Dura level-set evolution or export NIfTI files.
    """
    images = np.asarray(registration.recropped.images, dtype=np.float32)
    memberships = np.asarray(registration.post_recrop_initialization.memberships, dtype=np.float32)
    centroids = np.asarray(registration.post_recrop_initialization.centroids, dtype=np.float32)
    if images.ndim != 4 or images.shape[-1] != 1:
        raise ValueError("selected AtlasEM execution currently requires one T1 channel")
    if memberships.shape[:3] != images.shape[:3] or centroids.shape != (1, memberships.shape[-1]):
        raise ValueError("registered AtlasEM initialization does not match the recropped T1")
    descriptor = parse_atlas_em_descriptor(descriptor_path)
    if len(descriptor.structures) != memberships.shape[-1]:
        raise ValueError("descriptor class count does not match registered memberships")

    active = build_atlas_em_internal_mask(images)
    dimensions = tuple(int(value) for value in images.shape[:3])
    priors = transform_atlas_em_rigid_priors(
        registration.atlas.shape_priors,
        active,
        image_dimensions_xyz=dimensions,
        image_resolutions_mm=registration.recropped.signed_resolutions_mm,
        template_resolutions_mm=descriptor.shape_resolutions_mm,
        transform=registration.post_recrop_registration.transform,
        shape_slope=descriptor.shape_slope,
    )

    initial_deviations = compute_atlas_em_diagonal_covariances(
        images, memberships, centroids
    )
    covariance = np.zeros((memberships.shape[-1], 1, 1), dtype=np.float32)
    covariance[:, 0, 0] = initial_deviations[0]
    initialized_memberships, _initial_labels, initial_metrics = compute_atlas_em_atlas_prior_memberships(
        images,
        memberships,
        active,
        centroids,
        covariance,
        priors,
        smoothing=smoothing,
    )
    iterations = run_atlas_em_fixed_prior_iterations(
        images,
        initialized_memberships,
        active,
        centroids,
        covariance,
        priors,
        max_iterations=max_iterations,
        max_membership_change=max_membership_change,
        smoothing=smoothing,
    )
    segmentation = generate_atlas_em_segmentation(iterations.memberships)
    dura_inputs = generate_atlas_em_dura_inputs(
        iterations.memberships,
        segmentation,
        active,
        tuple(item.name for item in descriptor.structures),
    )
    return AtlasEMRegisteredExecutionResult(
        active,
        priors,
        initialized_memberships,
        np.asarray(initial_metrics, dtype=np.float32),
        iterations,
        segmentation,
        dura_inputs,
        tuple(item.name for item in descriptor.structures),
    )


def restore_atlas_em_native_nifti_products(
    original_t1: object,
    registration: AtlasEMTwoStageRegistrationResult,
    execution: AtlasEMRegisteredExecutionResult,
) -> AtlasEMNiftiProducts:
    """Restore the crop products to the input T1 grid and preserve RAS geometry."""
    if not isinstance(original_t1, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError("original_t1 must be a NIfTI image")
    t1 = normalize_image_to_ras(original_t1, source_name="AtlasEM output reference T1")
    if len(t1.shape) != 3 or tuple(t1.shape) != registration.recropped.original_dimensions_xyz:
        raise ValueError("original T1 must match the native grid recorded by AtlasEM recropping")
    crop = registration.recropped
    if execution.dura_inputs.shape[:3] != crop.images.shape[:3]:
        raise ValueError("AtlasEM Dura-input products must match the final crop grid")
    if execution.segmentation.shape != crop.images.shape[:3] or execution.mask.shape != crop.images.shape[:3]:
        raise ValueError("AtlasEM segmentation and mask must match the final crop grid")
    packed = np.concatenate(
        (
            np.asarray(execution.iterations.memberships, dtype=np.float32),
            np.asarray(execution.dura_inputs, dtype=np.float32),
            execution.mask[..., None].astype(np.float32),
            execution.segmentation[..., None].astype(np.float32),
        ),
        axis=-1,
    )
    native, _native_transform = uncrop_atlas_em_images(replace(crop, images=packed))
    classes = execution.iterations.memberships.shape[-1]
    if classes != len(execution.class_names) or native.shape[-1] != classes + 5:
        raise RuntimeError("native AtlasEM product assembly produced an unexpected channel count")

    def output(data: np.ndarray, name: str, dtype: np.dtype) -> nib.Nifti1Image:
        values = np.asarray(data, dtype=dtype)
        header = t1.header.copy()
        header.set_data_dtype(values.dtype)
        header.set_slope_inter(None, None)
        image = nib.Nifti1Image(values, t1.affine, header=header)
        _, qcode = t1.get_qform(coded=True)
        _, scode = t1.get_sform(coded=True)
        image.set_qform(t1.affine if qcode else None, code=int(qcode), strip_shears=True)
        image.set_sform(t1.affine if scode else None, code=int(scode))
        image.set_filename(name)
        return image

    return AtlasEMNiftiProducts(
        output(native[..., :classes], "atlasem_memberships.nii.gz", np.dtype(np.float32)),
        output(native[..., classes + 4], "atlasem_segmentation.nii.gz", np.dtype(np.int8)),
        output(native[..., classes + 3] > 0.0, "atlasem_internal_mask.nii.gz", np.dtype(np.uint8)),
        output(native[..., classes:classes + 3], "atlasem_dura_inputs.nii.gz", np.dtype(np.float32)),
    )


def save_atlas_em_native_nifti_products(
    products: AtlasEMNiftiProducts,
    output_directory: str | Path,
) -> dict[str, Path]:
    """Write the restored AtlasEM products as compressed RAS NIfTI files."""
    destination = Path(output_directory)
    names = {
        "memberships": "atlasem_memberships.nii.gz",
        "segmentation": "atlasem_segmentation.nii.gz",
        "internal_mask": "atlasem_internal_mask.nii.gz",
        "dura_inputs": "atlasem_dura_inputs.nii.gz",
    }
    return {
        key: save_nifti_ras(getattr(products, key), destination / filename)
        for key, filename in names.items()
    }
