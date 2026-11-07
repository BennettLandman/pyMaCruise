"""Compose the selected AtlasEM registration result with its fixed-prior EM tail."""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from collections.abc import Callable

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
from pymacruise.atlas_em_preprocessing import (
    crop_atlas_em_images,
    find_atlas_em_crop_bounds,
    normalize_atlas_em_images,
)
from pymacruise.atlas_em_registration import (
    ATLAS_EM_PLUGIN_DEFAULT_INITIAL_REGISTRATION_ITERATIONS,
    ATLAS_EM_PLUGIN_DEFAULT_LEVELS,
    ATLAS_EM_PLUGIN_DEFAULT_POST_RECROP_REGISTRATION_ITERATIONS,
    initialize_register_recrop_register_atlas_em_from_descriptor,
)
from pymacruise.image_io import normalize_image_to_ras, save_nifti_ras
from pymacruise.plus_atlas_em_dura import (
    AtlasEMPreLesionDuraResult,
    run_atlas_em_prelesion_dura,
)


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
    covariance_floor: float = 1.0e-6,
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
        covariance_floor=covariance_floor,
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


@dataclass(frozen=True)
class AtlasEMWorkflowResult:
    """Command result for descriptor AtlasEM execution and its Dura handoff."""

    registration: AtlasEMTwoStageRegistrationResult | None
    execution: AtlasEMRegisteredExecutionResult | None
    products: AtlasEMNiftiProducts
    prelesion_dura: AtlasEMPreLesionDuraResult
    output_paths: dict[str, Path]
    em_iterations_completed: int | None
    dura_outer_passes_completed: int


def run_atlas_em_descriptor_workflow(
    t1_image: nib.spatialimages.SpatialImage,
    brain_mask_image: nib.spatialimages.SpatialImage,
    descriptor_path: str | Path,
    output_directory: str | Path,
    *,
    modality: str = "T1_MPRAGE",
    background_ratio: float = 0.0,
    crop_border: int = 10,
    levels: int = ATLAS_EM_PLUGIN_DEFAULT_LEVELS,
    initial_registration_iterations: int = ATLAS_EM_PLUGIN_DEFAULT_INITIAL_REGISTRATION_ITERATIONS,
    post_recrop_registration_iterations: int = ATLAS_EM_PLUGIN_DEFAULT_POST_RECROP_REGISTRATION_ITERATIONS,
    max_em_iterations: int = 99,
    max_membership_change: float = 0.01,
    smoothing: float = 10.0,
    max_dura_iterations: int = 5,
    covariance_floor: float = 1.0e-6,
    precomputed_atlasem_products: AtlasEMNiftiProducts | None = None,
    on_atlasem_complete: Callable[[AtlasEMNiftiProducts, AtlasEMRegisteredExecutionResult | None], None] | None = None,
) -> AtlasEMWorkflowResult:
    """Run the currently supported one-channel descriptor AtlasEM workflow.

    This composes source-profiled intensity normalization, crop and initial
    memberships, both rigid-registration stages with recrop, fixed-prior EM,
    native-grid RAS restoration, and the selected pre-LesionTOADS Dura stage
    using the required T1-grid brain mask. It does not claim direct optimizer
    parity for the selected MIPAV plugin.
    """
    t1 = normalize_image_to_ras(t1_image, source_name="AtlasEM T1")
    if len(t1.shape) != 3:
        raise ValueError("AtlasEM currently requires a three-dimensional T1 image")
    descriptor = parse_atlas_em_descriptor(descriptor_path)
    if modality not in descriptor.intensity_atlas:
        raise ValueError(f"descriptor does not provide the requested modality {modality!r}")
    registration = None
    execution = None
    if precomputed_atlasem_products is None:
        native = np.asarray(t1.dataobj, dtype=np.float32)[..., None]
        normalized = normalize_atlas_em_images(native)
        lower, upper = find_atlas_em_crop_bounds(
            normalized.images,
            normalized.robust_minima,
            normalized.robust_maxima,
            background_ratio=background_ratio,
            normalized=True,
        )
        spacing = tuple(float(value) for value in nib.affines.voxel_sizes(t1.affine))
        crop = crop_atlas_em_images(
            normalized.images,
            lower,
            upper,
            border=crop_border,
            signed_resolutions_mm=spacing,
            transform=(0.0,) * 6,
        )
        registration = initialize_register_recrop_register_atlas_em_from_descriptor(
            crop,
            normalized,
            str(descriptor_path),
            (modality,),
            background_ratio=background_ratio,
            levels=levels,
            initial_registration_iterations=initial_registration_iterations,
            post_recrop_registration_iterations=post_recrop_registration_iterations,
        )
        execution = finish_atlas_em_registered_execution(
            registration,
            str(descriptor_path),
            max_iterations=max_em_iterations,
            max_membership_change=max_membership_change,
            smoothing=smoothing,
            covariance_floor=covariance_floor,
        )
        products = restore_atlas_em_native_nifti_products(t1, registration, execution)
        output_paths = save_atlas_em_native_nifti_products(products, output_directory)
    else:
        products = precomputed_atlasem_products
        output_paths = {
            "memberships": Path(output_directory) / "atlasem_memberships.nii.gz",
            "segmentation": Path(output_directory) / "atlasem_segmentation.nii.gz",
            "internal_mask": Path(output_directory) / "atlasem_internal_mask.nii.gz",
            "dura_inputs": Path(output_directory) / "atlasem_dura_inputs.nii.gz",
        }
    if on_atlasem_complete is not None:
        on_atlasem_complete(products, execution)
    prelesion_dura = run_atlas_em_prelesion_dura(
        t1,
        products.memberships,
        brain_mask_image,
        tuple(item.name for item in descriptor.structures),
        max_iterations=max_dura_iterations,
    )
    destination = Path(output_directory)
    dura_products = {
        "prelesion_dura_inputs": prelesion_dura.products,
        "prelesion_dura_segmentation": nib.Nifti1Image(
            np.asarray(prelesion_dura.segmentation, dtype=np.int8),
            t1.affine,
            header=t1.header.copy(),
        ),
        "dura_stripped_t1": prelesion_dura.dura.stripped_t1,
        "dura_corrected_white_matter": prelesion_dura.dura.corrected_white_matter,
        "dura_corrected_gray_matter": prelesion_dura.dura.corrected_gray_matter,
        "dura_brain_mask_level_set": prelesion_dura.dura.brain_mask_level_set,
        "dura_brain_mask": prelesion_dura.dura.brain_mask,
        "dura_wm_topology": prelesion_dura.dura.wm_topology,
        "dura_gm_topology": prelesion_dura.dura.gm_topology,
    }
    for name, image in dura_products.items():
        output_paths[name] = save_nifti_ras(image, destination / f"{name}.nii.gz")
    return AtlasEMWorkflowResult(
        registration,
        execution,
        products,
        prelesion_dura,
        output_paths,
        None if execution is None else execution.iterations.iterations,
        prelesion_dura.dura.evolution.completed_outer_passes,
    )
