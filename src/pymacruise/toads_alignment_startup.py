"""Compose the selected pre-registration TOADS alignment startup stages.

This is deliberately a bounded adapter around already source-checked pieces:
mask/normalize the native T1, find/crop image and topology support, initialize
intensity centroids, and compute the exact-membership array supplied to rigid
atlas registration. It does not perform atlas registration or TOADS evolution.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import nibabel as nib
import numpy as np
from numpy.typing import ArrayLike, NDArray

from .toads_crop import crop_toads_images, find_toads_crop_bounds, transform_toads_topology
from .toads_initialization import initial_exact_memberships, initial_intensity_centroids
from .toads_preprocessing import ToadsIntensityResult, ToadsRASGeometry, prepare_toads_intensity, toads_ras_geometry
from .toads_registration_sampling import (
    ToadsRigidRegistrationPyramidResult,
    register_toads_atlas_rigid_pyramid,
)
from .toads_atlas import (
    ToadsAtlas,
    ToadsAtlasRASData,
    load_toads_atlas,
    load_toads_atlas_ras_data,
    load_toads_intensity_priors,
    load_toads_lesion_priors,
    toads_atlas_shapes_in_structure_order,
)
from .toads_post_registration_template_refresh import (
    ToadsPostRegistrationTemplateRefresh,
    toads_post_registration_template_refresh,
)
from .toads_post_refresh_registration import (
    ToadsPostRefreshRegistration,
    register_toads_after_template_refresh,
)
from .toads_segmentation_initialization import (
    ToadsSegmentationInitialization,
    initialize_toads_segmentation_state,
)
from .toads_segmentation_constructor import (
    ToadsSegmentationConstructorState,
    initialize_toads_segmentation_constructor_state,
)
from .toads_initial_relation_state import toads_initialize_relation_state
from .toads_relations import ToadsRelationState
from .toads_segmentation_initial_variance import (
    ToadsInitialClassVariances,
    initialize_toads_class_variances,
)
from .toads_segmentation_initial_image_variance import (
    ToadsInitialImageVariance,
    initialize_toads_image_variance,
)
from .toads_segmentation_initial_similarity import (
    ToadsInitialCentroidSimilarity,
    initialize_toads_centroid_similarity,
)


@dataclass(frozen=True)
class ToadsAlignmentStartup:
    """Source-ordered state immediately before the first rigid atlas update."""

    intensity: ToadsIntensityResult
    geometry: ToadsRASGeometry
    input_orientation_codes: tuple[str, str, str]
    crop_bounds: tuple[int, int, int, int, int, int]
    cropped_images: NDArray[np.float32]
    cropped_topology: NDArray
    transform: NDArray[np.float32]
    template_resolution_mm: tuple[float, float, float]
    centroids: NDArray[np.float32]
    initial_memberships: NDArray[np.float32]
    stage_order: tuple[str, ...]


def prepare_toads_alignment_startup(
    t1_n4: nib.spatialimages.SpatialImage,
    brain_mask: nib.spatialimages.SpatialImage,
    topology_ras: ArrayLike,
    intensity_priors: ArrayLike,
    template_resolution_mm: tuple[float, float, float],
    *,
    padding: int = 10,
    background_label: int = 1,
    background_ratio: float = 0.0,
) -> ToadsAlignmentStartup:
    """Prepare one-channel base TOADS state up to rigid atlas registration.

    ``topology_ras`` is an atlas label array already loaded/reoriented into the
    selected RAS voxel-axis order; this function cannot infer orientation from
    an array without an affine. ``template_resolution_mm`` comes from the
    bundled atlas descriptor, which the selected Java loader uses in place of
    any conflicting NIfTI header spacing. The T1 and user mask are normalized
    to RAS by :func:`prepare_toads_intensity` and non-RAS inputs warn there.

    The selected wrapper initializes atlas alignment from exact memberships.
    This adapter stops at that registration input and does not claim an atlas
    fit, complete TOADS iteration, or segmentation equivalence.
    """
    input_orientation_codes = tuple(nib.aff2axcodes(t1_n4.affine))
    intensity = prepare_toads_intensity(t1_n4, brain_mask)
    geometry = toads_ras_geometry(intensity.image)
    values = np.asarray(geometry.image.dataobj, dtype=np.float32)
    topology = np.asarray(topology_ras)
    priors = np.asarray(intensity_priors, dtype=np.float32)
    template_resolution = np.asarray(template_resolution_mm, dtype=np.float32)

    if values.ndim != 3 or topology.ndim != 3 or topology.size == 0:
        raise ValueError("T1 and RAS atlas topology must be nonempty 3D arrays")
    if not np.issubdtype(topology.dtype, np.integer):
        raise ValueError("RAS atlas topology must contain integer structure labels")
    if priors.ndim != 2 or priors.shape[0] != 1 or priors.shape[1] < 1:
        raise ValueError("base TOADS intensity priors must have shape (1, classes)")
    if not np.isfinite(priors).all():
        raise ValueError("base TOADS intensity priors must be finite")
    if template_resolution.shape != (3,) or not np.isfinite(template_resolution).all() or np.any(template_resolution <= 0):
        raise ValueError("template_resolution_mm must contain three finite positive values")
    if not isinstance(padding, (int, np.integer)) or padding < 0:
        raise ValueError("padding must be a non-negative integer")

    image_channels = values[np.newaxis, ...]
    initial_transform = np.asarray(geometry.initial_transform, dtype=np.float32)
    robust_ranges = np.asarray(
        [[intensity.robust_minimum, intensity.robust_maximum]], dtype=np.float32
    )
    bounds = find_toads_crop_bounds(
        image_channels,
        topology,
        initial_transform,
        image_resolution=tuple(float(v) for v in np.abs(geometry.signed_image_resolutions_mm)),
        template_resolution=tuple(float(v) for v in template_resolution),
        orientation_signs=geometry.orientation_signs,
        robust_ranges=robust_ranges,
        background_ratio=background_ratio,
        background_label=background_label,
        normalized=True,
    )
    crop = crop_toads_images(
        image_channels,
        bounds,
        int(padding),
        initial_transform,
        image_resolution=tuple(float(v) for v in np.abs(geometry.signed_image_resolutions_mm)),
        orientation_signs=geometry.orientation_signs,
    )
    cropped_topology = transform_toads_topology(
        topology,
        bounds,
        int(padding),
        crop.transform,
        crop_shape=tuple(int(v) for v in crop.images.shape[1:]),
        image_resolution=tuple(float(v) for v in np.abs(geometry.signed_image_resolutions_mm)),
        template_resolution=tuple(float(v) for v in template_resolution),
        orientation_signs=geometry.orientation_signs,
        background_label=background_label,
    )
    centroids = initial_intensity_centroids(
        [intensity.robust_minimum], [intensity.robust_maximum], priors
    )
    memberships = initial_exact_memberships(crop.images, centroids, crop_border=int(padding))
    return ToadsAlignmentStartup(
        intensity=intensity,
        geometry=geometry,
        input_orientation_codes=input_orientation_codes,
        crop_bounds=bounds,
        cropped_images=crop.images,
        cropped_topology=cropped_topology,
        transform=crop.transform,
        template_resolution_mm=tuple(float(v) for v in template_resolution),
        centroids=centroids,
        initial_memberships=memberships,
        stage_order=(
            "normalizeImages",
            "findCroppingBoundaries",
            "cropImages",
            "initialCentroids",
            "initialMemberships",
        ),
    )


@dataclass(frozen=True)
class ToadsAlignmentRegistration:
    """Initial exact-membership atlas-pyramid registration result."""

    startup: ToadsAlignmentStartup
    setup: ToadsRigidRegistrationSetup
    registration: ToadsRigidRegistrationPyramidResult


@dataclass(frozen=True)
class ToadsAtlasAlignmentRun:
    """Loaded atlas state through initial segmentation initialization."""

    atlas: ToadsAtlas
    atlas_data: ToadsAtlasRASData
    alignment: ToadsAlignmentRegistration
    template_refresh: ToadsPostRegistrationTemplateRefresh
    post_refresh_registration: ToadsPostRefreshRegistration
    segmentation_initialization: ToadsSegmentationInitialization
    segmentation_constructor: ToadsSegmentationConstructorState
    relation_state: ToadsRelationState
    initial_variances: ToadsInitialClassVariances
    initial_image_variance: ToadsInitialImageVariance
    initial_centroid_similarity: ToadsInitialCentroidSimilarity
    modality: str


@dataclass(frozen=True)
class ToadsRigidRegistrationSetup:
    """Selected wrapper state immediately before rigid pyramid execution."""

    memberships: NDArray[np.float32]
    best_candidates: None
    initial_zero_transform: NDArray[np.float32]
    transform: NDArray[np.float32]
    transform_type: str
    approximate_class_count: int
    use_approximate_memberships: bool
    max_iterations: int
    levels: int
    subsample_factor: int
    image_dimensions_xyz: tuple[int, int, int]
    template_dimensions_xyz: tuple[int, int, int]
    mipav_image_orientation: int
    mipav_axis_orientation_codes: tuple[int, int, int]
    orientation_signs: tuple[float, float, float]
    signed_image_resolutions_mm: tuple[float, float, float]
    template_resolutions_mm: tuple[float, float, float]
    image_center_xyz: tuple[float, float, float]
    template_center_xyz: tuple[float, float, float]
    stage_order: tuple[str, ...]


def register_toads_alignment_startup(
    startup: ToadsAlignmentStartup,
    atlas_shape_priors: ArrayLike,
    *,
    levels: int,
    max_iterations: int = 50,
    energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
) -> ToadsAlignmentRegistration:
    """Run the selected exact-membership rigid pyramid from startup state.

    This connects the initialized memberships to the selected
    ``registerShapesPyramid`` numerical operation. As in the Java atlas
    method, class zero (MASK/background) is excluded from registration. Shape
    priors are on the template grid and use the same descriptor class order as
    the startup intensity priors. Atlas loading, the preceding
    ``updateShapeRegistration`` setup, and subsequent TOADS stages remain
    separate.
    """
    if not isinstance(startup, ToadsAlignmentStartup):
        raise TypeError("startup must be a ToadsAlignmentStartup")
    if not isinstance(levels, (int, np.integer)) or levels < 1:
        raise ValueError("levels must be a positive integer")
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations < 0:
        raise ValueError("max_iterations must be a non-negative integer")
    if startup.initial_memberships.dtype != np.float32:
        raise ValueError("startup memberships must retain the source float32 representation")
    shapes = np.asarray(atlas_shape_priors, dtype=np.float32)
    expected_classes = startup.initial_memberships.shape[-1]
    if shapes.ndim != 4 or shapes.shape[-1] != expected_classes:
        raise ValueError("atlas_shape_priors must be (template_x,template_y,template_z,classes) in startup class order")
    if not np.isfinite(shapes).all():
        raise ValueError("atlas_shape_priors must be finite")
    setup = ToadsRigidRegistrationSetup(
        memberships=startup.initial_memberships,
        best_candidates=None,
        initial_zero_transform=np.zeros(6, dtype=np.float32),
        transform=startup.transform,
        transform_type="rigid",
        approximate_class_count=int(startup.initial_memberships.shape[-1]),
        use_approximate_memberships=False,
        max_iterations=int(max_iterations),
        levels=int(levels),
        subsample_factor=2,
        image_dimensions_xyz=tuple(int(v) for v in startup.initial_memberships.shape[:3]),
        template_dimensions_xyz=tuple(int(v) for v in shapes.shape[:3]),
        mipav_image_orientation=startup.geometry.image_orientation,
        mipav_axis_orientation_codes=startup.geometry.axis_orientation_codes,
        orientation_signs=startup.geometry.orientation_signs,
        signed_image_resolutions_mm=tuple(
            float(value) for value in startup.geometry.signed_image_resolutions_mm
        ),
        template_resolutions_mm=tuple(float(value) for value in startup.template_resolution_mm),
        image_center_xyz=tuple(
            float(size) / 2.0 for size in startup.initial_memberships.shape[:3]
        ),
        template_center_xyz=tuple(float(size) / 2.0 for size in shapes.shape[:3]),
        stage_order=(
            "setImageInfo",
            "initShapeRegistration",
            "setTransform",
            "updateShapeRegistration",
            "registerShapesPyramid",
        ),
    )
    registered_classes = np.ones(expected_classes, dtype=bool)
    registered_classes[0] = False
    registration = register_toads_atlas_rigid_pyramid(
        setup.memberships,
        shapes,
        levels=setup.levels,
        image_resolutions_mm=setup.signed_image_resolutions_mm,
        template_resolutions_mm=setup.template_resolutions_mm,
        initial_transform=tuple(float(v) for v in setup.transform),
        registered_classes=registered_classes,
        max_iterations=setup.max_iterations,
        energy_tolerance=energy_tolerance,
        minimum_lambda=minimum_lambda,
    )
    return ToadsAlignmentRegistration(startup, setup, registration)


def register_toads_from_atlas_assets(
    t1_n4: nib.spatialimages.SpatialImage,
    brain_mask: nib.spatialimages.SpatialImage,
    atlas: ToadsAtlas | str | Path,
    *,
    levels: int = 4,
    max_iterations: int = 50,
    modality: str = "T1_MPRAGE",
    padding: int = 10,
    background_ratio: float = 0.0,
    refinement_iterations: int = 5,
    relation_limit: int = 3,
    energy_tolerance: float = 1e-6,
    minimum_lambda: float = 1e-3,
) -> ToadsAtlasAlignmentRun:
    """Wire descriptor assets through startup and initial rigid registration.

    This convenience composition loads the selected descriptor, RAS atlas
    arrays, and the requested intensity-prior profile, then runs the bounded
    one-channel startup and rigid-pyramid adapters. Registration levels and
    iteration limits remain explicit call arguments; this function does not
    infer wrapper parameters or claim TOADS iteration/workflow equivalence.
    """
    loaded_atlas = load_toads_atlas(atlas) if isinstance(atlas, (str, Path)) else atlas
    if not isinstance(loaded_atlas, ToadsAtlas):
        raise TypeError("atlas must be a ToadsAtlas or a descriptor path")
    atlas_data = load_toads_atlas_ras_data(loaded_atlas)
    priors = load_toads_intensity_priors(loaded_atlas.descriptor, [modality])
    lesion_priors = load_toads_lesion_priors(loaded_atlas.descriptor, [modality])
    mask_structure_indices = [
        index
        for index, structure in enumerate(loaded_atlas.structures)
        if structure.topology.lower() == "mask"
    ]
    if len(mask_structure_indices) != 1:
        raise ValueError("TOADS descriptor must identify exactly one mask/background structure")
    background_class_index = mask_structure_indices[0]
    background_label = loaded_atlas.structures[background_class_index].label
    startup = prepare_toads_alignment_startup(
        t1_n4,
        brain_mask,
        atlas_data.topology,
        priors,
        atlas_data.template_resolution_mm,
        padding=padding,
        background_ratio=background_ratio,
        background_label=background_label,
    )
    shape_priors = toads_atlas_shapes_in_structure_order(loaded_atlas, atlas_data)
    alignment = register_toads_alignment_startup(
        startup,
        shape_priors,
        levels=levels,
        max_iterations=max_iterations,
        energy_tolerance=energy_tolerance,
        minimum_lambda=minimum_lambda,
    )
    template_refresh = toads_post_registration_template_refresh(
        cropped_images=startup.cropped_images,
        native_shape_xyz=tuple(int(value) for value in startup.intensity.image.shape),
        previous_bounds=startup.crop_bounds,
        padding=padding,
        registered_crop_transform=alignment.registration.transform,
        topology_template=atlas_data.topology,
        image_resolution_mm=tuple(
            abs(value) for value in startup.geometry.signed_image_resolutions_mm
        ),
        template_resolution_mm=atlas_data.template_resolution_mm,
        orientation_signs=startup.geometry.orientation_signs,
        robust_ranges=np.asarray(
            [[startup.intensity.robust_minimum, startup.intensity.robust_maximum]],
            dtype=np.float32,
        ),
        background_ratio=background_ratio,
        background_label=background_label,
    )
    post_refresh_registration = register_toads_after_template_refresh(
        recropped_images=template_refresh.recropped_images,
        centroids=startup.centroids,
        atlas_shapes=shape_priors,
        initial_transform=template_refresh.recropped_transform,
        image_resolutions_mm=alignment.setup.signed_image_resolutions_mm,
        template_resolutions_mm=alignment.setup.template_resolutions_mm,
        background_class_index=background_class_index,
        crop_border=padding,
        max_iterations=refinement_iterations,
        energy_tolerance=energy_tolerance,
        minimum_lambda=minimum_lambda,
    )
    segmentation_constructor = initialize_toads_segmentation_constructor_state(
        template_refresh.transformed_topology,
        [structure.label for structure in loaded_atlas.structures],
    )
    segmentation_initialization = initialize_toads_segmentation_state(
        images=template_refresh.recropped_images,
        centroids=startup.centroids,
        robust_minimum=[startup.intensity.robust_minimum],
        robust_maximum=[startup.intensity.robust_maximum],
        lesion_prior_by_channel=lesion_priors,
        crop_border=padding,
    )
    relation_state = toads_initialize_relation_state(
        segmentation_constructor.segmentation,
        segmentation_initialization.memberships,
        relation_limit,
    )
    image_xyz_channels = np.moveaxis(template_refresh.recropped_images, 0, -1)
    initial_variances = initialize_toads_class_variances(
        images_xyz_channels=image_xyz_channels,
        memberships=segmentation_initialization.memberships,
        centroids=segmentation_initialization.centroids,
        segmentation=segmentation_constructor.segmentation,
        class_types=[structure.topology for structure in loaded_atlas.structures],
    )
    initial_image_variance = initialize_toads_image_variance(
        images_xyz_channels=image_xyz_channels,
        memberships=segmentation_initialization.memberships,
        centroids=segmentation_initialization.centroids,
        class_types=[structure.topology for structure in loaded_atlas.structures],
    )
    initial_centroid_similarity = initialize_toads_centroid_similarity(
        centroids=segmentation_initialization.centroids,
        class_standard_deviations=initial_variances.standard_deviations,
        image_variance_weights=initial_image_variance.image_variance_weights,
        class_types=[structure.topology for structure in loaded_atlas.structures],
    )
    return ToadsAtlasAlignmentRun(
        loaded_atlas,
        atlas_data,
        alignment,
        template_refresh,
        post_refresh_registration,
        segmentation_initialization,
        segmentation_constructor,
        relation_state,
        initial_variances,
        initial_image_variance,
        initial_centroid_similarity,
        modality,
    )
