"""Selected LesionTOADS post-registration uncrop and recrop continuation.

The preprocessing transform locates the second image crop. The independently
optimized atlas transform resamples topology and shapes on the new grid. Keep
these two transforms separate, as in ``AlgorithmLesionToads``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .lesion_atlas_priors import LesionAtlasTopology
from .lesion_initial_alignment_bounds import (
    LesionInitialAlignmentBounds,
    compute_lesion_initial_alignment_bounds,
)
from .lesion_initial_atlas_alignment import LesionLoadedAtlasInitialAlignment
from .lesion_initial_segmentation_state import (
    LesionInitialSegmentationState,
    initialize_lesion_segmentation_state,
)
from .lesion_centroid_handoff import (
    LesionToadsCentroidHandoffState,
    initialize_lesion_toads_segmentation_centroids,
)
from .image_io import normalize_image_to_ras
from .lesion_preprocessing_crop import (
    LesionPreprocessingCrop,
    crop_lesion_toads_preprocessing_inputs,
)
from .lesion_preprocessing_uncrop import uncrop_lesion_toads_preprocessing_inputs
from .lesion_registered_shapes import transform_registered_lesion_shape_priors
from .lesion_transformed_template import transform_lesion_toads_template
from .lesion_special_centroid_handoff import (
    LesionToadsSpecialCentroidState,
    initialize_lesion_toads_special_centroids,
)
from .lesion_toads_preprocessing import (
    LesionToadsCentroids,
    initialize_lesion_toads_centroids,
)


@dataclass(frozen=True)
class LesionRegisteredRecrop:
    """Image recrop and atlas-derived state on the post-registration grid."""

    alignment: LesionLoadedAtlasInitialAlignment
    topology: LesionAtlasTopology
    uncropped_images: NDArray[np.float32]
    uncropped_preprocessing_transform: NDArray[np.float32]
    preprocessing_crop: LesionPreprocessingCrop
    registration_bounds: LesionInitialAlignmentBounds
    transformed_template: NDArray[np.int8]
    transformed_shapes: NDArray[np.float32]
    constructor_state: LesionInitialSegmentationState
    class_centroid_state: LesionToadsCentroidHandoffState
    special_centroid_result: LesionToadsCentroids
    special_centroid_state: LesionToadsSpecialCentroidState
    stage_order: tuple[str, ...]
    excluded_wrapper_stages: tuple[str, ...]


def run_registered_lesion_recrop(
    alignment: LesionLoadedAtlasInitialAlignment,
    initial_crop: LesionPreprocessingCrop,
    topology: LesionAtlasTopology,
    original_dimensions_xyz: tuple[int, int, int],
    *,
    image_resolutions_mm: tuple[float, float, float],
    template_resolutions_mm: tuple[float, float, float],
    orientation_signs: tuple[float, float, float] = (1.0, 1.0, 1.0),
) -> LesionRegisteredRecrop:
    """Run Java's post-atlas-registration uncrop → recrop → atlas-grid update.

    ``initial_crop`` contains the normalized RAS channels and first crop
    geometry passed into the alignment. Its transform must match the
    preprocessing transform retained by that alignment. After uncropping,
    topology-aware recrop bounds are found using the restored preprocessing
    transform. Atlas topology and shape maps are then sampled using the
    independent transform optimized by atlas registration, on the newly
    cropped dimensions.

    The wrapper's normalization/initial transform estimation, centroid refresh,
    inhomogeneity correction, and segmentation memberships are outside this
    slice. Constructor initialization here means the source template-label
    constructor state only.
    """
    if not isinstance(alignment, LesionLoadedAtlasInitialAlignment):
        raise TypeError("alignment must be a descriptor-loaded LesionTOADS alignment")
    if not isinstance(initial_crop, LesionPreprocessingCrop):
        raise TypeError("initial_crop must be the selected preprocessing crop result")
    if not isinstance(topology, LesionAtlasTopology):
        raise TypeError("topology must be a RAS LesionAtlasTopology result")
    original_dimensions = tuple(int(size) for size in original_dimensions_xyz)
    if len(original_dimensions) != 3 or any(size < 2 for size in original_dimensions):
        raise ValueError("original_dimensions_xyz must contain three dimensions >= 2")
    if len(alignment.priors.structures) != 13:
        raise ValueError("selected LesionTOADS recrop requires the ordered 13-class atlas")
    ras_topology_image = normalize_image_to_ras(
        topology.image, source_name="LesionTOADS topology"
    )
    topology = LesionAtlasTopology(topology.path, ras_topology_image)
    topology_values = np.asarray(topology.image.dataobj)
    if topology_values.ndim != 3:
        raise ValueError("LesionTOADS topology must be three-dimensional")

    initial_transform = np.asarray(
        alignment.alignment.preprocessing_transform, dtype=np.float32
    )
    crop_transform = np.asarray(initial_crop.crop.transform, dtype=np.float32)
    if initial_transform.shape != (6,) or crop_transform.shape != (6,):
        raise ValueError("initial crop and alignment must carry six-value preprocessing transforms")
    if not np.array_equal(initial_transform.view(np.int32), crop_transform.view(np.int32)):
        raise ValueError("initial crop transform must match the alignment preprocessing transform")
    if initial_crop.crop.images.ndim != 4:
        raise ValueError("initial crop images must be channel-first [channel,x,y,z]")

    uncropped_images, uncropped_transform = uncrop_lesion_toads_preprocessing_inputs(
        initial_crop.crop.images,
        original_dimensions,
        initial_crop.crop.bounds,
        initial_crop.padding,
        initial_transform,
        image_resolution_mm=image_resolutions_mm,
        orientation_signs=orientation_signs,
    )
    recrop = crop_lesion_toads_preprocessing_inputs(
        uncropped_images,
        topology_values,
        uncropped_transform,
        image_resolution_mm=image_resolutions_mm,
        template_resolution_mm=template_resolutions_mm,
        orientation_signs=orientation_signs,
        padding=initial_crop.padding,
    )

    structures = alignment.priors.structures
    if alignment.shape_priors.structures != structures:
        raise ValueError("intensity, shape, and topology descriptor class orders must agree")
    atlas_shapes = tuple(
        alignment.shape_priors.shape_maps[item.name] for item in structures
    )
    recrop_dimensions = tuple(int(size) for size in recrop.crop.images.shape[1:])
    atlas_transform = np.asarray(
        alignment.alignment.registration.transform, dtype=np.float32
    )
    registration_bounds = compute_lesion_initial_alignment_bounds(
        alignment.alignment,
        atlas_shapes,
        image_dimensions_xyz=recrop_dimensions,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
    )
    if not np.array_equal(registration_bounds.transform.view(np.int32), atlas_transform.view(np.int32)):
        raise RuntimeError("recrop bounds changed the independently optimized atlas transform")

    background_label = structures[0].label
    transformed_template = transform_lesion_toads_template(
        topology_values,
        atlas_transform,
        recrop_dimensions,
        image_resolution_mm=image_resolutions_mm,
        template_resolution_mm=template_resolutions_mm,
        background_label=background_label,
    )
    transformed_shapes = transform_registered_lesion_shape_priors(
        registration_bounds,
        atlas_shapes,
        image_dimensions_xyz=recrop_dimensions,
        image_resolutions_mm=image_resolutions_mm,
        template_resolutions_mm=template_resolutions_mm,
    )
    constructor_state = initialize_lesion_segmentation_state(
        transformed_template,
        tuple(item.label for item in structures),
    )
    if alignment.alignment.initial_centroids.shape != (recrop.crop.images.shape[0], len(structures)):
        raise ValueError("retained class centroids must match recropped channel and atlas dimensions")
    class_centroid_state = LesionToadsCentroidHandoffState.empty(
        recrop.crop.images.shape[0], len(structures)
    )
    initialize_lesion_toads_segmentation_centroids(
        class_centroid_state, alignment.alignment.initial_centroids
    )
    special_centroids = initialize_lesion_toads_centroids(
        recrop.crop.images,
        alignment.priors.lesion_priors,
        alignment.priors.black_hole_priors,
    )
    special_centroid_state = LesionToadsSpecialCentroidState()
    initialize_lesion_toads_special_centroids(
        special_centroid_state,
        special_centroids.lesion_centroid,
        special_centroids.black_hole_centroid,
    )

    return LesionRegisteredRecrop(
        alignment=alignment,
        topology=topology,
        uncropped_images=uncropped_images,
        uncropped_preprocessing_transform=uncropped_transform,
        preprocessing_crop=recrop,
        registration_bounds=registration_bounds,
        transformed_template=transformed_template,
        transformed_shapes=transformed_shapes,
        constructor_state=constructor_state,
        class_centroid_state=class_centroid_state,
        special_centroid_result=special_centroids,
        special_centroid_state=special_centroid_state,
        stage_order=(
            *alignment.alignment.call_order,
            "atlas.registerShapesPyramid:optimizedAtlasTransform",
            "preprocess.uncropImages:preprocessingTransform",
            "preprocess.findCroppingBoundaries:preprocessingTransform",
            "preprocess.cropImages:preprocessingTransform",
            "atlas.setImageInfo:recroppedDimensions",
            "atlas.computeTransformedTemplate:optimizedAtlasTransform",
            "preprocess.updateTransformedTemplate",
            "MultipleToadSegmentation.constructor:recroppedGrid",
            "LesionToadSegmentation.initCentroids:retainedPreRegistrationResult",
            "preprocess.initialLesionCentroid:recroppedNormalizedImages",
            "preprocess.initialBlackHoleCentroid:recroppedNormalizedImages",
            "setLesionCentroid:reference",
            "setBlackHoleCentroid:reference",
        ),
        excluded_wrapper_stages=(
            "preprocess.normalizeImages",
            "preprocess.estimateInitialTransform",
            "preprocess.updateTransformedTemplate centroid refresh",
            "LesionToadSegmentation.computeMemberships",
            "initThinningLabels",
            "propagateThinning",
        ),
    )
