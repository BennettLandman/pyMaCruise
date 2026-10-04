"""Construct startup segmentation state from the transformed atlas template.

This ports the constructor initialization in VandyCRUISE's
``MultipleToadSegmentation``: each atlas template byte is mapped through
``templateLabel`` to its class index, and the same index initializes
segmentation, skeleton, and relation map. The constructor also seeds zero
ordering/availability and the lesion/outlier maps. This module intentionally
does not create the transformed template itself; rigid/nonrigid atlas
registration and interpolation are separate upstream stages.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_growth import LesionToadsThinningPass, propagate_lesion_toads_thinning
from .lesion_initial_atlas_alignment import LesionLoadedAtlasInitialAlignment
from .lesion_loaded_atlas_first_thinning import (
    LesionLoadedAtlasFirstThinning,
    run_loaded_lesion_atlas_first_thinning,
)


@dataclass(frozen=True)
class LesionInitialSegmentationState:
    """Arrays initialized from a transformed source template and constructor."""

    segmentation: NDArray[np.int8]
    skeleton: NDArray[np.int8]
    relation_map: NDArray[np.int8]
    available: NDArray[np.bool_]
    ordering: NDArray[np.float32]
    lesion_membership: NDArray[np.float32]
    outlier_membership: NDArray[np.float32]
    class_volume: NDArray[np.float32]


@dataclass(frozen=True)
class LesionLoadedTemplateFirstThinning:
    """Loaded-atlas first thinning plus its constructor-derived startup maps."""

    state: LesionInitialSegmentationState
    first_thinning: LesionLoadedAtlasFirstThinning
    thinning: LesionToadsThinningPass
    stage_order: tuple[str, ...]


def initialize_lesion_segmentation_state(
    transformed_template: ArrayLike,
    template_labels: Sequence[int],
) -> LesionInitialSegmentationState:
    """Map transformed template labels to constructor-initialized class state.

    ``transformed_template`` is the label image returned by the atlas's
    ``computeTransformedTemplate``. Its values are source template labels,
    whereas segmentation arrays store zero-based class indices. Every source
    voxel must match a declared label, mirroring the source's abort on an
    ``EMPTY`` class. The six helper maps follow the exact selected constructor
    initialization; actual memberships are supplied later by
    ``initialMemberships`` and therefore are not synthesized here.
    """
    template = np.asarray(transformed_template)
    if template.ndim != 3 or any(size == 0 for size in template.shape):
        raise ValueError("transformed_template must be a nonempty 3-D label image")
    if not np.issubdtype(template.dtype, np.integer):
        raise TypeError("transformed_template must contain integer atlas labels")
    labels = np.asarray(tuple(template_labels))
    if labels.ndim != 1 or labels.size < 2:
        raise ValueError("template_labels must contain at least two declared classes")
    if not np.issubdtype(labels.dtype, np.integer):
        raise TypeError("template_labels must be integers")
    if np.unique(labels).size != labels.size:
        raise ValueError("template_labels must be unique")
    if np.any(labels < -128) or np.any(labels > 127):
        raise ValueError("template_labels must fit the source signed-byte label range")
    if labels.size > 127:
        raise ValueError("class count must fit the source signed-byte class range")
    if np.any(template < -128) or np.any(template > 127):
        raise ValueError("transformed_template labels must fit the source signed-byte range")

    # Preserve the source's ordered per-class equality scan. It writes the
    # class index k, not the external atlas label, into the three label maps.
    segmentation = np.full(template.shape, -1, dtype=np.int8)
    volumes = np.zeros(labels.size, dtype=np.float32)
    for class_index, template_label in enumerate(labels):
        selected = template == template_label
        segmentation[selected] = np.int8(class_index)
        volumes[class_index] = np.float32(np.count_nonzero(selected))
    if np.any(segmentation == -1):
        missing = np.unique(template[segmentation == -1]).tolist()
        raise ValueError(f"transformed_template contains undeclared labels: {missing}")
    return LesionInitialSegmentationState(
        segmentation=segmentation,
        skeleton=segmentation.copy(),
        relation_map=segmentation.copy(),
        available=np.zeros((*template.shape, labels.size), dtype=np.bool_),
        ordering=np.zeros(template.shape, dtype=np.float32),
        # The selected LesionToadSegmentation constructor initializes this
        # auxiliary map to 0.0f. Descriptor lesion priors seed a later
        # membership calculation; they are not constructor map values.
        lesion_membership=np.zeros(template.shape, dtype=np.float32),
        outlier_membership=np.full(template.shape, np.float32(0.5), dtype=np.float32),
        class_volume=volumes,
    )


def run_loaded_template_first_thinning(
    alignment: LesionLoadedAtlasInitialAlignment,
    transformed_template: ArrayLike,
    transformed_shape: ArrayLike,
    black_hole_membership: ArrayLike,
    *,
    object_types: Sequence[str],
    lesion_class: int,
    volume_ratio: float,
    limit: float,
    use_lesion_max: bool = False,
) -> LesionLoadedTemplateFirstThinning:
    """Derive constructor maps from the transformed template and thin once.

    This removes explicit caller construction of segmentation, skeleton,
    availability, ordering, and initial lesion membership. Atlas registration,
    transformed shapes, and the source initial-membership maps remain inputs.
    """
    if not isinstance(alignment, LesionLoadedAtlasInitialAlignment):
        raise TypeError("alignment must be a LesionLoadedAtlasInitialAlignment")
    atlas_labels = tuple(item.label for item in alignment.priors.structures)
    state = initialize_lesion_segmentation_state(
        transformed_template, atlas_labels
    )
    if state.segmentation.shape != alignment.alignment.memberships.memberships.shape[:3]:
        raise ValueError("transformed_template must use the loaded atlas membership grid")
    if len(object_types) != len(atlas_labels):
        raise ValueError("object_types must match loaded atlas class order")
    first = run_loaded_lesion_atlas_first_thinning(
        alignment,
        state.segmentation,
        transformed_shape,
        state.lesion_membership,
        black_hole_membership,
        state.skeleton,
        state.available,
        state.ordering,
        object_types=object_types,
        lesion_class=lesion_class,
        volume_ratio=volume_ratio,
        limit=limit,
        use_lesion_max=use_lesion_max,
    )
    return LesionLoadedTemplateFirstThinning(
        state=state,
        first_thinning=first,
        thinning=first.thinning,
        stage_order=(
            "atlas.computeTransformedTemplate",
            "MultipleToadSegmentation.constructor:initializeClassState",
            "preprocess.initialMemberships",
            "atlas.registerShapesPyramid",
            "setGMflag",
            "initThinningLabels",
            "propagateThinning",
        ),
    )
