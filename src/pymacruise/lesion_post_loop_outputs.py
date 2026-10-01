"""Selected LesionTOADS post-loop output dispatch.

This composes the in-memory tail of ``AlgorithmLesionToads.runAlgorithm``:
snapshot label/subtype outputs, optionally create the membership-derived label
snapshot, and route the selected membership channels. JIST cleanup, image
creation, file I/O, progress/UI callbacks, and field export remain outside this
adapter. Input arrays use RAS-oriented ``(x, y, z)`` voxel order.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .lesion_cruise_inputs import lesion_toads_generate_cruise_inputs
from .lesion_dura_outputs import generate_lesion_dura_removal_outputs
from .lesion_state import LesionToadsFinalClassification
from .lesion_classification_update import update_lesion_classification_from_memberships


LesionOutputMode = Literal["class+member", "cruise inputs", "dura removal inputs"]


@dataclass(frozen=True)
class LesionToadsOutputProduct:
    """One ordered in-memory output array from the selected wrapper tail."""

    name: str
    values: NDArray


@dataclass(frozen=True)
class LesionToadsPostLoopOutputs:
    """Snapshots and mode-specific products selected after LesionTOADS."""

    classification: NDArray[np.int8]
    membership_classification: NDArray[np.int8] | None
    lesion_classification: NDArray[np.int16]
    result_mode: str
    products: tuple[LesionToadsOutputProduct, ...]
    source_order: tuple[str, ...]


def lesion_toads_post_loop_outputs(
    final_classification: LesionToadsFinalClassification,
    memberships: ArrayLike,
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
    segmentation: ArrayLike,
    template_labels: ArrayLike,
    atlas_names: tuple[str, ...] | list[str],
    lesion_class: int,
    *,
    result_mode: str,
    output_classification_from_membership: bool = False,
    include_lesion_classification: bool = False,
) -> LesionToadsPostLoopOutputs:
    """Assemble post-loop products without creating or writing image files.

    The selected source always exports classification and lesion subtype
    snapshots. If requested, it next updates the engine's classification from
    memberships and captures that separate result. It then routes one of the
    three selected result modes. ``include_lesion_classification`` is applied
    to the label-derived classification snapshot after mode dispatch, matching
    the wrapper's later in-memory encoding step; it does not alter the optional
    membership-derived snapshot. Unrecognized modes follow the selected
    wrapper's classification-only fallback. The incomplete legacy
    ``all_images`` branch is intentionally not represented as a supported
    product mode.
    """
    if not isinstance(result_mode, str) or not result_mode:
        raise ValueError("result_mode must be a nonempty string")
    if result_mode == "all_images":
        raise ValueError(
            "result_mode 'all_images' is unsupported: the selected source does not provide a complete in-memory export path"
        )
    if not isinstance(output_classification_from_membership, (bool, np.bool_)):
        raise TypeError("output_classification_from_membership must be boolean")
    if not isinstance(include_lesion_classification, (bool, np.bool_)):
        raise TypeError("include_lesion_classification must be boolean")
    labels = np.asarray(segmentation)
    mems = np.asarray(memberships, dtype=np.float32)
    lesion = np.asarray(lesion_membership, dtype=np.float32)
    black_hole = np.asarray(black_hole_membership, dtype=np.float32)
    names = tuple(atlas_names)
    label_output = np.asarray(final_classification.classification, dtype=np.int8).copy()
    subtype_output = np.asarray(final_classification.lesion_classification, dtype=np.int16).copy()
    if labels.ndim != 3 or label_output.shape != labels.shape or subtype_output.shape != labels.shape:
        raise ValueError("final classifications and segmentation must share one 3-D grid")
    if mems.ndim != 4 or mems.shape[:3] != labels.shape or mems.shape[-1] != len(names):
        raise ValueError("memberships and atlas_names must match segmentation")
    if lesion.shape != labels.shape or black_hole.shape != labels.shape:
        raise ValueError("lesion and black-hole maps must match segmentation")

    source_order = ["cleanUp", "exportClassification", "exportLesionClassification"]
    membership_classification: NDArray[np.int8] | None = None
    if output_classification_from_membership:
        membership_classification = label_output.copy()
        update_lesion_classification_from_memberships(
            membership_classification,
            template_labels,
            mems,
            lesion,
            lesion_class,
        )
        source_order.append("updateClassificationFromMemberships")
        source_order.append("getClassification")

    products: list[LesionToadsOutputProduct] = []
    routed_membership_product: LesionToadsOutputProduct | None = None
    if result_mode == "class+member":
        source_order.extend(("exportMemberships", "getLesionMemberships"))
        channel_first = np.moveaxis(mems, -1, 0).copy()
        stack = np.concatenate((channel_first, lesion[None, ...]), axis=0)
        routed_membership_product = LesionToadsOutputProduct("memberships_plus_lesion", stack)
    elif result_mode == "cruise inputs":
        source_order.append("generateCruiseInputs")
        channels = lesion_toads_generate_cruise_inputs(
            mems, lesion, black_hole, labels, names
        )
        routed_membership_product = LesionToadsOutputProduct("cruise_inputs", channels)
    elif result_mode == "dura removal inputs":
        source_order.append("generateDuraRemovalOutputs")
        channels = generate_lesion_dura_removal_outputs(
            mems, labels, names, lesion, black_hole
        ).as_channels()
        routed_membership_product = LesionToadsOutputProduct("dura_removal_inputs", channels)

    # Java's destination order is mode-specific: class+member exports its
    # membership stack first; CRUISE/dura modes compute their channel arrays
    # above but append them after classification products during image creation.
    if routed_membership_product is not None and result_mode == "class+member":
        products.append(routed_membership_product)

    products.append(LesionToadsOutputProduct("classification", label_output))
    if membership_classification is not None:
        products.append(
            LesionToadsOutputProduct("classification_from_membership", membership_classification)
        )
    products.append(LesionToadsOutputProduct("lesion_subtype", subtype_output))

    if routed_membership_product is not None and result_mode != "class+member":
        products.append(routed_membership_product)

    if include_lesion_classification:
        source_order.append("includeLesionnClassification")
        encoded = (
            label_output.astype(np.int16)
            - np.int16(15) * subtype_output.astype(np.int16)
        ).astype(np.int8)
        encoded[encoded == np.int8(26)] = np.int8(25)
        label_output[...] = encoded

    return LesionToadsPostLoopOutputs(
        classification=label_output,
        membership_classification=membership_classification,
        lesion_classification=subtype_output,
        result_mode=result_mode,
        products=tuple(products),
        source_order=tuple(source_order),
    )
