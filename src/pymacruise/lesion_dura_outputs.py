"""Selectable four-channel dura-removal products for LesionTOADS.

This is specifically ``LesionToadSegmentation.generateDuraRemovalOutputs``
from the selected VandyCRUISE source. It is separate from
``MultipleToadSegmentation.generateDuraRemovalOutputs``: the LesionTOADS
variant initializes white-matter fill with lesion plus black-hole membership.
It reproduces the source's ordered float32 loops, exact label strings,
last-matching-CSF overwrite, and the ``CerebellarGM)`` hard-mask typo. It does
not clip memberships or normalize the four outputs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


_CSF = ("SulcalCSF", "Sulcal-CSF")
_GM = (
    "CorticalGM",
    "Cerebrum-GM",
    "CerebralGM",
    "Cerebellum-GM",
    "CerebellarGM",
)
_WM_EXCLUDED = frozenset(
    (
        *_CSF,
        "CorticalGM",
        "Cerebrum-GM",
        "CerebralGM",
        "Cerebellum-GM",
        "CerebellarGM",
        "Background",
    )
)


@dataclass(frozen=True)
class LesionDuraRemovalOutputs:
    """Four scalar outputs in source order: CSF, GM, WM-fill, hard WM mask."""

    csf: NDArray[np.float32]
    gray_matter: NDArray[np.float32]
    white_matter_fill: NDArray[np.float32]
    white_matter_mask: NDArray[np.float32]

    def as_channels(self) -> NDArray[np.float32]:
        """Return ``[channel, x, y, z]`` outputs without conversion."""
        return np.stack(
            (self.csf, self.gray_matter, self.white_matter_fill, self.white_matter_mask),
            axis=0,
        )


def generate_lesion_dura_removal_outputs(
    memberships: ArrayLike,
    segmentation: ArrayLike,
    atlas_names: tuple[str, ...] | list[str],
    lesion_membership: ArrayLike,
    black_hole_membership: ArrayLike,
) -> LesionDuraRemovalOutputs:
    """Reproduce the LesionTOADS dura-removal output adapter.

    Args:
        memberships: Float-like ``[x,y,z,class]`` atlas/class memberships.
        segmentation: Integer ``[x,y,z]`` atlas class indices.
        atlas_names: Exact, case-sensitive atlas names indexed by segmentation.
        lesion_membership: Lesion membership map, added to WM-fill at init.
        black_hole_membership: Black-hole map, added after lesion membership.

    Each source accumulation occurs in atlas index order with a float32
    rounding after each addition. CSF aliases are assigned in order, so the
    final matching alias wins. GM aliases sum into channel 1. WM-fill starts
    as ``lesion + blackHole`` and adds memberships for all labels except the
    seven explicit exclusions. The hard WM channel marks those same included
    WM classes, plus voxels whose selected atlas name is ``Cerebellum-GM`` or
    the source's misspelled ``CerebellarGM)``. Correctly spelled
    ``CerebellarGM`` contributes to soft GM but does not trigger that special
    hard-mask branch. Results are unnormalized and unclipped.
    """
    mems = np.asarray(memberships)
    labels = np.asarray(segmentation)
    lesion = np.asarray(lesion_membership)
    black_hole = np.asarray(black_hole_membership)
    names = tuple(atlas_names)
    if mems.ndim != 4 or not names or mems.shape[-1] != len(names):
        raise ValueError("memberships must be [x,y,z,class] with one class per atlas name")
    if not np.issubdtype(mems.dtype, np.number):
        raise TypeError("memberships must be numeric")
    if labels.shape != mems.shape[:3] or labels.ndim != 3:
        raise ValueError("segmentation must match the membership [x,y,z] shape")
    if not np.issubdtype(labels.dtype, np.integer):
        raise TypeError("segmentation must contain integer atlas class indices")
    if np.any(labels < 0) or np.any(labels >= len(names)):
        raise ValueError("segmentation contains an atlas class index outside atlas_names")
    if lesion.shape != labels.shape or black_hole.shape != labels.shape:
        raise ValueError("lesion and black-hole maps must match the segmentation shape")
    if not np.issubdtype(lesion.dtype, np.number) or not np.issubdtype(black_hole.dtype, np.number):
        raise TypeError("lesion and black-hole maps must be numeric")

    mems32 = np.asarray(mems, dtype=np.float32)
    lesion32 = np.asarray(lesion, dtype=np.float32)
    black32 = np.asarray(black_hole, dtype=np.float32)
    csf = np.zeros(labels.shape, dtype=np.float32)
    gm = np.zeros(labels.shape, dtype=np.float32)
    wm_fill = np.zeros(labels.shape, dtype=np.float32)
    wm_mask = np.zeros(labels.shape, dtype=np.float32)

    for k, name in enumerate(names):
        if name in _CSF:
            csf[...] = mems32[..., k]

    # Match the source initialization: WM-fill receives lesion first, then
    # black-hole, with Java float addition rounding at this exact operation.
    wm_fill[...] = np.add(lesion32, black32, dtype=np.float32)
    selected_names = np.asarray(names, dtype=object)[labels]
    for k, name in enumerate(names):
        if name not in _WM_EXCLUDED:
            wm_fill[...] = np.add(wm_fill, mems32[..., k], dtype=np.float32)
            wm_mask[labels == k] = np.float32(1.0)

    cerebellar_selected = (selected_names == "Cerebellum-GM") | (
        selected_names == "CerebellarGM)"
    )
    for k, name in enumerate(names):
        if name in _GM:
            gm[...] = np.add(gm, mems32[..., k], dtype=np.float32)
            # This is deliberately outside the atlas-name test: each GM class
            # iteration repeats the same source selected-label predicate.
            wm_mask[cerebellar_selected] = np.float32(1.0)

    return LesionDuraRemovalOutputs(csf, gm, wm_fill, wm_mask)
