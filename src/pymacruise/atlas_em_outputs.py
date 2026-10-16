"""Selected AtlasEM result products consumed by the PLUS Dura stage."""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from pymacruise.mipav_morphology import mipav_fill_holes_25d


_WM_EXCLUSIONS = frozenset(
    {"SulcalCSF", "Sulcal-CSF", "CorticalGM", "Cerebrum-GM", "CerebralGM",
     "Cerebellum-GM", "CerebellarGM", "GM", "Background"}
)
_GM_NAMES = frozenset(
    {"CorticalGM", "Cerebrum-GM", "CerebralGM", "Cerebellum-GM", "GM", "CerebellarGM"}
)
_CEREBELLAR_NAMES = ("Cerebellum-GM", "CerebellarGM")


def generate_atlas_em_dura_inputs(
    memberships: ArrayLike,
    segmentation: ArrayLike,
    mask: ArrayLike,
    class_names: tuple[str, ...] | list[str],
) -> NDArray[np.float32]:
    """Port AtlasEM ``generateDuraRemovalOutputs`` for fixed maps and names.

    Returns ``(x, y, z, 3)`` float32 products in source order: GM membership,
    WM membership, and filled WM mask. The first two products are multiplied
    by 255 unconditionally, and the MIPAV 2.5D fill-holes method clears each
    in-slice border before filling holes.
    """
    mem = np.asarray(memberships, dtype=np.float32)
    labels = np.asarray(segmentation)
    active = np.asarray(mask, dtype=bool)
    names = tuple(class_names)
    if mem.ndim != 4 or any(n < 1 for n in mem.shape[:3]) or mem.shape[3] < 1:
        raise ValueError("memberships must have shape (x, y, z, classes)")
    spatial, classes = mem.shape[:3], mem.shape[3]
    if labels.shape != spatial or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("segmentation must be an integer grid matching memberships")
    if active.shape != spatial:
        raise ValueError("mask must match the membership grid")
    if len(names) != classes or any(not isinstance(name, str) for name in names):
        raise ValueError("class_names must contain one string per membership class")
    if np.any(labels < -1) or np.any(labels >= classes):
        raise ValueError("segmentation indices must be -1 or index a membership class")
    if not np.isfinite(mem).all():
        raise ValueError("memberships must be finite")

    products = np.zeros((*spatial, 3), dtype=np.float32)
    cerebellar_ids = tuple(names.index(name) if name in names else -1 for name in _CEREBELLAR_NAMES)
    for k, name in enumerate(names):
        if name not in _WM_EXCLUSIONS:
            for x in range(spatial[0]):
                for y in range(spatial[1]):
                    for z in range(spatial[2]):
                        if active[x, y, z]:
                            products[x, y, z, 1] = np.float32(products[x, y, z, 1] + mem[x, y, z, k])
                            if labels[x, y, z] == k:
                                products[x, y, z, 2] = np.float32(1.0)
    for k, name in enumerate(names):
        if name in _GM_NAMES:
            for x in range(spatial[0]):
                for y in range(spatial[1]):
                    for z in range(spatial[2]):
                        if active[x, y, z]:
                            products[x, y, z, 0] = np.float32(products[x, y, z, 0] + mem[x, y, z, k])
                            if labels[x, y, z] in cerebellar_ids:
                                products[x, y, z, 2] = np.float32(1.0)
    products[..., 2] = mipav_fill_holes_25d(products[..., 2]).astype(np.float32)
    products *= np.float32(255.0)
    return products
