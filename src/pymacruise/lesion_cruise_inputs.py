"""Selected LesionTOADS ``generateCruiseInputs`` adapter.

This reproduces the selected Java method's four channel layout and class-name
routing. Arrays are already RAS-oriented and use Java-style ``(x, y, z)`` voxel
order; this numerical helper does not inspect or reorient image affines.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray


_CSF_NAMES = frozenset({"SulcalCSF", "Sulcal-CSF"})
_GM_NAMES = frozenset({"CorticalGM", "Cerebrum-GM", "CerebralGM"})


def lesion_toads_generate_cruise_inputs(
    memberships: np.ndarray,
    lesion_membership: np.ndarray,
    black_hole_membership: np.ndarray,
    segmentation: np.ndarray,
    atlas_names: Sequence[str],
) -> NDArray[np.float32]:
    """Pack LesionTOADS maps into Java's ``[4][x][y][z]`` CRUISE inputs.

    Channel order is CSF, GM, WM membership, and hard-WM mask. The last atlas
    class matching a CSF or GM alias wins, as in the source's separate loops.
    WM begins with float32 ``lesion + black-hole`` and then adds memberships
    for each class except the exact CSF/GM aliases, names beginning with
    ``Cerebell``, ``Brainstem``, and ``Background``. The hard-WM channel is
    one where the segmentation label is any included class index.

    The selected Java method traverses its complete ``(x, y, z)`` grid and
    stores single-precision values. Inputs must therefore already share the
    same 3D shape, be ordered RAS in voxel space, and contain Java-equivalent
    float32 membership values. This function assumes (but cannot verify) RAS;
    callers loading NIfTI should use the repository's RAS image loader.
    """
    values = np.asarray(memberships, dtype=np.float32)
    lesion = np.asarray(lesion_membership, dtype=np.float32)
    black_hole = np.asarray(black_hole_membership, dtype=np.float32)
    labels = np.asarray(segmentation)
    names = tuple(atlas_names)

    if values.ndim != 4 or any(length == 0 for length in values.shape[:3]):
        raise ValueError("memberships must be a nonempty (x, y, z, classes) array")
    if values.shape[3] != len(names):
        raise ValueError("membership class axis must match atlas_names")
    spatial_shape = values.shape[:3]
    if lesion.shape != spatial_shape or black_hole.shape != spatial_shape:
        raise ValueError("lesion and black-hole maps must match membership spatial shape")
    if labels.shape != spatial_shape:
        raise ValueError("segmentation must match membership spatial shape")
    if not np.issubdtype(labels.dtype, np.integer):
        raise TypeError("segmentation must contain integer atlas indices")
    if any(not isinstance(name, str) for name in names):
        raise TypeError("atlas_names must contain only strings")
    if labels.size and (int(labels.min()) < 0 or int(labels.max()) >= len(names)):
        raise ValueError("segmentation contains an index outside atlas_names")

    output = np.zeros((4, *spatial_shape), dtype=np.float32)

    # Preserve source loop order: later matching aliases overwrite earlier
    # classes rather than being summed.
    for class_index, name in enumerate(names):
        if name in _CSF_NAMES:
            output[0] = values[..., class_index]
    for class_index, name in enumerate(names):
        if name in _GM_NAMES:
            output[1] = values[..., class_index]

    np.add(lesion, black_hole, out=output[2])
    included_class_indices: list[int] = []
    for class_index, name in enumerate(names):
        if (
            name not in _CSF_NAMES
            and name not in _GM_NAMES
            and not name.startswith("Cerebell")
            and name != "Brainstem"
            and name != "Background"
        ):
            np.add(output[2], values[..., class_index], out=output[2])
            included_class_indices.append(class_index)
    if included_class_indices:
        output[3] = np.isin(labels, np.asarray(included_class_indices)).astype(np.float32)

    return output
