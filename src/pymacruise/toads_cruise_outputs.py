"""Selected VandyCRUISE adapter for ``generateCruiseOutputs``.

The adapter preserves the four returned float32 channels in Java order:
CSF, GM, WM fill, and assigned-label WM mask. It is array based; NIfTI
coordinate handling belongs to the caller. The source method does not perform
normalization or apply a spatial transform.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pymacruise.toads_crop import uncrop_toads_buffer


_CSF = frozenset(("SulcalCSF", "Sulcal-CSF"))
_GM = frozenset(("CorticalGM", "Cerebrum-GM", "CerebralGM"))
_WM_EXCLUDED = _CSF | _GM | {"Brainstem", "Background"}


@dataclass(frozen=True)
class ToadsCruiseProducts:
    """Four cropped-volume cruise products in Java's channel order."""

    csf: np.ndarray
    gray_matter: np.ndarray
    white_matter_fill: np.ndarray
    white_matter_mask: np.ndarray

    def as_channels(self) -> np.ndarray:
        return np.stack(
            (self.csf, self.gray_matter, self.white_matter_fill, self.white_matter_mask),
            axis=0,
        )


def generate_toads_cruise_outputs(
    memberships: np.ndarray,
    segmentation: np.ndarray,
    atlas_names: tuple[str, ...] | list[str],
    *,
    use_approx: bool,
    best: np.ndarray | None = None,
) -> ToadsCruiseProducts:
    """Adapt the selected ``MultipleToadSegmentation.generateCruiseOutputs``.

    Parameters follow the Java class arrays: ``memberships[x,y,z,class]`` and
    ``segmentation[x,y,z]`` contains class indices. In approximate mode,
    ``best[x,y,z,candidate]`` contains candidate class indices and membership
    values are read from the candidate slot, as Java does; exact mode reads the
    atlas class slot. The selected Java loops process atlas classes in order,
    then approximate slots in order. Matching CSF/GM aliases overwrite earlier
    matches, while eligible WM probabilities add in float32 operation order.

    The fourth output depends only on the assigned segmentation label for
    eligible WM classes. Approximation does not alter that label test. No
    normalization or clipping is performed. Class-name comparisons and the
    ``Cerebell`` WM exclusion prefix are case sensitive, preserving near misses.
    """
    mem = np.asarray(memberships)
    labels = np.asarray(segmentation)
    names = tuple(atlas_names)
    if mem.ndim != 4 or labels.ndim != 3 or mem.shape[:3] != labels.shape:
        raise ValueError("memberships must be [x,y,z,class] and match segmentation")
    if not names or mem.shape[3] != len(names):
        raise ValueError("membership class axis must match a nonempty atlas name list")
    if not np.issubdtype(mem.dtype, np.number) or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("memberships must be numeric and segmentation class indices must be integers")
    if np.any(labels < 0) or np.any(labels >= len(names)):
        raise ValueError("segmentation contains an atlas class index outside atlas_names")
    if not isinstance(use_approx, (bool, np.bool_)):
        raise ValueError("use_approx must be boolean")
    if use_approx:
        if best is None:
            raise ValueError("best candidate class indices are required in approximate mode")
        candidates = np.asarray(best)
        if candidates.ndim != 4 or candidates.shape[:3] != labels.shape:
            raise ValueError("best must have shape [x,y,z,candidate]")
        if candidates.shape[3] < 1 or not np.issubdtype(candidates.dtype, np.integer):
            raise ValueError("best must contain one or more integer candidate class indices")
        if np.any(candidates < 0) or np.any(candidates >= len(names)):
            raise ValueError("best contains an atlas class index outside atlas_names")
    else:
        candidates = None

    mem32 = np.asarray(mem, dtype=np.float32)
    shape = labels.shape
    csf = np.zeros(shape, dtype=np.float32)
    gm = np.zeros(shape, dtype=np.float32)
    wm_fill = np.zeros(shape, dtype=np.float32)
    wm_mask = np.zeros(shape, dtype=np.float32)

    def assign_matched_class(output: np.ndarray, class_index: int) -> None:
        if candidates is None:
            output[...] = mem32[..., class_index]
            return
        # Java loops candidate slots inside each atlas class: the last slot
        # containing this class overwrites earlier slots for this same class.
        for slot in range(candidates.shape[3]):
            matches = candidates[..., slot] == class_index
            output[matches] = mem32[..., slot][matches]

    for class_index, name in enumerate(names):
        if name in _CSF:
            assign_matched_class(csf, class_index)

    for class_index, name in enumerate(names):
        if name in _GM:
            assign_matched_class(gm, class_index)

    for class_index, name in enumerate(names):
        if name in _WM_EXCLUDED or name.startswith("Cerebell"):
            continue
        if candidates is None:
            wm_fill[...] = np.add(wm_fill, mem32[..., class_index], dtype=np.float32)
        else:
            for slot in range(candidates.shape[3]):
                matches = candidates[..., slot] == class_index
                additions = np.where(matches, mem32[..., slot], np.float32(0.0))
                wm_fill[...] = np.add(wm_fill, additions, dtype=np.float32)
        wm_mask[labels == class_index] = np.float32(1.0)

    return ToadsCruiseProducts(csf, gm, wm_fill, wm_mask)


def pack_toads_cruise_outputs(
    products: ToadsCruiseProducts,
    original_shape: tuple[int, int, int],
    bounds: tuple[int, int, int, int, int, int],
    padding: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Restore the wrapper's four scalar FLOAT outputs as separate buffers.

    This corresponds only to the ``AlgorithmMultipleToads`` cruise-output
    wrapper calls to ``uncropAndBuffer(memberships[k], 0.0f)``. It returns raw
    pre-scale buffers; the subsequent image scaling and MIPAV image writing
    are outside this method slice.
    """
    return tuple(
        uncrop_toads_buffer(channel, original_shape, bounds, padding, 0.0)
        for channel in products.as_channels()
    )
