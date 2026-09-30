"""Selected VandyCRUISE TOADS dura-removal membership output adapter.

This module follows ``MultipleToadSegmentation.generateDuraRemovalOutputs``
and the corresponding four-volume output branch in ``AlgorithmMultipleToads``.
It intentionally preserves the selected Java typo that recognizes
``CerebellarGM)`` (with a trailing parenthesis) when constructing the hard
mask. The correct ``CerebellarGM`` spelling contributes to soft GM, but does
not set the hard mask through that branch. See the fixture provenance beside
the reference data and the selected-method Java harness.

Arrays use the public Python ``[x, y, z, ...]`` convention. The companion
packing function delegates crop placement and x-fastest, channel-major buffer
packing to the independently validated ``uncrop_toads_buffer`` adapter.
``adapt_toads_wrapper_outputs`` also models the selected TOADS08 wrapper's
unconditional multiplication of the four restored float products by 255, while
retaining the raw restored buffers as a separate stage.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pymacruise.toads_crop import uncrop_toads_buffer


_CSF_NAMES = frozenset(("SulcalCSF", "Sulcal-CSF"))
_GM_NAMES = frozenset(
    (
        "CorticalGM",
        "Cerebrum-GM",
        "CerebralGM",
        "Cerebellum-GM",
        "CerebellarGM",
    )
)
_WM_EXCLUDED_NAMES = _CSF_NAMES | _GM_NAMES | {"Background"}


@dataclass(frozen=True)
class ToadsDuraRemovalProducts:
    """Four crop-grid scalar outputs in Java's returned order."""

    csf: np.ndarray
    gray_matter: np.ndarray
    white_matter_fill: np.ndarray
    white_matter_mask: np.ndarray

    def as_channels(self) -> np.ndarray:
        """Stack as ``[output, x, y, z]`` without changing float32 values."""
        return np.stack(
            (self.csf, self.gray_matter, self.white_matter_fill, self.white_matter_mask),
            axis=0,
        )


@dataclass(frozen=True)
class ToadsWrapperOutputBuffers:
    """Hard labels plus raw and always-scaled restored output buffers."""

    classification: np.ndarray
    raw_float_buffers: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    wrapper_float_buffers: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]


def generate_toads_dura_removal_outputs(
    memberships: np.ndarray,
    segmentation: np.ndarray,
    atlas_names: tuple[str, ...] | list[str],
) -> ToadsDuraRemovalProducts:
    """Build source-faithful CSF, GM, WM-fill and WM-mask TOADS products.

    ``memberships`` is ``[x,y,z,class]`` and ``segmentation`` is ``[x,y,z]``;
    class indices in the latter refer to ``atlas_names``. The selected Java
    code iterates atlas classes in order and adds float values into float
    output buffers, so each addition is rounded to float32 here as well.

    Source behavior intentionally retained:

    * Exact class-name comparisons are case-sensitive.
    * CSF output is assigned for each matching class, so the last matching
      CSF class wins if both aliases occur.
    * GM/WM probabilities are accumulated in atlas order.
    * WM hard-mask membership comes from selected WM labels.
    * The WM exclusion checks the correct ``CerebellarGM`` spelling, so an
      atlas label literally named ``CerebellarGM)`` falls through into WM
      fill/mask.
    * The GM pass adds hard-mask voxels for ``Cerebellum-GM`` and the typo
      ``CerebellarGM)``. A label named ``CerebellarGM`` does not match that
      hard-mask test, even though its probability is in the GM output.
    """
    mems = np.asarray(memberships)
    labels = np.asarray(segmentation)
    names = tuple(atlas_names)
    if mems.ndim != 4 or labels.ndim != 3 or mems.shape[:3] != labels.shape:
        raise ValueError("memberships must be [x,y,z,class] and match [x,y,z] segmentation")
    if mems.shape[3] != len(names) or not names:
        raise ValueError("membership class axis must match a nonempty atlas name list")
    if not np.issubdtype(mems.dtype, np.number) or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("memberships must be numeric and segmentation must contain integer class indices")
    if np.any(labels < 0) or np.any(labels >= len(names)):
        raise ValueError("segmentation contains an atlas class index outside atlas_names")

    shape = labels.shape
    # Java's Mems arrays are float[][][]; source membership values are float.
    csf = np.zeros(shape, dtype=np.float32)
    gm = np.zeros(shape, dtype=np.float32)
    wm_fill = np.zeros(shape, dtype=np.float32)
    wm_mask = np.zeros(shape, dtype=np.float32)
    mems_f32 = np.asarray(mems, dtype=np.float32)

    for class_index, name in enumerate(names):
        if name in _CSF_NAMES:
            csf[...] = mems_f32[..., class_index]

    for class_index, name in enumerate(names):
        if name not in _WM_EXCLUDED_NAMES:
            wm_fill[...] = np.add(wm_fill, mems_f32[..., class_index], dtype=np.float32)
            wm_mask[labels == class_index] = np.float32(1.0)

    selected_labels = np.asarray(names, dtype=object)[labels]
    for class_index, name in enumerate(names):
        if name in _GM_NAMES:
            gm[...] = np.add(gm, mems_f32[..., class_index], dtype=np.float32)
            # Keep the selected source's misspelled CerebellarGM literal.
            cerebellar_hard = (selected_labels == "Cerebellum-GM") | (
                selected_labels == "CerebellarGM)"
            )
            wm_mask[cerebellar_hard] = np.float32(1.0)

    return ToadsDuraRemovalProducts(csf, gm, wm_fill, wm_mask)


def pack_toads_dura_removal_outputs(
    products: ToadsDuraRemovalProducts,
    original_shape: tuple[int, int, int],
    bounds: tuple[int, int, int, int, int, int],
    padding: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Apply the selected four-image adapter's zero-background uncrop.

    This mirrors ``AlgorithmMultipleToads``: create four scalar FLOAT images,
    call ``uncropAndBuffer(memberships[k], 0.0f)`` for k=0..3, and import each
    x-fastest buffer into its own image. Returned arrays are linear float32
    buffers, in order CSF, GM, WM fill, WM mask. NIfTI metadata/writing and the
    rest of the MIPAV ``ModelImage`` lifecycle are not included.
    """
    channels = products.as_channels()
    return (
        uncrop_toads_buffer(channels[0], original_shape, bounds, padding, 0.0),
        uncrop_toads_buffer(channels[1], original_shape, bounds, padding, 0.0),
        uncrop_toads_buffer(channels[2], original_shape, bounds, padding, 0.0),
        uncrop_toads_buffer(channels[3], original_shape, bounds, padding, 0.0),
    )


def adapt_toads_wrapper_outputs(
    products: ToadsDuraRemovalProducts,
    classification: np.ndarray,
    original_shape: tuple[int, int, int],
    bounds: tuple[int, int, int, int, int, int],
    padding: int,
) -> ToadsWrapperOutputBuffers:
    """Restore wrapper outputs and apply the selected unconditional 255 scale.

    ``AlgorithmMultipleToads`` first restores the hard classification to a
    byte buffer with background value 1. ``MedicAlgorithmToads08`` then wraps
    that image and unconditionally multiplies the four float output images by
    ``255.0f``; its former scale-data parameter and guard are commented out in
    the selected source. Both unscaled restored buffers and final wrapper
    buffers are returned so the source stages remain separately inspectable.
    Label values use Java signed-byte narrowing.
    """
    labels = np.asarray(classification)
    if labels.shape != products.csf.shape or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("classification must be an integer crop-grid array matching TOADS products")
    # Java byte casts retain the low eight bits and expose them as signed.
    labels_i64 = labels.astype(np.int64, copy=False)
    java_bytes = ((labels_i64 + 128) % 256 - 128).astype(np.int8)
    class_buffer = uncrop_toads_buffer(java_bytes, original_shape, bounds, padding, np.int8(1))
    raw = pack_toads_dura_removal_outputs(products, original_shape, bounds, padding)
    scaled = tuple(np.multiply(buffer, np.float32(255.0), dtype=np.float32) for buffer in raw)
    return ToadsWrapperOutputBuffers(class_buffer, raw, scaled)
