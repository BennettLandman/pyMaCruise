"""Source-mapped MaCRUISE multi-atlas segmentation refinement.

Ports the selected V3 ``run_SegmentationRefine`` thresholds and intended
neighborhood-mode rules. The MATLAB helper accidentally assigns the row index
returned by ``max`` instead of the corresponding label value; this module
returns the actual modal label and records counts/masks for every rule.
Source reference: ``MAsegRefineBySurf.m`` and ``run_SegmentationRefine.m``,
MaCRUISE SVN r9360.

Published-method context: Huo et al., *Consistent cortical reconstruction
and multi-atlas brain segmentation* (NeuroImage, 2016), §2.4, describes
Cortical Consistent Segmentation Editing (CCSE): compare labels with the
inner/outer surface level sets and consistency distances. This function
implements the selected V3 MATLAB helper's three label rules and neighborhood
fill order; it is not a complete reproduction of every CCSE detail in §2.4.
In particular, the V3 wrapper's 0.2 mm thresholds, label groups, bounded
5×5×5 propagation, and historical ``maxind`` behavior follow source code.
The paper's reported 0.5 mm coefficients are evaluation settings, not defaults
for this source-selected wrapper.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import nibabel as nib
import numpy as np
from nibabel.affines import voxel_sizes

from pymacruise.image_io import normalize_image_to_ras


@dataclass(frozen=True)
class SegmentationRefinementResult:
    """Refined integer label image and per-rule audit masks/counts."""

    image: nib.Nifti1Image
    candidate_masks: Mapping[str, np.ndarray]
    changed_masks: Mapping[str, np.ndarray]
    changed_voxels: Mapping[str, int]
    unresolved_voxels: Mapping[str, int]
    sweeps: Mapping[str, int]
    threshold_voxels: Mapping[str, float]


def _prepare_refinement_inputs(
    segmentation: nib.spatialimages.SpatialImage,
    inner_levelset: nib.spatialimages.SpatialImage,
    outer_levelset: nib.spatialimages.SpatialImage,
) -> tuple[nib.spatialimages.SpatialImage, np.ndarray, np.ndarray, np.ndarray]:
    images = []
    for image, name in (
        (segmentation, "MaSeg segmentation"),
        (inner_levelset, "inner level set"),
        (outer_levelset, "outer level set"),
    ):
        if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
            raise ValueError(f"{name} must be a NIfTI image")
        if len(image.shape) != 3:
            raise ValueError(f"{name} must be a 3D image")
        images.append(normalize_image_to_ras(image, source_name=name))
    reference = images[0]
    if any(
        image.shape != reference.shape
        or not np.allclose(image.affine, reference.affine, rtol=0, atol=1e-4)
        for image in images[1:]
    ):
        raise ValueError("Segmentation and level sets must share the same RAS grid")

    label_values = np.asarray(reference.dataobj)
    if not np.issubdtype(label_values.dtype, np.integer):
        if not np.isfinite(label_values).all() or not np.equal(label_values, np.rint(label_values)).all():
            raise ValueError("MaSeg segmentation must contain finite integer labels")
        label_values = label_values.astype(np.int32)
    inner = np.asarray(images[1].dataobj, dtype=np.float32)
    outer = np.asarray(images[2].dataobj, dtype=np.float32)
    if not np.isfinite(inner).all() or not np.isfinite(outer).all():
        raise ValueError("Segmentation level sets must contain finite values")
    return reference, label_values.copy(), inner, outer


def _fill_by_local_mode(
    labels: np.ndarray,
    pending: np.ndarray,
    candidate_rule: str,
    *,
    radius: int,
    max_sweeps: int,
    source_faithful: bool,
) -> tuple[np.ndarray, int, int, int]:
    """Sequentially fill pending voxels from a clipped local modal label.

    The source visits MATLAB column-major linear indices and processes updates
    in place, so later voxels in a sweep can use earlier assignments. Ties are
    resolved to the numerically lowest candidate label, as MATLAB ``tabulate``
    sorts its distinct values and ``max`` selects the first maximum.
    """
    pending = np.asarray(pending, dtype=bool).copy()
    initial_count = int(pending.sum())
    sweep_count = 0
    changed = 0
    for sweep in range(max_sweeps):
        indices = np.flatnonzero(pending.ravel(order="F"))
        if indices.size == 0:
            break
        sweep_count = sweep + 1
        changed_this_sweep = 0
        for linear_index in indices:
            coordinate = np.unravel_index(int(linear_index), labels.shape, order="F")
            bounds = tuple(
                slice(max(0, value - radius), min(size, value + radius + 1))
                for value, size in zip(coordinate, labels.shape, strict=True)
            )
            neighborhood = labels[bounds]
            if candidate_rule == "white_matter":
                candidate_mask = (neighborhood == 44) | (neighborhood == 45)
            elif candidate_rule == "cortical":
                candidate_mask = neighborhood >= 100
            else:
                raise ValueError(f"Unknown refinement candidate rule: {candidate_rule}")
            values = neighborhood[candidate_mask]
            if values.size == 0:
                continue
            modes, counts = np.unique(values, return_counts=True)
            mode_index = int(np.flatnonzero(counts == counts.max())[0])
            if source_faithful:
                # MATLAB assigns the 1-based tabulate row index (`maxind`)
                # instead of the label in column 1 (`foundtable(maxind,1)`).
                mode = mode_index + 1
            else:
                mode = int(modes[mode_index])
            if mode != labels[coordinate]:
                labels[coordinate] = mode
                pending[coordinate] = False
                changed += 1
                changed_this_sweep += 1
        if changed_this_sweep == 0:
            break
    return labels, changed, int(pending.sum()), sweep_count


def refine_segmentation_by_surfaces(
    segmentation: nib.spatialimages.SpatialImage,
    inner_levelset: nib.spatialimages.SpatialImage,
    outer_levelset: nib.spatialimages.SpatialImage,
    *,
    inner_coefficient_mm: float = 0.2,
    outer_coefficient_mm: float = 0.2,
    inner_central_coefficient_mm: float = 0.2,
    outer_central_coefficient_mm: float = 0.2,
    deal_with_central: bool = True,
    neighborhood_radius: int = 2,
    max_sweeps: int = 50,
    source_faithful: bool = True,
) -> SegmentationRefinementResult:
    """Refine MaSeg labels around the CRUISE inner/outer surfaces.

    ``source_faithful=True`` is the default and preserves the historical
    ``maxind`` row-index assignment in MATLAB. Set it to ``False`` to assign
    the actual modal label; that intentional correction remains opt-in pending
    project-owner sign-off. The source converts the four millimetre coefficients to voxel values by
    dividing by the largest segmentation voxel spacing. It changes labels
    ``>=100`` on the inner side to local-mode WM labels 44/45, clears remaining
    labels ``>=100`` beyond the outer surface, and optionally fills central
    labels 0/44/45 from the local mode of neighboring ``>=100`` cortical
    labels. Each rule and any unresolved candidate voxels are returned.
    """
    coefficients = (
        inner_coefficient_mm,
        outer_coefficient_mm,
        inner_central_coefficient_mm,
        outer_central_coefficient_mm,
    )
    if any(not np.isfinite(value) or value < 0 for value in coefficients):
        raise ValueError("Refinement coefficients must be finite and nonnegative")
    if neighborhood_radius < 0 or max_sweeps < 1:
        raise ValueError("Neighborhood radius must be nonnegative and max_sweeps positive")

    seg_image, original_labels, inner, outer = _prepare_refinement_inputs(
        segmentation, inner_levelset, outer_levelset
    )
    labels = original_labels.copy()
    max_spacing = float(np.max(voxel_sizes(seg_image.affine)))
    if not np.isfinite(max_spacing) or max_spacing <= 0:
        raise ValueError("Segmentation grid has invalid voxel spacing")
    thresholds = {
        "inner_wm": inner_coefficient_mm / max_spacing,
        "outer_gm": outer_coefficient_mm / max_spacing,
        "central_wm": inner_central_coefficient_mm / max_spacing,
        "outer_central": outer_central_coefficient_mm / max_spacing,
        "central_gm": outer_central_coefficient_mm / max_spacing,
    }

    candidates = {
        "inner_wm": (labels >= 100) & (inner < -thresholds["inner_wm"]),
    }
    masks: dict[str, np.ndarray] = {name: mask.copy() for name, mask in candidates.items()}
    changed_masks: dict[str, np.ndarray] = {}
    changed_voxels: dict[str, int] = {}
    unresolved: dict[str, int] = {}
    sweeps: dict[str, int] = {}

    before = labels.copy()
    labels, changed, left, n_sweeps = _fill_by_local_mode(
        labels,
        candidates["inner_wm"],
        "white_matter",
        radius=neighborhood_radius,
        max_sweeps=max_sweeps,
        source_faithful=source_faithful,
    )
    changed_masks["inner_wm"] = labels != before
    changed_voxels["inner_wm"] = changed
    unresolved["inner_wm"] = left
    sweeps["inner_wm"] = n_sweeps

    candidates["outer_clear"] = (labels >= 100) & (outer > thresholds["outer_gm"])
    masks["outer_clear"] = candidates["outer_clear"].copy()
    changed_masks["outer_clear"] = candidates["outer_clear"].copy()
    changed_voxels["outer_clear"] = int(candidates["outer_clear"].sum())
    labels[candidates["outer_clear"]] = 0
    unresolved["outer_clear"] = 0
    sweeps["outer_clear"] = 1

    if deal_with_central:
        candidates["central_wm"] = ((labels == 44) | (labels == 45)) & (
            inner > thresholds["central_wm"]
        )
        candidates["central_gm"] = (
            (labels == 0)
            & (outer < -thresholds["outer_central"])
            & (inner > thresholds["central_wm"])
        )
        central_pending = candidates["central_wm"] | candidates["central_gm"]
        masks["central_wm"] = candidates["central_wm"].copy()
        masks["central_gm"] = candidates["central_gm"].copy()
        before = labels.copy()
        labels, changed, left, n_sweeps = _fill_by_local_mode(
            labels,
            central_pending,
            "cortical",
            radius=neighborhood_radius,
            max_sweeps=max_sweeps,
            source_faithful=source_faithful,
        )
        changed_masks["central_fill"] = labels != before
        changed_voxels["central_fill"] = changed
        unresolved["central_fill"] = left
        sweeps["central_fill"] = n_sweeps
    else:
        candidates["central_wm"] = np.zeros(labels.shape, dtype=bool)
        candidates["central_gm"] = np.zeros(labels.shape, dtype=bool)
        masks["central_wm"] = candidates["central_wm"].copy()
        masks["central_gm"] = candidates["central_gm"].copy()
        changed_masks["central_fill"] = np.zeros(labels.shape, dtype=bool)
        changed_voxels["central_fill"] = 0
        unresolved["central_fill"] = 0
        sweeps["central_fill"] = 0

    header = seg_image.header.copy()
    header.set_data_dtype(original_labels.dtype)
    header.set_slope_inter(None, None)
    output = nib.Nifti1Image(labels.astype(original_labels.dtype, copy=False), seg_image.affine, header)
    qform, qcode = seg_image.get_qform(coded=True)
    _, scode = seg_image.get_sform(coded=True)
    output.set_qform(qform if qcode else None, code=int(qcode), strip_shears=True)
    output.set_sform(seg_image.affine, code=int(scode) if scode else 1)
    return SegmentationRefinementResult(
        output,
        masks,
        changed_masks,
        changed_voxels,
        unresolved,
        sweeps,
        thresholds,
    )
