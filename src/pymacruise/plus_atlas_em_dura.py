"""Compose AtlasEM products with the selected pre-LesionTOADS Dura stage."""
from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Sequence

import nibabel as nib
import numpy as np
from numpy.typing import NDArray

from pymacruise.atlas_em_outputs import (
    generate_atlas_em_dura_inputs,
    generate_atlas_em_segmentation,
)
from pymacruise.image_io import normalize_image_to_ras
from pymacruise.toads_dura import DuraRemovalResult, run_dura_removal


@dataclass(frozen=True)
class AtlasEMPreLesionDuraResult:
    products: nib.Nifti1Image
    segmentation: NDArray[np.int8]
    dura: DuraRemovalResult


def run_atlas_em_prelesion_dura(
    original_t1: nib.spatialimages.SpatialImage,
    memberships: nib.spatialimages.SpatialImage,
    brain_mask: nib.spatialimages.SpatialImage,
    class_names: Sequence[str],
    *,
    max_iterations: int = 5,
    record_passes: bool = False,
) -> AtlasEMPreLesionDuraResult:
    """Run the selected AtlasEM-output → pre-LesionTOADS Dura handoff.

    The class-axis membership NIfTI, T1, and supplied brain mask must share
    the same physical grid (they are axis-reoriented to RAS without
    interpolation first). The source AtlasEM adapter generates the GM, WM,
    and WM-mask products; the selected RemoveDura wrapper then consumes them
    with target/membership level 0.3 and V3 WM/GM topology settings 0.05/0.02.
    Atlas registration, AtlasEM iteration, and descriptor parsing happen
    upstream and are not part of this adapter.
    """
    for name, image in (("original T1", original_t1), ("AtlasEM memberships", memberships), ("brain mask", brain_mask)):
        if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
            raise ValueError(f"{name} must be a NIfTI image")
    t1 = normalize_image_to_ras(original_t1, source_name="pre-LesionTOADS original T1")
    mem_image = normalize_image_to_ras(memberships, source_name="AtlasEM memberships")
    mask_image = normalize_image_to_ras(brain_mask, source_name="AtlasEM brain mask")
    t1_values = np.asanyarray(t1.dataobj)
    mem_values = np.asanyarray(mem_image.dataobj)
    mask_values = np.asanyarray(mask_image.dataobj)
    if t1_values.ndim != 3 or mem_values.ndim != 4 or mask_values.ndim != 3:
        raise ValueError("T1 and brain mask must be 3D; AtlasEM memberships must be 4D with classes last")
    if tuple(mem_values.shape[:3]) != t1.shape or mask_values.shape != t1.shape:
        raise ValueError("T1, AtlasEM memberships, and brain mask shapes must match on the RAS grid")
    for name, image in (("AtlasEM memberships", mem_image), ("brain mask", mask_image)):
        if not np.allclose(image.affine, t1.affine, rtol=0, atol=1e-4):
            raise ValueError(f"{name} must match the original T1 physical RAS grid")
    if not np.issubdtype(t1_values.dtype, np.number) or not np.isfinite(t1_values).all():
        raise ValueError("original T1 must contain finite numeric voxels")
    if not np.issubdtype(mem_values.dtype, np.number) or not np.isfinite(mem_values).all():
        raise ValueError("AtlasEM memberships must contain finite numeric voxels")
    if mask_values.dtype != np.bool_ and not np.isin(mask_values, (0, 1)).all():
        raise ValueError("brain mask must be binary")
    names = tuple(class_names)
    if len(names) != mem_values.shape[3] or any(not isinstance(name, str) for name in names):
        raise ValueError("class_names must provide one string per AtlasEM membership class")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations < 0:
        raise ValueError("max_iterations must be a nonnegative integer")

    mask = mask_values.astype(bool, copy=False)
    segmentation = generate_atlas_em_segmentation(mem_values)
    products = generate_atlas_em_dura_inputs(mem_values, segmentation, mask, names)
    header = t1.header.copy()
    header.set_data_dtype(np.float32)
    product_image = nib.Nifti1Image(products, t1.affine, header=header)
    _, qcode = t1.get_qform(coded=True)
    _, scode = t1.get_sform(coded=True)
    product_image.set_qform(t1.affine if qcode else None, code=int(qcode), strip_shears=True)
    product_image.set_sform(t1.affine if scode else None, code=int(scode))
    product_image.set_filename("atlasem_dura_inputs.nii.gz")

    dura = run_dura_removal(
        t1,
        nib.Nifti1Image(products[..., 1], t1.affine, header=header),
        nib.Nifti1Image(products[..., 0], t1.affine, header=header),
        nib.Nifti1Image(products[..., 2], t1.affine, header=header),
        membership_level=0.3,
        wm_topology_regularization=0.05,
        gm_topology_regularization=0.02,
        curvature_force=0.5,
        alignment_force=1.0,
        pressure_force=0.5,
        alignment_alpha=0.5,
        max_iterations=max_iterations,
        mask_translation=2.0,
        level_set_translation=0.0,
        record_passes=record_passes,
    )
    return AtlasEMPreLesionDuraResult(product_image, segmentation, dura)
