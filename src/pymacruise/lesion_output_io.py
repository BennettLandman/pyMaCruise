"""Create and save selected LesionTOADS output images on the native T1 grid."""

from __future__ import annotations

from pathlib import Path
import re

import nibabel as nib
import numpy as np
from numpy.typing import ArrayLike

from .image_io import normalize_image_to_ras, save_nifti_ras
from .lesion_post_loop_outputs import LesionToadsPostLoopOutputs


_CHANNEL_NAMES = ("csf", "gm", "wmFill", "wmMask")
_PRODUCT_FILENAMES = {
    "classification": "_lToads_class.nii.gz",
    "classification_from_membership": "_lToads_class_max_mem.nii.gz",
    "lesion_subtype": "_lToads_lesions_hard.nii.gz",
    "memberships_plus_lesion": "_lToads_mems.nii.gz",
    "inhomogeneity_field": "_lToads_field.nii.gz",
}


def _make_image(
    reference: nib.spatialimages.SpatialImage,
    values: np.ndarray,
    *,
    dtype: np.dtype,
) -> nib.spatialimages.SpatialImage:
    data = np.asarray(values).astype(dtype, copy=False)
    if isinstance(reference, nib.Nifti2Image):
        image_type = nib.Nifti2Image
    else:
        image_type = nib.Nifti1Image
    header = reference.header.copy()
    header.set_data_dtype(data.dtype)
    header.set_slope_inter(None, None)
    output = image_type(data, reference.affine, header=header)
    qform, qcode = reference.get_qform(coded=True)
    _, scode = reference.get_sform(coded=True)
    if qcode:
        output.set_qform(qform, code=int(qcode), strip_shears=True)
    else:
        output.set_qform(None, code=0)
    output.set_sform(reference.affine, code=int(scode) if scode else 1)
    return output


def lesion_toads_output_images(
    outputs: LesionToadsPostLoopOutputs,
    t1_reference: nib.spatialimages.SpatialImage,
    *,
    inhomogeneity_field: ArrayLike | None = None,
) -> dict[str, nib.spatialimages.SpatialImage]:
    """Build typed RAS NIfTI images for already-uncropped output arrays.

    Arrays from the selected output adapter are already in RAS voxel order and
    on the original T1 grid. The reference image supplies only physical
    geometry and non-spatial header metadata. UBYTE output preserves Java's
    signed-byte bit pattern while representing it as NIfTI uint8. Supply an
    already-restored field only when both source field-output options are on;
    its 3-D/4-D shape selects the source `_field`/`_fields` output form.
    """
    if not isinstance(outputs, LesionToadsPostLoopOutputs):
        raise TypeError("outputs must be LesionToadsPostLoopOutputs")
    if not isinstance(t1_reference, (nib.Nifti1Image, nib.Nifti2Image)):
        raise TypeError("t1_reference must be a NIfTI image")
    reference = normalize_image_to_ras(t1_reference, source_name="T1 output reference")
    shape = tuple(reference.shape[:3])
    if len(reference.shape) != 3:
        raise ValueError("T1 output reference must be a 3-D image")

    images: dict[str, nib.spatialimages.SpatialImage] = {}
    for product in outputs.products:
        values = np.asarray(product.values)
        if product.name in {"classification", "classification_from_membership"}:
            if values.shape != shape:
                raise ValueError(f"{product.name} must match the T1 reference grid")
            images[product.name] = _make_image(reference, values, dtype=np.dtype(np.uint8))
        elif product.name == "lesion_subtype":
            if values.shape != shape:
                raise ValueError("lesion_subtype must match the T1 reference grid")
            images[product.name] = _make_image(reference, values, dtype=np.dtype(np.int16))
        elif product.name == "memberships_plus_lesion":
            if values.ndim != 4 or values.shape[:3] != shape:
                raise ValueError("memberships_plus_lesion must be shaped [x, y, z, channels] on the T1 grid")
            images[product.name] = _make_image(reference, values, dtype=np.dtype(np.float32))
        elif product.name in {"cruise_inputs", "dura_removal_inputs"}:
            if values.shape != (4, *shape):
                raise ValueError(f"{product.name} must be shaped [4, x, y, z] on the T1 grid")
            for index, channel_name in enumerate(_CHANNEL_NAMES):
                images[channel_name] = _make_image(
                    reference, values[index], dtype=np.dtype(np.float32)
                )
        else:
            raise ValueError(f"unsupported LesionTOADS output product: {product.name}")

    if inhomogeneity_field is not None:
        field = np.asarray(inhomogeneity_field, dtype=np.float32)
        if field.ndim == 3 and field.shape == shape:
            images["inhomogeneity_field"] = _make_image(
                reference, field, dtype=np.dtype(np.float32)
            )
        elif field.ndim == 4 and field.shape[:3] == shape and field.shape[3] > 0:
            images["inhomogeneity_field"] = _make_image(
                reference, field, dtype=np.dtype(np.float32)
            )
        else:
            raise ValueError(
                "inhomogeneity_field must be [x, y, z] or [x, y, z, channels] on the T1 grid"
            )

    expected = {"classification", "lesion_subtype"}
    if outputs.membership_classification is not None:
        expected.add("classification_from_membership")
    if outputs.result_mode == "class+member":
        expected.add("memberships_plus_lesion")
    elif outputs.result_mode in {"cruise inputs", "dura removal inputs"}:
        expected.update(_CHANNEL_NAMES)
    if inhomogeneity_field is not None:
        expected.add("inhomogeneity_field")
    if set(images) != expected:
        raise ValueError(
            f"output products do not match selected mode: expected {sorted(expected)}, got {sorted(images)}"
        )
    return images


def _safe_stem(input_name: str | Path) -> str:
    name = Path(input_name).name
    if name.lower().endswith(".nii.gz"):
        name = name[:-7]
    else:
        name = Path(name).stem
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._-")
    if not name:
        raise ValueError("input_name must contain a usable filename stem")
    return name


def save_lesion_toads_outputs(
    outputs: LesionToadsPostLoopOutputs,
    t1_reference: nib.spatialimages.SpatialImage,
    output_directory: str | Path,
    *,
    input_name: str | Path,
    inhomogeneity_field: ArrayLike | None = None,
) -> dict[str, Path]:
    """Save selected output products as compressed NIfTI files in RAS.

    Existing paths are rejected to avoid silently overwriting prior results.
    Product names and types follow the selected VandyCRUISE wrapper branch.
    An optional field must already be restored to the original T1 grid.
    """
    images = lesion_toads_output_images(
        outputs, t1_reference, inhomogeneity_field=inhomogeneity_field
    )
    stem = _safe_stem(input_name)
    directory = Path(output_directory)
    suffixes = dict(_PRODUCT_FILENAMES)
    for name, suffix in zip(_CHANNEL_NAMES, ("_lToads_csf.nii.gz", "_lToads_gm.nii.gz", "_lToads_wmFill.nii.gz", "_lToads_wmMask.nii.gz"), strict=True):
        suffixes[name] = suffix
    destinations = {}
    for name, image in images.items():
        suffix = suffixes[name]
        if name == "inhomogeneity_field" and image.ndim == 4:
            suffix = "_lToads_fields.nii.gz"
        destinations[name] = directory / f"{stem}{suffix}"
    collisions = [str(path) for path in destinations.values() if path.exists()]
    if collisions:
        raise FileExistsError(f"refusing to overwrite existing output(s): {', '.join(collisions)}")
    for name, image in images.items():
        save_nifti_ras(image, destinations[name])
    return destinations
