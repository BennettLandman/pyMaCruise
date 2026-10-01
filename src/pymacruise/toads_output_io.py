"""Create and save selected base-TOADS output images on the native T1 grid."""

from __future__ import annotations

from pathlib import Path
import re
from typing import Literal

import nibabel as nib
import numpy as np
from numpy.typing import ArrayLike

from .image_io import normalize_image_to_ras, nifti_ras_on_reference, save_nifti_ras
from .toads_cruise_outputs import ToadsCruiseProducts


ToadsOutputMode = Literal[
    "hard_segmentation", "full_segmentation", "Cruise_input", "DuraRemoval_input", "thalamus_seg"
]
_CHANNELS = (
    ("csf", "_csf.nii.gz"),
    ("gray_matter", "_gm.nii.gz"),
    ("white_matter_fill", "_wmfill.nii.gz"),
    ("white_matter_mask", "_wmmask.nii.gz"),
)


def toads_output_images(
    t1_reference: nib.spatialimages.SpatialImage,
    classification: ArrayLike,
    *,
    result_mode: ToadsOutputMode,
    memberships: ArrayLike | None = None,
    cruise_products: ToadsCruiseProducts | None = None,
    thalamus: ArrayLike | None = None,
    inhomogeneity_field: ArrayLike | None = None,
) -> dict[str, nib.spatialimages.SpatialImage]:
    """Build typed NIfTI outputs from already-restored RAS arrays.

    Only the selected wrapper's primary output modes are exposed. The debug
    ``all_images`` mode is deliberately excluded. Arrays must already be on
    the full, RAS T1 grid; crop restoration belongs to the upstream adapter.
    """
    if result_mode not in {
        "hard_segmentation", "full_segmentation", "Cruise_input",
        "DuraRemoval_input", "thalamus_seg",
    }:
        raise ValueError(f"unsupported base TOADS output mode: {result_mode!r}")
    if memberships is not None and result_mode != "full_segmentation":
        raise ValueError("memberships are only exported by full_segmentation")
    if cruise_products is not None and result_mode not in {"Cruise_input", "DuraRemoval_input"}:
        raise ValueError("cruise_products require Cruise_input or DuraRemoval_input")
    if thalamus is not None and result_mode != "thalamus_seg":
        raise ValueError("thalamus is only exported by thalamus_seg")
    if inhomogeneity_field is not None and result_mode != "full_segmentation":
        raise ValueError("inhomogeneity_field is only exported by full_segmentation")
    if not isinstance(t1_reference, (nib.Nifti1Image, nib.Nifti2Image)):
        raise TypeError("t1_reference must be a NIfTI image")
    reference = normalize_image_to_ras(t1_reference, source_name="T1 output reference")
    if len(reference.shape) != 3:
        raise ValueError("T1 output reference must be a 3-D image")
    shape = tuple(reference.shape)
    labels = np.asarray(classification)
    if labels.shape != shape or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("classification must be an integer map matching the T1 reference grid")

    images: dict[str, nib.spatialimages.SpatialImage] = {
        "classification": nifti_ras_on_reference(reference, labels, dtype=np.uint8)
    }
    if result_mode == "full_segmentation":
        if memberships is None:
            raise ValueError("full_segmentation requires memberships")
        mem = np.asarray(memberships)
        if mem.ndim != 4 or mem.shape[:3] != shape or mem.shape[3] < 1:
            raise ValueError("memberships must have shape [x, y, z, classes] on the T1 grid")
        images["memberships"] = nifti_ras_on_reference(reference, mem, dtype=np.float32)
    elif result_mode in {"Cruise_input", "DuraRemoval_input"}:
        if not isinstance(cruise_products, ToadsCruiseProducts):
            raise ValueError(f"{result_mode} requires ToadsCruiseProducts")
        channels = cruise_products.as_channels()
        if channels.shape != (4, *shape):
            raise ValueError("CRUISE products must contain four maps matching the T1 grid")
        for index, (name, _) in enumerate(_CHANNELS):
            images[name] = nifti_ras_on_reference(
                reference, channels[index], dtype=np.float32
            )
    elif result_mode == "thalamus_seg":
        if thalamus is None:
            raise ValueError("thalamus_seg requires a thalamus-membership map")
        values = np.asarray(thalamus)
        if values.shape == (1, *shape):
            values = values[0]
        if values.shape != shape:
            raise ValueError("thalamus map must match the T1 reference grid")
        images["thalamus"] = nifti_ras_on_reference(
            reference, values, dtype=np.float32
        )

    if inhomogeneity_field is not None:
        field = np.asarray(inhomogeneity_field, dtype=np.float32)
        if field.shape == shape:
            images["inhomogeneity_field"] = nifti_ras_on_reference(
                reference, field, dtype=np.float32
            )
        elif field.ndim == 4 and field.shape[:3] == shape and field.shape[3] > 0:
            images["inhomogeneity_field"] = nifti_ras_on_reference(
                reference, field, dtype=np.float32
            )
        else:
            raise ValueError("inhomogeneity_field must be [x,y,z] or [x,y,z,channels]")
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


def save_toads_outputs(
    t1_reference: nib.spatialimages.SpatialImage,
    classification: ArrayLike,
    output_directory: str | Path,
    *,
    input_name: str | Path,
    result_mode: ToadsOutputMode,
    memberships: ArrayLike | None = None,
    cruise_products: ToadsCruiseProducts | None = None,
    thalamus: ArrayLike | None = None,
    inhomogeneity_field: ArrayLike | None = None,
) -> dict[str, Path]:
    """Save selected base-TOADS output images as compressed RAS NIfTI files."""
    images = toads_output_images(
        t1_reference,
        classification,
        result_mode=result_mode,
        memberships=memberships,
        cruise_products=cruise_products,
        thalamus=thalamus,
        inhomogeneity_field=inhomogeneity_field,
    )
    stem = _safe_stem(input_name)
    suffixes = {
        "classification": "_Toads_seg.nii.gz",
        "memberships": "_Toads_mems.nii.gz",
        "thalamus": "_thalamus.nii.gz",
        "inhomogeneity_field": "_Toads_field.nii.gz",
    }
    suffixes.update(dict(_CHANNELS))
    destinations: dict[str, Path] = {}
    for name, image in images.items():
        suffix = suffixes[name]
        if name == "inhomogeneity_field" and image.ndim == 4:
            suffix = "_Toads_fields.nii.gz"
        destinations[name] = Path(output_directory) / f"{stem}{suffix}"
    collisions = [str(path) for path in destinations.values() if path.exists()]
    if collisions:
        raise FileExistsError(f"refusing to overwrite existing output(s): {', '.join(collisions)}")
    for name, image in images.items():
        save_nifti_ras(image, destinations[name])
    return destinations
