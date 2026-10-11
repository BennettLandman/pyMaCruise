"""Input validation for the selected MaCRUISE workflow profiles."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import normalize_image_to_ras


class WorkflowProfile(str, Enum):
    """Supported MaCRUISE runner profiles."""

    BASE = "base"
    PLUS = "plus"


class WorkflowInputError(ValueError):
    """Raised when workflow inputs violate the documented contract."""


def has_finite_integer_values(data: np.ndarray) -> bool:
    """Return whether an array stores finite integer-valued labels.

    NIfTI label maps are sometimes stored as float types even when every voxel
    is an exact integer (for example, source MaSeg exports). Validate label
    values rather than requiring a particular storage dtype.
    """
    if np.issubdtype(data.dtype, np.integer):
        return True
    return bool(np.isfinite(data).all() and np.equal(data, np.floor(data)).all())


@dataclass(frozen=True)
class ValidatedInputs:
    """Loaded, RAS-oriented images and original orientation metadata."""

    profile: WorkflowProfile
    t1: nib.spatialimages.SpatialImage
    maseg: nib.spatialimages.SpatialImage
    brain_mask: nib.spatialimages.SpatialImage
    flair: nib.spatialimages.SpatialImage | None
    input_orientations: dict[str, tuple[str, str, str]]


def _validate_image(
    image: nib.spatialimages.SpatialImage,
    *,
    name: str,
    labels: bool = False,
) -> nib.spatialimages.SpatialImage:
    if len(image.shape) != 3:
        raise WorkflowInputError(f"{name} must be a 3D single-volume image; got {image.shape}")
    normalized = normalize_image_to_ras(image, source_name=name)
    data = np.asanyarray(normalized.dataobj)
    if not np.isfinite(data).all():
        raise WorkflowInputError(f"{name} contains NaN or infinite voxel values")
    if labels and not has_finite_integer_values(data):
        raise WorkflowInputError(f"MaSeg labels must have integer-valued voxels: {name}")
    return normalized


def validate_workflow_inputs(
    profile: WorkflowProfile | str,
    t1: nib.spatialimages.SpatialImage,
    maseg: nib.spatialimages.SpatialImage,
    flair: nib.spatialimages.SpatialImage | None = None,
    *,
    brain_mask: nib.spatialimages.SpatialImage,
) -> ValidatedInputs:
    """Validate required images and return RAS-normalized inputs.

    MaSeg and a binary brain mask must share the T1 grid. FLAIR is checked as a
    valid 3D intensity image but may use another grid because explicit rigid
    registration is a later MaCRUISEplus stage.
    """
    try:
        selected_profile = WorkflowProfile(profile)
    except ValueError as exc:
        raise WorkflowInputError(f"Unsupported workflow profile: {profile!r}") from exc
    if selected_profile is WorkflowProfile.PLUS and flair is None:
        raise WorkflowInputError("The plus profile requires a FLAIR image")
    if selected_profile is WorkflowProfile.BASE and flair is not None:
        raise WorkflowInputError("The base profile does not accept a FLAIR input")

    original_orientations = {
        "t1": nib.aff2axcodes(t1.affine),
        "maseg": nib.aff2axcodes(maseg.affine),
        "brain_mask": nib.aff2axcodes(brain_mask.affine),
    }
    normalized_t1 = _validate_image(t1, name="T1")
    normalized_maseg = _validate_image(maseg, name="MaSeg", labels=True)
    normalized_brain_mask = _validate_image(brain_mask, name="brain mask", labels=True)
    if normalized_t1.shape != normalized_maseg.shape:
        raise WorkflowInputError(
            f"T1 and MaSeg shapes differ after orientation normalization: "
            f"{normalized_t1.shape} versus {normalized_maseg.shape}"
        )
    if not np.allclose(normalized_t1.affine, normalized_maseg.affine, rtol=0, atol=1e-4):
        raise WorkflowInputError("T1 and MaSeg are not aligned on the same physical grid")
    mask_data = np.asanyarray(normalized_brain_mask.dataobj)
    if not np.isin(mask_data, (0, 1, 255)).all():
        raise WorkflowInputError("Brain mask must be binary with values 0/1 or 0/255")
    if not np.any(mask_data):
        raise WorkflowInputError("Brain mask must contain at least one foreground voxel")
    if normalized_t1.shape != normalized_brain_mask.shape:
        raise WorkflowInputError("T1 and brain mask shapes differ after RAS normalization")
    if not np.allclose(
        normalized_t1.affine, normalized_brain_mask.affine, rtol=0, atol=1e-4
    ):
        raise WorkflowInputError("Brain mask and T1 are not aligned on the same physical grid")

    normalized_flair = None
    if flair is not None:
        original_orientations["flair"] = nib.aff2axcodes(flair.affine)
        normalized_flair = _validate_image(flair, name="FLAIR")

    return ValidatedInputs(
        profile=selected_profile,
        t1=normalized_t1,
        maseg=normalized_maseg,
        brain_mask=normalized_brain_mask,
        flair=normalized_flair,
        input_orientations=original_orientations,
    )


def validate_workflow_files(
    profile: WorkflowProfile | str,
    t1_path: str | Path,
    maseg_path: str | Path,
    flair_path: str | Path | None = None,
    *,
    brain_mask_path: str | Path,
) -> ValidatedInputs:
    """Read and validate compressed NIfTI inputs with path-specific errors."""
    selected_profile = WorkflowProfile(profile)
    paths = {
        "T1": Path(t1_path),
        "MaSeg": Path(maseg_path),
        "brain mask": Path(brain_mask_path),
    }
    if flair_path is not None:
        paths["FLAIR"] = Path(flair_path)
    images: dict[str, nib.spatialimages.SpatialImage] = {}
    for name, path in paths.items():
        if not path.is_file():
            raise WorkflowInputError(f"Missing {name} input file: {path}")
        if not path.name.lower().endswith(".nii.gz"):
            raise WorkflowInputError(f"{name} input must use .nii.gz: {path}")
        try:
            images[name] = nib.load(str(path))
        except (OSError, nib.filebasedimages.ImageFileError) as exc:
            raise WorkflowInputError(f"Cannot load {name} image {path}: {exc}") from exc
    return validate_workflow_inputs(
        selected_profile,
        images["T1"],
        images["MaSeg"],
        images.get("FLAIR"),
        brain_mask=images["brain mask"],
    )
