"""Python equivalent of the selected workflow's `prepare_raw` stage."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import nibabel as nib

from pymacruise.contracts import validate_workflow_inputs
from pymacruise.image_io import save_nifti_ras


@dataclass(frozen=True)
class RawPreparationResult:
    """Prepared T1 and MaSeg files plus original input axis codes."""

    t1_path: Path
    maseg_path: Path
    input_orientations: dict[str, tuple[str, str, str]]


def _load_nifti_input(path: Path, role: str) -> nib.spatialimages.SpatialImage:
    if not (path.name.lower().endswith(".nii") or path.name.lower().endswith(".nii.gz")):
        raise ValueError(
            f"{role} must be a NIfTI image (.nii or .nii.gz); convert other formats first: {path}"
        )
    if not path.is_file():
        raise FileNotFoundError(f"Missing {role} input file: {path}")
    try:
        image = nib.load(str(path))
    except (OSError, nib.filebasedimages.ImageFileError) as exc:
        raise ValueError(f"Cannot load {role} NIfTI image {path}: {exc}") from exc
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError(f"{role} input is not a NIfTI image: {path}")
    return image


def prepare_raw(
    t1_input: str | Path,
    maseg_input: str | Path,
    output_dir: str | Path,
    *,
    overwrite: bool = False,
) -> RawPreparationResult:
    """Validate and persist `T1.nii.gz` and `Seg.nii.gz` in RAS.

    Unlike the MATLAB helper, existing outputs are never silently reused.
    Uncompressed `.nii` input is accepted and normalized to `.nii.gz`; other
    input formats must first be converted to NIfTI by an explicit converter.
    """
    t1_source = Path(t1_input)
    maseg_source = Path(maseg_input)
    destination = Path(output_dir)
    t1_output = destination / "T1.nii.gz"
    maseg_output = destination / "Seg.nii.gz"
    existing = [path for path in (t1_output, maseg_output) if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(
            "Prepared output already exists; choose another directory or pass overwrite=True: "
            + ", ".join(str(path) for path in existing)
        )

    t1 = _load_nifti_input(t1_source, "T1")
    maseg = _load_nifti_input(maseg_source, "MaSeg")
    validated = validate_workflow_inputs("base", t1, maseg)
    destination.mkdir(parents=True, exist_ok=True)
    save_nifti_ras(validated.t1, t1_output)
    save_nifti_ras(validated.maseg, maseg_output)
    return RawPreparationResult(
        t1_path=t1_output,
        maseg_path=maseg_output,
        input_orientations=validated.input_orientations,
    )
