"""Parse numeric priors from a source-format LesionTOADS atlas descriptor.

The numeric-prior parser reads structures, intensity rows, lesion priors, and
black-hole priors without loading atlas images. A separate focused loader can
resolve the topology NIfTI; it does not inspect or load shape maps. Prior
vectors retain the structure order declared in the descriptor.

Paper context: the selected numeric rows supply statistical intensity and
lesion/black-hole prior inputs for the Lesion-TOADS model described by Shiee et
al. (2010). This parser does not realize the paper's topology/shape atlas or
registration/segmentation steps. See ``docs/papers/lesion-toads.md`` for the
precise mapping and source-version caveat.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import nibabel as nib
import numpy as np

from .image_io import load_nifti_ras


class LesionAtlasDescriptorError(ValueError):
    """Raised when LesionTOADS descriptor metadata or atlas input is invalid."""


@dataclass(frozen=True)
class LesionAtlasStructure:
    name: str
    label: int
    topology: str


@dataclass(frozen=True)
class LesionAtlasPriors:
    """Selected modality priors, with each intensity row in source order."""

    structures: tuple[LesionAtlasStructure, ...]
    modalities: tuple[str, ...]
    intensity_priors: np.ndarray
    lesion_priors: np.ndarray
    black_hole_priors: np.ndarray


@dataclass(frozen=True)
class LesionAtlasTopology:
    """Topology NIfTI referenced by a descriptor; no shape maps are loaded."""

    path: Path
    image: nib.spatialimages.SpatialImage


@dataclass(frozen=True)
class LesionAtlasShapePriors:
    """Descriptor-ordered shape maps as RAS float32 voxel arrays."""

    structures: tuple[LesionAtlasStructure, ...]
    shape_maps: dict[str, np.ndarray]
    files: dict[str, Path]


_MODALITY_ALIASES = {
    "T1": "T1_SPGR",
    "T1_SPGR": "T1_SPGR",
    "T2": "T2",
    "FLAIR": "FLAIR",
    "MPRAGE": "T1_MPRAGE",
    "T1_MPRAGE": "T1_MPRAGE",
    "MPRAGE_RAW": "T1_RAW",
    "T1_RAW": "T1_RAW",
    "PD": "PD",
    "PDFSE": "PDFSE",
    "DIR": "DIR",
}

_SECTION_HEADINGS = (
    "Topology Atlas",
    "Shape Atlas",
    "Intensity Atlas",
    "Intensity Variance",
    "Lesions",
    "Black Holes",
    "Slope",
    "Modality weights",
    "Registered Shapes",
    "Optimized factors",
)


def load_lesion_atlas_priors(
    descriptor_path: str | Path,
    modalities: Sequence[str],
) -> LesionAtlasPriors:
    """Load selected modality intensity, lesion, and black-hole priors.

    ``intensity_priors[i]`` maps directly to the ``i``th declared structure.
    The requested modality order is retained in all returned arrays. Names and
    aliases follow ``DemonToadDeformableAtlas.modalityId``. Requested modality
    rows and both scalar priors must be present in the descriptor; missing
    values are rejected rather than silently replaced with the Java class's
    zero/label fallbacks.
    """
    path = Path(descriptor_path)
    if not path.is_file():
        raise LesionAtlasDescriptorError(f"Missing LesionTOADS atlas descriptor: {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].strip() != "Structure Atlas File (edit at your own risks)":
        raise LesionAtlasDescriptorError(f"Not a LesionTOADS atlas descriptor: {path}")

    requested = tuple(modalities)
    if not requested:
        raise LesionAtlasDescriptorError("At least one modality must be selected")
    canonical = tuple(_canonical_modality(item) for item in requested)
    if len(set(canonical)) != len(canonical):
        raise LesionAtlasDescriptorError("Selected modalities must be unique")

    structures = _parse_structures(lines)
    intensity_rows = _parse_intensity_rows(lines, len(structures))
    lesion_values = _parse_modality_values(lines, "Lesions")
    black_hole_values = _parse_modality_values(lines, "Black Holes")

    selected_intensity: list[np.ndarray] = []
    selected_lesion: list[np.float32] = []
    selected_black_hole: list[np.float32] = []
    for modality_name, modality_key in zip(requested, canonical, strict=True):
        if modality_key not in intensity_rows:
            raise LesionAtlasDescriptorError(
                f"No intensity row for selected modality {modality_name!r}"
            )
        if modality_key not in lesion_values:
            raise LesionAtlasDescriptorError(
                f"No lesion prior for selected modality {modality_name!r}"
            )
        if modality_key not in black_hole_values:
            raise LesionAtlasDescriptorError(
                f"No black-hole prior for selected modality {modality_name!r}"
            )
        selected_intensity.append(intensity_rows[modality_key])
        selected_lesion.append(lesion_values[modality_key])
        selected_black_hole.append(black_hole_values[modality_key])

    intensity_array = np.stack(selected_intensity).astype(np.float32, copy=False)
    lesion_array = np.asarray(selected_lesion, dtype=np.float32)
    black_hole_array = np.asarray(selected_black_hole, dtype=np.float32)
    for array in (intensity_array, lesion_array, black_hole_array):
        array.setflags(write=False)
    return LesionAtlasPriors(
        structures=structures,
        modalities=requested,
        intensity_priors=intensity_array,
        lesion_priors=lesion_array,
        black_hole_priors=black_hole_array,
    )


def load_lesion_atlas_topology(
    descriptor_path: str | Path,
) -> LesionAtlasTopology:
    """Load just the descriptor's topology atlas as an RAS ``.nii.gz`` image.

    Shape maps and the numeric priors remain untouched. The NIfTI is loaded via
    the shared image I/O path, which warns and axis-reorients non-RAS inputs.
    The descriptor reference must be relative to its directory and remain
    within that atlas directory.
    """
    descriptor = Path(descriptor_path).resolve()
    if not descriptor.is_file():
        raise LesionAtlasDescriptorError(
            f"Missing LesionTOADS atlas descriptor: {descriptor}"
        )
    lines = descriptor.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].strip() != "Structure Atlas File (edit at your own risks)":
        raise LesionAtlasDescriptorError(
            f"Not a LesionTOADS atlas descriptor: {descriptor}"
        )
    section_index = next(
        (i for i, line in enumerate(lines) if line.strip().startswith("Topology Atlas:")),
        None,
    )
    if section_index is None:
        raise LesionAtlasDescriptorError("Missing Topology Atlas section")
    reference = None
    declared_dimensions = None
    for line in lines[section_index + 1 :]:
        stripped = line.strip()
        if not stripped:
            continue
        if _starts_section(stripped):
            break
        key, separator, value = stripped.partition(":")
        if not separator:
            continue
        if key.strip() == "File":
            reference = value.strip()
        elif key.strip() == "Dimensions":
            try:
                parsed_dimensions = tuple(int(item) for item in value.split())
            except ValueError as exc:
                raise LesionAtlasDescriptorError(
                    f"Invalid topology atlas dimensions: {value}"
                ) from exc
            if len(parsed_dimensions) != 3 or any(size <= 0 for size in parsed_dimensions):
                raise LesionAtlasDescriptorError(
                    f"Topology atlas dimensions must be three positive integers: {value}"
                )
            declared_dimensions = parsed_dimensions
    if not reference:
        raise LesionAtlasDescriptorError("Topology Atlas section has no File reference")
    relative = Path(reference)
    if relative.is_absolute() or ".." in relative.parts:
        raise LesionAtlasDescriptorError(
            f"Unsafe topology atlas path in descriptor: {reference}"
        )
    if relative.suffixes[-2:] != [".nii", ".gz"]:
        raise LesionAtlasDescriptorError(
            f"Topology atlas must use compressed .nii.gz: {reference}"
        )
    atlas_path = (descriptor.parent / relative).resolve()
    if not atlas_path.is_relative_to(descriptor.parent):
        raise LesionAtlasDescriptorError(
            f"Topology atlas path escapes descriptor directory: {reference}"
        )
    try:
        image = load_nifti_ras(atlas_path)
    except (OSError, ValueError) as exc:
        raise LesionAtlasDescriptorError(
            f"Cannot load LesionTOADS topology atlas {atlas_path}: {exc}"
        ) from exc
    if len(image.shape) != 3:
        raise LesionAtlasDescriptorError(
            f"LesionTOADS topology atlas must be 3D, got {image.shape}: {atlas_path}"
        )
    if declared_dimensions is not None and tuple(image.shape) != declared_dimensions:
        raise LesionAtlasDescriptorError(
            "LesionTOADS topology atlas RAS shape "
            f"{image.shape} does not match descriptor dimensions {declared_dimensions}: "
            f"{atlas_path}"
        )
    return LesionAtlasTopology(atlas_path, image)


def load_lesion_atlas_shape_priors(
    descriptor_path: str | Path,
) -> LesionAtlasShapePriors:
    """Load all shape maps from the bundled descriptor's NIfTI mapping.

    Arrays follow the main Structures section order and are normalized using
    shared RAS-aware NIfTI I/O. Loading is intentionally separate from affine
    registration and transformed-prior sampling. This handles the bundled
    converted NIfTI schema, not the selected Java loader's raw float format.
    """
    descriptor = Path(descriptor_path).resolve()
    lines = _read_descriptor(descriptor)
    structures = _parse_structures(lines)
    dimensions = _parse_shape_dimensions(lines)
    asset_paths = lesion_atlas_shape_asset_paths(descriptor)
    arrays: dict[str, np.ndarray] = {}
    paths: dict[str, Path] = {}
    loaded: dict[Path, np.ndarray] = {}
    for structure in structures:
        atlas_path = asset_paths[structure.name]
        if atlas_path in loaded:
            array = loaded[atlas_path]
        else:
            try:
                image = load_nifti_ras(atlas_path)
            except (OSError, ValueError) as exc:
                raise LesionAtlasDescriptorError(
                    f"Cannot load LesionTOADS shape atlas for {structure.name} "
                    f"({atlas_path}): {exc}"
                ) from exc
            if len(image.shape) != 3:
                raise LesionAtlasDescriptorError(
                    f"LesionTOADS shape atlas for {structure.name} must be 3D, "
                    f"got {image.shape}: {atlas_path}"
                )
            if tuple(image.shape) != dimensions:
                raise LesionAtlasDescriptorError(
                    f"LesionTOADS shape atlas for {structure.name} has RAS shape "
                    f"{image.shape}; descriptor declares {dimensions}: {atlas_path}"
                )
            array = np.asarray(image.dataobj, dtype=np.float32)
            if not np.isfinite(array).all():
                raise LesionAtlasDescriptorError(
                    f"LesionTOADS shape atlas for {structure.name} contains non-finite values"
                )
            array.setflags(write=False)
            loaded[atlas_path] = array
        arrays[structure.name] = array
        paths[structure.name] = atlas_path
    return LesionAtlasShapePriors(structures, arrays, paths)


def lesion_atlas_shape_asset_paths(descriptor_path: str | Path) -> dict[str, Path]:
    """Resolve descriptor-owned shape images without loading their voxel data."""
    descriptor = Path(descriptor_path).resolve()
    lines = _read_descriptor(descriptor)
    structures = _parse_structures(lines)
    mapping = _parse_shape_files(lines, structures)
    paths: dict[str, Path] = {}
    for structure in structures:
        reference = mapping[structure.name]
        relative = Path(reference)
        if relative.is_absolute() or ".." in relative.parts:
            raise LesionAtlasDescriptorError(
                f"Unsafe shape atlas path for {structure.name}: {reference}"
            )
        if relative.suffixes[-2:] != [".nii", ".gz"]:
            raise LesionAtlasDescriptorError(
                f"Shape atlas must use compressed .nii.gz for {structure.name}: {reference}"
            )
        atlas_path = (descriptor.parent / relative).resolve()
        if not atlas_path.is_relative_to(descriptor.parent):
            raise LesionAtlasDescriptorError(
                f"Shape atlas path escapes descriptor directory for {structure.name}: {reference}"
            )
        if not atlas_path.is_file() or atlas_path.is_symlink():
            raise LesionAtlasDescriptorError(
                f"Cannot load LesionTOADS shape atlas for {structure.name}: "
                f"missing regular file {atlas_path}"
            )
        paths[structure.name] = atlas_path
    return paths


def _canonical_modality(name: str) -> str:
    try:
        return _MODALITY_ALIASES[name]
    except KeyError as exc:
        raise LesionAtlasDescriptorError(
            f"Unsupported LesionTOADS modality {name!r}"
        ) from exc


def _read_descriptor(path: Path) -> list[str]:
    if not path.is_file():
        raise LesionAtlasDescriptorError(
            f"Missing LesionTOADS atlas descriptor: {path}"
        )
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].strip() != "Structure Atlas File (edit at your own risks)":
        raise LesionAtlasDescriptorError(f"Not a LesionTOADS atlas descriptor: {path}")
    return lines


def _parse_structures(lines: list[str]) -> tuple[LesionAtlasStructure, ...]:
    header_index = _find_section(lines, "Structures")
    _, separator, value = lines[header_index].partition(":")
    if not separator:
        raise LesionAtlasDescriptorError("Malformed Structures declaration")
    try:
        expected = int(value.strip())
    except ValueError as exc:
        raise LesionAtlasDescriptorError("Invalid Structures count") from exc
    if expected <= 0:
        raise LesionAtlasDescriptorError("Structures count must be positive")

    rows: list[LesionAtlasStructure] = []
    for line in lines[header_index + 1 :]:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("Topology Atlas"):
            break
        parts = stripped.split()
        if len(parts) != 3:
            raise LesionAtlasDescriptorError(f"Malformed structure row: {line!r}")
        try:
            rows.append(LesionAtlasStructure(parts[0], int(parts[1]), parts[2]))
        except ValueError as exc:
            raise LesionAtlasDescriptorError(f"Malformed structure row: {line!r}") from exc
    if len(rows) != expected:
        raise LesionAtlasDescriptorError(
            f"Descriptor declares {expected} structures but lists {len(rows)}"
        )
    if len({row.name for row in rows}) != expected or len({row.label for row in rows}) != expected:
        raise LesionAtlasDescriptorError("Structure names and labels must be unique")
    return tuple(rows)


def _parse_shape_files(
    lines: list[str], structures: tuple[LesionAtlasStructure, ...]
) -> dict[str, str]:
    """Read the supplied NIfTI descriptor's Structure/File pairs."""
    header_index = _find_section(lines, "Shape Atlas")
    files: dict[str, str] = {}
    current: str | None = None
    structure_names = {row.name for row in structures}
    for line in lines[header_index + 1 :]:
        stripped = line.strip()
        if not stripped:
            continue
        if _starts_section(stripped):
            break
        key, separator, value = stripped.partition(":")
        if not separator:
            raise LesionAtlasDescriptorError(f"Malformed Shape Atlas row: {line!r}")
        key, value = key.strip(), value.strip()
        if key in {"Dimensions", "Resolutions", "Orientations"}:
            # These describe the original atlas grid; this loader preserves
            # each NIfTI's own RAS geometry and does not resample to it.
            continue
        if key == "Structure":
            if current is not None:
                raise LesionAtlasDescriptorError(
                    f"Shape Atlas has no File reference for {current}"
                )
            if not value or value not in structure_names:
                raise LesionAtlasDescriptorError(
                    f"Shape Atlas references invalid structure {value!r}"
                )
            if value in files:
                raise LesionAtlasDescriptorError(
                    f"Duplicate Shape Atlas structure {value!r}"
                )
            current = value
        elif key == "File":
            if current is None:
                raise LesionAtlasDescriptorError(
                    "Shape Atlas File reference has no preceding Structure"
                )
            if not value:
                raise LesionAtlasDescriptorError(
                    f"Shape Atlas has an empty File reference for {current}"
                )
            files[current] = value
            current = None
        else:
            raise LesionAtlasDescriptorError(f"Unexpected Shape Atlas field {key!r}")
    if current is not None:
        raise LesionAtlasDescriptorError(f"Shape Atlas has no File reference for {current}")
    missing = sorted(row.name for row in structures if row.name not in files)
    if missing:
        raise LesionAtlasDescriptorError(
            f"Shape Atlas is missing File mappings for: {', '.join(missing)}"
        )
    return files


def _parse_shape_dimensions(lines: list[str]) -> tuple[int, int, int]:
    """Read the global voxel dimensions declared by the Shape Atlas section."""
    header_index = _find_section(lines, "Shape Atlas")
    dimensions: tuple[int, int, int] | None = None
    for line in lines[header_index + 1 :]:
        stripped = line.strip()
        if not stripped:
            continue
        if _starts_section(stripped):
            break
        key, separator, value = stripped.partition(":")
        if not separator or key.strip() != "Dimensions":
            continue
        if dimensions is not None:
            raise LesionAtlasDescriptorError("Duplicate Shape Atlas Dimensions declaration")
        try:
            parsed = tuple(int(token) for token in value.split())
        except ValueError as exc:
            raise LesionAtlasDescriptorError(
                "Shape Atlas Dimensions must contain three positive integers"
            ) from exc
        if len(parsed) != 3 or any(size <= 0 for size in parsed):
            raise LesionAtlasDescriptorError(
                "Shape Atlas Dimensions must contain three positive integers"
            )
        dimensions = parsed
    if dimensions is None:
        raise LesionAtlasDescriptorError("Shape Atlas has no Dimensions declaration")
    return dimensions


def _parse_intensity_rows(
    lines: list[str], structure_count: int
) -> dict[str, np.ndarray]:
    header_index = _find_section(lines, "Intensity Atlas")
    rows: dict[str, np.ndarray] = {}
    for line in lines[header_index + 1 :]:
        stripped = line.strip()
        if not stripped:
            continue
        if _starts_section(stripped):
            break
        parts = stripped.split()
        if not parts or parts[0] not in _MODALITY_ALIASES:
            continue
        modality = _canonical_modality(parts[0])
        if modality in rows:
            raise LesionAtlasDescriptorError(f"Duplicate intensity row for {modality}")
        if len(parts) != structure_count + 1:
            raise LesionAtlasDescriptorError(
                f"Intensity row for {parts[0]} has {len(parts) - 1} values; "
                f"expected {structure_count}"
            )
        rows[modality] = _finite_float32(parts[1:], f"intensity row {parts[0]}")
    return rows


def _parse_modality_values(lines: list[str], heading: str) -> dict[str, np.float32]:
    header_index = _find_section(lines, heading)
    values: dict[str, np.float32] = {}
    for line in lines[header_index + 1 :]:
        stripped = line.strip()
        if not stripped:
            continue
        if _starts_section(stripped):
            break
        tokens = stripped.split()
        if len(tokens) % 2:
            raise LesionAtlasDescriptorError(
                f"Malformed {heading} modality/value row: {line!r}"
            )
        for offset in range(0, len(tokens), 2):
            name, value = tokens[offset : offset + 2]
            if name not in _MODALITY_ALIASES:
                continue
            modality = _canonical_modality(name)
            if modality in values:
                raise LesionAtlasDescriptorError(f"Duplicate {heading} value for {modality}")
            values[modality] = _finite_float32([value], f"{heading} value {name}")[0]
    return values


def _finite_float32(tokens: Sequence[str], name: str) -> np.ndarray:
    try:
        values = np.asarray([np.float32(token) for token in tokens], dtype=np.float32)
    except (TypeError, ValueError, OverflowError) as exc:
        raise LesionAtlasDescriptorError(f"Invalid numeric value in {name}") from exc
    if not np.isfinite(values).all():
        raise LesionAtlasDescriptorError(f"Non-finite numeric value in {name}")
    return values


def _find_section(lines: list[str], name: str) -> int:
    for index, line in enumerate(lines):
        if line.strip().startswith(f"{name}:"):
            return index
    raise LesionAtlasDescriptorError(f"Missing {name} section")


def _starts_section(line: str) -> bool:
    return any(line.startswith(f"{heading}:") for heading in _SECTION_HEADINGS)
