"""Parse numeric priors from a source-format LesionTOADS atlas descriptor.

This parser reads only the descriptor's structures, intensity rows, lesion
priors, and black-hole priors. It does not inspect or load topology/shape image
arrays. Prior vectors retain the structure order declared in the descriptor.

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

import numpy as np


class LesionAtlasDescriptorError(ValueError):
    """Raised when a LesionTOADS text descriptor has invalid prior metadata."""


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


def _canonical_modality(name: str) -> str:
    try:
        return _MODALITY_ALIASES[name]
    except KeyError as exc:
        raise LesionAtlasDescriptorError(
            f"Unsupported LesionTOADS modality {name!r}"
        ) from exc


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
