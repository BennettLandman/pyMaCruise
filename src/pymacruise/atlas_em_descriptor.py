"""Metadata reader for the selected VandyCRUISE AtlasEM text descriptor.

This reads descriptor metadata and inventories paths only. The Java class loads
the referenced topology bytes and little-endian shape-float files; this module
does not load those payloads or claim AtlasEM execution parity.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


class AtlasEMDescriptorError(ValueError):
    """Raised when selected AtlasEM descriptor metadata is malformed."""


@dataclass(frozen=True)
class AtlasEMStructure:
    name: str
    label: int
    topology: str


@dataclass(frozen=True)
class AtlasEMDescriptor:
    path: Path
    structures: tuple[AtlasEMStructure, ...]
    topology_file: Path
    topology_dimensions_xyz: tuple[int, int, int]
    topology_resolutions_mm: tuple[float, float, float]
    shape_samples: int
    shape_dimensions_xyz: tuple[int, int, int]
    shape_resolutions_mm: tuple[float, float, float]
    orientations: tuple[int, int, int, int]
    shape_files: dict[str, Path]
    registered_shapes: tuple[bool, ...]
    referenced_files: tuple[Path, ...]
    missing_files: tuple[Path, ...]


def _fields(line: str, tag: str, count: int, line_number: int) -> list[str]:
    parts = line.split("\t")
    if len(parts) != count + 1 or parts[0].strip() != tag:
        raise AtlasEMDescriptorError(
            f"Line {line_number}: expected {tag} and {count} tab-separated values"
        )
    return parts[1:]


def _integer_fields(line: str, tag: str, count: int, line_number: int) -> tuple[int, ...]:
    try:
        values = tuple(int(value) for value in _fields(line, tag, count, line_number))
    except ValueError as exc:
        raise AtlasEMDescriptorError(f"Line {line_number}: invalid {tag} integer") from exc
    if any(value <= 0 for value in values):
        raise AtlasEMDescriptorError(f"Line {line_number}: {tag} values must be positive")
    return values


def _float_fields(line: str, tag: str, count: int, line_number: int) -> tuple[float, ...]:
    import math

    try:
        values = tuple(float(value) for value in _fields(line, tag, count, line_number))
    except ValueError as exc:
        raise AtlasEMDescriptorError(f"Line {line_number}: invalid {tag} number") from exc
    if not all(math.isfinite(value) and value > 0.0 for value in values):
        raise AtlasEMDescriptorError(f"Line {line_number}: {tag} values must be finite and positive")
    return values


def _asset_path(descriptor: Path, value: str, line_number: int) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise AtlasEMDescriptorError(f"Line {line_number}: unsafe atlas payload path {value!r}")
    resolved = (descriptor.parent / relative).resolve()
    if not resolved.is_relative_to(descriptor.parent):
        raise AtlasEMDescriptorError(f"Line {line_number}: atlas payload escapes descriptor directory")
    return resolved


def parse_atlas_em_descriptor(path: str | Path) -> AtlasEMDescriptor:
    """Parse source class order, spatial metadata, and referenced payload paths.

    Java's loader treats the topology dimensions as the dimensions for every
    shape payload. The per-structure metadata text following each ``File`` row
    is descriptive and is not used to determine the raw image allocation.
    Missing payloads are returned for inventory/reporting; use
    :func:`require_atlas_em_payloads` before attempting to load an atlas.
    """
    descriptor = Path(path).expanduser().resolve()
    try:
        lines = descriptor.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise AtlasEMDescriptorError(f"Cannot read AtlasEM descriptor {descriptor}: {exc}") from exc
    if not lines or lines[0] != "Structure Atlas File (edit at your own risks)":
        raise AtlasEMDescriptorError(f"Unsupported AtlasEM descriptor header: {descriptor}")

    structure_at = next((i for i, row in enumerate(lines) if row.startswith("Structures:")), None)
    topology_at = next((i for i, row in enumerate(lines) if row.startswith("Topology Atlas")), None)
    shape_at = next((i for i, row in enumerate(lines) if row.startswith("Shape Atlas:")), None)
    if structure_at is None or topology_at is None or shape_at is None:
        raise AtlasEMDescriptorError("Descriptor must contain Structures, Topology Atlas, and Shape Atlas sections")
    try:
        declared = int(_fields(lines[structure_at], "Structures:", 1, structure_at + 1)[0])
        shape_samples = int(_fields(lines[shape_at], "Shape Atlas:", 1, shape_at + 1)[0])
    except ValueError as exc:
        raise AtlasEMDescriptorError("Invalid structure count or shape sample count") from exc
    if declared <= 0 or shape_samples <= 0:
        raise AtlasEMDescriptorError("Structure and shape sample counts must be positive")

    structures = []
    for index in range(structure_at + 1, topology_at):
        row = lines[index].strip()
        if not row:
            continue
        fields = row.split("\t")
        if len(fields) != 3:
            raise AtlasEMDescriptorError(f"Line {index + 1}: malformed structure row")
        try:
            structures.append(AtlasEMStructure(fields[0], int(fields[1]), fields[2]))
        except ValueError as exc:
            raise AtlasEMDescriptorError(f"Line {index + 1}: malformed structure label") from exc
    if len(structures) != declared:
        raise AtlasEMDescriptorError(f"Descriptor declares {declared} structures but lists {len(structures)}")
    if len({s.name for s in structures}) != declared or len({s.label for s in structures}) != declared:
        raise AtlasEMDescriptorError("Structure names and labels must be unique")

    topology_fields = {}
    for tag, count in (("File:", 1), ("Dimensions:", 3), ("Resolutions:", 3)):
        line_at = next((i for i in range(topology_at + 1, min(topology_at + 5, len(lines)))
                        if lines[i].startswith(tag)), None)
        if line_at is None:
            raise AtlasEMDescriptorError(f"Topology Atlas section is missing {tag}")
        topology_fields[tag] = (line_at, _fields(lines[line_at], tag, count, line_at + 1))
    topology_at_file, topology_file_values = topology_fields["File:"]
    topology_file = _asset_path(descriptor, topology_file_values[0], topology_at_file + 1)
    dims_at, dims_values = topology_fields["Dimensions:"]
    topology_dims = _integer_fields(lines[dims_at], "Dimensions:", 3, dims_at + 1)
    res_at, _ = topology_fields["Resolutions:"]
    topology_res = _float_fields(lines[res_at], "Resolutions:", 3, res_at + 1)

    # The selected source reads four contiguous Shape Atlas metadata rows.
    shape_dims = _integer_fields(lines[shape_at + 1], "Dimensions:", 3, shape_at + 2)
    shape_res = _float_fields(lines[shape_at + 2], "Resolutions:", 3, shape_at + 3)
    orientations = _integer_fields(lines[shape_at + 3], "Orientations:", 4, shape_at + 4)

    shape_files: dict[str, Path] = {}
    referenced = [topology_file]
    names = {structure.name for structure in structures}
    for index, row in enumerate(lines[shape_at + 4 :], start=shape_at + 5):
        if not row.startswith("Structure:"):
            continue
        fields = row.split("\t")
        if len(fields) != 2:
            raise AtlasEMDescriptorError(f"Line {index}: malformed Shape Atlas structure entry")
        name = fields[1]
        if name not in names:
            raise AtlasEMDescriptorError(f"Line {index}: shape entry references unknown structure {name!r}")
        file_at = next((j for j in range(index, len(lines))
                        if lines[j].startswith("File:") or lines[j].startswith("Structure:")
                        or lines[j].startswith("Intensity Atlas:")
                        or lines[j].startswith("Registered Shapes")), None)
        if file_at is None or not lines[file_at].startswith("File:"):
            raise AtlasEMDescriptorError(f"Line {index}: shape {name!r} has no File row")
        file_values = _fields(lines[file_at], "File:", 1, file_at + 1)
        if name in shape_files:
            raise AtlasEMDescriptorError(f"Duplicate Shape Atlas entry for {name!r}")
        shape_files[name] = _asset_path(descriptor, file_values[0], file_at + 1)
        referenced.append(shape_files[name])

    if not shape_files:
        raise AtlasEMDescriptorError("Shape Atlas section contains no structure files")

    # Java defaults every class to registered except class 0. If the optional
    # section exists, it replaces the complete boolean vector from its value row.
    registered = [True] * declared
    registered[0] = False
    registered_at = next((i for i, row in enumerate(lines) if row.startswith("Registered Shapes")), None)
    if registered_at is not None:
        if registered_at + 2 >= len(lines):
            raise AtlasEMDescriptorError("Registered Shapes section lacks index and value rows")
        values = lines[registered_at + 2].split("\t")
        if len(values) != declared:
            raise AtlasEMDescriptorError("Registered Shapes row must contain one value per class")
        try:
            bits = [int(value) for value in values]
        except ValueError as exc:
            raise AtlasEMDescriptorError("Registered Shapes values must be 0 or 1") from exc
        # Source uses MipavUtil.getInt(value) == 1; every other integer is false.
        registered = [value == 1 for value in bits]

    unique_files = tuple(dict.fromkeys(referenced))
    missing = tuple(item for item in unique_files if not item.is_file())
    return AtlasEMDescriptor(
        descriptor, tuple(structures), topology_file, topology_dims, topology_res,
        shape_samples, shape_dims, shape_res, orientations, shape_files,
        tuple(registered), unique_files, missing,
    )


def require_atlas_em_payloads(descriptor: AtlasEMDescriptor) -> None:
    """Fail before execution when any listed topology/shape payload is absent."""
    if descriptor.missing_files:
        missing = ", ".join(str(item) for item in descriptor.missing_files)
        raise FileNotFoundError(f"AtlasEM descriptor payloads are missing: {missing}")
