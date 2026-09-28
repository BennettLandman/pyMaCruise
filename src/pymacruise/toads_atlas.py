"""Read and preflight the selected source TOADS atlas descriptor and assets."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import nibabel as nib
import numpy as np


class ToadsAtlasError(ValueError):
    """Raised when a TOADS atlas descriptor or one of its images is invalid."""


@dataclass(frozen=True)
class ToadsStructure:
    name: str
    label: int
    topology: str


@dataclass(frozen=True)
class ToadsAtlasImage:
    role: str
    path: Path
    shape: tuple[int, int, int]
    voxel_sizes_mm: tuple[float, float, float]
    descriptor_voxel_sizes_mm: tuple[float, float, float]


@dataclass(frozen=True)
class ToadsAtlas:
    descriptor: Path
    structures: tuple[ToadsStructure, ...]
    topology: ToadsAtlasImage
    shapes: dict[str, ToadsAtlasImage]
    shape_samples: int
    shape_dimensions: tuple[int, int, int]
    shape_voxel_sizes_mm: tuple[float, float, float]
    orientations: tuple[int, int, int, int]
    geometry_warnings: tuple[str, ...]


def _tag_values(line: str, tag: str, expected: int) -> tuple[float, ...]:
    key, separator, value = line.partition(":")
    if not separator or key.strip() != tag:
        raise ToadsAtlasError(f"Expected {tag!r} metadata in atlas descriptor; got {line!r}")
    try:
        result = tuple(float(part) for part in value.split())
    except ValueError as exc:
        raise ToadsAtlasError(f"Invalid {tag} values in atlas descriptor: {line!r}") from exc
    if len(result) != expected or not np.isfinite(result).all():
        raise ToadsAtlasError(f"Expected {expected} finite {tag} values; got {line!r}")
    return result


def _metadata(lines: list[str], start: int, tag: str, count: int) -> tuple[float, ...]:
    for line in lines[start : start + 6]:
        if line.strip().startswith(f"{tag}:"):
            return _tag_values(line.strip(), tag, count)
    raise ToadsAtlasError(f"Missing {tag} metadata near descriptor line {start + 1}")


def _resolve_image(
    descriptor: Path,
    relative_path: str,
    role: str,
    descriptor_voxel_sizes_mm: tuple[float, ...],
) -> ToadsAtlasImage:
    relative = Path(relative_path)
    if relative.is_absolute() or ".." in relative.parts:
        raise ToadsAtlasError(f"Unsafe atlas image path for {role}: {relative_path}")
    if relative.suffixes[-2:] != [".nii", ".gz"]:
        raise ToadsAtlasError(f"TOADS atlas image must use .nii.gz: {relative_path}")
    path = (descriptor.parent / relative).resolve()
    if not path.is_file():
        raise ToadsAtlasError(f"Missing {role} atlas image: {path}")
    try:
        image = nib.load(str(path))
    except (OSError, nib.filebasedimages.ImageFileError) as exc:
        raise ToadsAtlasError(f"Cannot load {role} atlas image {path}: {exc}") from exc
    if len(image.shape) != 3:
        raise ToadsAtlasError(f"{role} atlas image must be 3D, got {image.shape}: {path}")
    if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
        raise ToadsAtlasError(f"{role} atlas image is not stored in RAS voxel orientation: {path}")
    zooms = tuple(float(value) for value in nib.affines.voxel_sizes(image.affine))
    return ToadsAtlasImage(
        role,
        path,
        tuple(int(value) for value in image.shape),
        zooms,
        tuple(float(value) for value in descriptor_voxel_sizes_mm),
    )


def load_toads_atlas(descriptor_path: str | Path) -> ToadsAtlas:
    """Resolve and validate atlas maps listed by a source-format TOADS file.

    The returned NIfTI images are inspected lazily; voxel payloads are not
    materialized. This performs asset preflight only and does not run TOADS.
    The descriptor's grid spacing is kept separately from each NIfTI affine:
    the selected source topology map declares 0.9375 mm while its source MIPAV
    XML and normalized NIfTI declare 1 mm. That conflict is surfaced as a
    geometry warning instead of silently changing either value.
    """
    descriptor = Path(descriptor_path).resolve()
    if not descriptor.is_file():
        raise ToadsAtlasError(f"Missing TOADS atlas descriptor: {descriptor}")
    lines = descriptor.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].strip() != "Structure Atlas File (edit at your own risks)":
        raise ToadsAtlasError(f"Not a supported source TOADS atlas descriptor: {descriptor}")

    structures_index = next(
        (index for index, line in enumerate(lines) if line.startswith("Structures:")), None
    )
    topology_index = next(
        (index for index, line in enumerate(lines) if line.startswith("Topology Atlas")), None
    )
    shape_index = next(
        (index for index, line in enumerate(lines) if line.startswith("Shape Atlas:")), None
    )
    if None in (structures_index, topology_index, shape_index):
        raise ToadsAtlasError("TOADS descriptor lacks structures, topology, or shape sections")

    count = int(_tag_values(lines[structures_index], "Structures", 1)[0])
    structures: list[ToadsStructure] = []
    for line in lines[structures_index + 1 : topology_index]:
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) != 3:
            raise ToadsAtlasError(f"Malformed TOADS structure row: {line!r}")
        try:
            structures.append(ToadsStructure(parts[0], int(parts[1]), parts[2]))
        except ValueError as exc:
            raise ToadsAtlasError(f"Malformed TOADS structure row: {line!r}") from exc
    if len(structures) != count:
        raise ToadsAtlasError(f"Descriptor declares {count} structures but lists {len(structures)}")
    if len({item.name for item in structures}) != count or len({item.label for item in structures}) != count:
        raise ToadsAtlasError("TOADS structure names and labels must be unique")

    topology_file = _metadata_file(lines, topology_index)
    topology_dimensions = tuple(
        int(value) for value in _metadata(lines, topology_index, "Dimensions", 3)
    )
    topology_resolution = _metadata(lines, topology_index, "Resolutions", 3)
    topology = _resolve_image(
        descriptor, topology_file, "topology", topology_resolution
    )
    geometry_warnings = list(
        _match_descriptor_image(topology, topology_dimensions, topology_resolution)
    )

    samples = int(_tag_values(lines[shape_index], "Shape Atlas", 1)[0])
    dimensions = tuple(int(value) for value in _metadata(lines, shape_index, "Dimensions", 3))
    resolutions = _metadata(lines, shape_index, "Resolutions", 3)
    orientations = tuple(
        int(value) for value in _metadata(lines, shape_index, "Orientations", 4)
    )
    shapes: dict[str, ToadsAtlasImage] = {}
    structure_names = {item.name for item in structures}
    for index, line in enumerate(lines[shape_index + 1 :]):
        if not line.startswith("Structure:"):
            continue
        name = line.partition(":")[2].strip()
        if name not in structure_names:
            raise ToadsAtlasError(f"Shape entry references unknown structure {name!r}")
        metadata_start = shape_index + 1 + index
        path = _metadata_file(lines, metadata_start)
        image_dims = tuple(
            int(value) for value in _metadata(lines, metadata_start, "Dimensions", 3)
        )
        image_resolution = _metadata(lines, metadata_start, "Resolutions", 3)
        file_image = _resolve_image(descriptor, path, f"shape {name}", image_resolution)
        geometry_warnings.extend(
            _match_descriptor_image(file_image, image_dims, image_resolution)
        )
        if name in shapes:
            raise ToadsAtlasError(f"Duplicate shape entry for {name!r}")
        shapes[name] = file_image
    if not shapes:
        raise ToadsAtlasError("TOADS descriptor contains no shape image entries")

    return ToadsAtlas(
        descriptor,
        tuple(structures),
        topology,
        shapes,
        samples,
        dimensions,
        tuple(float(value) for value in resolutions),
        orientations,
        tuple(geometry_warnings),
    )


def _metadata_file(lines: list[str], start: int) -> str:
    for line in lines[start : start + 6]:
        if line.strip().startswith("File:"):
            return line.partition(":")[2].strip()
    raise ToadsAtlasError(f"Missing atlas image file metadata near descriptor line {start + 1}")


def _match_descriptor_image(
    image: ToadsAtlasImage,
    dimensions: tuple[int, int, int],
    resolutions: tuple[float, ...],
) -> tuple[str, ...]:
    if image.shape != dimensions:
        raise ToadsAtlasError(
            f"{image.role} image shape {image.shape} disagrees with descriptor {dimensions}"
        )
    if np.allclose(image.voxel_sizes_mm, resolutions, rtol=0, atol=1e-4):
        return ()
    return (
        f"{image.role} NIfTI voxel sizes {image.voxel_sizes_mm} disagree with "
        f"TOADS descriptor resolutions {tuple(float(value) for value in resolutions)}; "
        "preserve both metadata sources until the historical loader behavior is verified",
    )
