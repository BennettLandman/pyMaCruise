"""Read and preflight the selected source TOADS atlas descriptor and assets."""

from __future__ import annotations

import json
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


@dataclass(frozen=True)
class ToadsAtlasData:
    """Atlas arrays restored to the original descriptor/Java index order."""

    topology: np.ndarray
    shapes: dict[str, np.ndarray]
    source_axis_codes: dict[str, tuple[str, str, str]]


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
    The descriptor's grid spacing is kept separately from each NIfTI affine.
    The selected Java loader reads descriptor resolutions and raw voxel
    payloads, not the image-header affine, so descriptor spacing controls the
    source algorithm. The topology NIfTI's 1 mm zoom (versus 0.9375 mm in the
    descriptor) is retained as a visible conversion/header warning; values
    are neither resampled nor relabeled as having the descriptor geometry.
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


def load_toads_intensity_priors(
    descriptor_path: str | Path,
    modalities: tuple[str, ...] | list[str],
) -> np.ndarray:
    """Return selected descriptor intensity priors as ``(channels, classes)``.

    The parser follows ``MultipleToadDeformableAtlas.modalityId`` and its
    ``Intensity Atlas`` table: aliases map to the same modality identifier,
    class values remain in descriptor structure order, and the caller's
    modality ordering determines the output channel order. The selected base
    workflow has explicit ``T1_MPRAGE`` priors. Source fallback behavior for
    missing modality profiles is not silently synthesized here; unsupported
    or absent profiles raise an actionable error. The source class is
    `toads-cruise/src/edu/jhmi/rad/medic/methods/MultipleToadDeformableAtlas.java`
    (June 2006; Pierre-Louis Bazin).
    """
    descriptor = Path(descriptor_path).resolve()
    if not descriptor.is_file():
        raise ToadsAtlasError(f"Missing TOADS atlas descriptor: {descriptor}")
    if not modalities:
        raise ToadsAtlasError("At least one image modality is required")
    try:
        lines = descriptor.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise ToadsAtlasError(f"Cannot read TOADS atlas descriptor {descriptor}: {exc}") from exc

    structures_index = next(
        (index for index, line in enumerate(lines) if line.startswith("Structures:")), None
    )
    intensity_index = next(
        (index for index, line in enumerate(lines) if line.startswith("Intensity Atlas:")), None
    )
    if structures_index is None or intensity_index is None:
        raise ToadsAtlasError("TOADS descriptor lacks Structures or Intensity Atlas sections")
    class_count = int(_tag_values(lines[structures_index], "Structures", 1)[0])

    profiles: dict[str, np.ndarray] = {}
    for line in lines[intensity_index + 1 :]:
        tokens = line.split()
        if not tokens:
            continue
        modality = _MODALITY_ALIASES.get(tokens[0])
        if modality is None:
            break
        if modality in profiles:
            raise ToadsAtlasError(f"Duplicate TOADS intensity profile for {modality}")
        if len(tokens) != class_count + 1:
            raise ToadsAtlasError(
                f"TOADS intensity profile {tokens[0]} has {len(tokens) - 1} values; "
                f"expected {class_count}"
            )
        try:
            values = np.asarray([float(value) for value in tokens[1:]], dtype=np.float32)
        except ValueError as exc:
            raise ToadsAtlasError(
                f"Invalid TOADS intensity prior in profile {tokens[0]}"
            ) from exc
        if not np.isfinite(values).all():
            raise ToadsAtlasError(f"Non-finite TOADS intensity prior in profile {tokens[0]}")
        profiles[modality] = values

    rows: list[np.ndarray] = []
    for requested in modalities:
        modality = _MODALITY_ALIASES.get(requested)
        if modality is None:
            raise ToadsAtlasError(f"Unsupported TOADS image modality: {requested}")
        try:
            rows.append(profiles[modality])
        except KeyError as exc:
            raise ToadsAtlasError(
                f"TOADS descriptor has no explicit intensity profile for {requested}"
            ) from exc
    return np.stack(rows, axis=0)


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
        "the selected Java loader uses descriptor spacing and reads raw voxel payloads "
        "without the image-header affine",
    )


def load_toads_atlas_data(
    atlas: ToadsAtlas,
    conversion_audit_path: str | Path | None = None,
) -> ToadsAtlasData:
    """Load atlas arrays in the historical Java index order.

    The checked-in NIfTIs are RAS-oriented, while the legacy Java code indexes
    headerless data in the original axis order. The conversion audit records
    the pre-RAS axis codes and is required to undo only that axis permutation
    and flipping. No interpolation, resampling, or physical-coordinate
    transformation is performed.
    """
    audit_path = (
        Path(conversion_audit_path)
        if conversion_audit_path is not None
        else atlas.descriptor.parent / "CONVERSION_AUDIT.json"
    )
    try:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ToadsAtlasError(f"Cannot read atlas conversion audit {audit_path}: {exc}") from exc
    rows = {
        str(row["nifti_output"]): row
        for row in audit.get("images", [])
        if "nifti_output" in row
    }
    loaded: dict[str, np.ndarray] = {}
    source_codes: dict[str, tuple[str, str, str]] = {}

    def read_image(item: ToadsAtlasImage) -> np.ndarray:
        try:
            relative = item.path.relative_to(atlas.descriptor.parent).as_posix()
        except ValueError as exc:
            raise ToadsAtlasError(
                f"Atlas image is outside descriptor directory: {item.path}"
            ) from exc
        row = rows.get(relative)
        if row is None:
            raise ToadsAtlasError(f"No conversion-audit row for atlas image {relative}")
        try:
            codes = tuple(str(code) for code in row["input_axis_codes"])
            if len(codes) != 3:
                raise ValueError("expected three axis codes")
            source_orientation = nib.orientations.axcodes2ornt(codes)
            image = nib.load(str(item.path))
            ras_orientation = nib.orientations.io_orientation(image.affine)
            transform = nib.orientations.ornt_transform(ras_orientation, source_orientation)
            array = np.asanyarray(image.dataobj)
            source_order = nib.orientations.apply_orientation(array, transform)
        except (KeyError, OSError, ValueError, nib.filebasedimages.ImageFileError) as exc:
            raise ToadsAtlasError(
                f"Cannot restore source index order for atlas image {relative}: {exc}"
            ) from exc
        expected_shape = tuple(int(value) for value in row["input_shape"])
        if source_order.shape != expected_shape or source_order.shape != item.shape:
            raise ToadsAtlasError(
                f"Restored {item.role} shape {source_order.shape} disagrees with "
                f"audit {expected_shape} or descriptor image {item.shape}"
            )
        source_codes[relative] = codes  # type: ignore[assignment]
        return np.array(source_order, copy=True)

    loaded["topology"] = read_image(atlas.topology)
    for name, image in atlas.shapes.items():
        loaded[name] = read_image(image)
    return ToadsAtlasData(loaded.pop("topology"), loaded, source_codes)
