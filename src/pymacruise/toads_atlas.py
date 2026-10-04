"""Read and preflight the selected source TOADS atlas descriptor and assets.

The statistical shape-prior and topology-template roles of these assets follow
the atlas framework in Bazin & Pham (2008). Parsing, file ordering, and stored
metadata are governed by the selected VandyCRUISE descriptor/source contract;
loading the assets does not itself perform TOADS segmentation. See
``docs/papers/toads.md``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import normalize_image_to_ras


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


@dataclass(frozen=True)
class ToadsAtlasRASData:
    """Atlas voxel arrays in canonical RAS axis order for Python algorithms."""

    topology: np.ndarray
    shapes: dict[str, np.ndarray]
    source_axis_codes: dict[str, tuple[str, str, str]]
    template_resolution_mm: tuple[float, float, float]


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
    # Preflight records the RAS array dimensions/spacing that the Python path
    # will use. The data-loading API performs the user-visible warning and
    # returns the actual reoriented voxels; this metadata-only operation keeps
    # atlas resolution separate from algorithm execution.
    image = nib.as_closest_canonical(image)
    if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
        raise ToadsAtlasError(f"{role} atlas image cannot be represented in RAS voxel orientation: {path}")
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


def load_toads_intensity_variance_factors(
    descriptor_path: str | Path,
    modalities: tuple[str, ...] | list[str],
) -> np.ndarray:
    """Return source TOADS intensity-variance factors as ``(channels, classes)``.

    This reads the distinct ``Intensity Variance Factor`` descriptor section
    used by ``MultipleToadDeformableAtlas.getIntensityVariancePriors``. Every
    Java variance row is initialized to float32 ones before parsing; recognized
    modality profiles replace those defaults, and unknown requested modality
    names retain the all-ones fallback. Aliases and requested row order follow
    :meth:`MultipleToadDeformableAtlas.modalityId`.
    """
    descriptor = Path(descriptor_path).resolve()
    if not descriptor.is_file():
        raise ToadsAtlasError(f"Missing TOADS atlas descriptor: {descriptor}")
    try:
        lines = descriptor.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise ToadsAtlasError(f"Cannot read TOADS atlas descriptor {descriptor}: {exc}") from exc

    structures_index = next(
        (index for index, line in enumerate(lines) if line.startswith("Structures:")), None
    )
    if structures_index is None:
        raise ToadsAtlasError("TOADS descriptor lacks Structures section")
    raw_class_count = _tag_values(lines[structures_index], "Structures", 1)[0]
    class_count = int(raw_class_count)
    if class_count < 1 or raw_class_count != class_count:
        raise ToadsAtlasError(f"Invalid TOADS structure count: {class_count}")

    # Java allocates [INTENSITY][classes] and explicitly fills it with 1.0f
    # while loading Structures (MultipleToadDeformableAtlas.java:939-942).
    modality_ids = tuple(dict.fromkeys(_MODALITY_ALIASES.values()))
    defaults = np.ones((len(modality_ids), class_count), dtype=np.float32)
    profiles = {
        modality: defaults[index].copy() for index, modality in enumerate(modality_ids)
    }
    variance_index = next(
        (index for index, line in enumerate(lines) if line.startswith("Intensity Variance Factor:")),
        None,
    )
    if variance_index is not None:
        for line in lines[variance_index + 1 :]:
            fields = line.split()
            if not fields:
                continue
            modality = _MODALITY_ALIASES.get(fields[0])
            # The Java parser stops at the first non-modality section/row.
            if modality is None:
                break
            if len(fields) != class_count + 1:
                raise ToadsAtlasError(
                    f"TOADS intensity variance profile {fields[0]} has {len(fields) - 1} values; "
                    f"expected {class_count}"
                )
            try:
                values = np.asarray(fields[1:], dtype=np.float32)
            except ValueError as exc:
                raise ToadsAtlasError(
                    f"Invalid TOADS intensity variance factor in profile {fields[0]}"
                ) from exc
            if not np.isfinite(values).all():
                raise ToadsAtlasError(
                    f"Non-finite TOADS intensity variance factor in profile {fields[0]}"
                )
            profiles[modality] = values

    result = np.ones((len(modalities), class_count), dtype=np.float32)
    for channel, requested in enumerate(modalities):
        modality = _MODALITY_ALIASES.get(requested)
        if modality is not None:
            result[channel] = profiles[modality]
    return result


def load_toads_optimized_factors(
    descriptor_path: str | Path,
    modalities: tuple[str, ...] | list[str],
) -> np.ndarray:
    """Return source optimized CSF/GM/WM exponents for each modality.

    ``MultipleToadDeformableAtlas`` initializes all three exponents to 1.0,
    then reads the optional ``Optimized factors`` descriptor section. The
    selected ``AlgorithmMultipleToads`` wrapper passes those values to its
    CSF/GM/WM ``PowerTable`` instances. Missing rows preserve the Java ones.
    """
    descriptor = Path(descriptor_path).resolve()
    if not descriptor.is_file():
        raise ToadsAtlasError(f"Missing TOADS atlas descriptor: {descriptor}")
    try:
        lines = descriptor.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise ToadsAtlasError(f"Cannot read TOADS atlas descriptor {descriptor}: {exc}") from exc

    modality_ids = tuple(dict.fromkeys(_MODALITY_ALIASES.values()))
    profiles = {modality: np.ones(3, dtype=np.float32) for modality in modality_ids}
    section = next(
        (index for index, line in enumerate(lines) if line.startswith("Optimized factors")),
        None,
    )
    if section is not None:
        for line in lines[section + 1 :]:
            fields = line.split()
            if not fields:
                continue
            modality = _MODALITY_ALIASES.get(fields[0])
            if modality is None:
                break
            if len(fields) != 4:
                raise ToadsAtlasError(
                    f"TOADS optimized factor row {fields[0]} has {len(fields) - 1} values; expected 3"
                )
            try:
                values = np.asarray(fields[1:], dtype=np.float32)
            except ValueError as exc:
                raise ToadsAtlasError(
                    f"Invalid TOADS optimized factors in profile {fields[0]}"
                ) from exc
            if not np.isfinite(values).all() or np.any(values <= 0):
                raise ToadsAtlasError(
                    f"TOADS optimized factors in profile {fields[0]} must be finite and positive"
                )
            profiles[modality] = values

    rows: list[np.ndarray] = []
    for requested in modalities:
        modality = _MODALITY_ALIASES.get(requested)
        if modality is None:
            raise ToadsAtlasError(f"Unsupported TOADS image modality: {requested}")
        rows.append(profiles[modality])
    return np.stack(rows, axis=0)


def load_toads_lesion_priors(
    descriptor_path: str | Path,
    modalities: tuple[str, ...] | list[str],
) -> np.ndarray:
    """Return the selected atlas lesion prior per requested image channel.

    This follows ``MultipleToadAtlas.getLesionPriors``: absent lesion metadata
    yields zero for every channel, and a listed but unsupported modality also
    yields zero. Values are in caller modality order.
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

    prior_by_modality: dict[str, np.float32] = {}
    lesion_index = next((i for i, line in enumerate(lines) if line.startswith("Lesions")), None)
    if lesion_index is not None:
        if lesion_index + 1 >= len(lines):
            raise ToadsAtlasError("TOADS Lesions section is missing its modality values")
        tokens = lines[lesion_index + 1].split()
        if len(tokens) % 2:
            raise ToadsAtlasError("TOADS Lesions values must contain modality/value pairs")
        for modality, raw_value in zip(tokens[::2], tokens[1::2], strict=True):
            canonical = _MODALITY_ALIASES.get(modality)
            try:
                value = np.float32(float(raw_value))
            except ValueError as exc:
                raise ToadsAtlasError(f"Invalid lesion prior for modality {modality}") from exc
            if not np.isfinite(value):
                raise ToadsAtlasError(f"Non-finite lesion prior for modality {modality}")
            if canonical is not None:
                prior_by_modality[canonical] = value

    result: list[np.float32] = []
    for requested in modalities:
        canonical = _MODALITY_ALIASES.get(requested)
        if canonical is None:
            raise ToadsAtlasError(f"Unsupported TOADS image modality: {requested}")
        result.append(prior_by_modality.get(canonical, np.float32(0.0)))
    return np.asarray(result, dtype=np.float32)


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
    if conversion_audit_path is not None:
        audit_path = Path(conversion_audit_path)
    else:
        candidates = (atlas.descriptor.parent, *atlas.descriptor.parent.parents)
        audit_path = next(
            (
                parent / "CONVERSION_AUDIT.json"
                for parent in candidates
                if (parent / "CONVERSION_AUDIT.json").is_file()
            ),
            atlas.descriptor.parent / "CONVERSION_AUDIT.json",
        )
    audit_available = audit_path.exists()
    audit_root = audit_path.parent if audit_available else atlas.descriptor.parent
    if audit_available:
        try:
            audit = json.loads(audit_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ToadsAtlasError(f"Cannot read atlas conversion audit {audit_path}: {exc}") from exc
    else:
        # Standalone synthetic or user-created NIfTI atlases have no prior
        # raw payload to restore; their native array order is their source order.
        audit = {}
    rows = {
        str(row["nifti_output"]): row
        for row in audit.get("images", [])
        if "nifti_output" in row
    }
    loaded: dict[str, np.ndarray] = {}
    source_codes: dict[str, tuple[str, str, str]] = {}

    def read_image(item: ToadsAtlasImage) -> np.ndarray:
        try:
            relative = item.path.relative_to(audit_root).as_posix()
        except ValueError as exc:
            raise ToadsAtlasError(
                f"Atlas image is outside descriptor directory: {item.path}"
            ) from exc
        row = rows.get(relative)
        if row is None and audit_available:
            raise ToadsAtlasError(f"No conversion-audit row for atlas image {relative}")
        try:
            image = nib.load(str(item.path))
            codes = (
                tuple(str(code) for code in row["input_axis_codes"])
                if row is not None
                else tuple(str(code) for code in nib.aff2axcodes(image.affine))
            )
            if len(codes) != 3:
                raise ValueError("expected three axis codes")
            source_orientation = nib.orientations.axcodes2ornt(codes)
            if tuple(int(value) for value in nib.as_closest_canonical(image).shape) != item.shape:
                raise ValueError("canonical image shape disagrees with atlas metadata")
            ras_orientation = nib.orientations.io_orientation(image.affine)
            transform = nib.orientations.ornt_transform(ras_orientation, source_orientation)
            array = np.asanyarray(image.dataobj)
            source_order = nib.orientations.apply_orientation(array, transform)
        except (KeyError, OSError, ValueError, nib.filebasedimages.ImageFileError) as exc:
            raise ToadsAtlasError(
                f"Cannot restore source index order for atlas image {relative}: {exc}"
            ) from exc
        expected_shape = (
            tuple(int(value) for value in row["input_shape"])
            if row is not None
            else tuple(int(value) for value in image.shape)
        )
        if source_order.shape != expected_shape:
            raise ToadsAtlasError(
                f"Restored {item.role} shape {source_order.shape} disagrees with "
                f"conversion audit {expected_shape}"
            )
        source_codes[relative] = codes  # type: ignore[assignment]
        return np.array(source_order, copy=True)

    loaded["topology"] = read_image(atlas.topology)
    for name, image in atlas.shapes.items():
        loaded[name] = read_image(image)
    return ToadsAtlasData(loaded.pop("topology"), loaded, source_codes)


def load_toads_atlas_ras_data(atlas: ToadsAtlas) -> ToadsAtlasRASData:
    """Load topology and shape maps in RAS array order for Python processing.

    Unlike :func:`load_toads_atlas_data`, this returns the compressed NIfTI
    voxel order after axis-only RAS normalization. It is intended for the
    Python registration/segmentation APIs, whose image arrays are RAS ordered;
    it must not be substituted into Java-source-index comparison harnesses.
    Descriptor spacing is returned separately because the selected atlas
    loader uses those values, even when NIfTI header zooms differ.
    """
    source_codes: dict[str, tuple[str, str, str]] = {}

    def read_ras(item: ToadsAtlasImage) -> np.ndarray:
        try:
            relative = item.path.relative_to(atlas.descriptor.parent).as_posix()
            image = nib.load(str(item.path))
        except (ValueError, OSError, nib.filebasedimages.ImageFileError) as exc:
            raise ToadsAtlasError(f"Cannot load atlas image {item.path}: {exc}") from exc
        source_orientation = tuple(nib.aff2axcodes(image.affine))
        if len(image.shape) != 3:
            raise ToadsAtlasError(f"Atlas image must be 3D: {item.path}")
        try:
            normalized = normalize_image_to_ras(image, source_name=f"TOADS atlas {item.role}")
        except ValueError as exc:
            raise ToadsAtlasError(f"Cannot normalize atlas image {item.path}: {exc}") from exc
        if tuple(nib.aff2axcodes(normalized.affine)) != ("R", "A", "S"):
            raise ToadsAtlasError(f"Atlas image could not be normalized to RAS: {item.path}")
        values = np.asanyarray(normalized.dataobj)
        if values.ndim != 3 or not np.isfinite(values).all():
            raise ToadsAtlasError(f"Atlas image contains invalid voxel data: {item.path}")
        if tuple(int(v) for v in values.shape) != item.shape:
            raise ToadsAtlasError(
                f"RAS atlas shape {values.shape} disagrees with preflight shape {item.shape}: {item.path}"
            )
        source_codes[relative] = source_orientation  # type: ignore[assignment]
        return np.array(values, copy=True)

    topology = read_ras(atlas.topology)
    shapes = {name: read_ras(image) for name, image in atlas.shapes.items()}
    for name, values in shapes.items():
        if values.shape != topology.shape:
            raise ToadsAtlasError(
                f"Shape map {name!r} has grid {values.shape}, expected topology grid {topology.shape}"
            )
    if not np.issubdtype(topology.dtype, np.integer):
        raise ToadsAtlasError("TOADS topology image must contain integer labels")
    return ToadsAtlasRASData(
        topology,
        shapes,
        source_codes,
        atlas.shape_voxel_sizes_mm,
    )


def toads_atlas_shapes_in_structure_order(
    atlas: ToadsAtlas,
    atlas_data: ToadsAtlasData | ToadsAtlasRASData,
) -> np.ndarray:
    """Stack shape maps as ``(x,y,z,classes)`` in descriptor order.

    The returned voxel order is the order used by the supplied atlas data.
    The source TOADS alignment path uses the audited Java index order; callers
    may also explicitly pass canonical RAS data for independent RAS workflows.
    """
    if not isinstance(atlas_data, (ToadsAtlasData, ToadsAtlasRASData)):
        raise TypeError("atlas_data must be ToadsAtlasData or ToadsAtlasRASData")
    names = tuple(structure.name for structure in atlas.structures)
    missing = [name for name in names if name not in atlas_data.shapes]
    if missing:
        raise ToadsAtlasError(f"RAS atlas data is missing structure shape(s): {missing}")
    if not names:
        raise ToadsAtlasError("Atlas descriptor has no structures")
    return np.stack(
        tuple(np.asarray(atlas_data.shapes[name], dtype=np.float32) for name in names),
        axis=-1,
    )
