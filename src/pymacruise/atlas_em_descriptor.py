"""Read selected AtlasEM descriptor metadata and compressed NIfTI payloads.

The source Java loader consumes raw topology bytes and little-endian shape
floats. PyMaCRUISE's payload path requires NIfTI and keeps its RAS storage
contract; this does not claim byte-identical raw-file loading or AtlasEM parity.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import normalize_image_to_ras


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
    intensity_samples: int
    intensity_atlas: dict[str, tuple[float, ...]]
    intensity_variance: dict[str, tuple[float, ...]]
    lesion_models: dict[str, float]
    black_hole_models: dict[str, float]
    modality_weights: dict[str, tuple[float, ...]]
    optimized_factors: dict[str, tuple[float, ...]]
    shape_slope: float
    referenced_files: tuple[Path, ...]
    missing_files: tuple[Path, ...]


@dataclass(frozen=True)
class AtlasEMNiftiPayloads:
    """Selected AtlasEM atlas voxels loaded from compressed NIfTI assets."""

    topology: np.ndarray
    shapes: dict[str, np.ndarray]
    source_axis_codes: dict[str, tuple[str, str, str]]


@dataclass(frozen=True)
class AtlasEMRegistrationAtlasInputs:
    """Class-ordered atlas arrays ready for the rigid-registration adapter."""

    topology: np.ndarray
    shape_priors: np.ndarray
    class_names: tuple[str, ...]
    class_labels: tuple[int, ...]
    registered_classes: np.ndarray


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


_MODALITY_ALIASES = {
    "T1": "T1_SPGR", "T1_SPGR": "T1_SPGR", "T2": "T2", "FLAIR": "FLAIR",
    "MPRAGE": "T1_MPRAGE", "T1_MPRAGE": "T1_MPRAGE",
    "MPRAGE_RAW": "T1_RAW", "T1_RAW": "T1_RAW", "PD": "PD",
    "PDFSE": "PDFSE", "DIR": "DIR",
}
_MODALITIES = tuple(dict.fromkeys(_MODALITY_ALIASES.values()))
_SECTION_PREFIXES = (
    "Structures:", "Topology Atlas", "Shape Atlas:", "Intensity Atlas:",
    "Intensity Variance", "Lesions", "Black Holes", "Slope", "Modality weights",
    "Registered Shapes", "Optimized factors",
)


def _profile_rows(lines: list[str], header_at: int, width: int) -> dict[str, tuple[float, ...]]:
    """Read modality rows as the selected loader's float profiles."""
    profiles = {}
    for index in range(header_at + 1, len(lines)):
        row = lines[index]
        if any(row.startswith(prefix) for prefix in _SECTION_PREFIXES):
            break
        fields = row.split("\t")
        if not fields or fields[0] not in _MODALITY_ALIASES:
            continue
        if len(fields) != width + 1:
            raise AtlasEMDescriptorError(
                f"Line {index + 1}: {fields[0]} profile must have {width} values"
            )
        try:
            values = tuple(float(value) for value in fields[1:])
        except ValueError as exc:
            raise AtlasEMDescriptorError(f"Line {index + 1}: invalid modality profile value") from exc
        import math

        if not all(math.isfinite(value) for value in values):
            raise AtlasEMDescriptorError(f"Line {index + 1}: modality profile values must be finite")
        # Java assigns by modality ID, so a later alias row replaces the earlier one.
        profiles[_MODALITY_ALIASES[fields[0]]] = values
    return profiles


def _scalar_models(lines: list[str], header_at: int) -> dict[str, float]:
    models = {}
    if header_at + 1 >= len(lines):
        raise AtlasEMDescriptorError(f"Line {header_at + 1}: missing modality scalar row")
    fields = lines[header_at + 1].split("\t")
    if len(fields) % 2:
        raise AtlasEMDescriptorError(f"Line {header_at + 2}: expected modality/value pairs")
    for offset in range(0, len(fields), 2):
        modality = _MODALITY_ALIASES.get(fields[offset])
        if modality is None:
            continue
        try:
            value = float(fields[offset + 1])
        except ValueError as exc:
            raise AtlasEMDescriptorError(f"Line {header_at + 2}: invalid modality scalar") from exc
        import math

        if not math.isfinite(value):
            raise AtlasEMDescriptorError(f"Line {header_at + 2}: modality scalar must be finite")
        models[modality] = value
    return models


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

    intensity_at = next((i for i, row in enumerate(lines) if row.startswith("Intensity Atlas:")), None)
    try:
        intensity_samples = (
            int(_fields(lines[intensity_at], "Intensity Atlas:", 1, intensity_at + 1)[0])
            if intensity_at is not None else 0
        )
    except ValueError as exc:
        raise AtlasEMDescriptorError("Invalid intensity sample count") from exc
    if intensity_samples < 0:
        raise AtlasEMDescriptorError("Intensity sample count must not be negative")

    intensity_atlas = {name: (0.0,) * declared for name in _MODALITIES}
    intensity_variance = {name: (1.0,) * declared for name in _MODALITIES}
    modality_weights = {name: (1.0, 1.0) for name in _MODALITIES}
    optimized_factors = {name: (1.0,) * 4 for name in _MODALITIES}
    lesion_models = {name: 0.0 for name in _MODALITIES}
    black_hole_models = {name: 0.0 for name in _MODALITIES}
    for index, row in enumerate(lines):
        if row.startswith("Intensity Atlas:"):
            intensity_atlas.update(_profile_rows(lines, index, declared))
        elif row.startswith("Intensity Variance"):
            intensity_variance.update(_profile_rows(lines, index, declared))
        elif row.startswith("Modality weights"):
            modality_weights.update(_profile_rows(lines, index, 2))
        elif row.startswith("Optimized factors"):
            optimized_factors.update(_profile_rows(lines, index, 4))
    lesions_at = next((i for i, row in enumerate(lines) if row.startswith("Lesions")), None)
    if lesions_at is not None:
        lesion_models.update(_scalar_models(lines, lesions_at))
    black_holes_at = next((i for i, row in enumerate(lines) if row.startswith("Black Holes")), None)
    if black_holes_at is not None:
        black_hole_models.update(_scalar_models(lines, black_holes_at))
    slope = 2.0
    slope_at = next((i for i, row in enumerate(lines) if row.startswith("Slope")), None)
    if slope_at is not None:
        fields = lines[slope_at].split("\t")
        if len(fields) != 2:
            raise AtlasEMDescriptorError(f"Line {slope_at + 1}: Slope requires one value")
        try:
            slope = float(fields[1])
        except ValueError as exc:
            raise AtlasEMDescriptorError(f"Line {slope_at + 1}: invalid Slope value") from exc
        import math

        if not math.isfinite(slope):
            raise AtlasEMDescriptorError(f"Line {slope_at + 1}: Slope must be finite")

    unique_files = tuple(dict.fromkeys(referenced))
    missing = tuple(item for item in unique_files if not item.is_file())
    return AtlasEMDescriptor(
        descriptor, tuple(structures), topology_file, topology_dims, topology_res,
        shape_samples, shape_dims, shape_res, orientations, shape_files,
        tuple(registered), intensity_samples, intensity_atlas, intensity_variance,
        lesion_models, black_hole_models, modality_weights, optimized_factors, slope,
        unique_files, missing,
    )


def require_atlas_em_payloads(descriptor: AtlasEMDescriptor) -> None:
    """Fail before execution when any listed topology/shape payload is absent."""
    if descriptor.missing_files:
        missing = ", ".join(str(item) for item in descriptor.missing_files)
        raise FileNotFoundError(f"AtlasEM descriptor payloads are missing: {missing}")


def load_atlas_em_nifti_payloads(descriptor: AtlasEMDescriptor) -> AtlasEMNiftiPayloads:
    """Load topology and per-structure maps from project-standard `.nii.gz` assets.

    Java consumes raw byte and little-endian float arrays. This Python path
    reads NIfTI headers, warns and reorients non-RAS maps, and returns arrays in
    RAS order. Descriptor dimensions/resolutions remain the source algorithm's
    registration geometry; header spacing discrepancies are reported instead
    of silently replacing descriptor values.
    """
    require_atlas_em_payloads(descriptor)

    def load_one(path: Path, role: str, expected_dimensions: tuple[int, int, int]):
        if path.suffixes[-2:] != [".nii", ".gz"]:
            raise AtlasEMDescriptorError(f"{role} payload must be stored as .nii.gz: {path}")
        try:
            source = nib.load(str(path))
        except (OSError, nib.filebasedimages.ImageFileError) as exc:
            raise AtlasEMDescriptorError(f"Cannot load {role} NIfTI {path}: {exc}") from exc
        if len(source.shape) != 3:
            raise AtlasEMDescriptorError(f"{role} payload must be a 3D scalar image: {path}")
        axis_codes = tuple(str(value) for value in nib.aff2axcodes(source.affine))
        ras = normalize_image_to_ras(source, source_name=f"AtlasEM {role}")
        if tuple(ras.shape) != expected_dimensions:
            raise AtlasEMDescriptorError(
                f"{role} payload RAS dimensions {ras.shape} disagree with descriptor {expected_dimensions}"
            )
        header_spacing = tuple(float(value) for value in nib.affines.voxel_sizes(ras.affine))
        reference_spacing = (
            descriptor.topology_resolutions_mm if role == "topology"
            else descriptor.shape_resolutions_mm
        )
        if not np.allclose(header_spacing, reference_spacing, rtol=0.0, atol=1e-4):
            warnings.warn(
                f"AtlasEM {role} NIfTI voxel sizes {header_spacing} mm differ from descriptor "
                f"resolutions {reference_spacing} mm; descriptor geometry is retained.",
                UserWarning,
                stacklevel=2,
            )
        values = np.asanyarray(ras.dataobj)
        if not np.isfinite(values).all():
            raise AtlasEMDescriptorError(f"{role} payload contains non-finite voxels: {path}")
        return values, axis_codes, np.asarray(ras.affine, dtype=np.float64)

    topology, topology_codes, reference_affine = load_one(
        descriptor.topology_file, "topology", descriptor.topology_dimensions_xyz
    )
    if not np.issubdtype(topology.dtype, np.integer):
        if not np.equal(topology, np.floor(topology)).all():
            raise AtlasEMDescriptorError("topology NIfTI must contain integer labels")
        topology = topology.astype(np.int32)
    shapes = {}
    axis_codes = {"topology": topology_codes}
    for name, path in descriptor.shape_files.items():
        values, codes, affine = load_one(path, f"shape {name}", descriptor.topology_dimensions_xyz)
        if not np.allclose(affine, reference_affine, rtol=0.0, atol=1e-4):
            raise AtlasEMDescriptorError(
                f"shape {name} NIfTI grid does not match the topology NIfTI grid"
            )
        shapes[name] = np.asarray(values, dtype=np.float32)
        axis_codes[f"shape:{name}"] = codes
    return AtlasEMNiftiPayloads(topology, shapes, axis_codes)


def load_atlas_em_registration_inputs(
    descriptor: AtlasEMDescriptor,
) -> AtlasEMRegistrationAtlasInputs:
    """Load source-ordered RAS atlas arrays for the rigid optimizer interface.

    The shared-transform `registerShapesPyramid` source indexes every shape
    array by descriptor class ID. Require exactly one payload for each class
    rather than silently shifting or fabricating missing maps.
    """
    payloads = load_atlas_em_nifti_payloads(descriptor)
    names = tuple(item.name for item in descriptor.structures)
    missing = tuple(name for name in names if name not in payloads.shapes)
    unexpected = tuple(name for name in payloads.shapes if name not in names)
    if missing or unexpected:
        raise AtlasEMDescriptorError(
            f"registration requires one shape prior per descriptor class; missing={missing}, "
            f"unexpected={unexpected}"
        )
    ordered = np.stack([payloads.shapes[name] for name in names], axis=-1)
    return AtlasEMRegistrationAtlasInputs(
        payloads.topology,
        ordered,
        names,
        tuple(item.label for item in descriptor.structures),
        np.asarray(descriptor.registered_shapes, dtype=bool),
    )
