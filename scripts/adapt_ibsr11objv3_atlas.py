#!/usr/bin/env python3
"""Convert the selected MIPAV raw AtlasEM package to Python RAS NIfTI inputs.

The output is a Python-loader descriptor and a set of derived NIfTI payloads;
it does not modify or replace the source atlas. The source descriptor's axis
codes and spacings define the array reorientation and output geometry.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import nibabel as nib
import numpy as np
from nibabel.orientations import apply_orientation, axcodes2ornt, ornt_transform


AXIS_CODES = {
    1: (0, -1),  # Right to Left
    2: (0, 1),   # Left to Right
    3: (1, 1),   # Posterior to Anterior
    4: (1, -1),  # Anterior to Posterior
    5: (2, 1),   # Inferior to Superior
    6: (2, -1),  # Superior to Inferior
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _field(line: str, name: str) -> str:
    key, separator, value = line.partition(":")
    if not separator or key.strip() != name:
        raise ValueError(f"Expected {name}: row, got {line!r}")
    return value.strip()


def _triple(value: str, cast, name: str):
    try:
        result = tuple(cast(item) for item in value.split())
    except ValueError as exc:
        raise ValueError(f"Invalid {name}: {value!r}") from exc
    if len(result) != 3:
        raise ValueError(f"{name} must have three values: {value!r}")
    return result


def _source_assets(descriptor: Path):
    lines = descriptor.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].strip() != "Structure Atlas File (edit at your own risks)":
        raise ValueError(f"Not a LesionTOADS structure atlas descriptor: {descriptor}")
    struct_index = next((i for i, line in enumerate(lines) if line.startswith("Structures:")), None)
    if struct_index is None:
        raise ValueError("Missing Structures section")
    count = int(_field(lines[struct_index], "Structures"))
    structures = []
    for line in lines[struct_index + 1:struct_index + 1 + count]:
        fields = line.split()
        if len(fields) != 3:
            raise ValueError(f"Malformed structure row: {line!r}")
        structures.append((fields[0], int(fields[1]), fields[2]))
    if len(structures) != count or len({row[0] for row in structures}) != count:
        raise ValueError("Structure rows do not match the declared unique structure count")

    top_index = next((i for i, line in enumerate(lines) if line.startswith("Topology Atlas:")), None)
    shape_index = next((i for i, line in enumerate(lines) if line.startswith("Shape Atlas:")), None)
    intensity_index = next((i for i, line in enumerate(lines) if line.startswith("Intensity Atlas:")), None)
    if None in (top_index, shape_index, intensity_index) or not top_index < shape_index < intensity_index:
        raise ValueError("Expected Topology Atlas, Shape Atlas, and Intensity Atlas sections in order")

    def section_fields(start: int, end: int) -> dict[str, str]:
        fields = {}
        for line in lines[start + 1:end]:
            if not line.strip():
                continue
            key, separator, value = line.partition(":")
            if separator:
                fields[key.strip()] = value.strip()
        return fields

    topology_fields = section_fields(top_index, shape_index)
    top_ref = topology_fields["File"]
    dimensions = _triple(topology_fields["Dimensions"], int, "topology dimensions")
    topology_spacing = _triple(topology_fields["Resolutions"], float, "topology resolutions")

    shape_header = _field(lines[shape_index], "Shape Atlas")
    shape_fields = section_fields(shape_index, intensity_index)
    shape_dimensions = _triple(shape_fields["Dimensions"], int, "shape dimensions")
    spacing = _triple(shape_fields["Resolutions"], float, "shape resolutions")
    if dimensions != shape_dimensions or topology_spacing != spacing:
        raise ValueError("Topology and shape atlas dimensions/resolutions must match")
    orientation_line = next((line for line in lines[shape_index + 1:intensity_index]
                             if line.strip().startswith("Orientations:")), None)
    if orientation_line is None:
        raise ValueError("Shape Atlas section is missing Orientations")
    try:
        orientation_values = tuple(int(item) for item in _field(orientation_line, "Orientations").split())
    except ValueError as exc:
        raise ValueError("Invalid MIPAV orientation codes") from exc
    if len(orientation_values) != 4:
        raise ValueError(f"MIPAV Orientations must have plane and three axis codes: {orientation_values}")
    axis_codes = orientation_values[1:]
    orientation = [AXIS_CODES.get(code) for code in axis_codes]
    if any(code is None for code in orientation) or len(set(code[0] for code in orientation) if all(orientation) else ()) != 3:
        raise ValueError(f"Unsupported or non-orthogonal MIPAV orientation: {orientation_values}")
    source_ornt = np.asarray([[axis, direction] for axis, direction in orientation], dtype=float)
    transform = ornt_transform(source_ornt, axcodes2ornt(("R", "A", "S")))

    shape_records: dict[str, tuple[str, tuple[int, int, int], tuple[float, float, float]]] = {}
    current: str | None = None
    i = shape_index + 1
    while i < intensity_index:
        line = lines[i].strip()
        if line.startswith("Structure:"):
            current = _field(line, "Structure")
            if current not in {row[0] for row in structures} or current in shape_records:
                raise ValueError(f"Invalid or duplicate shape structure: {current!r}")
            if i + 4 >= intensity_index:
                raise ValueError(f"Incomplete shape record for {current}")
            file_ref = _field(lines[i + 1].strip(), "File")
            file_dims = _triple(_field(lines[i + 2].strip(), "Dimensions"), int, "shape map dimensions")
            file_spacing = _triple(_field(lines[i + 3].strip(), "Resolutions"), float, "shape map resolutions")
            if file_dims != dimensions or file_spacing != spacing:
                raise ValueError(f"Per-map grid differs from global grid for {current}")
            if not lines[i + 4].strip().startswith("Center:"):
                raise ValueError(f"Missing Center row for {current}")
            shape_records[current] = (file_ref, file_dims, file_spacing)
            i += 5
        else:
            i += 1
    if set(shape_records) != {row[0] for row in structures}:
        raise ValueError("Shape maps do not cover every declared structure in source order")

    def resolve(ref: str) -> Path:
        relative = Path(ref)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe source asset path: {ref}")
        candidates = [descriptor.parent / relative]
        if relative.parts and relative.parts[0] == descriptor.parent.name:
            candidates.append(descriptor.parent.parent / relative)
        for candidate in candidates:
            if candidate.is_file() and not candidate.is_symlink():
                return candidate.resolve()
        raise FileNotFoundError(f"Cannot resolve atlas payload {ref!r} relative to {descriptor.parent}")

    return {
        "lines": lines,
        "structures": structures,
        "topology": resolve(top_ref),
        "shape_records": {name: resolve(record[0]) for name, record in shape_records.items()},
        "dimensions": dimensions,
        "spacing": spacing,
        "orientation_values": orientation_values,
        "transform": transform,
        "shape_header": shape_header,
        "intensity_index": intensity_index,
        "structure_index": struct_index,
    }


def _read_raw(path: Path, dimensions: tuple[int, int, int], dtype: str) -> np.ndarray:
    expected = int(np.prod(dimensions))
    item_size = np.dtype(dtype).itemsize
    if path.stat().st_size != expected * item_size:
        raise ValueError(f"Raw payload size mismatch for {path}: expected {expected * item_size}, got {path.stat().st_size}")
    # MaCRUISE Java reads x-fastest: buffer[x + Nx*y + Nx*Ny*z]. Convert that
    # byte order to NumPy's explicit (x,y,z) coordinates before orientation.
    raw = np.fromfile(path, dtype=dtype)
    return raw.reshape((dimensions[2], dimensions[1], dimensions[0])).transpose(2, 1, 0)


def _nifti_ras(data_xyz: np.ndarray, transform: np.ndarray, spacing_xyz: tuple[float, float, float], path: Path) -> dict:
    ras = np.asarray(apply_orientation(data_xyz, transform))
    ras_spacing = tuple(float(spacing_xyz[int(axis)]) for axis in transform[:, 0])
    affine = np.eye(4, dtype=np.float64)
    affine[:3, :3] = np.diag(ras_spacing)
    # Source AtlasEM coordinates use dimensions/2 as the atlas center.
    affine[:3, 3] = -np.asarray(ras_spacing) * np.asarray(ras.shape, dtype=np.float64) / 2.0
    image = nib.Nifti1Image(ras, affine)
    image.header.set_xyzt_units("mm")
    nib.save(image, str(path))
    check = nib.load(str(path))
    if check.shape != ras.shape or not np.array_equal(np.asarray(check.dataobj), ras):
        raise RuntimeError(f"NIfTI round-trip changed atlas values or dimensions: {path}")
    return {"shape_ras": list(ras.shape), "spacing_ras_mm": list(ras_spacing), "sha256": sha256(path)}


def adapt(descriptor_path: Path, output_dir: Path) -> Path:
    descriptor = descriptor_path.expanduser().resolve()
    if not descriptor.is_file() or descriptor.is_symlink():
        raise FileNotFoundError(f"Missing regular source descriptor: {descriptor}")
    output_dir = output_dir.expanduser().resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Output directory must be empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    source = _source_assets(descriptor)
    dimensions = source["dimensions"]
    spacing = source["spacing"]
    transform = source["transform"]
    ras_shape = tuple(np.asarray(apply_orientation(np.empty(dimensions, dtype=np.uint8), transform)).shape)
    ras_spacing = tuple(float(spacing[int(axis)]) for axis in transform[:, 0])
    outputs: dict[str, dict] = {}

    topology_name = "topology.nii.gz"
    top_array = _read_raw(source["topology"], dimensions, "i1")
    outputs[topology_name] = _nifti_ras(top_array, transform, spacing, output_dir / topology_name)
    shape_names: dict[str, str] = {}
    for index, (name, _, _) in enumerate(source["structures"]):
        filename = f"shape-{index:02d}-{re.sub(r'[^A-Za-z0-9_.-]+', '-', name)}.nii.gz"
        shape_array = _read_raw(source["shape_records"][name], dimensions, "<f4")
        if not np.isfinite(shape_array).all():
            raise ValueError(f"Non-finite shape-prior values in {source['shape_records'][name]}")
        outputs[filename] = _nifti_ras(shape_array, transform, spacing, output_dir / filename)
        shape_names[name] = filename

    lines = source["lines"]
    struct_index = source["structure_index"]
    adapted_lines = [lines[0].strip(), f"Structures:\t{len(source['structures'])}"]
    adapted_lines.extend("\t".join((name, str(label), topology))
                         for name, label, topology in source["structures"])
    adapted_lines += [
        "", "Topology Atlas:", f"File:\t{topology_name}",
        "Dimensions:\t" + "\t".join(map(str, ras_shape)),
        "Resolutions:\t" + "\t".join(f"{value:.9g}" for value in ras_spacing),
        "", f"Shape Atlas:\t{source['shape_header']}",
        "Dimensions:\t" + "\t".join(map(str, ras_shape)),
        "Resolutions:\t" + "\t".join(f"{value:.9g}" for value in ras_spacing),
        # The AtlasEM adapter stores canonical RAS voxel arrays. These legacy
        # orientation fields are parsed but not used by Python registration.
        "Orientations:\t1\t1\t1\t1",
    ]
    for name, _, _ in source["structures"]:
        adapted_lines.extend((f"Structure:\t{name}", f"File:\t{shape_names[name]}"))
    # Preserve intensity means, lesion priors, and subsequent descriptor
    # options verbatim while removing source-grid-only metadata from map rows.
    intensity_index = source["intensity_index"]
    tail = lines[intensity_index:]
    # The Python parser expects Lesions and optional Black Holes; this source
    # has Lesions but no Black Holes, so use explicit zero defaults as the Java
    # loader does for an absent black-hole row.
    if not any(line.strip().startswith("Black Holes:") for line in tail):
        tail += ["", "Black Holes:", "T1\t0\tT2\t0\tFLAIR\t0\tMPRAGE\t0\tPD\t0"]
    adapted_lines.extend(["", *[line.rstrip() for line in tail]])
    adapted_descriptor = output_dir / "LesionToads-11-obj-atlasV3-atlasem-python.txt"
    adapted_descriptor.write_text("\n".join(adapted_lines) + "\n", encoding="utf-8")
    manifest = {
        "schema_version": 1,
        "source_descriptor": str(descriptor),
        "source_descriptor_sha256": sha256(descriptor),
        "source_dimensions_xyz": list(dimensions),
        "source_spacing_mm": list(spacing),
        "source_orientation_codes": list(source["orientation_values"]),
        "raw_indexing": "selected DemonToadDeformableAtlas x-fastest buffer[x+Nx*y+Nx*Ny*z]",
        "orientation_transform_to_ras": transform.tolist(),
        "ras_dimensions_xyz": list(ras_shape),
        "ras_spacing_mm": list(ras_spacing),
        "structure_order": [{"name": name, "label": label, "topology": topology}
                            for name, label, topology in source["structures"]],
        "source_payloads": {
            "topology": {"path": str(source["topology"]), "sha256": sha256(source["topology"])},
            **{name: {"path": str(path), "sha256": sha256(path)}
               for name, path in source["shape_records"].items()},
        },
        "outputs": {
            "descriptor": {"path": str(adapted_descriptor), "sha256": sha256(adapted_descriptor)},
            "payloads": outputs,
        },
        "limitations": [
            "Python-adapted descriptor and canonical RAS arrays; not a Java raw-atlas descriptor.",
            "NIfTI affine centers atlas index dimensions/2 at RAS origin to match the selected centered AtlasEM transform convention.",
            "Source MIPAV XML headers are not used because their 1 mm unknown-orientation metadata conflicts with the selected descriptor's 0.9375 mm geometry.",
        ],
    }
    manifest_path = output_dir / "adaptation-manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return adapted_descriptor


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--descriptor", type=Path, required=True, help="Source ibsr11objV3 descriptor")
    parser.add_argument("--output-dir", type=Path, required=True, help="New or empty output directory")
    args = parser.parse_args(argv)
    try:
        result = adapt(args.descriptor, args.output_dir)
    except (OSError, ValueError, KeyError) as exc:
        print(f"Atlas adaptation failed: {exc}", file=sys.stderr)
        return 2
    print(result)
    print(result.parent / "adaptation-manifest.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
