"""GIFTI geometry and vertex measurements, using NiBabel without a renderer.

Coordinates are RAS+ millimetres. A declared coordinate transform is applied
once on import and written as identity on export. Unknown space stays unknown;
export never asserts scanner or template registration that was not supplied.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from typing import Sequence

import nibabel as nib
import numpy as np

from pymacruise.surfaces import SurfaceMesh, transform_surface_ras

POINTSET = 1008
TRIANGLE = 1009
LABEL = 1002
SHAPE = 2005


def _load(path: str | Path) -> nib.gifti.GiftiImage:
    image = nib.load(str(path))
    if not isinstance(image, nib.gifti.GiftiImage):
        raise ValueError(f"Expected a GIFTI file: {path}")
    return image


def _destination(path: str | Path) -> Path:
    destination = Path(path)
    if destination.suffix.lower() != ".gii":
        raise ValueError("GIFTI outputs must use .gii (for example .surf.gii or .shape.gii)")
    return destination


def _with_point_data(mesh: SurfaceMesh, image: nib.gifti.GiftiImage) -> SurfaceMesh:
    values = dict(mesh.point_data)
    metadata = dict(mesh.point_data_metadata)
    intents = dict(mesh.point_data_intents)
    for index, array in enumerate(image.darrays):
        if array.intent in (POINTSET, TRIANGLE):
            continue
        if array.intent == 2002:
            raise ValueError("Sparse GIFTI NODE_INDEX data is not supported")
        name = array.meta.get("Name") or f"array_{index}"
        if name in values:
            raise ValueError(f"Duplicate GIFTI point-data name: {name!r}")
        data = np.asarray(array.data)
        if data.ndim not in (1, 2) or data.shape[0] != len(mesh.vertices):
            raise ValueError(f"GIFTI point data {name!r} must have one row per vertex")
        if array.intent == LABEL and (
            not np.issubdtype(data.dtype, np.integer)
        ):
            raise ValueError("GIFTI label arrays must have integer values")
        values[name] = data
        metadata[name] = dict(array.meta)
        intents[name] = int(array.intent)
    table = mesh.label_table
    if image.labeltable.labels:
        if table is not None and table.labels and table.to_xml() != image.labeltable.to_xml():
            raise ValueError("Conflicting GIFTI label tables")
        table = deepcopy(image.labeltable)
    return replace(mesh, point_data=values, point_data_metadata=metadata,
                   point_data_intents=intents, label_table=table)


def read_gifti(path: str | Path) -> SurfaceMesh:
    """Read one triangulated surface and its embedded vertex information.

    ASCII and base64 (including compressed) arrays are handled by NiBabel.
    Nonidentity transforms require a known target space. Reflection reverses
    triangle winding; vertex information remains attached by vertex index.
    """
    image = _load(path)
    points = image.get_arrays_from_intent(POINTSET)
    triangles = image.get_arrays_from_intent(TRIANGLE)
    if len(points) != 1 or len(triangles) != 1:
        raise ValueError("GIFTI surface requires exactly one POINTSET and one TRIANGLE array")
    if not np.issubdtype(triangles[0].data.dtype, np.integer):
        raise ValueError("GIFTI triangle indices must be integers")
    meta = dict(image.meta)
    # Pointset metadata includes anatomical structure and surface type.
    # Keep these entries on the mesh as well as file-level information.
    meta.update(dict(points[0].meta))
    if meta.get("CoordinateFrame", "RAS").upper().replace("+", "") != "RAS":
        raise ValueError("GIFTI coordinates must be RAS+ millimetres")
    if meta.get("Units", "mm") != "mm":
        raise ValueError("GIFTI coordinates must be RAS+ millimetres")
    mesh = SurfaceMesh(points[0].data, triangles[0].data, metadata=meta)
    coords = points[0].coordsys
    if coords is not None:
        matrix = np.asarray(coords.xform, dtype=np.float64)
        identity = np.array_equal(matrix, np.eye(4))
        if not identity and coords.xformspace == 0:
            raise ValueError("Nonidentity GIFTI coordinate transform requires a known target space")
        mesh = transform_surface_ras(mesh, matrix)
        mesh = replace(mesh, gifti_dataspace=int(coords.xformspace or coords.dataspace))
    return _with_point_data(mesh, image)


def load_gifti_point_data(mesh: SurfaceMesh, path: str | Path) -> SurfaceMesh:
    """Attach a dense .shape/.func/.label.gii file in the same vertex order.

    Matching counts cannot prove vertex correspondence; the caller must use
    measurements generated for this mesh. Sparse NODE_INDEX files are rejected.
    """
    image = _load(path)
    if not image.darrays or any(a.intent in (POINTSET, TRIANGLE, 2002) for a in image.darrays):
        raise ValueError("Expected dense GIFTI vertex data without geometry or NODE_INDEX")
    metadata = dict(image.meta)
    metadata.update(mesh.metadata)
    return _with_point_data(replace(mesh, metadata=metadata), image)


def _data_arrays(mesh: SurfaceMesh, names: Sequence[str] | None = None) -> list[nib.gifti.GiftiDataArray]:
    arrays = []
    for name in mesh.point_data if names is None else names:
        data = mesh.point_data[name]
        intent = mesh.point_data_intents.get(name, SHAPE if data.ndim == 1 else 0)
        if intent in (POINTSET, TRIANGLE, 2002):
            raise ValueError(f"Invalid GIFTI vertex-data intent for {name!r}")
        if intent == LABEL:
            if np.any(data != np.floor(data)) or np.any(data < np.iinfo(np.int32).min) or np.any(data > np.iinfo(np.int32).max):
                raise ValueError(f"GIFTI labels {name!r} must fit int32")
            encoded = data.astype(np.int32)
        else:
            encoded = _float32(data, name)
        meta = dict(mesh.point_data_metadata.get(name, {}))
        meta["Name"] = name
        arrays.append(nib.gifti.GiftiDataArray(encoded, intent=intent,
                      meta=nib.gifti.GiftiMetaData(meta)))
    return arrays


def _float32(data: np.ndarray, name: str) -> np.ndarray:
    if np.any(np.abs(data) > np.finfo(np.float32).max):
        raise ValueError(f"GIFTI {name!r} values exceed the float32 range")
    return np.asarray(data, dtype=np.float32)


def _save(mesh: SurfaceMesh, arrays: list[nib.gifti.GiftiDataArray], path: Path) -> Path:
    meta = dict(mesh.metadata)
    meta.update(CoordinateFrame="RAS", Units="mm")
    image = nib.gifti.GiftiImage(darrays=arrays, meta=nib.gifti.GiftiMetaData(meta),
                               labeltable=deepcopy(mesh.label_table))
    path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(image, str(path))
    return path


def write_gifti(mesh: SurfaceMesh, path: str | Path) -> Path:
    """Write geometry plus named vertex arrays using float32/int32 GIFTI.

    Float32 export is interoperable with standard GIFTI readers; it can round
    binary64 coordinates and measurements. VTK remains the lossless geometry
    format for comparisons requiring binary64 coordinates.
    """
    destination = _destination(path)
    if np.max(mesh.faces) > np.iinfo(np.int32).max:
        raise ValueError("GIFTI triangle indices exceed the int32 range")
    coords = nib.gifti.GiftiCoordSystem(mesh.gifti_dataspace, mesh.gifti_dataspace, np.eye(4))
    arrays = [
        nib.gifti.GiftiDataArray(_float32(mesh.vertices, "coordinates"), intent=POINTSET, coordsys=coords,
                                meta=nib.gifti.GiftiMetaData(dict(mesh.metadata))),
        nib.gifti.GiftiDataArray(mesh.faces.astype(np.int32), intent=TRIANGLE),
        *_data_arrays(mesh),
    ]
    return _save(mesh, arrays, destination)


def write_gifti_point_data(
    mesh: SurfaceMesh, path: str | Path, *, names: Sequence[str] | None = None,
) -> Path:
    """Write measurements without geometry, optionally selecting field names."""
    destination = _destination(path)
    arrays = _data_arrays(mesh, names)
    if not arrays:
        raise ValueError("GIFTI vertex-data output requires at least one field")
    return _save(mesh, arrays, destination)
