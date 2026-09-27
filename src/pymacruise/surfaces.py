"""Scientific triangular surfaces in RAS+ millimetres.

Level-set extraction follows MaCRUISE's ``GetMeshSurfaceFromLevelSet``
contract (zero level, connectivity 18, non-inclusive exact zeros, and the
source's near-zero perturbation). Marching Cubes is supplied by scikit-image;
vertices are transformed through the complete NIfTI affine before export.

Source references: ``GetMeshSurfaceFromLevelSet.java``, CVS r1.4; JIST
``IsoSurfaceOnGrid.java``, CVS r1.3, author Blake Lucas. The Python triangulator
is an explicit, documented implementation substitution; it does not copy
JIST source code.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import warnings
from typing import Mapping

import nibabel as nib
import numpy as np
from skimage.measure import marching_cubes

from pymacruise.image_io import normalize_image_to_ras


@dataclass(frozen=True)
class SurfaceMesh:
    """A triangular mesh whose coordinates are RAS+ millimetres."""

    vertices: np.ndarray
    faces: np.ndarray
    coordinate_frame: str = "RAS"
    units: str = "mm"
    point_data: Mapping[str, np.ndarray] = field(default_factory=dict)

    def __post_init__(self) -> None:
        vertices = np.asarray(self.vertices, dtype=np.float64)
        faces = np.asarray(self.faces, dtype=np.int64)
        if vertices.ndim != 2 or vertices.shape[1] != 3 or vertices.shape[0] < 4:
            raise ValueError("Surface vertices must have shape (N, 3), N >= 4")
        if faces.ndim != 2 or faces.shape[1] != 3 or faces.shape[0] < 4:
            raise ValueError("Surface faces must have shape (M, 3), M >= 4")
        if not np.isfinite(vertices).all():
            raise ValueError("Surface vertices contain NaN or infinite coordinates")
        if np.any(faces < 0) or np.any(faces >= len(vertices)):
            raise ValueError("Surface face index is outside the vertex array")
        if np.any((faces[:, 0] == faces[:, 1]) | (faces[:, 1] == faces[:, 2]) | (faces[:, 0] == faces[:, 2])):
            raise ValueError("Surface contains a face with repeated vertex indices")
        if self.coordinate_frame != "RAS" or self.units != "mm":
            raise ValueError("SurfaceMesh coordinates must be RAS+ millimetres")
        point_data: dict[str, np.ndarray] = {}
        for name, values in self.point_data.items():
            array = np.asarray(values, dtype=np.float64)
            if not name or array.ndim not in (1, 2) or array.shape[0] != len(vertices):
                raise ValueError(f"Point data {name!r} must have one value row per vertex")
            if array.ndim == 2 and array.shape[1] < 1:
                raise ValueError(f"Point data {name!r} has no components")
            if not np.isfinite(array).all():
                raise ValueError(f"Point data {name!r} contains NaN or infinite values")
            point_data[name] = array
        object.__setattr__(self, "vertices", vertices)
        object.__setattr__(self, "faces", faces)
        object.__setattr__(self, "point_data", point_data)


def _signed_mesh_volume(vertices: np.ndarray, faces: np.ndarray) -> float:
    triangles = vertices[faces]
    return float(
        np.einsum(
            "ij,ij->i",
            triangles[:, 0],
            np.cross(triangles[:, 1], triangles[:, 2]),
        ).sum()
        / 6.0
    )


def extract_levelset_surface(
    image: nib.spatialimages.SpatialImage,
    *,
    level: float = 0.0,
) -> SurfaceMesh:
    """Extract a closed triangular surface from a 3D scalar NIfTI image.

    This ports the active JIST defaults: iso-level 0, digital connectivity 18,
    and ``inclusive=false``. Source ``solveOriginal`` moves values within
    0.1 of the iso-level to +/-0.1 (exact zeros to -0.01) before extraction.
    """
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError("Level-set input must be a NIfTI image")
    if len(image.shape) != 3:
        raise ValueError("Level-set input must be a 3D scalar image")
    if not np.isfinite(level):
        raise ValueError("Iso-level must be finite")
    ras_image = normalize_image_to_ras(image, source_name="level-set image")
    values = np.asarray(ras_image.dataobj, dtype=np.float32)
    if not np.isfinite(values).all():
        raise ValueError("Level-set image contains NaN or infinite values")

    shifted = values - np.float32(level)
    shifted = shifted.copy()
    shifted[(shifted > np.float32(-0.1)) & (shifted <= 0)] = np.float32(-0.1)
    shifted[(shifted < np.float32(0.1)) & (shifted > 0)] = np.float32(0.1)
    if float(shifted.min()) >= 0 or float(shifted.max()) <= 0:
        raise ValueError("Level-set does not cross the requested iso-level")

    try:
        index_vertices, faces, _normals, _values = marching_cubes(
            shifted,
            level=0.0,
            gradient_direction="ascent",
            allow_degenerate=False,
            method="lewiner",
        )
    except (RuntimeError, ValueError) as exc:
        raise ValueError(f"Could not extract level-set surface: {exc}") from exc

    affine = ras_image.affine
    vertices = index_vertices @ affine[:3, :3].T + affine[:3, 3]
    faces = np.asarray(faces, dtype=np.int64)
    if _signed_mesh_volume(vertices, faces) < 0:
        faces = faces[:, ::-1]
    return SurfaceMesh(vertices, faces)


def convert_surface_to_ras(
    mesh: SurfaceMesh,
    *,
    input_frame: str,
) -> SurfaceMesh:
    """Convert an explicitly identified RAS+ or LPS surface to RAS+ mm."""
    frame = input_frame.upper().replace("+", "")
    if frame in {"RAS", "RASMM"}:
        return SurfaceMesh(mesh.vertices.copy(), mesh.faces.copy(), point_data=mesh.point_data)
    if frame not in {"LPS", "LPSMM"}:
        raise ValueError("Surface coordinate frame must be RAS or LPS")
    warnings.warn(
        "Surface input coordinates are LPS; converting vertices and winding to RAS.",
        UserWarning,
        stacklevel=2,
    )
    vertices = mesh.vertices.copy()
    vertices[:, :2] *= -1
    return SurfaceMesh(
        vertices,
        mesh.faces[:, ::-1].copy(),
        point_data=mesh.point_data,
    )


def surface_vertex_normals(mesh: SurfaceMesh) -> np.ndarray:
    """Compute outward area-weighted unit normals from the triangle winding."""
    triangles = mesh.vertices[mesh.faces]
    face_normals = np.cross(
        triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    )
    normals = np.zeros_like(mesh.vertices)
    for corner in range(3):
        np.add.at(normals, mesh.faces[:, corner], face_normals)
    lengths = np.linalg.norm(normals, axis=1)
    if np.any(lengths <= np.finfo(float).eps):
        raise ValueError("Surface has vertices without a defined normal")
    return normals / lengths[:, None]


def transform_surface_ras(mesh: SurfaceMesh, transform: np.ndarray) -> SurfaceMesh:
    """Apply a full-space RAS millimetre affine to surface points.

    ``transform`` maps source RAS+ millimetres to target RAS+ millimetres.
    Triangle winding is reversed if the affine contains a reflection. Scalar
    point-data fields (including thickness) are preserved without resampling.
    """
    matrix = np.asarray(transform, dtype=np.float64)
    if matrix.shape != (4, 4) or not np.isfinite(matrix).all():
        raise ValueError("Surface transform must be a finite 4x4 affine")
    if not np.allclose(matrix[3], [0, 0, 0, 1], rtol=0, atol=1e-10):
        raise ValueError("Surface transform must be affine in homogeneous coordinates")
    if abs(np.linalg.det(matrix[:3, :3])) < 1e-12:
        raise ValueError("Surface transform has a singular spatial component")
    vertices = mesh.vertices @ matrix[:3, :3].T + matrix[:3, 3]
    faces = mesh.faces[:, ::-1] if np.linalg.det(matrix[:3, :3]) < 0 else mesh.faces
    return SurfaceMesh(vertices, faces.copy(), point_data=mesh.point_data)


def write_legacy_vtk(mesh: SurfaceMesh, path: str | Path) -> Path:
    """Write a scientific mesh as an ASCII legacy VTK POLYDATA file."""
    destination = Path(path)
    if destination.suffix.lower() != ".vtk":
        raise ValueError("Scientific surface outputs must use the .vtk extension")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="ascii", newline="\n") as stream:
        stream.write("# vtk DataFile Version 3.0\n")
        stream.write("PyMaCRUISE coordinate_frame=RAS units=mm\n")
        stream.write("ASCII\nDATASET POLYDATA\n")
        stream.write(f"POINTS {len(mesh.vertices)} double\n")
        np.savetxt(stream, mesh.vertices, fmt="%.10g %.10g %.10g")
        stream.write(f"POLYGONS {len(mesh.faces)} {4 * len(mesh.faces)}\n")
        polygons = np.column_stack((np.full(len(mesh.faces), 3), mesh.faces))
        np.savetxt(stream, polygons, fmt="%d %d %d %d")
        if mesh.point_data:
            stream.write(f"POINT_DATA {len(mesh.vertices)}\n")
            stream.write(f"FIELD PointData {len(mesh.point_data)}\n")
            for name, data in mesh.point_data.items():
                values = data.reshape(len(mesh.vertices), -1)
                stream.write(f"{name} {values.shape[1]} {len(mesh.vertices)} double\n")
                np.savetxt(stream, values, fmt=" ".join(["%.10g"] * values.shape[1]))
    return destination


def read_legacy_vtk(
    path: str | Path,
    *,
    input_frame: str | None = None,
) -> SurfaceMesh:
    """Read ASCII legacy VTK triangular POLYDATA and normalize frame to RAS.

    A frame may be passed explicitly. Files produced here declare RAS in their
    title. Older untagged files require ``input_frame`` rather than silently
    guessing whether their coordinates are RAS or LPS.
    """
    source = Path(path)
    lines = source.read_text(encoding="ascii").splitlines()
    if len(lines) < 5 or not lines[0].startswith("# vtk DataFile"):
        raise ValueError(f"Not a legacy VTK file: {source}")
    if lines[2].strip().upper() != "ASCII":
        raise ValueError("Only ASCII legacy VTK surface files are supported")
    title = lines[1].upper()
    declared_frame = "LPS" if "COORDINATE_FRAME=LPS" in title else (
        "RAS" if "COORDINATE_FRAME=RAS" in title else None
    )
    if input_frame is not None and declared_frame is not None:
        if input_frame.upper().replace("+", "") != declared_frame:
            raise ValueError("Explicit surface frame conflicts with VTK frame metadata")
    frame = input_frame or declared_frame
    if frame is None:
        raise ValueError("Untagged VTK surface requires input_frame='RAS' or 'LPS'")

    tokens = "\n".join(lines[3:]).split()
    try:
        point_header = tokens.index("POINTS")
        n_vertices = int(tokens[point_header + 1])
        point_start = point_header + 3
        n_point_values = n_vertices * 3
        vertices = np.asarray(tokens[point_start : point_start + n_point_values], dtype=float).reshape(-1, 3)
        polygon_header = tokens.index("POLYGONS", point_start + n_point_values)
        n_faces = int(tokens[polygon_header + 1])
        polygon_tokens = tokens[polygon_header + 3 :]
        faces: list[list[int]] = []
        cursor = 0
        for _ in range(n_faces):
            count = int(polygon_tokens[cursor])
            if count != 3:
                raise ValueError("Only triangular VTK POLYDATA faces are supported")
            faces.append([int(value) for value in polygon_tokens[cursor + 1 : cursor + 4]])
            cursor += 4
        point_data: dict[str, np.ndarray] = {}
        trailing = polygon_header + 3 + cursor
        if "POINT_DATA" in tokens[trailing:]:
            point_header = tokens.index("POINT_DATA", trailing)
            n_point_data = int(tokens[point_header + 1])
            if n_point_data != n_vertices:
                raise ValueError("VTK POINT_DATA count does not match the vertex count")
            data_cursor = point_header + 2
            while data_cursor < len(tokens):
                kind = tokens[data_cursor]
                if kind == "FIELD":
                    n_arrays = int(tokens[data_cursor + 2])
                    data_cursor += 3
                    for _ in range(n_arrays):
                        name = tokens[data_cursor]
                        n_components = int(tokens[data_cursor + 1])
                        n_tuples = int(tokens[data_cursor + 2])
                        data_cursor += 4
                        if n_tuples != n_vertices:
                            raise ValueError(f"VTK point-data field {name!r} has the wrong tuple count")
                        count = n_components * n_tuples
                        field_values = np.asarray(
                            tokens[data_cursor : data_cursor + count], dtype=np.float64
                        ).reshape(n_tuples, n_components)
                        point_data[name] = (
                            field_values[:, 0] if n_components == 1 else field_values
                        )
                        data_cursor += count
                elif kind == "SCALARS":
                    name = tokens[data_cursor + 1]
                    n_components = int(tokens[data_cursor + 3]) if tokens[data_cursor + 3] != "LOOKUP_TABLE" else 1
                    data_cursor += 4 if n_components > 1 else 3
                    if tokens[data_cursor] != "LOOKUP_TABLE":
                        raise ValueError("Malformed VTK SCALARS lookup table")
                    data_cursor += 2
                    count = n_vertices * n_components
                    scalar_values = np.asarray(
                        tokens[data_cursor : data_cursor + count], dtype=np.float64
                    ).reshape(n_vertices, n_components)
                    point_data[name] = (
                        scalar_values[:, 0] if n_components == 1 else scalar_values
                    )
                    data_cursor += count
                else:
                    raise ValueError(f"Unsupported VTK point-data section {kind!r}")
        mesh = SurfaceMesh(vertices, np.asarray(faces, dtype=np.int64), point_data=point_data)
    except (ValueError, IndexError) as exc:
        raise ValueError(f"Malformed triangular legacy VTK surface: {source}") from exc
    return convert_surface_to_ras(mesh, input_frame=frame)
