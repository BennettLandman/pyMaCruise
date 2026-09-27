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

from dataclasses import dataclass
from pathlib import Path
import warnings

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
        object.__setattr__(self, "vertices", vertices)
        object.__setattr__(self, "faces", faces)


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
        return SurfaceMesh(mesh.vertices.copy(), mesh.faces.copy())
    if frame not in {"LPS", "LPSMM"}:
        raise ValueError("Surface coordinate frame must be RAS or LPS")
    warnings.warn(
        "Surface input coordinates are LPS; converting vertices and winding to RAS.",
        UserWarning,
        stacklevel=2,
    )
    vertices = mesh.vertices.copy()
    vertices[:, :2] *= -1
    return SurfaceMesh(vertices, mesh.faces[:, ::-1].copy())


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
        mesh = SurfaceMesh(vertices, np.asarray(faces, dtype=np.int64))
    except (ValueError, IndexError) as exc:
        raise ValueError(f"Malformed triangular legacy VTK surface: {source}") from exc
    return convert_surface_to_ras(mesh, input_frame=frame)
