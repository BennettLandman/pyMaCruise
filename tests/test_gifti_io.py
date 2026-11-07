from dataclasses import replace
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.cli import main
from pymacruise.gifti_io import (
    read_gifti, write_gifti, load_gifti_point_data, write_gifti_point_data,
)
from pymacruise.surfaces import SurfaceMesh, read_surface, write_surface, transform_surface_ras


@pytest.fixture
def mesh():
    return SurfaceMesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]),
        point_data={"thickness_mm": np.array([1.1, 2.2, 3.3, 4.4]),
                    "directions": np.arange(12).reshape(4, 3)},
        metadata={"AnatomicalStructurePrimary": "CortexLeft", "GeometricType": "Anatomical"},
        point_data_metadata={"thickness_mm": {"Description": "Example thickness", "Units": "mm"}},
    )


def test_geometry_and_information_roundtrip(mesh, tmp_path):
    path = write_surface(mesh, tmp_path / "cortex.surf.gii")
    # Check the public format with an independent NiBabel load.
    image = nib.load(path)
    assert image.get_arrays_from_intent("pointset")[0].data.dtype == np.float32
    assert image.get_arrays_from_intent("triangle")[0].data.dtype == np.int32
    assert image.get_arrays_from_intent("pointset")[0].coordsys.dataspace == 0
    restored = read_surface(path)
    np.testing.assert_array_equal(restored.faces, mesh.faces)
    np.testing.assert_allclose(restored.vertices, mesh.vertices)
    for name in mesh.point_data:
        np.testing.assert_allclose(restored.point_data[name], mesh.point_data[name], rtol=1e-6)
    assert restored.metadata["AnatomicalStructurePrimary"] == "CortexLeft"
    assert restored.point_data_metadata["thickness_mm"]["Units"] == "mm"
    assert restored.point_data_intents["thickness_mm"] == 2005
    assert transform_surface_ras(restored, np.eye(4)).metadata == restored.metadata


@pytest.mark.parametrize("encoding", ["ASCII", "Base64Binary", "GZipBase64Binary"])
def test_external_file_transform_is_applied_once(mesh, tmp_path, encoding):
    matrix = np.diag([-2., 3., 4., 1.])
    matrix[:3, 3] = [7, 8, 9]
    points = nib.gifti.GiftiDataArray(mesh.vertices.astype(np.float32), intent="pointset",
        encoding=encoding, coordsys=nib.gifti.GiftiCoordSystem(0, 4, matrix))
    triangles = nib.gifti.GiftiDataArray(mesh.faces.astype(np.int32), intent="triangle", encoding=encoding)
    source = tmp_path / "external.surf.gii"
    # Triangles deliberately precede points: intent, not array order, defines geometry.
    nib.save(nib.gifti.GiftiImage(darrays=[triangles, points]), source)
    restored = read_gifti(source)
    np.testing.assert_allclose(restored.vertices, mesh.vertices @ matrix[:3, :3].T + matrix[:3, 3])
    np.testing.assert_array_equal(restored.faces, mesh.faces[:, ::-1])
    assert restored.gifti_dataspace == 4
    assert transform_surface_ras(restored, matrix).gifti_dataspace == 0
    again = read_gifti(write_gifti(restored, tmp_path / "again.surf.gii"))
    np.testing.assert_array_equal(again.vertices, restored.vertices)
    assert again.gifti_dataspace == 4


def test_separate_vertex_data_and_labels(mesh, tmp_path):
    table = nib.gifti.GiftiLabelTable()
    label = nib.gifti.GiftiLabel(1, red=1, green=0, blue=0, alpha=1)
    label.label = "Cortex"
    table.labels.append(label)
    labeled = replace(mesh, point_data={"region": np.ones(4)},
                      point_data_intents={"region": 1002}, label_table=table)
    data_path = write_gifti_point_data(labeled, tmp_path / "region.label.gii")
    assert nib.load(data_path).darrays[0].data.dtype == np.int32
    restored = load_gifti_point_data(mesh, data_path)
    assert restored.label_table.get_labels_as_dict() == {1: "Cortex"}
    assert restored.point_data_intents["region"] == 1002
    again = read_gifti(write_gifti(restored, tmp_path / "labeled.surf.gii"))
    assert again.label_table.labels[0].rgba == (1, 0, 0, 1)
    with pytest.raises(ValueError, match="Duplicate"):
        load_gifti_point_data(restored, data_path)
    with pytest.raises(ValueError, match="one POINTSET"):
        read_gifti(data_path)


def test_malformed_or_ambiguous_files_are_rejected(mesh, tmp_path):
    path = write_gifti(mesh, tmp_path / "surface.gii")
    image = nib.load(path)
    image.darrays[-1] = nib.gifti.GiftiDataArray(np.zeros((3, 3), np.float32),
        meta=nib.gifti.GiftiMetaData({"Name": "directions"}))
    nib.save(image, path)
    with pytest.raises(ValueError, match="one row per vertex"):
        read_gifti(path)
    image = nib.load(write_gifti(mesh, path))
    image.darrays[0].coordsys.xform[0, 3] = 10
    nib.save(image, path)
    with pytest.raises(ValueError, match="known target space"):
        read_gifti(path)
    image.darrays[0].coordsys.xformspace = 1
    image.darrays[0].coordsys.xform[3, 3] = 0
    nib.save(image, path)
    with pytest.raises(ValueError, match="homogeneous"):
        read_gifti(path)
    with pytest.raises(ValueError, match=".vtk or .gii"):
        write_surface(mesh, tmp_path / "mesh.obj")


def test_sparse_data_rejected_and_vertex_counts_checked(mesh, tmp_path):
    path = tmp_path / "data.shape.gii"
    for intent, data, message in [
        (2002, np.arange(4, dtype=np.int32), "NODE_INDEX"),
        (2005, np.zeros(3, np.float32), "one row per vertex"),
    ]:
        nib.save(nib.gifti.GiftiImage(darrays=[nib.gifti.GiftiDataArray(data, intent=intent)]), path)
        with pytest.raises(ValueError, match=message):
            load_gifti_point_data(mesh, path)


def test_cli_conversion_data_export_and_info(mesh, tmp_path, capsys):
    vtk = write_surface(mesh, tmp_path / "source.vtk")
    gifti = tmp_path / "surface.surf.gii"
    data = tmp_path / "metrics.shape.gii"
    provenance = tmp_path / "provenance.json"
    assert main(["surface-convert", "--input", str(vtk), "--output", str(gifti),
                 "--data-output", str(data), "--provenance", str(provenance)]) == 0
    assert json.loads(provenance.read_text())["outputs"]["point_data"] == str(data)
    assert len(nib.load(data).darrays) == 2
    assert main(["surface-convert", "--input", str(gifti), "--output", str(tmp_path / "back.vtk"),
                 "--provenance", str(provenance)]) == 0
    capsys.readouterr()
    assert main(["surface-info", "--input", str(gifti)]) == 0
    info = json.loads(capsys.readouterr().out)
    assert info["vertex_count"] == 4
    assert info["point_data"]["thickness_mm"]["shape"] == [4]


def test_cli_extract_and_thickness_gifti(tmp_path):
    fixtures = Path(__file__).parent / "data/synthetic/surfaces/thickness"
    for name in ("inner", "outer"):
        assert main(["surface-extract", "--level-set", str(fixtures / f"{name}.nii.gz"),
                     "--output", str(tmp_path / f"{name}.surf.gii")]) == 0
    output = tmp_path / "thickness.surf.gii"
    assert main(["thickness-surface", "--inner", str(tmp_path / "inner.surf.gii"),
                 "--outer", str(tmp_path / "outer.surf.gii"), "--inner-output", str(output)]) == 0
    measured = read_gifti(output)
    assert "mean_thickness_mm" in measured.point_data
    assert np.isfinite(measured.point_data["mean_thickness_mm"]).all()


def test_open_triangle_surface_supported(tmp_path):
    mesh = SurfaceMesh(np.eye(3), np.array([[0, 1, 2]]))
    restored = read_gifti(write_gifti(mesh, tmp_path / "triangle.surf.gii"))
    np.testing.assert_array_equal(restored.faces, mesh.faces)


@pytest.mark.parametrize("format", ["vtk", "gii"])
def test_base_surface_thickness_output_format(tmp_path, format):
    fixtures = Path(__file__).parent / "data/synthetic/surfaces/thickness"
    output = tmp_path / "products"
    provenance = tmp_path / "run.json"
    assert main(["base-surface-thickness", "--inner-levelset", str(fixtures / "inner.nii.gz"),
                 "--central-levelset", str(fixtures / "central.nii.gz"),
                 "--outer-levelset", str(fixtures / "outer.nii.gz"),
                 "--output-dir", str(output), "--surface-format", format,
                 "--provenance", str(provenance)]) == 0
    surface = output / f"inner_surface_thickness.{format}"
    assert "mean_thickness_mm" in read_surface(surface).point_data
    assert json.loads(provenance.read_text())["outputs"]["inner_surface_with_thickness"] == str(surface)
