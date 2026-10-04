from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_atlas import (
    ToadsAtlas,
    ToadsAtlasError,
    ToadsAtlasImage,
    ToadsStructure,
    load_toads_atlas,
    load_toads_atlas_data,
    load_toads_atlas_ras_data,
    toads_atlas_shapes_in_structure_order,
)


ROOT = Path(__file__).parents[1]
ATLAS_ROOT = ROOT / "data" / "atlases"


def test_bundled_atlases_are_hashed_ras_nifti_and_source_matched() -> None:
    audit = json.loads((ATLAS_ROOT / "CONVERSION_AUDIT.json").read_text())
    checksum_entries = {
        line.split("  ", 1)[1]: line.split("  ", 1)[0]
        for line in (ATLAS_ROOT / "SHA256SUMS").read_text().splitlines()
    }
    assert len(audit["images"]) == 35
    assert len(list(ATLAS_ROOT.rglob("*.nii.gz"))) == 35
    assert not list(ATLAS_ROOT.rglob("*.raw"))
    assert not list(ATLAS_ROOT.rglob("*.img"))
    assert not list(ATLAS_ROOT.rglob("*.hdr"))

    for row in audit["images"]:
        image_path = ATLAS_ROOT / row["nifti_output"]
        relative_path = image_path.relative_to(ATLAS_ROOT).as_posix()
        image = nib.load(str(image_path))
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert row["output_axis_codes"] == ["R", "A", "S"]
        assert row["source_voxel_values_exact_after_ras_reorientation"] is True
        assert checksum_entries[relative_path] == hashlib.sha256(image_path.read_bytes()).hexdigest()


def test_atlas_descriptors_reference_only_compressed_nifti() -> None:
    for descriptor in ATLAS_ROOT.rglob("*.txt"):
        content = descriptor.read_text()
        assert ".raw" not in content
        assert ".img" not in content
        assert not content or ".nii.gz" in content or descriptor.name.startswith("oasis-")

    toads_descriptor = ATLAS_ROOT / "cruise-atlas-10obj-toads2008.txt"
    assert toads_descriptor.is_file()
    for line in toads_descriptor.read_text().splitlines():
        if line.startswith("File:"):
            referenced = line.split(":", 1)[1].strip()
            assert referenced.endswith(".nii.gz")
            assert (ATLAS_ROOT / referenced).is_file(), referenced


def test_selected_toads_2008_descriptor_resolves_all_ras_atlas_images() -> None:
    atlas = load_toads_atlas(ATLAS_ROOT / "cruise-atlas-10obj-toads2008.txt")

    assert len(atlas.structures) == 11
    assert {item.label for item in atlas.structures} == {1, 5, 6, 14, 15, 16, 17, 18, 22, 24, 25}
    assert atlas.topology.shape == (172, 216, 164)
    assert atlas.topology.voxel_sizes_mm == (1.0, 1.0, 1.0)
    assert atlas.topology.descriptor_voxel_sizes_mm == (0.9375, 0.9375, 0.9375)
    assert len(atlas.shapes) == 11
    assert atlas.shape_samples == 8
    assert atlas.shape_dimensions == (172, 216, 164)
    assert atlas.shape_voxel_sizes_mm == (0.9375, 0.9375, 0.9375)
    assert all(
        image.descriptor_voxel_sizes_mm == (0.9375, 0.9375, 0.9375)
        for image in atlas.shapes.values()
    )
    assert atlas.orientations == (1, 1, 6, 3)
    assert len(atlas.geometry_warnings) == 1
    assert "topology NIfTI voxel sizes" in atlas.geometry_warnings[0]
    assert "descriptor resolutions" in atlas.geometry_warnings[0]
    assert "Java loader uses descriptor spacing" in atlas.geometry_warnings[0]
    assert "without the image-header affine" in atlas.geometry_warnings[0]
    assert all(image.path.suffixes[-2:] == [".nii", ".gz"] for image in atlas.shapes.values())


def test_toads_descriptor_rejects_missing_atlas_images(tmp_path) -> None:
    descriptor = ATLAS_ROOT / "cruise-atlas-10obj-toads2008.txt"
    broken = tmp_path / descriptor.name
    broken.write_text(descriptor.read_text().replace(".nii.gz", "-absent.nii.gz"))

    with pytest.raises(ToadsAtlasError, match="Missing topology atlas image"):
        load_toads_atlas(broken)


@pytest.mark.parametrize("source_codes", [("L", "P", "S"), ("L", "I", "A")])
def test_toads_atlas_arrays_restore_original_java_index_order(tmp_path, source_codes) -> None:
    shape = (2, 3, 4)
    source = np.arange(np.prod(shape), dtype=np.int16).reshape(shape)
    source_orientation = nib.orientations.axcodes2ornt(source_codes)
    ras_orientation = nib.orientations.axcodes2ornt(("R", "A", "S"))
    to_ras = nib.orientations.ornt_transform(source_orientation, ras_orientation)
    ras_array = nib.orientations.apply_orientation(source, to_ras)
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    affine[0, 3] = -(shape[0] - 1)
    affine[1, 3] = -(shape[1] - 1)
    topology_path = tmp_path / "topology.nii.gz"
    shape_path = tmp_path / "shape.nii.gz"
    nib.save(nib.Nifti1Image(ras_array, affine), topology_path)
    nib.save(nib.Nifti1Image(ras_array, affine), shape_path)
    descriptor = tmp_path / "atlas.txt"
    descriptor.write_text("synthetic descriptor")
    audit_path = tmp_path / "CONVERSION_AUDIT.json"
    audit_path.write_text(
        json.dumps(
            {
                "images": [
                    {
                        "nifti_output": name,
                        "input_axis_codes": list(source_codes),
                        "input_shape": list(shape),
                    }
                    for name in ("topology.nii.gz", "shape.nii.gz")
                ]
            }
        )
    )
    ras_shape = tuple(int(value) for value in ras_array.shape)
    topology = ToadsAtlasImage("topology", topology_path, ras_shape, (1.0,) * 3, (1.0,) * 3)
    shape_image = ToadsAtlasImage("shape Test", shape_path, ras_shape, (1.0,) * 3, (1.0,) * 3)
    atlas = ToadsAtlas(
        descriptor,
        (),
        topology,
        {"Test": shape_image},
        1,
        shape,
        (1.0,) * 3,
        (1, 1, 1, 1),
        (),
    )

    data = load_toads_atlas_data(atlas, audit_path)

    np.testing.assert_array_equal(data.topology, source)
    np.testing.assert_array_equal(data.shapes["Test"], source)
    assert data.source_axis_codes["topology.nii.gz"] == source_codes


def test_toads_atlas_python_arrays_reorient_to_ras_and_follow_descriptor_order(tmp_path) -> None:
    source = np.arange(2 * 3 * 4, dtype=np.int16).reshape((2, 3, 4))
    affine = np.diag((-1.0, -1.0, 1.0, 1.0))
    affine[:3, 3] = (1.0, 2.0, 0.0)
    topology_path = tmp_path / "topology.nii.gz"
    background_path = tmp_path / "background.nii.gz"
    foo_path = tmp_path / "foo.nii.gz"
    for path, values in (
        (topology_path, source),
        (background_path, source + 100),
        (foo_path, source + 200),
    ):
        nib.save(nib.Nifti1Image(values, affine), path)
    atlas = ToadsAtlas(
        tmp_path / "atlas.txt",
        (ToadsStructure("Foo", 2, "object"), ToadsStructure("Background", 1, "mask")),
        ToadsAtlasImage("topology.nii.gz", topology_path, source.shape, (1.0,) * 3, (0.8,) * 3),
        {
            "Background": ToadsAtlasImage("background", background_path, source.shape, (1.0,) * 3, (0.8,) * 3),
            "Foo": ToadsAtlasImage("foo", foo_path, source.shape, (1.0,) * 3, (0.8,) * 3),
        },
        1,
        source.shape,
        (0.8, 0.8, 0.8),
        (1, 1, 1, 1),
        (),
    )

    with pytest.warns(UserWarning, match="reorienting the image to RAS") as captured:
        ras = load_toads_atlas_ras_data(atlas)

    expected = source[::-1, ::-1, :]
    assert len(captured) == 3
    assert ras.source_axis_codes["topology.nii.gz"] == ("L", "P", "S")
    assert ras.template_resolution_mm == (0.8, 0.8, 0.8)
    np.testing.assert_array_equal(ras.topology, expected)
    np.testing.assert_array_equal(ras.shapes["Background"], expected + 100)
    np.testing.assert_array_equal(ras.shapes["Foo"], expected + 200)
    stacked = toads_atlas_shapes_in_structure_order(atlas, ras)
    assert stacked.shape == (*source.shape, 2)
    np.testing.assert_array_equal(stacked[..., 0], (expected + 200).astype(np.float32))
    np.testing.assert_array_equal(stacked[..., 1], (expected + 100).astype(np.float32))
