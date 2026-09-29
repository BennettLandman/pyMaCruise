from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.lesion_atlas_priors import (
    LesionAtlasDescriptorError,
    load_lesion_atlas_topology,
    load_lesion_atlas_shape_priors,
    load_lesion_atlas_priors,
)


ROOT = Path(__file__).parents[1]
ATLAS_ROOT = ROOT / "data" / "atlases" / "lesionToads-atlas-2012"
SYNTHETIC = Path(__file__).parent / "data" / "synthetic" / "lesion_atlas_priors"


def test_minimal_lesion_atlas_preserves_structure_order_and_modality_order() -> None:
    priors = load_lesion_atlas_priors(
        SYNTHETIC / "minimal-lesion-atlas.txt", ["FLAIR", "MPRAGE"]
    )

    assert [structure.name for structure in priors.structures] == [
        "Background",
        "Lesion-WM",
        "Cortex",
    ]
    assert priors.modalities == ("FLAIR", "MPRAGE")
    np.testing.assert_array_equal(
        priors.intensity_priors,
        np.asarray([[0.1, 1.2, 0.3], [0.0, 0.4, 0.8]], dtype=np.float32),
    )
    np.testing.assert_array_equal(
        priors.lesion_priors, np.asarray([1.2, 0.4], dtype=np.float32)
    )
    np.testing.assert_array_equal(
        priors.black_hole_priors, np.asarray([0.6, 0.1], dtype=np.float32)
    )
    assert not priors.intensity_priors.flags.writeable


def test_bundled_2012_descriptor_matches_selected_t1_flair_priors() -> None:
    priors = load_lesion_atlas_priors(
        ATLAS_ROOT / "cruise-atlas-12obj-lesiontoads2012.txt",
        ["T1_MPRAGE", "FLAIR"],
    )

    assert [structure.name for structure in priors.structures] == [
        "Background",
        "Sulcal-CSF",
        "Ventricles",
        "Cerebellum-GM",
        "Cerebrum-GM",
        "Caudate",
        "Thalamus",
        "Putamen",
        "Brainstem",
        "Cerebellum-WM",
        "OuterVentricular-WM",
        "InterVentricular-WM",
        "Frontal-WM",
    ]
    np.testing.assert_array_equal(
        priors.intensity_priors[0],
        np.asarray(
            [0.0, 0.1, 0.1, 0.7, 0.7, 0.6, 0.6, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0],
            dtype=np.float32,
        ),
    )
    np.testing.assert_array_equal(
        priors.intensity_priors[1],
        np.asarray(
            [0.0, 0.5, 0.5, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7, 0.7, 0.7, 1.2, 0.75],
            dtype=np.float32,
        ),
    )
    np.testing.assert_array_equal(
        priors.lesion_priors, np.asarray([0.4, 1.2], dtype=np.float32)
    )
    np.testing.assert_array_equal(
        priors.black_hole_priors, np.asarray([0.1, 0.6], dtype=np.float32)
    )


def test_lesion_atlas_parser_rejects_incomplete_selected_modality(tmp_path) -> None:
    source = (SYNTHETIC / "minimal-lesion-atlas.txt").read_text()
    descriptor = tmp_path / "incomplete.txt"
    descriptor.write_text(source.replace("FLAIR\t0.1\t1.2\t0.3\n", ""))

    with pytest.raises(LesionAtlasDescriptorError, match="No intensity row.*FLAIR"):
        load_lesion_atlas_priors(descriptor, ["T1_MPRAGE", "FLAIR"])


def test_lesion_atlas_parser_rejects_row_with_wrong_structure_count(tmp_path) -> None:
    source = (SYNTHETIC / "minimal-lesion-atlas.txt").read_text()
    descriptor = tmp_path / "wrong-count.txt"
    descriptor.write_text(source.replace("T1_MPRAGE\t0.0\t0.4\t0.8", "T1_MPRAGE\t0.0\t0.4"))

    with pytest.raises(LesionAtlasDescriptorError, match="expected 3"):
        load_lesion_atlas_priors(descriptor, ["T1_MPRAGE"])


def test_bundled_lesion_topology_is_loaded_as_ras_nifti_only() -> None:
    topology = load_lesion_atlas_topology(
        ATLAS_ROOT / "cruise-atlas-12obj-lesiontoads2012.txt"
    )

    assert topology.path.name == "brain-atlas-cruise-12obj-lesions-topology_calc.nii.gz"
    assert topology.path.is_file()
    assert topology.image.shape == (172, 216, 164)
    assert nib.aff2axcodes(topology.image.affine) == ("R", "A", "S")


def test_lesion_topology_loader_reorients_non_ras_and_does_not_load_shape_maps(
    tmp_path,
) -> None:
    # The unreferenced shape file intentionally does not exist: this step only
    # resolves and loads the topology atlas.
    values = np.arange(24, dtype=np.int16).reshape((2, 3, 4))
    source_affine = np.diag([-1.0, -1.0, 1.0, 1.0])
    source_affine[0, 3] = 1.0
    source_affine[1, 3] = 2.0
    atlas_path = tmp_path / "topology.nii.gz"
    nib.save(nib.Nifti1Image(values, source_affine), atlas_path)
    descriptor = tmp_path / "atlas.txt"
    descriptor.write_text(
        "Structure Atlas File (edit at your own risks)\n"
        "Topology Atlas:\n"
        "File:\ttopology.nii.gz\n"
        "Shape Atlas:\t8\n"
        "File:\tmissing-shape.nii.gz\n"
    )

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        topology = load_lesion_atlas_topology(descriptor)

    expected = nib.orientations.apply_orientation(
        values,
        nib.orientations.ornt_transform(
            nib.orientations.axcodes2ornt(("L", "P", "S")),
            nib.orientations.axcodes2ornt(("R", "A", "S")),
        ),
    )
    assert nib.aff2axcodes(topology.image.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(np.asanyarray(topology.image.dataobj), expected)


@pytest.mark.parametrize(
    "reference, message",
    [
        ("topology.nii", "compressed .nii.gz"),
        ("../outside.nii.gz", "Unsafe topology atlas path"),
    ],
)
def test_lesion_topology_loader_rejects_nonbundled_path_forms(
    tmp_path, reference, message
) -> None:
    descriptor = tmp_path / "atlas.txt"
    descriptor.write_text(
        "Structure Atlas File (edit at your own risks)\n"
        "Topology Atlas:\n"
        f"File:\t{reference}\n"
    )

    with pytest.raises(LesionAtlasDescriptorError, match=message):
        load_lesion_atlas_topology(descriptor)


def test_lesion_topology_loader_rejects_non_3d_image(tmp_path) -> None:
    atlas_path = tmp_path / "topology.nii.gz"
    nib.save(nib.Nifti1Image(np.zeros((2, 3, 4, 2), dtype=np.int16), np.eye(4)), atlas_path)
    descriptor = tmp_path / "atlas.txt"
    descriptor.write_text(
        "Structure Atlas File (edit at your own risks)\n"
        "Topology Atlas:\n"
        "File:\ttopology.nii.gz\n"
    )

    with pytest.raises(LesionAtlasDescriptorError, match="must be 3D"):
        load_lesion_atlas_topology(descriptor)


def _write_shape_descriptor(
    directory: Path, structures, shape_rows, dimensions=(1, 1, 2)
) -> Path:
    rows = [
        "Structure Atlas File (edit at your own risks)",
        f"Structures:\t{len(structures)}",
        *(f"{name}\t{index}\tobj" for index, name in enumerate(structures, start=1)),
        "Topology Atlas:",
        "File:\ttopology.nii.gz",
        "Shape Atlas:\t8",
        "Dimensions:\t" + "\t".join(str(size) for size in dimensions),
        "Resolutions:\t1\t1\t1",
        "Orientations:\t1\t1\t6\t3",
    ]
    for name, reference in shape_rows:
        rows.extend((f"Structure:\t{name}", f"File:\t{reference}"))
    rows.extend(("Intensity Atlas: 0", "T1_MPRAGE\t0.0"))
    descriptor = directory / "atlas.txt"
    descriptor.write_text("\n".join(rows) + "\n")
    return descriptor


def test_shape_loader_uses_structure_order_and_exact_float32_values(tmp_path) -> None:
    structures = ("Alpha", "Beta", "Gamma")
    values = {
        "Alpha": np.asarray([[[0.125, 2.5]]], dtype=np.float64),
        "Beta": np.asarray([[[3.25, -4.5]]], dtype=np.float64),
        "Gamma": np.asarray([[[8.0, 9.75]]], dtype=np.float64),
    }
    for name, data in values.items():
        nib.save(nib.Nifti1Image(data, np.eye(4)), tmp_path / f"{name}.nii.gz")
    descriptor = _write_shape_descriptor(
        tmp_path,
        structures,
        [("Gamma", "Gamma.nii.gz"), ("Alpha", "Alpha.nii.gz"), ("Beta", "Beta.nii.gz")],
    )

    bundle = load_lesion_atlas_shape_priors(descriptor)

    assert tuple(bundle.shape_maps) == structures
    assert tuple(bundle.files) == structures
    for name in structures:
        assert bundle.shape_maps[name].dtype == np.float32
        np.testing.assert_array_equal(bundle.shape_maps[name], values[name].astype(np.float32))
        assert not bundle.shape_maps[name].flags.writeable
    assert bundle.files["Alpha"].name == "Alpha.nii.gz"


def test_bundled_shape_atlas_loads_every_declared_map_in_structure_order() -> None:
    descriptor = ATLAS_ROOT / "cruise-atlas-12obj-lesiontoads2012.txt"

    bundle = load_lesion_atlas_shape_priors(descriptor)

    assert tuple(bundle.shape_maps) == tuple(
        structure.name for structure in load_lesion_atlas_priors(descriptor, ["T1_MPRAGE"]).structures
    )
    assert len(bundle.files) == 13
    assert bundle.files["OuterVentricular-WM"].name == "atlas_Cerebrum-WM.nii.gz"
    assert all(array.dtype == np.float32 for array in bundle.shape_maps.values())
    assert all(array.shape == (172, 216, 164) for array in bundle.shape_maps.values())


def test_shape_loader_warns_and_reorients_non_ras_map(tmp_path) -> None:
    values = np.arange(24, dtype=np.float32).reshape((2, 3, 4))
    affine = np.diag([-1.0, -1.0, 1.0, 1.0])
    affine[:2, 3] = (1.0, 2.0)
    nib.save(nib.Nifti1Image(values, affine), tmp_path / "A.nii.gz")
    descriptor = _write_shape_descriptor(
        tmp_path,
        ("A",),
        [("A", "A.nii.gz")],
        dimensions=(2, 3, 4),
    )

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        bundle = load_lesion_atlas_shape_priors(descriptor)

    expected = nib.orientations.apply_orientation(
        values,
        nib.orientations.ornt_transform(
            nib.orientations.axcodes2ornt(("L", "P", "S")),
            nib.orientations.axcodes2ornt(("R", "A", "S")),
        ),
    )
    np.testing.assert_array_equal(bundle.shape_maps["A"], expected)


def test_shape_loader_reuses_declared_cerebrum_wm_map_for_outer_alias(tmp_path) -> None:
    values = np.asarray([[[0.0, 0.25], [0.5, 1.0]]], dtype=np.float32)
    nib.save(nib.Nifti1Image(values, np.eye(4)), tmp_path / "cerebrum.nii.gz")
    descriptor = _write_shape_descriptor(
        tmp_path,
        ("Cerebrum-WM", "OuterVentricular-WM"),
        [
            ("Cerebrum-WM", "cerebrum.nii.gz"),
            ("OuterVentricular-WM", "cerebrum.nii.gz"),
        ],
        dimensions=(1, 2, 2),
    )

    bundle = load_lesion_atlas_shape_priors(descriptor)

    assert bundle.shape_maps["Cerebrum-WM"] is bundle.shape_maps["OuterVentricular-WM"]
    np.testing.assert_array_equal(bundle.shape_maps["OuterVentricular-WM"], values)
    assert bundle.files["OuterVentricular-WM"] == bundle.files["Cerebrum-WM"]


@pytest.mark.parametrize(
    ("reference", "message"),
    [
        ("missing.nii.gz", "Cannot load LesionTOADS shape atlas"),
        ("../outside.nii.gz", "Unsafe shape atlas path"),
        ("map.nii", "compressed .nii.gz"),
    ],
)
def test_shape_loader_rejects_missing_or_unsafe_paths(tmp_path, reference, message) -> None:
    descriptor = _write_shape_descriptor(
        tmp_path, ("A",), [("A", reference)]
    )
    with pytest.raises(LesionAtlasDescriptorError, match=message):
        load_lesion_atlas_shape_priors(descriptor)


def test_shape_loader_rejects_malformed_and_incomplete_file_mapping(tmp_path) -> None:
    descriptor = _write_shape_descriptor(tmp_path, ("A",), [("A", "A.nii.gz")])
    text = descriptor.read_text().replace("File:\tA.nii.gz\n", "")
    descriptor.write_text(text)
    with pytest.raises(LesionAtlasDescriptorError, match="no File reference"):
        load_lesion_atlas_shape_priors(descriptor)

    descriptor = _write_shape_descriptor(tmp_path, ("A", "B"), [("A", "A.nii.gz")])
    with pytest.raises(LesionAtlasDescriptorError, match="missing File mappings.*B"):
        load_lesion_atlas_shape_priors(descriptor)


def test_shape_loader_does_not_transform_or_sample_loaded_priors(tmp_path) -> None:
    values = np.asarray([[[0.0, 1.0]]], dtype=np.float32)
    nib.save(nib.Nifti1Image(values, np.eye(4)), tmp_path / "A.nii.gz")
    descriptor = _write_shape_descriptor(
        tmp_path, ("A",), [("A", "A.nii.gz")]
    )

    bundle = load_lesion_atlas_shape_priors(descriptor)

    # A descriptor file is loaded byte-for-voxel into its RAS grid; there is
    # no transform argument or hidden resampling in this loader stage.
    np.testing.assert_array_equal(bundle.shape_maps["A"], values)
    assert bundle.shape_maps["A"].shape == (1, 1, 2)


def test_shape_loader_rejects_shape_that_differs_from_descriptor_dimensions(tmp_path) -> None:
    values = np.zeros((1, 1, 2), dtype=np.float32)
    nib.save(nib.Nifti1Image(values, np.eye(4)), tmp_path / "A.nii.gz")
    descriptor = _write_shape_descriptor(
        tmp_path,
        ("A",),
        [("A", "A.nii.gz")],
        dimensions=(1, 1, 3),
    )

    with pytest.raises(LesionAtlasDescriptorError, match=r"descriptor declares \(1, 1, 3\)"):
        load_lesion_atlas_shape_priors(descriptor)


def test_shape_loader_rejects_non_3d_map_before_dimension_comparison(tmp_path) -> None:
    values = np.zeros((1, 1, 2, 2), dtype=np.float32)
    nib.save(nib.Nifti1Image(values, np.eye(4)), tmp_path / "A.nii.gz")
    descriptor = _write_shape_descriptor(
        tmp_path, ("A",), [("A", "A.nii.gz")]
    )

    with pytest.raises(LesionAtlasDescriptorError, match="must be 3D"):
        load_lesion_atlas_shape_priors(descriptor)
