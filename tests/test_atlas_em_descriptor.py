from __future__ import annotations

from pathlib import Path

import pytest
import nibabel as nib
import numpy as np

from pymacruise.atlas_em_descriptor import (
    AtlasEMDescriptorError,
    load_atlas_em_registration_inputs,
    load_atlas_em_nifti_payloads,
    parse_atlas_em_descriptor,
    require_atlas_em_payloads,
)


def _descriptor_text(topology_file: str = "topology.img", shape_file: str = "shape.img") -> str:
    return "\n".join(
        (
            "Structure Atlas File (edit at your own risks)",
            "",
            "Structures:\t2",
            "Background\t1\tmask",
            "Cortex\t2\tclosed",
            "",
            "Topology Atlas:",
            f"File:\t{topology_file}",
            "Dimensions:\t5\t6\t7",
            "Resolutions:\t0.8\t0.9\t1.2",
            "",
            "Shape Atlas:\t4",
            "Dimensions:\t5\t6\t7",
            "Resolutions:\t0.8\t0.9\t1.2",
            "Orientations:\t1\t1\t6\t3",
            "Structure:\tBackground",
            f"File:\t{shape_file}",
            "Dimensions:\t5\t6\t7",
            "Resolutions:\t0.8\t0.9\t1.2",
            "Center:\t0\t0\t0",
            "Structure:\tCortex",
            "File:\tcortex.img",
            "Dimensions:\t5\t6\t7",
            "Resolutions:\t0.8\t0.9\t1.2",
            "Center:\t0\t0\t0",
            "",
            "Registered Shapes (0/1: true/false)",
            "0\t1",
            "0\t1",
            "Intensity Atlas:\t12",
            "T1\t0.1\t0.9",
            "Intensity Variance:",
            "T1_MPRAGE\t1.5\t2.5",
            "Lesions:",
            "T1\t0.2\tFLAIR\t0.3",
            "Black Holes:",
            "MPRAGE\t0.05",
            "Slope:\t3.5",
            "Modality weights:",
            "T1_SPGR\t0.4\t0.6",
            "Optimized factors:",
            "MPRAGE_RAW\t1.1\t1.2\t1.3\t1.4",
        )
    ) + "\n"


def test_descriptor_parser_matches_selected_class_order_metadata_and_asset_inventory(tmp_path: Path) -> None:
    path = tmp_path / "atlas.txt"
    path.write_text(_descriptor_text(), encoding="utf-8")

    parsed = parse_atlas_em_descriptor(path)

    assert [(item.name, item.label, item.topology) for item in parsed.structures] == [
        ("Background", 1, "mask"),
        ("Cortex", 2, "closed"),
    ]
    assert parsed.topology_dimensions_xyz == (5, 6, 7)
    assert parsed.topology_resolutions_mm == (0.8, 0.9, 1.2)
    assert parsed.shape_samples == 4
    assert parsed.shape_dimensions_xyz == (5, 6, 7)
    assert parsed.orientations == (1, 1, 6, 3)
    assert parsed.registered_shapes == (False, True)
    assert parsed.intensity_samples == 12
    assert parsed.intensity_atlas["T1_SPGR"] == (0.1, 0.9)
    assert parsed.intensity_variance["T1_MPRAGE"] == (1.5, 2.5)
    assert parsed.lesion_models["T1_SPGR"] == 0.2
    assert parsed.lesion_models["FLAIR"] == 0.3
    assert parsed.black_hole_models["T1_MPRAGE"] == 0.05
    assert parsed.modality_weights["T1_SPGR"] == (0.4, 0.6)
    assert parsed.optimized_factors["T1_RAW"] == (1.1, 1.2, 1.3, 1.4)
    assert parsed.shape_slope == 3.5
    assert parsed.intensity_atlas["T2"] == (0.0, 0.0)
    assert parsed.intensity_variance["T2"] == (1.0, 1.0)
    assert parsed.modality_weights["T2"] == (1.0, 1.0)
    assert parsed.optimized_factors["T2"] == (1.0, 1.0, 1.0, 1.0)
    assert [item.name for item in parsed.referenced_files] == [
        "topology.img", "shape.img", "cortex.img"
    ]
    assert len(parsed.missing_files) == 3
    with pytest.raises(FileNotFoundError, match="payloads are missing"):
        require_atlas_em_payloads(parsed)

    for item in parsed.referenced_files:
        item.touch()
    refreshed = parse_atlas_em_descriptor(path)
    assert refreshed.missing_files == ()
    require_atlas_em_payloads(refreshed)


def test_descriptor_parser_preserves_java_default_registered_shape_vector(tmp_path: Path) -> None:
    path = tmp_path / "atlas.txt"
    text = _descriptor_text().split("\nRegistered Shapes")[0]
    path.write_text(text + "\n", encoding="utf-8")

    parsed = parse_atlas_em_descriptor(path)

    assert parsed.registered_shapes == (False, True)


def test_registered_shape_values_follow_java_equals_one_semantics(tmp_path: Path) -> None:
    path = tmp_path / "atlas.txt"
    path.write_text(_descriptor_text().replace("0\t1\n0\t1", "0\t1\n2\t1"), encoding="utf-8")

    parsed = parse_atlas_em_descriptor(path)

    assert parsed.registered_shapes == (False, True)


def test_nifti_payloads_reorient_to_ras_and_keep_descriptor_geometry(tmp_path: Path) -> None:
    raw = np.arange(5 * 6 * 7, dtype=np.int16).reshape(5, 6, 7)
    topology = raw % 3
    first_shape = raw.astype(np.float32) / raw.size
    second_shape = (raw.astype(np.float32) + 2.0) / raw.size
    las_affine = np.diag((-0.8, 0.9, 1.2, 1.0))
    nib.save(nib.Nifti1Image(topology, las_affine), tmp_path / "topology.nii.gz")
    nib.save(nib.Nifti1Image(first_shape, las_affine), tmp_path / "background.nii.gz")
    nib.save(nib.Nifti1Image(second_shape, las_affine), tmp_path / "cortex.nii.gz")
    path = tmp_path / "atlas.txt"
    path.write_text(
        _descriptor_text("topology.nii.gz", "background.nii.gz").replace(
            "cortex.img", "cortex.nii.gz"
        ),
        encoding="utf-8",
    )
    descriptor = parse_atlas_em_descriptor(path)

    with pytest.warns(UserWarning, match="reorienting the image to RAS") as captured:
        loaded = load_atlas_em_nifti_payloads(descriptor)

    assert len(captured) == 3
    assert loaded.source_axis_codes == {
        "topology": ("L", "A", "S"),
        "shape:Background": ("L", "A", "S"),
        "shape:Cortex": ("L", "A", "S"),
    }
    np.testing.assert_array_equal(loaded.topology, np.flip(topology, axis=0))
    np.testing.assert_array_equal(loaded.shapes["Background"], np.flip(first_shape, axis=0))
    np.testing.assert_array_equal(loaded.shapes["Cortex"], np.flip(second_shape, axis=0))
    assert loaded.shapes["Background"].dtype == np.float32
    assert loaded.topology.shape == (5, 6, 7)
    with pytest.warns(UserWarning):
        registration_inputs = load_atlas_em_registration_inputs(descriptor)
    assert registration_inputs.class_names == ("Background", "Cortex")
    assert registration_inputs.class_labels == (1, 2)
    np.testing.assert_array_equal(registration_inputs.registered_classes, [False, True])
    np.testing.assert_array_equal(registration_inputs.shape_priors[..., 0], loaded.shapes["Background"])
    np.testing.assert_array_equal(registration_inputs.shape_priors[..., 1], loaded.shapes["Cortex"])


def test_nifti_payload_loader_rejects_shape_prior_on_another_grid(tmp_path: Path) -> None:
    data = np.ones((5, 6, 7), dtype=np.float32)
    affine = np.diag((0.8, 0.9, 1.2, 1.0))
    nib.save(nib.Nifti1Image(data.astype(np.uint8), affine), tmp_path / "topology.nii.gz")
    nib.save(nib.Nifti1Image(data, affine), tmp_path / "background.nii.gz")
    shifted = affine.copy()
    shifted[0, 3] = 1.0
    nib.save(nib.Nifti1Image(data, shifted), tmp_path / "cortex.nii.gz")
    path = tmp_path / "atlas.txt"
    path.write_text(
        _descriptor_text("topology.nii.gz", "background.nii.gz").replace(
            "cortex.img", "cortex.nii.gz"
        ),
        encoding="utf-8",
    )

    with pytest.raises(AtlasEMDescriptorError, match="grid does not match"):
        load_atlas_em_nifti_payloads(parse_atlas_em_descriptor(path))


@pytest.mark.parametrize(
    ("topology_values", "message"),
    (
        (np.full((5, 6, 7), np.nan, dtype=np.float32), "non-finite"),
        (np.full((5, 6, 7), 0.5, dtype=np.float32), "integer labels"),
    ),
)
def test_nifti_payload_loader_rejects_invalid_topology_values(
    tmp_path: Path, topology_values: np.ndarray, message: str
) -> None:
    affine = np.diag((0.8, 0.9, 1.2, 1.0))
    shape = np.ones((5, 6, 7), dtype=np.float32)
    nib.save(nib.Nifti1Image(topology_values, affine), tmp_path / "topology.nii.gz")
    nib.save(nib.Nifti1Image(shape, affine), tmp_path / "background.nii.gz")
    nib.save(nib.Nifti1Image(shape, affine), tmp_path / "cortex.nii.gz")
    path = tmp_path / "atlas.txt"
    path.write_text(
        _descriptor_text("topology.nii.gz", "background.nii.gz").replace(
            "cortex.img", "cortex.nii.gz"
        ),
        encoding="utf-8",
    )

    with pytest.raises(AtlasEMDescriptorError, match=message):
        load_atlas_em_nifti_payloads(parse_atlas_em_descriptor(path))


@pytest.mark.parametrize(
    "replacement",
    (
        ("Cortex\t2\tclosed", "Cortex\t2"),
        ("File:\ttopology.img", "File:\t../topology.img"),
        ("0\t1\n0\t1", "0\t1\n0\tinvalid"),
    ),
)
def test_descriptor_parser_rejects_malformed_or_unsafe_metadata(tmp_path: Path, replacement) -> None:
    path = tmp_path / "atlas.txt"
    path.write_text(_descriptor_text().replace(*replacement), encoding="utf-8")

    with pytest.raises(AtlasEMDescriptorError):
        parse_atlas_em_descriptor(path)
