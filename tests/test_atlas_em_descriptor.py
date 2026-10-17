from __future__ import annotations

from pathlib import Path

import pytest

from pymacruise.atlas_em_descriptor import (
    AtlasEMDescriptorError,
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
