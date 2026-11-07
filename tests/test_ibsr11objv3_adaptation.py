from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import apply_orientation

from pymacruise.atlas_em_descriptor import (
    load_atlas_em_registration_inputs,
    parse_atlas_em_descriptor,
    require_atlas_em_payloads,
)


SCRIPT = Path(__file__).parents[1] / "scripts/adapt_ibsr11objv3_atlas.py"
SPEC = importlib.util.spec_from_file_location("adapt_ibsr11objv3_atlas", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _source_package(root: Path) -> tuple[Path, dict[str, np.ndarray], np.ndarray]:
    atlas = root / "ibsr11objV3"
    atlas.mkdir()
    dims = (2, 3, 4)
    topology = np.arange(np.prod(dims), dtype=np.int8).reshape(dims) - 12
    maps = {
        "Background": np.arange(np.prod(dims), dtype=np.float32).reshape(dims),
        "Thalamus": np.arange(np.prod(dims), dtype=np.float32).reshape(dims) + 100,
    }

    def write_raw(name: str, array: np.ndarray) -> None:
        # The selected Java loader indexes x + Nx*y + Nx*Ny*z.
        array.transpose(2, 1, 0).tofile(atlas / name)

    write_raw("topology.raw", topology)
    for name, values in maps.items():
        write_raw(f"{name}.raw", values.astype("<f4"))

    rows = [
        "Structure Atlas File (edit at your own risks)",
        "Structures:\t2",
        "Background\t1\toutmask",
        "Thalamus\t17\tobj",
        "",
        "Topology Atlas:",
        "File:\tibsr11objV3/topology.raw",
        "Dimensions:\t2\t3\t4",
        "Resolutions:\t1\t2\t3",
        "",
        "Shape Atlas:\t8",
        "Dimensions:\t2\t3\t4",
        "Resolutions:\t1\t2\t3",
        "Orientations:\t1\t1\t6\t3",
    ]
    for name in maps:
        rows.extend((
            f"Structure:\t{name}", f"File:\tibsr11objV3/{name}.raw",
            "Dimensions:\t2\t3\t4", "Resolutions:\t1\t2\t3",
            "Center:\t0\t0\t0",
        ))
    rows += [
        "", "Intensity Atlas:\t0", "T1\t0.0\t0.2", "FLAIR\t0.0\t0.4",
        "", "Shape Priors:", "Overlap\t1.0\t1.0",
        "", "Lesions:", "T1\t0.4\tFLAIR\t1.2",
    ]
    descriptor = atlas / "atlas.txt"
    descriptor.write_text("\n".join(rows) + "\n")
    return descriptor, maps, topology


def test_selected_raw_atlas_adapts_to_ras_nifti_and_loader_schema(tmp_path: Path) -> None:
    descriptor, source_maps, source_topology = _source_package(tmp_path)
    output = tmp_path / "converted"
    adapted_descriptor = MODULE.adapt(descriptor, output)

    source_ornt = np.asarray([[0, -1], [2, -1], [1, 1]], dtype=float)
    transform = MODULE.ornt_transform(source_ornt, MODULE.axcodes2ornt(("R", "A", "S")))
    expected_shape = apply_orientation(next(iter(source_maps.values())), transform)
    loaded = load_atlas_em_registration_inputs(parse_atlas_em_descriptor(adapted_descriptor))
    for name, values in source_maps.items():
        class_index = loaded.class_names.index(name)
        np.testing.assert_array_equal(loaded.shape_priors[..., class_index], apply_orientation(values, transform))
    assert loaded.class_names == ("Background", "Thalamus")
    assert loaded.shape_priors.shape == expected_shape.shape + (2,) == (2, 4, 3, 2)
    np.testing.assert_array_equal(loaded.topology, apply_orientation(source_topology, transform))
    assert loaded.class_labels == (1, 17)
    descriptor_meta = parse_atlas_em_descriptor(adapted_descriptor)
    require_atlas_em_payloads(descriptor_meta)
    assert descriptor_meta.shape_resolutions_mm == (1.0, 3.0, 2.0)
    assert descriptor_meta.intensity_atlas["T1_SPGR"] == (0.0, 0.2)
    assert descriptor_meta.intensity_atlas["FLAIR"] == (0.0, 0.4)

    manifest = json.loads((output / "adaptation-manifest.json").read_text())
    assert manifest["source_dimensions_xyz"] == [2, 3, 4]
    assert manifest["ras_dimensions_xyz"] == [2, 4, 3]
    assert manifest["structure_order"][1]["label"] == 17


def test_selected_raw_atlas_converter_refuses_nonempty_output_directory(tmp_path: Path) -> None:
    descriptor, _, _ = _source_package(tmp_path)
    output = tmp_path / "occupied"
    output.mkdir()
    (output / "keep.txt").write_text("keep")
    try:
        MODULE.adapt(descriptor, output)
    except FileExistsError:
        pass
    else:
        raise AssertionError("nonempty output directory was not rejected")
