from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.ace import ace_enhance
from pymacruise.cli import build_parser
from pymacruise.gvf import generalized_gvf_image
from pymacruise.surfaces import (
    extract_levelset_surface,
    read_legacy_vtk,
)
from pymacruise.tgdm import tgdm_nested_images
from pymacruise.thickness import closest_point_thickness, levelset_subtraction_thickness
from pymacruise.toads_dura import run_dura_removal


REPO = Path(__file__).resolve().parents[1]
SYNTHETIC = REPO / "tests" / "data" / "synthetic"


def _run_cli(*arguments: str) -> str:
    completed = subprocess.run(
        [sys.executable, "-m", "pymacruise.cli", *map(str, arguments)],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, (
        f"CLI failed ({completed.returncode}): {completed.args}\n"
        f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    )
    return completed.stdout


def test_standalone_commands_are_listed_but_full_toads_is_not() -> None:
    parser = build_parser()
    command_names = set(parser._subparsers._group_actions[0].choices)
    assert {
        "ace", "gvf", "tgdm", "dura", "surface-extract",
        "thickness-levelset", "thickness-surface",
    } <= command_names
    assert "toads" not in command_names
    assert "lesiontoads" not in command_names


def test_cli_standalone_commands_match_public_apis_on_committed_synthetic_data(tmp_path) -> None:
    # ACE's 41^3 folded-cortex inputs exercise more than a concentric toy; zero
    # diffusion iterations keeps this CLI contract test quick and deterministic.
    folded = SYNTHETIC / "folded_cortex"
    gm = nib.load(folded / "gray_matter.nii.gz")
    wm = nib.load(folded / "white_matter.nii.gz")
    atlas_csf = nib.load(folded / "multi_atlas_csf.nii.gz")
    ace_expected = ace_enhance(
        gm, wm, multi_atlas_csf=atlas_csf, diffusion_iterations=0
    )
    ace_paths = {
        "enhanced_gray_matter": tmp_path / "ace-enhanced.nii.gz",
        "skeleton": tmp_path / "ace-skeleton.nii.gz",
        "thinned_skeleton": tmp_path / "ace-thinned.nii.gz",
    }
    _run_cli(
        "ace", "--gray", folded / "gray_matter.nii.gz",
        "--white", folded / "white_matter.nii.gz",
        "--atlas-csf", folded / "multi_atlas_csf.nii.gz",
        "--enhanced", ace_paths["enhanced_gray_matter"],
        "--skeleton", ace_paths["skeleton"],
        "--thinned-skeleton", ace_paths["thinned_skeleton"],
        "--diffusion-iterations", "0",
    )
    for name, path in ace_paths.items():
        written = nib.load(path)
        expected = getattr(ace_expected, name)
        assert nib.aff2axcodes(written.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(np.asarray(written.dataobj), np.asarray(expected.dataobj))

    # GVF uses the committed plus-stage scalar phantom; compare every vector
    # component against a direct call with the identical solver settings.
    gvf_input_path = SYNTHETIC / "plus_gvf" / "scalar.nii.gz"
    gvf_input = nib.load(gvf_input_path)
    gvf_parameters = {
        "lambda_value": 0.8,
        "max_iterations": 200,
        "solver_tolerance": 1e-5,
    }
    gvf_expected = generalized_gvf_image(gvf_input, **gvf_parameters)
    gvf_output = tmp_path / "gvf.nii.gz"
    gvf_provenance = tmp_path / "gvf-run.json"
    _run_cli("gvf", "--input", gvf_input_path, "--output", gvf_output,
             "--lambda", "0.8", "--max-iterations", "200",
             "--solver-tolerance", "1e-5", "--provenance", gvf_provenance)
    gvf_written = nib.load(gvf_output)
    assert gvf_written.shape == (*gvf_expected.vector_field.shape[:3], 3)
    assert nib.aff2axcodes(gvf_written.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(
        np.asarray(gvf_written.dataobj), np.asarray(gvf_expected.vector_field.dataobj)
    )
    provenance = json.loads(gvf_provenance.read_text(encoding="utf-8"))
    assert provenance["schema_version"] == 1
    assert provenance["command"] == "gvf"
    assert provenance["software"]["name"] == "pymacruise"
    assert provenance["software"]["version"]
    assert provenance["spatial_convention"] == "RAS"
    assert provenance["inputs"] == {"scalar_image": str(gvf_input_path.resolve())}
    assert provenance["parameters"] == gvf_parameters
    assert provenance["outputs"] == {"vector_field": str(gvf_output.resolve())}

    # Nested TGDM compares all six output products to its public API on the
    # committed source-mapped 25^3 synthetic phantom.
    tgdm_dir = SYNTHETIC / "tgdm"
    tgdm_inputs = [nib.load(tgdm_dir / f"{name}.nii.gz") for name in (
        "initial_membership", "white_matter", "gray_matter", "ace_gray_matter"
    )]
    tgdm_skeleton_path = tgdm_dir / "ace_skeleton.nii.gz"
    tgdm_expected = tgdm_nested_images(
        *tgdm_inputs, ace_skeleton=nib.load(tgdm_skeleton_path)
    )
    tgdm_output_dir = tmp_path / "tgdm"
    _run_cli(
        "tgdm", "--initial", tgdm_dir / "initial_membership.nii.gz",
        "--white", tgdm_dir / "white_matter.nii.gz",
        "--gray", tgdm_dir / "gray_matter.nii.gz",
        "--ace-gray", tgdm_dir / "ace_gray_matter.nii.gz",
        "--ace-skeleton", tgdm_skeleton_path,
        "--output-dir", tgdm_output_dir,
    )
    for name in (
        "initial_level_set", "tissue_sign", "initial_inner", "inner", "central", "outer"
    ):
        written = nib.load(tgdm_output_dir / f"{name}.nii.gz")
        assert nib.aff2axcodes(written.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(
            np.asarray(written.dataobj), np.asarray(getattr(tgdm_expected, name).dataobj)
        )

    # Dura is exercised with the committed 32x28x24 tissue/T1 phantom and zero
    # evolution passes; compare every public output volume, including topology.
    dura_inputs = [nib.load(SYNTHETIC / f"{name}.nii.gz") for name in (
        "T1", "WM_mask", "GM_mask", "WM_mask"
    )]
    dura_expected = run_dura_removal(*dura_inputs, max_iterations=0)
    dura_output_dir = tmp_path / "dura"
    _run_cli(
        "dura", "--t1", SYNTHETIC / "T1.nii.gz",
        "--white", SYNTHETIC / "WM_mask.nii.gz",
        "--gray", SYNTHETIC / "GM_mask.nii.gz",
        "--white-mask", SYNTHETIC / "WM_mask.nii.gz",
        "--output-dir", dura_output_dir, "--max-iterations", "0",
    )
    for name in (
        "stripped_t1", "corrected_white_matter", "corrected_gray_matter",
        "brain_mask_level_set", "brain_mask", "wm_topology", "gm_topology",
    ):
        written = nib.load(dura_output_dir / f"{name}.nii.gz")
        assert nib.aff2axcodes(written.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(
            np.asarray(written.dataobj), np.asarray(getattr(dura_expected, name).dataobj)
        )


def test_surface_and_thickness_cli_match_public_apis_on_committed_phantoms(tmp_path) -> None:
    surfaces_dir = SYNTHETIC / "plus_surfaces"
    inner_path = surfaces_dir / "inner.nii.gz"
    outer_path = surfaces_dir / "outer.nii.gz"

    # Surface extraction is run in a fresh interpreter; compare RAS-mm vertices
    # and triangles against the direct public extractor on the same level set.
    direct_inner_mesh = extract_levelset_surface(nib.load(inner_path))
    direct_outer_mesh = extract_levelset_surface(nib.load(outer_path))
    inner_vtk = tmp_path / "inner.vtk"
    outer_vtk = tmp_path / "outer.vtk"
    _run_cli("surface-extract", "--level-set", inner_path, "--output", inner_vtk)
    _run_cli("surface-extract", "--level-set", outer_path, "--output", outer_vtk)
    cli_inner_mesh = read_legacy_vtk(inner_vtk)
    cli_outer_mesh = read_legacy_vtk(outer_vtk)
    assert cli_inner_mesh.coordinate_frame == "RAS"
    # Legacy VTK text is written at ten significant digits.
    np.testing.assert_allclose(cli_inner_mesh.vertices, direct_inner_mesh.vertices, rtol=0, atol=1e-8)
    np.testing.assert_array_equal(cli_inner_mesh.faces, direct_inner_mesh.faces)

    # The volume thickness CLI returns the same array and printed summary as
    # the public level-set subtraction API.
    inner_image = nib.load(inner_path)
    outer_image = nib.load(outer_path)
    thickness_expected = levelset_subtraction_thickness(inner_image, outer_image)
    thickness_path = tmp_path / "thickness.nii.gz"
    stdout = _run_cli(
        "thickness-levelset", "--inner", inner_path, "--outer", outer_path,
        "--output", thickness_path,
    )
    written_volume = nib.load(thickness_path)
    np.testing.assert_array_equal(
        np.asarray(written_volume.dataobj), np.asarray(thickness_expected.image.dataobj)
    )
    assert f"crossing_voxels={thickness_expected.crossing_voxel_count}" in stdout
    assert nib.aff2axcodes(written_volume.affine) == ("R", "A", "S")

    # Surface thickness consumes VTK files produced from the committed phantom.
    # Compare all pointwise output fields against the direct surface API.
    expected_surface_thickness = closest_point_thickness(direct_inner_mesh, direct_outer_mesh)
    measured_inner_path = tmp_path / "inner-measured.vtk"
    stdout = _run_cli(
        "thickness-surface", "--inner", inner_vtk, "--outer", outer_vtk,
        "--inner-output", measured_inner_path,
    )
    measured_inner = read_legacy_vtk(measured_inner_path)
    # The VTK writer stores RAS coordinates at ten significant digits. The
    # source-faithful JIST query converts those coordinates to Point3f float32,
    # so the direct mesh and CLI round trip can differ by a few 1e-6 mm.
    for field_name, expected_values in expected_surface_thickness.inner_surface.point_data.items():
        np.testing.assert_allclose(
            measured_inner.point_data[field_name], expected_values, rtol=0, atol=5e-6
        )
    assert f"nonzero_vertices={expected_surface_thickness.nonzero_vertex_count}" in stdout
