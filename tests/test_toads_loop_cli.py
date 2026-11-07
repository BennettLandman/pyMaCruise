from __future__ import annotations

import json
from pathlib import Path
import shutil

import nibabel as nib
import numpy as np
import pytest

from pymacruise.cli import build_parser, main
from pymacruise.toads_dura import run_dura_removal


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "toads_registered_first_iteration_cycle"
ATLAS = ROOT / "data" / "synthetic" / "toads_registered_field_refit" / "atlas-11.txt"


def test_toads_loop_defaults_match_selected_v3_registration_schedule() -> None:
    args = build_parser().parse_args(
        [
            "toads-loop",
            "--t1-n4",
            "t1.nii.gz",
            "--brain-mask",
            "mask.nii.gz",
            "--output",
            "classification.nii.gz",
        ]
    )
    assert args.registration_levels == 4
    assert args.registration_iterations == 50
    assert args.post_refresh_registration_iterations == 5
    assert args.registration_max_iterations == 5


def test_toads_loop_cli_writes_terminal_classification_in_native_ras(tmp_path, capsys) -> None:
    output = tmp_path / "classification.nii.gz"
    dura_output_dir = tmp_path / "dura-inputs"
    result = main(
        [
            "toads-loop",
            "--t1-n4",
            str(FIXTURE / "t1.nii.gz"),
            "--brain-mask",
            str(FIXTURE / "brain-mask.nii.gz"),
            "--atlas-descriptor",
            str(ATLAS),
            "--output",
            str(output),
            "--dura-output-dir",
            str(dura_output_dir),
            "--registration-levels",
            "1",
            "--registration-iterations",
            "1",
            "--post-refresh-registration-iterations",
            "1",
            "--registration-max-iterations",
            "1",
            "--padding",
            "2",
            "--maximum-iterations",
            "2",
            "--relation-limit",
            "1",
        ]
    )
    assert result == 0
    summary = capsys.readouterr().out
    assert "iterations=2" in summary
    assert "iteration_limit_reached=true" in summary
    assert output.is_file()
    saved = nib.load(str(output))
    reference = nib.load(str(FIXTURE / "t1.nii.gz"))
    assert saved.shape == reference.shape
    assert nib.aff2axcodes(saved.affine) == ("R", "A", "S")
    assert saved.get_data_dtype().name == "uint8"
    dura_outputs = sorted(dura_output_dir.glob("*.nii.gz"))
    assert len(dura_outputs) == 5  # Classification plus the four float inputs.
    for path in dura_outputs:
        image = nib.load(str(path))
        assert image.shape == reference.shape
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        if path.name.endswith("_Toads_seg.nii.gz"):
            assert image.get_data_dtype().name == "uint8"
        else:
            assert image.get_data_dtype().name == "float32"


def test_toads_loop_source_axis_compatibility_writes_all_products_in_ras(tmp_path, capsys) -> None:
    t1 = nib.load(FIXTURE / "t1.nii.gz")
    mask = nib.load(FIXTURE / "brain-mask.nii.gz")
    flip = np.diag((-1.0, -1.0, 1.0, 1.0))
    flip[0, 3] = t1.shape[0] - 1
    flip[1, 3] = t1.shape[1] - 1
    t1_path = tmp_path / "t1-lps.nii.gz"
    mask_path = tmp_path / "mask-lps.nii.gz"
    nib.save(
        nib.Nifti1Image(np.asarray(t1.dataobj)[::-1, ::-1, :], t1.affine @ flip),
        t1_path,
    )
    nib.save(
        nib.Nifti1Image(np.asarray(mask.dataobj)[::-1, ::-1, :], mask.affine @ flip),
        mask_path,
    )
    output = tmp_path / "classification-RAS.nii.gz"
    dura_output_dir = tmp_path / "dura-RAS"

    with pytest.warns(UserWarning, match="reorienting the image to RAS"):
        result = main([
            "toads-loop", "--t1-n4", str(t1_path), "--brain-mask", str(mask_path),
            "--atlas-descriptor", str(ATLAS), "--output", str(output),
            "--dura-output-dir", str(dura_output_dir), "--registration-levels", "1",
            "--registration-iterations", "1", "--post-refresh-registration-iterations", "1",
            "--registration-max-iterations", "1", "--padding", "2",
            "--maximum-iterations", "2", "--relation-limit", "1",
            "--preserve-input-axis-order",
        ])

    assert result == 0
    assert "(RAS)" in capsys.readouterr().out
    outputs = [output, *dura_output_dir.glob("*.nii.gz")]
    assert len(outputs) == 6
    for path in outputs:
        image = nib.load(path)
        assert image.shape == t1.shape
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")

    # BASE supplies a RAS working T1. Explicit LPS processing must reproduce
    # the established path that receives an already-LPS source image.
    forced_output = tmp_path / 'forced-classification.nii.gz'
    forced_dura = tmp_path / 'forced-dura'
    with pytest.warns(UserWarning, match='reorienting the image to RAS'):
        assert main([
            'toads-loop', '--t1-n4', str(FIXTURE / 't1.nii.gz'),
            '--brain-mask', str(FIXTURE / 'brain-mask.nii.gz'),
            '--atlas-descriptor', str(ATLAS), '--output', str(forced_output),
            '--dura-output-dir', str(forced_dura), '--registration-levels', '1',
            '--registration-iterations', '1', '--post-refresh-registration-iterations', '1',
            '--registration-max-iterations', '1', '--padding', '2',
            '--maximum-iterations', '2', '--relation-limit', '1', '--source-axis-order', 'LPS',
        ]) == 0
    np.testing.assert_array_equal(nib.load(output).dataobj, nib.load(forced_output).dataobj)
    for suffix in ('_Toads_seg.nii.gz', '_csf.nii.gz', '_gm.nii.gz', '_wmfill.nii.gz', '_wmmask.nii.gz'):
        before = next(p for p in dura_output_dir.glob('*.nii.gz') if p.name.endswith(suffix))
        after = next(p for p in forced_dura.glob('*.nii.gz') if p.name.endswith(suffix))
        np.testing.assert_array_equal(nib.load(before).dataobj, nib.load(after).dataobj)
        np.testing.assert_allclose(nib.load(after).affine, t1.affine, atol=1e-6)


def test_toads_loop_cli_runs_supported_multilevel_pyramid_after_nan_trial(tmp_path, capsys) -> None:
    output = tmp_path / "classification.nii.gz"
    result = main(
        [
            "toads-loop",
            "--t1-n4",
            str(FIXTURE / "t1.nii.gz"),
            "--brain-mask",
            str(FIXTURE / "brain-mask.nii.gz"),
            "--atlas-descriptor",
            str(ATLAS),
            "--output",
            str(output),
            "--registration-levels",
            "3",
            "--registration-iterations",
            "2",
            "--registration-max-iterations",
            "2",
            "--maximum-iterations",
            "2",
        ]
    )
    assert result == 0
    summary = capsys.readouterr().out
    assert "iterations=2" in summary
    assert "energy_threshold_reached=false" in summary
    assert "iteration_limit_reached=true" in summary
    image = nib.load(str(output))
    values = np.asanyarray(image.dataobj)
    assert values.dtype == np.uint8
    assert image.shape == (13, 13, 13)
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert np.isfinite(values).all()
    assert set(np.unique(values)) == set(range(1, 12))


def test_toads_loop_help_exposes_required_native_inputs(capsys) -> None:
    try:
        main(["toads-loop", "--help"])
    except SystemExit as exc:
        assert exc.code == 0
    help_text = capsys.readouterr().out
    assert "--t1-n4" in help_text
    assert "--brain-mask" in help_text
    assert "--output" in help_text


def test_toads_wrapper_outputs_feed_the_dura_cli_and_match_api(tmp_path) -> None:
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    for asset in ATLAS.parent.glob("*.nii.gz"):
        shutil.copy2(asset, atlas_dir / asset.name)
    descriptor_text = ATLAS.read_text(encoding="utf-8")
    dura_names = {
        "Csf1": "SulcalCSF",
        "Csf2": "Ventricles",
        "Gm1": "CerebellarGM",
        "Gm2": "CerebralGM",
        "Gm3": "Caudate",
        "Gm4": "Thalamus",
        "Gm5": "Putamen",
        "Wm1": "Brainstem",
        "Wm2": "CerebellarWM",
        "Wm3": "CerebralWM",
    }
    for source_name, output_name in dura_names.items():
        descriptor_text = descriptor_text.replace(source_name, output_name)
    descriptor = atlas_dir / "atlas-dura.txt"
    descriptor.write_text(descriptor_text, encoding="utf-8")

    t1 = FIXTURE / "t1.nii.gz"
    mask = FIXTURE / "brain-mask.nii.gz"
    toads_output = tmp_path / "toads" / "classification.nii.gz"
    toads_dura_inputs = tmp_path / "toads-dura-inputs"
    assert main(
        [
            "toads-loop", "--t1-n4", str(t1), "--brain-mask", str(mask),
            "--atlas-descriptor", str(descriptor), "--output", str(toads_output),
            "--dura-output-dir", str(toads_dura_inputs),
            "--registration-levels", "1", "--registration-iterations", "1",
            "--registration-max-iterations", "1", "--padding", "2",
            "--maximum-iterations", "2", "--relation-limit", "1",
        ]
    ) == 0

    white = nib.load(str(toads_dura_inputs / "t1_wmfill.nii.gz"))
    gray = nib.load(str(toads_dura_inputs / "t1_gm.nii.gz"))
    white_mask = nib.load(str(toads_dura_inputs / "t1_wmmask.nii.gz"))
    expected = run_dura_removal(
        nib.load(str(t1)), white, gray, white_mask, max_iterations=1,
        membership_level=0.3,
    )
    dura_output_dir = tmp_path / "dura-output"
    assert main(
        [
            "dura", "--t1", str(t1), "--white", str(toads_dura_inputs / "t1_wmfill.nii.gz"),
            "--gray", str(toads_dura_inputs / "t1_gm.nii.gz"),
            "--white-mask", str(toads_dura_inputs / "t1_wmmask.nii.gz"),
            "--output-dir", str(dura_output_dir), "--max-iterations", "1",
        ]
    ) == 0
    for name in (
        "stripped_t1", "corrected_white_matter", "corrected_gray_matter",
        "brain_mask_level_set", "brain_mask", "wm_topology", "gm_topology",
    ):
        written = nib.load(str(dura_output_dir / f"{name}.nii.gz"))
        reference = getattr(expected, name)
        assert written.shape == nib.load(str(t1)).shape
        assert nib.aff2axcodes(written.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(np.asarray(written.dataobj), np.asarray(reference.dataobj))


def test_dura_stripped_t1_and_mask_feed_second_toads_cli_pass(tmp_path) -> None:
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    for asset in ATLAS.parent.glob("*.nii.gz"):
        shutil.copy2(asset, atlas_dir / asset.name)
    descriptor_text = ATLAS.read_text(encoding="utf-8")
    dura_names = {
        "Csf1": "SulcalCSF",
        "Csf2": "Ventricles",
        "Gm1": "CerebellarGM",
        "Gm2": "CerebralGM",
        "Gm3": "Caudate",
        "Gm4": "Thalamus",
        "Gm5": "Putamen",
        "Wm1": "Brainstem",
        "Wm2": "CerebellarWM",
        "Wm3": "CerebralWM",
    }
    for source_name, output_name in dura_names.items():
        descriptor_text = descriptor_text.replace(source_name, output_name)
    descriptor = atlas_dir / "atlas-dura.txt"
    descriptor.write_text(descriptor_text, encoding="utf-8")

    t1 = FIXTURE / "t1.nii.gz"
    mask = FIXTURE / "brain-mask.nii.gz"
    t1_n4 = tmp_path / "t1-n4.nii.gz"
    n4_provenance = tmp_path / "n4.json"
    assert main([
        "n4", "--input", str(t1), "--output", str(t1_n4),
        "--provenance", str(n4_provenance),
    ]) == 0
    assert t1_n4.is_file() and n4_provenance.is_file()
    first_dura_inputs = tmp_path / "first-toads-dura-inputs"
    first_pass = tmp_path / "first-pass.nii.gz"
    common_options = [
        "--atlas-descriptor", str(descriptor),
        "--registration-levels", "1",
        "--registration-iterations", "1",
        "--registration-max-iterations", "1",
        "--padding", "2",
        "--maximum-iterations", "2",
        "--relation-limit", "1",
    ]
    assert main([
        "toads-loop", "--t1-n4", str(t1_n4), "--brain-mask", str(mask),
        "--output", str(first_pass), "--dura-output-dir", str(first_dura_inputs),
        *common_options,
    ]) == 0

    dura_output_dir = tmp_path / "dura"
    assert main([
        "dura", "--t1", str(t1_n4),
        "--white", str(first_dura_inputs / "t1-n4_wmfill.nii.gz"),
        "--gray", str(first_dura_inputs / "t1-n4_gm.nii.gz"),
        "--white-mask", str(first_dura_inputs / "t1-n4_wmmask.nii.gz"),
        "--output-dir", str(dura_output_dir), "--max-iterations", "1",
    ]) == 0

    second_pass = tmp_path / "second-pass.nii.gz"
    second_dura_inputs = tmp_path / "second-toads-dura-inputs"
    assert main([
        "toads-loop",
        "--t1-n4", str(dura_output_dir / "stripped_t1.nii.gz"),
        "--brain-mask", str(dura_output_dir / "brain_mask.nii.gz"),
        "--output", str(second_pass), "--dura-output-dir", str(second_dura_inputs),
        *common_options,
    ]) == 0

    reference = nib.load(str(t1_n4))
    for path in (first_pass, second_pass, second_dura_inputs / "stripped_t1_gm.nii.gz"):
        assert path.is_file()
    second_labels = nib.load(str(second_pass))
    assert second_labels.shape == reference.shape
    assert nib.aff2axcodes(second_labels.affine) == ("R", "A", "S")
    assert second_labels.get_data_dtype().name == "uint8"
    second_gm = nib.load(str(second_dura_inputs / "stripped_t1_gm.nii.gz"))
    assert second_gm.shape == reference.shape
    assert nib.aff2axcodes(second_gm.affine) == ("R", "A", "S")
    assert second_gm.get_data_dtype().name == "float32"

    # Continue this same first-pass -> Dura -> second-pass run through the
    # selected base CRUISE core, rather than validating the second-pass seam
    # in isolation. The small fixture makes the default algorithm path usable
    # in a focused test while retaining the source crop and level-set stages.
    core_output_dir = tmp_path / "cruise-core"
    assert main([
        "toads-cruise-core",
        "--gm", str(second_dura_inputs / "stripped_t1_gm.nii.gz"),
        "--wm-fill", str(second_dura_inputs / "stripped_t1_wmfill.nii.gz"),
        "--wm-mask", str(second_dura_inputs / "stripped_t1_wmmask.nii.gz"),
        "--edge-width", "2",
        "--output-dir", str(core_output_dir),
    ]) == 0
    for filename in ("inner.nii.gz", "central.nii.gz", "outer.nii.gz"):
        image = nib.load(str(core_output_dir / filename))
        assert image.shape == reference.shape
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    for filename in ("inner.nii.gz", "central.nii.gz", "outer.nii.gz"):
        image = nib.load(str(core_output_dir / filename))
        values = np.asarray(image.dataobj)
        assert image.get_data_dtype().name == "float32"
        assert np.isfinite(values).all()
    enhanced_gm = nib.load(str(core_output_dir / "ace_enhanced_gm.nii.gz"))
    gvf = nib.load(str(core_output_dir / "gvf.nii.gz"))
    assert enhanced_gm.shape == gvf.shape[:3]
    assert gvf.shape == (*enhanced_gm.shape, 3)
    assert nib.aff2axcodes(enhanced_gm.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(gvf.affine) == ("R", "A", "S")

    # Exercise the optional MaSeg branch on the same second-pass products.
    # Place all eight selected labels on positive GM voxels so the provenance
    # proves that cleanup happened before crop selection and core evolution.
    gm_values = np.asarray(second_gm.dataobj)
    selected_voxels = np.argwhere(gm_values > 0)
    assert selected_voxels.shape[0] >= 8
    maseg_values = np.zeros(second_gm.shape, dtype=np.uint8)
    excluded_labels = (35, 38, 39, 40, 41, 71, 72, 73)
    for coordinate, label in zip(selected_voxels[:8], excluded_labels, strict=True):
        maseg_values[tuple(coordinate)] = label
    maseg_path = tmp_path / "synthetic-maseg.nii.gz"
    nib.save(nib.Nifti1Image(maseg_values, second_gm.affine), str(maseg_path))
    cleanup_output_dir = tmp_path / "cruise-core-cleanup"
    cleanup_provenance = tmp_path / "cruise-core-cleanup.json"
    assert main([
        "toads-cruise-core",
        "--gm", str(second_dura_inputs / "stripped_t1_gm.nii.gz"),
        "--wm-fill", str(second_dura_inputs / "stripped_t1_wmfill.nii.gz"),
        "--wm-mask", str(second_dura_inputs / "stripped_t1_wmmask.nii.gz"),
        "--segmentation", str(maseg_path), "--remove-cerebellum",
        "--edge-width", "2", "--output-dir", str(cleanup_output_dir),
        "--provenance", str(cleanup_provenance),
    ]) == 0
    cleanup_parameters = json.loads(cleanup_provenance.read_text())["parameters"]
    assert cleanup_parameters["remove_cerebellum"] is True
    assert cleanup_parameters["removed_voxel_counts"] == [8, 8, 8]
    assert cleanup_parameters["stage_order"][0] == "removeBrainstemAndCerebellum"
    for filename in ("inner.nii.gz", "central.nii.gz", "outer.nii.gz"):
        cleaned_levelset = nib.load(str(cleanup_output_dir / filename))
        assert cleaned_levelset.shape == reference.shape
        assert nib.aff2axcodes(cleaned_levelset.affine) == ("R", "A", "S")
        assert np.isfinite(np.asarray(cleaned_levelset.dataobj)).all()
