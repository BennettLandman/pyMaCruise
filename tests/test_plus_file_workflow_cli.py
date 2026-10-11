from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pymacruise.cli as cli
import pymacruise.plus_file_workflow as plus_file_workflow


def test_plus_run_help_identifies_brain_tissue_fusion_mask(capsys) -> None:
    try:
        cli.main(["plus-run", "--help"])
    except SystemExit as exc:
        assert exc.code == 0
    else:
        raise AssertionError("argparse --help should exit")

    output = capsys.readouterr().out
    assert "Atlas_BrainTissueMask" in output
    assert "not the WM-only Atlas_Mask" in output


def test_plus_run_cli_routes_all_file_inputs_and_options(tmp_path: Path, monkeypatch, capsys) -> None:
    captured = {}

    def fake_run(**kwargs):
        captured.update(kwargs)
        return SimpleNamespace(
            manifest_path=kwargs["output_dir"] / "run.json",
            outputs={"inner_surface": kwargs["output_dir"] / "inner_surface.vtk"},
            resumed_stages=("transfer_crop",),
        )

    monkeypatch.setattr(plus_file_workflow, "run_plus_file_workflow", fake_run)
    arguments = [
        "plus-run",
        "--t1-reference", str(tmp_path / "t1.nii.gz"),
        "--segmentation", str(tmp_path / "seg.nii.gz"),
        "--dura-csf", str(tmp_path / "csf.nii.gz"),
        "--dura-gm", str(tmp_path / "gm.nii.gz"),
        "--dura-wm-fill", str(tmp_path / "wmfill.nii.gz"),
        "--dura-wm-mask", str(tmp_path / "wmmask.nii.gz"),
        "--multi-atlas-gm", str(tmp_path / "atlasgm.nii.gz"),
        "--multi-atlas-wm", str(tmp_path / "atlaswm.nii.gz"),
        "--multi-atlas-mask", str(tmp_path / "atlasmask.nii.gz"),
        "--brain-mask", str(tmp_path / "brain.nii.gz"),
        "--output-dir", str(tmp_path / "run"),
        "--resume", "--crop-edge-width", "3", "--gvf-lambda", "0.15",
        "--wm-iso-level", "0.4", "--no-refine-topology-by-maseg",
    ]
    assert cli.main(arguments) == 0
    assert captured["resume"] is True
    assert captured["crop_edge_width"] == 3
    assert captured["options"].gvf_lambda == 0.15
    assert captured["options"].wm_iso_level == 0.4
    assert captured["options"].refine_topology_by_maseg is False
    assert captured["t1_reference_path"] == Path(arguments[2])
    assert captured["dura_wm_mask_path"] == tmp_path / "wmmask.nii.gz"
    out = capsys.readouterr().out
    assert "Verified stages skipped: transfer_crop" in out
    assert "inner_surface:" in out
