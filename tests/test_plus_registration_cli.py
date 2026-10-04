from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pymacruise.cli as cli
import pymacruise.plus_workflow as plus_workflow


def test_plus_registration_cli_routes_settings_and_resume(tmp_path: Path, monkeypatch, capsys) -> None:
    captured = {}

    def fake_checkpoint(t1, flair, output_dir, *, config, resume):
        captured.update(t1=t1, flair=flair, output_dir=output_dir, config=config, resume=resume)
        return SimpleNamespace(
            resumed=resume,
            warped_flair_path=output_dir / "flair_in_t1.nii.gz",
            manifest_path=output_dir / "flair-registration.json",
        )

    monkeypatch.setattr(plus_workflow, "run_plus_flair_registration_checkpoint", fake_checkpoint)
    output_dir = tmp_path / "registration"
    assert cli.main([
        "plus-register-flair",
        "--t1", str(tmp_path / "t1.nii.gz"),
        "--flair", str(tmp_path / "flair.nii.gz"),
        "--output-dir", str(output_dir),
        "--resume",
        "--metric", "meansquares",
        "--metric-bins", "24",
        "--sampling-rate", "0.75",
        "--iterations", "70", "40", "0",
        "--shrink-factors", "4", "2", "1",
        "--smoothing-sigmas", "2", "1", "0",
    ]) == 0

    config = captured["config"]
    assert captured["resume"] is True
    assert captured["t1"] == tmp_path / "t1.nii.gz"
    assert captured["flair"] == tmp_path / "flair.nii.gz"
    assert captured["output_dir"] == output_dir
    assert config.transform_type == "Rigid"
    assert config.metric == "meansquares"
    assert config.metric_bins == 24
    assert config.sampling_rate == 0.75
    assert config.iterations == (70, 40, 0)
    assert config.shrink_factors == (4, 2, 1)
    assert config.smoothing_sigmas == (2.0, 1.0, 0.0)
    assert "Verified and resumed" in capsys.readouterr().out
