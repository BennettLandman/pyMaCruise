from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from pymacruise.base_workflow import BaseRunOptions
from pymacruise.cli import build_parser, main


def test_base_cli_and_api_default_to_full_selected_source_profile(monkeypatch, tmp_path) -> None:
    required = ["base-run", "--t1", "T1.nii.gz", "--brain-mask", "mask.nii.gz",
                "--maseg", "MaSeg.nii.gz", "--atlas-descriptor", "atlas.txt",
                "--output-dir", str(tmp_path)]
    expected = {
        "registration_levels": 4, "registration_iterations": 50,
        "post_refresh_registration_iterations": 5, "registration_max_iterations": 5,
        "toads_iterations": 50, "dura_iterations": 5, "ace_diffusion_iterations": 25,
        "toads_axis_order": "LPS",
        "wm_iso_level": .5, "ace_use_updated_csf_in_ace": 2,
        "n4_estimation_mask": "all",
    }
    defaults = BaseRunOptions()
    parsed = build_parser().parse_args(required)
    for key, value in expected.items():
        assert getattr(defaults, key) == getattr(parsed, key) == value
    assert defaults.effective_tgdm_iterations() == (6, 4, 7)
    seen = {}

    def capture(**kwargs):
        seen.update(kwargs)
        return Path(tmp_path) / "run.json"

    monkeypatch.setattr("pymacruise.base_workflow.run_base_workflow", capture)
    assert main(required) == 0
    assert asdict(seen["options"]) == asdict(defaults)
    assert main([*required, "--n4-estimation-mask", "legacy-otsu"]) == 0
    assert seen["options"].n4_estimation_mask == "legacy-otsu"
    assert main([*required, "--tgdm-iterations", "1", "--tgdm-inner-iterations", "2"]) == 0
    assert seen["options"].effective_tgdm_iterations() == (1, 2, 1)


def test_base_tgdm_iteration_options_keep_legacy_fallback_and_allow_independent_counts() -> None:
    legacy = BaseRunOptions(tgdm_iterations=3)
    assert legacy.effective_tgdm_iterations() == (3, 3, 3)

    selected = BaseRunOptions(
        tgdm_iterations=2,
        tgdm_initial_iterations=6,
        tgdm_inner_iterations=4,
        tgdm_central_iterations=7,
    )
    assert selected.effective_tgdm_iterations() == (6, 4, 7)

    partial = BaseRunOptions(tgdm_iterations=2, tgdm_inner_iterations=5)
    assert partial.effective_tgdm_iterations() == (2, 5, 2)
    assert BaseRunOptions(tgdm_inner_iterations=5).effective_tgdm_iterations() == (6, 5, 7)
