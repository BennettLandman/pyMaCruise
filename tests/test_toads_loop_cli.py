from __future__ import annotations

from pathlib import Path

import nibabel as nib

from pymacruise.cli import main


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "toads_registered_first_iteration_cycle"
ATLAS = ROOT / "data" / "synthetic" / "toads_registered_field_refit" / "atlas-11.txt"


def test_toads_loop_cli_writes_terminal_classification_in_native_ras(tmp_path) -> None:
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
            "1",
            "--registration-iterations",
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
    assert output.is_file()
    saved = nib.load(str(output))
    reference = nib.load(str(FIXTURE / "t1.nii.gz"))
    assert saved.shape == reference.shape
    assert nib.aff2axcodes(saved.affine) == ("R", "A", "S")
    assert saved.get_data_dtype().name == "uint8"


def test_toads_loop_help_exposes_required_native_inputs(capsys) -> None:
    try:
        main(["toads-loop", "--help"])
    except SystemExit as exc:
        assert exc.code == 0
    help_text = capsys.readouterr().out
    assert "--t1-n4" in help_text
    assert "--brain-mask" in help_text
    assert "--output" in help_text
