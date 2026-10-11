from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts/compare_fine_grid_lesion_dura_java.py"


def test_iteration_six_dura_outputs_match_selected_java_exactly(tmp_path: Path) -> None:
    if not shutil.which("javac") or not shutil.which("java"):
        pytest.skip("selected Java reference requires a JDK")
    report_path = tmp_path / "dura-java-comparison.json"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--report", str(report_path)],
        check=True,
        capture_output=True,
        text=True,
        cwd=REPO,
    )
    report = json.loads(report_path.read_text())
    assert report["completed_iteration"] == 6
    assert report["synthetic_only"] is True
    assert report["java_source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert report["all_four_exact"] is True
    assert set(report["comparisons"]) == {"csf", "gm", "wm_fill", "wm_mask"}
    for comparison in report["comparisons"].values():
        assert comparison["exact_float32_match"] is True
        assert comparison["different_voxels"] == 0
