from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.plus_segmentation_refinement import run_plus_segmentation_refinement_checkpoint

FIXTURE = Path(__file__).parent / "data/synthetic/matlab_refinement_reference"


def test_plus_final_refinement_is_hash_checked_and_repairs_corrupt_output(tmp_path: Path) -> None:
    shape = nib.load(str(FIXTURE / "segmentation.nii.gz")).shape
    segmentation = FIXTURE / "segmentation.nii.gz"
    inner_path = FIXTURE / "inner.nii.gz"
    outer_path = FIXTURE / "outer.nii.gz"
    expected = nib.load(str(FIXTURE / "matlab-refined.nii.gz"))

    output_dir = tmp_path / "refinement"
    first = run_plus_segmentation_refinement_checkpoint(
        segmentation, inner_path, outer_path, output_dir
    )
    assert not first.resumed
    result = nib.load(str(first.output_path))
    assert result.shape == shape
    assert nib.aff2axcodes(result.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(result.dataobj, expected.dataobj)
    manifest = json.loads(first.manifest_path.read_text())
    assert manifest["status"] == "complete"
    assert set(manifest["implementation_signature"]) == {
        "plus_segmentation_refinement.py", "segmentation_refine.py", "image_io.py"
    }
    assert set(manifest["rule_metrics"]) == {
        "inner_wm", "outer_clear", "central_wm", "central_gm"
    }

    second = run_plus_segmentation_refinement_checkpoint(
        segmentation, inner_path, outer_path, output_dir, resume=True
    )
    assert second.resumed
    second.output_path.write_bytes(b"corrupt")
    repaired = run_plus_segmentation_refinement_checkpoint(
        segmentation, inner_path, outer_path, output_dir, resume=True
    )
    assert not repaired.resumed
    assert nib.load(str(repaired.output_path)).shape == shape

    shifted_affine = nib.load(str(outer_path)).affine.copy()
    shifted_affine[0, 3] += 1.0
    mismatched = tmp_path / "mismatched-outer.nii.gz"
    original_outer = nib.load(str(outer_path))
    nib.save(nib.Nifti1Image(np.asarray(original_outer.dataobj), shifted_affine), mismatched)
    failure_dir = tmp_path / "failed-refinement"
    with pytest.raises(ValueError, match="same RAS grid"):
        run_plus_segmentation_refinement_checkpoint(
            segmentation, inner_path, mismatched, failure_dir
        )
    failed = json.loads((failure_dir / "plus-segmentation-refinement.json").read_text())
    assert failed["status"] == "failed"
    assert failed["error"].startswith("ValueError:")


def test_plus_final_refinement_rejects_changed_inputs_and_implementation_on_resume(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import pymacruise.plus_segmentation_refinement as workflow

    output_dir = tmp_path / "refinement"
    first = run_plus_segmentation_refinement_checkpoint(
        FIXTURE / "segmentation.nii.gz", FIXTURE / "inner.nii.gz",
        FIXTURE / "outer.nii.gz", output_dir,
    )
    changed = nib.load(str(FIXTURE / "outer.nii.gz"))
    values = np.asarray(changed.dataobj).copy()
    values.flat[0] += 0.25
    changed_path = tmp_path / "changed-outer.nii.gz"
    nib.save(nib.Nifti1Image(values, changed.affine), changed_path)
    with pytest.raises(ValueError, match="inputs or settings differ"):
        run_plus_segmentation_refinement_checkpoint(
            FIXTURE / "segmentation.nii.gz", FIXTURE / "inner.nii.gz",
            changed_path, output_dir, resume=True,
        )

    original_hash = workflow._sha256
    implementation = Path(workflow.__file__).with_name("segmentation_refine.py").resolve()
    monkeypatch.setattr(
        workflow, "_sha256",
        lambda path: ("0" * 64 if Path(path).resolve() == implementation else original_hash(path)),
    )
    with pytest.raises(ValueError, match="inputs or settings differ"):
        run_plus_segmentation_refinement_checkpoint(
            FIXTURE / "segmentation.nii.gz", FIXTURE / "inner.nii.gz",
            FIXTURE / "outer.nii.gz", output_dir, resume=True,
        )

