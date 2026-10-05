from __future__ import annotations

import json
import hashlib
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.cli import main
from pymacruise.lesion_dura_outputs import generate_lesion_dura_removal_outputs


ROOT = Path(__file__).parent
SYNTHETIC = ROOT / "data/synthetic/plus_lesion_atlas_init"
ATLAS = SYNTHETIC / "mini-lesion-toads-atlas.txt"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    modalities = nib.load(str(SYNTHETIC / "synthetic_t1_flair.nii.gz"))
    values = np.asarray(modalities.dataobj, dtype=np.float32)
    t1 = tmp_path / "t1.nii.gz"
    flair = tmp_path / "flair-shifted-grid.nii.gz"
    mask = tmp_path / "brain-mask.nii.gz"
    nib.save(nib.Nifti1Image(values[..., 0], modalities.affine), t1)
    shifted = modalities.affine.copy()
    shifted[0, 3] += 0.25
    nib.save(nib.Nifti1Image(values[..., 1], shifted), flair)
    topology = np.asarray(nib.load(str(SYNTHETIC / "synthetic_topology.nii.gz")).dataobj)
    nib.save(nib.Nifti1Image((topology > 0).astype(np.uint8), modalities.affine), mask)
    return t1, flair, mask


def _command(t1: Path, flair: Path, mask: Path, output: Path, *extra: str) -> list[str]:
    return [
        "plus-full-run", "--t1", str(t1), "--flair", str(flair),
        "--brain-mask", str(mask), "--atlas-descriptor", str(ATLAS),
        "--output-dir", str(output),
        "--registration-iterations", "2", "--shrink-factors", "1",
        "--smoothing-sigmas", "0", "--registration-sampling-rate", "1",
        "--n4-shrink-factor", "1", "--n4-iterations", "2",
        "--n4-spline-mesh", "1", "1", "1", "--crop-padding", "1",
        "--init-align-iterations", "0", "--pyramid-levels", "2",
        "--lesion-iterations", "1", "--first-thinning-limit", "1000",
        "--loop-convergence-threshold", "0", "--n-iter-max", "2",
        "--loop-thinning-limit", "1000", "--finite-zero-variance", *extra,
    ]


def test_plus_full_run_is_restartable_from_raw_inputs_to_terminal_outputs(tmp_path: Path) -> None:
    t1, flair, mask = _inputs(tmp_path)
    output = tmp_path / "plus-full"
    command = _command(t1, flair, mask, output)
    assert main(command) == 0

    manifest_path = output / "plus-full-workflow.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["status"] == "complete"
    assert manifest["terminal"] is True
    assert manifest["stage_order"] == [
        "preprocessing", "crop", "initial_alignment", "registered_recrop",
        "lesion_workflow", "dura_outputs",
    ]
    for record in manifest["stages"].values():
        child = Path(record["manifest"]["path"])
        assert child.is_file()
        assert record["manifest"]["sha256"]
        child_manifest = json.loads(child.read_text())
        assert child_manifest["status"] == "complete"
        for output_record in child_manifest.get("outputs", {}).values():
            path = Path(output_record["path"])
            assert path.is_file()
            assert output_record["sha256"]

    prep = json.loads(Path(manifest["stages"]["preprocessing"]["manifest"]["path"]).read_text())
    assert prep["channel_order"] == ["T1", "FLAIR"]
    loop_path = Path(manifest["stages"]["lesion_workflow"]["manifest"]["path"])
    loop = json.loads(loop_path.read_text())
    registered_loop = json.loads(Path(loop["stages"]["registered_loop"]["manifest"]["path"]).read_text())
    classification_path = Path(registered_loop["outputs"]["classification"]["path"])
    subtype_path = Path(registered_loop["outputs"]["lesion_subtype"]["path"])
    recrop = json.loads(Path(manifest["stages"]["registered_recrop"]["manifest"]["path"]).read_text())
    reference = nib.load(recrop["outputs"]["modalities"]["path"])
    for path in (classification_path, subtype_path):
        image = nib.load(str(path))
        assert path.name.endswith(".nii.gz")
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert image.shape == reference.shape[:3]
        np.testing.assert_allclose(image.affine, reference.affine, rtol=0, atol=1e-5)
    classification = np.asarray(nib.load(str(classification_path)).dataobj)
    atlas_labels = {
        int(fields[1])
        for line in ATLAS.read_text().splitlines()
        if len(fields := line.split()) >= 3 and fields[1].isdigit()
    }
    assert set(np.unique(classification)).issubset(atlas_labels)

    loop_outputs = registered_loop["outputs"]
    packed = np.asarray(nib.load(loop_outputs["memberships"]["path"]).dataobj, dtype=np.float32)
    class_indices = np.asarray(nib.load(loop_outputs["segmentation"]["path"]).dataobj)
    recrop_summary = json.loads(Path(recrop["outputs"]["summary"]["path"]).read_text())
    expected_dura = generate_lesion_dura_removal_outputs(
        packed[..., :13], class_indices, tuple(recrop_summary["class_names"]),
        packed[..., 13], packed[..., 14],
    )
    expected_channels = {
        "csf": expected_dura.csf,
        "gm": expected_dura.gray_matter,
        "wm_fill": expected_dura.white_matter_fill,
        "wm_mask": expected_dura.white_matter_mask,
    }
    dura_manifest_path = Path(manifest["stages"]["dura_outputs"]["manifest"]["path"])
    dura_manifest = json.loads(dura_manifest_path.read_text())
    assert dura_manifest["channel_order"] == ["csf", "gm", "wm_fill", "wm_mask"]
    for name, expected in expected_channels.items():
        image = nib.load(dura_manifest["outputs"][name]["path"])
        assert image.get_data_dtype() == np.dtype(np.float32)
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        np.testing.assert_allclose(image.affine, reference.affine, rtol=0, atol=1e-5)
        np.testing.assert_array_equal(np.asarray(image.dataobj), expected)

    assert main(_command(t1, flair, mask, output, "--resume")) == 0
    resumed = json.loads(manifest_path.read_text())
    assert resumed["status"] == "complete"
    assert all(stage["resumed"] for stage in resumed["stages"].values())

    prior_manifests = {
        name: record["manifest"]["sha256"]
        for name, record in resumed["stages"].items()
        if name != "lesion_workflow"
    }
    subtype_path.write_bytes(b"corrupt synthetic checkpoint output")
    assert main(_command(t1, flair, mask, output, "--resume")) == 0
    repaired_top = json.loads(manifest_path.read_text())
    repaired_loop_workflow = json.loads(
        Path(repaired_top["stages"]["lesion_workflow"]["manifest"]["path"]).read_text()
    )
    repaired_loop = json.loads(
        Path(repaired_loop_workflow["stages"]["registered_loop"]["manifest"]["path"]).read_text()
    )
    assert repaired_loop["status"] == "complete"
    repaired_subtype = Path(repaired_loop["outputs"]["lesion_subtype"]["path"])
    assert repaired_subtype.stat().st_size > 100
    assert nib.load(str(repaired_subtype)).shape == reference.shape[:3]
    assert _sha256(repaired_subtype) == repaired_loop["outputs"]["lesion_subtype"]["sha256"]
    assert all(
        repaired_top["stages"][name]["manifest"]["sha256"] == digest
        for name, digest in prior_manifests.items()
    )

    repaired_dura_manifest = Path(
        repaired_top["stages"]["dura_outputs"]["manifest"]["path"]
    )
    repaired_dura = json.loads(repaired_dura_manifest.read_text())
    prior_stage_hashes = {
        name: record["manifest"]["sha256"]
        for name, record in repaired_top["stages"].items()
        if name in {"preprocessing", "crop", "initial_alignment", "registered_recrop"}
    }
    prior_registered_loop = Path(
        json.loads(Path(repaired_top["stages"]["lesion_workflow"]["manifest"]["path"]).read_text())
        ["stages"]["registered_loop"]["manifest"]["path"]
    )
    prior_registered_loop_hash = _sha256(prior_registered_loop)
    dura_tail = Path(repaired_dura["outputs"]["wm_fill"]["path"])
    dura_tail.write_bytes(b"corrupt Dura output")
    assert main(_command(t1, flair, mask, output, "--resume")) == 0
    repaired_dura_top = json.loads(manifest_path.read_text())
    repaired_dura = json.loads(
        Path(repaired_dura_top["stages"]["dura_outputs"]["manifest"]["path"]).read_text()
    )
    repaired_dura_tail = Path(repaired_dura["outputs"]["wm_fill"]["path"])
    assert nib.load(str(repaired_dura_tail)).shape == reference.shape[:3]
    assert _sha256(repaired_dura_tail) == repaired_dura["outputs"]["wm_fill"]["sha256"]
    assert all(
        repaired_dura_top["stages"][name]["manifest"]["sha256"] == digest
        for name, digest in prior_stage_hashes.items()
    )
    repaired_loop_workflow = json.loads(
        Path(repaired_dura_top["stages"]["lesion_workflow"]["manifest"]["path"]).read_text()
    )
    assert _sha256(Path(
        repaired_loop_workflow["stages"]["registered_loop"]["manifest"]["path"]
    )) == prior_registered_loop_hash

    with pytest.raises(ValueError, match="inputs or settings differ"):
        main(_command(t1, flair, mask, output, "--resume", "--crop-padding", "2"))

    wrong_grid_mask = tmp_path / "wrong-grid-mask.nii.gz"
    t1_image = nib.load(str(t1))
    nib.save(
        nib.Nifti1Image(np.ones((10, 11, 11), dtype=np.uint8), t1_image.affine),
        wrong_grid_mask,
    )
    with pytest.raises(ValueError, match="brain mask must match the T1 reference grid"):
        main(_command(t1, flair, wrong_grid_mask, tmp_path / "rejected-input"))
