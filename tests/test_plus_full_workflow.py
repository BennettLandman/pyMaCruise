from __future__ import annotations

import json
import hashlib
from pathlib import Path
from types import SimpleNamespace

import nibabel as nib
import numpy as np
import pytest

from pymacruise.cli import build_parser, main
from pymacruise.image_io import is_axis_aligned_ras_grid
from pymacruise.lesion_dura_outputs import generate_lesion_dura_removal_outputs
from pymacruise.spatial_transfer import transfer_to_reference_space
from pymacruise import workflow
from pymacruise import plus_full_workflow as plus_full_workflow_module


ROOT = Path(__file__).parent
SYNTHETIC = ROOT / "data/synthetic/plus_lesion_atlas_init"
FOLDED = ROOT / "data/synthetic/plus_pre_tgdm_composition"
ATLAS = SYNTHETIC / "mini-lesion-toads-atlas.txt"
LARGER_GRID = ROOT / "data/synthetic/plus_larger_grid"
ATLAS_EM_HANDOFF = ROOT / "data/synthetic/plus_atlasem_handoff/atlasem.txt"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    modalities = nib.load(str(FOLDED / "multi_atlas_gm.nii.gz"))
    gm = np.asarray(modalities.dataobj) > 0
    wm = np.asarray(nib.load(str(FOLDED / "multi_atlas_wm.nii.gz")).dataobj) > 0
    brain = np.asarray(nib.load(str(FOLDED / "multi_atlas_brain_mask.nii.gz")).dataobj) > 0
    values = (20.0 + 70.0 * gm + 110.0 * wm).astype(np.float32)
    t1 = tmp_path / "t1.nii.gz"
    flair = tmp_path / "flair-shifted-grid.nii.gz"
    mask = tmp_path / "brain-mask.nii.gz"
    maseg = tmp_path / "maseg.nii.gz"
    nib.save(nib.Nifti1Image(values, modalities.affine), t1)
    shifted = modalities.affine.copy()
    shifted[0, 3] += 0.25
    nib.save(nib.Nifti1Image((values * 0.7 + 10).astype(np.float32), shifted), flair)
    nib.save(nib.Nifti1Image(brain.astype(np.uint8), modalities.affine), mask)
    labels = np.zeros(modalities.shape, dtype=np.uint16)
    labels[gm] = 100
    labels[wm] = 4
    nib.save(nib.Nifti1Image(labels, modalities.affine), maseg)
    return t1, flair, mask, maseg


def test_plus_full_workflow_regrids_oblique_t1_mask_and_maseg_by_default(
    tmp_path: Path, monkeypatch
) -> None:
    t1, flair, mask, maseg = _inputs(tmp_path / "inputs")
    original = nib.load(str(t1))
    affine = original.affine.copy()
    theta = np.deg2rad(8.0)
    rotation = np.array(
        [[np.cos(theta), -np.sin(theta), 0.0],
         [np.sin(theta), np.cos(theta), 0.0],
         [0.0, 0.0, 1.0]]
    )
    center = nib.affines.apply_affine(affine, (np.asarray(original.shape) - 1) / 2)
    affine[:3, :3] = rotation @ affine[:3, :3]
    affine[:3, 3] = center - affine[:3, :3] @ ((np.asarray(original.shape) - 1) / 2)
    for path in (t1, mask, maseg):
        image = nib.load(str(path))
        nib.save(nib.Nifti1Image(np.asarray(image.dataobj), affine), path)
    flair_image = nib.load(str(flair))
    flair_affine = flair_image.affine.copy()
    flair_center_index = (np.asarray(flair_image.shape) - 1) / 2
    flair_center = nib.affines.apply_affine(flair_affine, flair_center_index)
    flair_affine[:3, :3] = rotation @ flair_affine[:3, :3]
    flair_affine[:3, 3] = flair_center - flair_affine[:3, :3] @ flair_center_index
    nib.save(nib.Nifti1Image(np.asarray(flair_image.dataobj), flair_affine), flair)

    def stop_after_preprocessing_inputs(*args, **kwargs):
        raise RuntimeError("captured regridded inputs")

    monkeypatch.setattr(
        plus_full_workflow_module,
        "run_plus_lesion_preprocessing_checkpoint",
        stop_after_preprocessing_inputs,
    )
    no_regrid_output = tmp_path / "no-regrid-plus"
    with pytest.raises(RuntimeError, match="captured regridded inputs"):
        plus_full_workflow_module.run_plus_full_workflow(
            t1, flair, mask, maseg, ATLAS, no_regrid_output,
            regrid_inputs=False,
        )
    no_regrid_manifest = json.loads(
        (no_regrid_output / "plus-full-workflow.json").read_text()
    )
    assert no_regrid_manifest["settings"]["regrid_inputs"] is False
    assert no_regrid_manifest["derived_resources"]["regridded_inputs"] == []
    oriented = no_regrid_manifest["derived_resources"]["orientation_only_inputs"]
    assert {item["role"] for item in oriented} == {"t1", "flair", "brain_mask", "maseg"}
    oriented_images = {item["role"]: nib.load(item["path"]) for item in oriented}
    assert all(is_axis_aligned_ras_grid(image) for image in oriented_images.values())
    assert len({(oriented_images[role].shape, tuple(oriented_images[role].affine.ravel()))
                for role in ("t1", "brain_mask", "maseg")}) == 1
    assert oriented_images["flair"].shape == nib.load(str(flair)).shape

    output = tmp_path / "regrid-plus"
    with pytest.raises(RuntimeError, match="captured regridded inputs"):
        plus_full_workflow_module.run_plus_full_workflow(
            t1, flair, mask, maseg, ATLAS, output,
        )

    manifest = json.loads((output / "plus-full-workflow.json").read_text())
    resources = manifest["derived_resources"]["regridded_inputs"]
    assert {item["role"] for item in resources} == {"t1", "brain_mask", "maseg"}
    assert all(_sha256(Path(item["path"])) == item["sha256"] for item in resources)
    derived = {item["role"]: Path(item["path"]) for item in resources}
    images = {role: nib.load(str(path)) for role, path in derived.items()}
    assert all(nib.aff2axcodes(image.affine) == ("R", "A", "S") for image in images.values())
    assert len({(image.shape, tuple(image.affine.ravel())) for image in images.values()}) == 1
    assert set(np.unique(np.asarray(images["brain_mask"].dataobj))) <= {0, 1}
    assert set(np.unique(np.asarray(images["maseg"].dataobj))).issubset(
        set(np.unique(np.asarray(nib.load(str(maseg)).dataobj)))
    )


def test_plus_full_workflow_regrid_cli_and_config_defaults() -> None:
    parser = build_parser()
    required = [
        "plus-full-run", "--t1", "t1.nii.gz", "--flair", "flair.nii.gz",
        "--brain-mask", "mask.nii.gz", "--maseg", "maseg.nii.gz",
        "--atlas-descriptor", str(ATLAS), "--output-dir", "out",
    ]
    assert parser.parse_args(required).regrid_inputs is True
    assert parser.parse_args(required).dura_transfer_mode == "physical-affine"
    assert parser.parse_args(required).atlasem_descriptor is None
    assert parser.parse_args([*required, "--atlasem-descriptor", "atlasem.txt"]).atlasem_descriptor == Path("atlasem.txt")
    assert parser.parse_args([*required, "--no-regrid-inputs"]).regrid_inputs is False
    assert parser.parse_args([*required, "--dura-transfer-mode", "jist-unchanged"]).dura_transfer_mode == "jist-unchanged"
    assert workflow._validate_run_options(
        workflow.WorkflowProfile.PLUS, {"plus": {"regrid_inputs": False}}
    )["plus"]["regrid_inputs"] is False
    assert workflow._validate_run_options(
        workflow.WorkflowProfile.PLUS, {"plus": {"dura_transfer_mode": "jist-unchanged"}}
    )["plus"]["dura_transfer_mode"] == "jist-unchanged"
    with pytest.raises(ValueError, match="dura_transfer_mode"):
        workflow._validate_run_options(
            workflow.WorkflowProfile.PLUS, {"plus": {"dura_transfer_mode": "nearestish"}}
        )


def _command(t1: Path, flair: Path, mask: Path, maseg: Path, output: Path, *extra: str) -> list[str]:
    return [
        "plus-full-run", "--t1", str(t1), "--flair", str(flair),
        "--brain-mask", str(mask), "--maseg", str(maseg), "--atlas-descriptor", str(ATLAS),
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


def test_plus_full_workflow_routes_brain_tissue_mask_to_surface_stage(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    t1, flair, mask, maseg = _inputs(tmp_path / "inputs")

    def artifact(name: str) -> Path:
        path = tmp_path / "stub-artifacts" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name, encoding="utf-8")
        return path

    def result(name: str, **attributes: object) -> SimpleNamespace:
        attributes.setdefault("manifest_path", artifact(f"{name}.json"))
        return SimpleNamespace(**attributes)

    monkeypatch.setattr(plus_full_workflow_module, "_verify_stage", lambda *_args: {})
    monkeypatch.setattr(
        plus_full_workflow_module,
        "run_plus_lesion_preprocessing_checkpoint",
        lambda *_args, **_kwargs: result("preprocessing", modalities_path=artifact("modalities.nii.gz")),
    )
    monkeypatch.setattr(
        plus_full_workflow_module,
        "run_plus_lesion_crop_checkpoint",
        lambda *_args, **_kwargs: result(
            "crop", normalized_path=artifact("normalized.nii.gz"),
            cropped_path=artifact("cropped.nii.gz"),
        ),
    )
    monkeypatch.setattr(
        plus_full_workflow_module,
        "run_plus_lesion_initial_atlas_alignment_checkpoint",
        lambda *_args, **_kwargs: result("alignment"),
    )
    monkeypatch.setattr(
        plus_full_workflow_module,
        "run_plus_lesion_registered_recrop_checkpoint",
        lambda *_args, **_kwargs: result("recrop"),
    )
    monkeypatch.setattr(
        plus_full_workflow_module,
        "run_plus_lesion_loop_workflow",
        lambda *_args, **_kwargs: result(
            "lesion-loop", registered_loop_manifest=artifact("registered-loop.json"),
            classification_path=None, lesion_subtype_path=None,
        ),
    )
    monkeypatch.setattr(
        plus_full_workflow_module,
        "run_plus_lesion_dura_output_checkpoint",
        lambda *_args, **_kwargs: result("dura", outputs={"gm": artifact("dura-gm.nii.gz")}),
    )
    native_outputs = {
        name: artifact(f"native-{name}.nii.gz")
        for name in ("csf", "gm", "wm_fill", "wm_mask")
    }
    captured_transfer: dict[str, object] = {}

    def capture_native_transfer(*args: object, **kwargs: object) -> SimpleNamespace:
        captured_transfer.update(kwargs)
        return result("native-dura", outputs=native_outputs)

    monkeypatch.setattr(
        plus_full_workflow_module,
        "run_plus_lesion_native_dura_transfer",
        capture_native_transfer,
    )
    cortical_outputs = {
        name: artifact(f"{name}.nii.gz")
        for name in ("Atlas_GM", "Atlas_WM", "Atlas_Mask", "Atlas_BrainTissueMask")
    }
    monkeypatch.setattr(
        plus_full_workflow_module,
        "run_plus_maseg_cortical_checkpoint",
        lambda *_args, **_kwargs: result("cortical", outputs=cortical_outputs),
    )
    captured: dict[str, object] = {}

    def capture_plus_surface(**kwargs: object) -> SimpleNamespace:
        captured.update(kwargs)
        return result("plus-surface", outputs={})

    monkeypatch.setattr(plus_full_workflow_module, "run_plus_file_workflow", capture_plus_surface)

    plus_full_workflow_module.run_plus_full_workflow(
        t1, flair, mask, maseg, ATLAS, tmp_path / "full-workflow",
        dura_transfer_mode="jist-unchanged",
    )

    assert captured["multi_atlas_mask_path"] == cortical_outputs["Atlas_BrainTissueMask"]
    assert captured["multi_atlas_mask_path"] != cortical_outputs["Atlas_Mask"]
    assert captured_transfer["transfer_mode"] == "jist-unchanged"
    assert captured["dura_transfer_mode"] == "jist-unchanged"


def test_plus_full_workflow_inserts_atlasem_dura_before_lesion_crop(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    t1, flair, mask, maseg = _inputs(tmp_path / "inputs")
    atlasem_descriptor = tmp_path / "atlasem-descriptor.txt"
    atlasem_descriptor.write_text("mocked synthetic descriptor\n", encoding="utf-8")

    def artifact(name: str) -> Path:
        path = tmp_path / "atlasem-stage-stubs" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name, encoding="utf-8")
        return path

    preprocessing_manifest = artifact("preprocessing.json")
    atlasem_manifest = artifact("atlasem-dura.json")
    prepared_t1 = artifact("t1-brain.nii.gz")
    prepared_flair = artifact("flair-brain.nii.gz")
    atlasem_modalities = artifact("dura-conditioned-modalities.nii.gz")
    captured: dict[str, Path] = {}
    monkeypatch.setattr(plus_full_workflow_module, "_verify_stage", lambda *_args: {"outputs": {}})
    monkeypatch.setattr(
        plus_full_workflow_module,
        "run_plus_lesion_preprocessing_checkpoint",
        lambda *_args, **_kwargs: SimpleNamespace(
            manifest_path=preprocessing_manifest,
            t1_brain_path=prepared_t1,
            flair_brain_path=prepared_flair,
            modalities_path=artifact("original-modalities.nii.gz"),
            resumed=False,
        ),
    )

    def atlasem_stage(t1_path, flair_path, mask_path, descriptor_path, output_dir, **kwargs):
        assert t1_path == prepared_t1
        assert flair_path == prepared_flair
        assert mask_path == mask.resolve()
        assert descriptor_path == atlasem_descriptor.resolve()
        return SimpleNamespace(
            manifest_path=atlasem_manifest,
            modalities_path=atlasem_modalities,
            resumed=False,
        )

    monkeypatch.setattr(
        plus_full_workflow_module,
        "run_plus_atlas_em_dura_handoff_checkpoint",
        atlasem_stage,
    )

    def stop_at_crop(modalities_path, *_args, **_kwargs):
        captured["modalities"] = Path(modalities_path)
        raise RuntimeError("crop reached after AtlasEM/Dura")

    monkeypatch.setattr(
        plus_full_workflow_module, "run_plus_lesion_crop_checkpoint", stop_at_crop
    )
    output = tmp_path / "plus-with-atlasem"
    with pytest.raises(RuntimeError, match="crop reached after AtlasEM/Dura"):
        plus_full_workflow_module.run_plus_full_workflow(
            t1, flair, mask, maseg, ATLAS, output,
            atlasem_descriptor_path=atlasem_descriptor,
            regrid_inputs=False,
        )
    assert captured["modalities"] == atlasem_modalities
    manifest = json.loads((output / "plus-full-workflow.json").read_text())
    assert manifest["stage_order"][:3] == ["preprocessing", "atlasem_dura", "crop"]
    assert manifest["stages"]["preprocessing"]["status"] == "complete"
    assert manifest["stages"]["atlasem_dura"]["status"] == "complete"


def test_plus_full_run_is_restartable_from_raw_inputs_to_terminal_outputs(tmp_path: Path) -> None:
    t1, flair, mask, maseg = _inputs(tmp_path)
    output = tmp_path / "plus-full"
    command = _command(t1, flair, mask, maseg, output)
    assert main(command) == 0

    manifest_path = output / "plus-full-workflow.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["status"] == "complete"
    assert manifest["terminal"] is True
    assert manifest["stage_order"] == [
        "preprocessing", "crop", "initial_alignment", "registered_recrop",
        "lesion_workflow", "dura_outputs", "native_dura_transfer",
        "maseg_cortical_separation", "plus_surface_workflow",
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
        "csf": np.multiply(expected_dura.csf, np.float32(255.0), dtype=np.float32),
        "gm": np.multiply(expected_dura.gray_matter, np.float32(255.0), dtype=np.float32),
        "wm_fill": np.multiply(expected_dura.white_matter_fill, np.float32(255.0), dtype=np.float32),
        "wm_mask": np.multiply(expected_dura.white_matter_mask, np.float32(255.0), dtype=np.float32),
    }
    dura_manifest_path = Path(manifest["stages"]["dura_outputs"]["manifest"]["path"])
    dura_manifest = json.loads(dura_manifest_path.read_text())
    assert dura_manifest["channel_order"] == ["csf", "gm", "wm_fill", "wm_mask"]
    assert dura_manifest["wrapper_output_scale_factor"] == 255.0
    assert dura_manifest["source_order"][-1] == "MedicAlgorithmLesionToads.scaleFloatValue(255)"
    for name, expected in expected_channels.items():
        image = nib.load(dura_manifest["outputs"][name]["path"])
        assert image.get_data_dtype() == np.dtype(np.float32)
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        np.testing.assert_allclose(image.affine, reference.affine, rtol=0, atol=1e-5)
        np.testing.assert_array_equal(np.asarray(image.dataobj), expected)

    assert main([*command, "--resume"]) == 0
    resumed = json.loads(manifest_path.read_text())
    assert resumed["status"] == "complete"
    assert all(stage["resumed"] for stage in resumed["stages"].values())

    prior_manifests = {
        name: record["manifest"]["sha256"]
        for name, record in resumed["stages"].items()
        if name != "lesion_workflow"
    }
    subtype_path.write_bytes(b"corrupt synthetic checkpoint output")
    assert main(_command(t1, flair, mask, maseg, output, "--resume")) == 0
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
    assert main(_command(t1, flair, mask, maseg, output, "--resume")) == 0
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

    native_manifest = json.loads(
        Path(repaired_dura_top["stages"]["native_dura_transfer"]["manifest"]["path"]).read_text()
    )
    assert native_manifest["source_order"] == [
        "transfer CSF auxiliary", "transfer GM", "transfer WM-fill", "transfer WM-mask"
    ]
    native_reference = nib.load(str(t1))
    for name in ("gm", "wm_fill", "wm_mask"):
        source = nib.load(repaired_dura["outputs"][name]["path"])
        expected_native = transfer_to_reference_space(source, native_reference, cval=0.0)
        actual_native_path = Path(native_manifest["outputs"][name]["path"])
        actual_native = nib.load(str(actual_native_path))
        assert actual_native.shape == native_reference.shape
        assert actual_native.get_data_dtype() == np.dtype(np.float32)
        assert nib.aff2axcodes(actual_native.affine) == ("R", "A", "S")
        np.testing.assert_allclose(actual_native.affine, native_reference.affine, rtol=0, atol=1e-5)
        np.testing.assert_array_equal(
            np.asarray(actual_native.dataobj), np.asarray(expected_native.dataobj)
        )

    expected_surface_products = {
        "inner_level_set", "central_level_set", "outer_level_set", "thickness_volume",
        "inner_surface", "central_surface", "outer_surface", "thickness_metrics",
    }
    top_outputs = repaired_dura_top["final_outputs"]
    assert {
        name.removeprefix("plus_")
        for name in top_outputs
        if name.startswith("plus_")
    } == expected_surface_products
    for name in expected_surface_products - {"thickness_metrics"}:
        record = top_outputs[f"plus_{name}"]
        artifact = Path(record["path"])
        assert _sha256(artifact) == record["sha256"]
        if artifact.suffix == ".gz":
            image = nib.load(str(artifact))
            assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
            assert np.isfinite(np.asarray(image.dataobj)).all()
        else:
            vtk = artifact.read_text()
            assert "DATASET POLYDATA" in vtk and "POINTS " in vtk
    metrics_path = Path(top_outputs["plus_thickness_metrics"]["path"])
    metrics = json.loads(metrics_path.read_text())
    assert metrics["spatial_convention"] == "RAS"
    assert metrics["surface_units"] == "mm"
    assert np.isfinite(metrics["surface"]["mean_mm"])
    # The source-correct Atlas_BrainTissueMask now preserves tissue through
    # fusion, so this workflow fixture has actual inner/outer sign crossings.
    # Validate that the reported volume statistics are present and finite;
    # anatomical ribbon-quality acceptance is measured by the dedicated audit.
    assert metrics["volume"]["crossing_voxel_count"] > 0
    assert metrics["volume"]["available"] is True
    assert np.isfinite(metrics["volume"]["mean_mm"])

    native_tail = Path(native_manifest["outputs"]["gm"]["path"])
    native_tail.write_bytes(b"corrupt native transfer")
    prior_dura_manifest_hash = _sha256(repaired_dura_manifest)
    assert main(_command(t1, flair, mask, maseg, output, "--resume")) == 0
    repaired_native_top = json.loads(manifest_path.read_text())
    repaired_native_manifest_path = Path(
        repaired_native_top["stages"]["native_dura_transfer"]["manifest"]["path"]
    )
    repaired_native_manifest = json.loads(repaired_native_manifest_path.read_text())
    repaired_native_tail = Path(repaired_native_manifest["outputs"]["gm"]["path"])
    assert nib.load(str(repaired_native_tail)).shape == native_reference.shape
    assert _sha256(repaired_native_tail) == repaired_native_manifest["outputs"]["gm"]["sha256"]
    assert _sha256(repaired_dura_manifest) == prior_dura_manifest_hash

    with pytest.raises(ValueError, match="inputs or settings differ"):
        main(_command(t1, flair, mask, maseg, output, "--resume", "--crop-padding", "2"))

    wrong_grid_mask = tmp_path / "wrong-grid-mask.nii.gz"
    t1_image = nib.load(str(t1))
    nib.save(
        nib.Nifti1Image(np.ones((10, 11, 11), dtype=np.uint8), t1_image.affine),
        wrong_grid_mask,
    )
    with pytest.raises(ValueError, match="brain mask must match the T1 reference grid"):
        main(_command(t1, flair, wrong_grid_mask, maseg, tmp_path / "rejected-input"))


def test_plus_full_run_with_atlasem_handoff_completes_and_resumes_all_ten_stages(
    tmp_path: Path,
) -> None:
    t1, flair, mask, maseg = _inputs(tmp_path / "inputs")
    output = tmp_path / "plus-full-atlasem"
    command = [
        *_command(t1, flair, mask, maseg, output),
        "--atlasem-descriptor", str(ATLAS_EM_HANDOFF),
    ]
    assert main(command) == 0

    manifest_path = output / "plus-full-workflow.json"
    manifest = json.loads(manifest_path.read_text())
    expected_order = [
        "preprocessing", "atlasem_dura", "crop", "initial_alignment",
        "registered_recrop", "lesion_workflow", "dura_outputs",
        "native_dura_transfer", "maseg_cortical_separation",
        "plus_surface_workflow",
    ]
    assert manifest["status"] == "complete"
    assert manifest["terminal"] is True
    assert manifest["stage_order"] == expected_order
    assert manifest["stages"]["atlasem_dura"]["status"] == "complete"
    assert manifest["stages"]["atlasem_dura"]["manifest"]["sha256"]
    for stage in expected_order:
        child_path = Path(manifest["stages"][stage]["manifest"]["path"])
        assert child_path.is_file()
        child = json.loads(child_path.read_text())
        assert child["status"] == "complete"

    def output_nifti_paths(value: object) -> set[Path]:
        found: set[Path] = set()
        if isinstance(value, dict):
            path = value.get("path")
            if isinstance(path, str) and path.lower().endswith(".nii.gz"):
                found.add(Path(path))
            for nested in value.values():
                found.update(output_nifti_paths(nested))
        elif isinstance(value, list):
            for nested in value:
                found.update(output_nifti_paths(nested))
        return found

    products: set[Path] = set()
    for stage in expected_order:
        child = json.loads(
            Path(manifest["stages"][stage]["manifest"]["path"]).read_text()
        )
        products.update(output_nifti_paths(child.get("outputs", {})))
    assert products
    for path in products:
        assert path.is_file()
        image = nib.load(str(path))
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S"), path
        assert np.isfinite(np.asarray(image.dataobj)).all(), path

    assert main([*command, "--resume"]) == 0
    resumed = json.loads(manifest_path.read_text())
    assert resumed["status"] == "complete" and resumed["terminal"] is True
    assert all(stage["resumed"] for stage in resumed["stages"].values())

def test_plus_full_run_completes_with_jist_unchanged_dura_transfer(tmp_path: Path) -> None:
    t1, flair, mask, maseg = _inputs(tmp_path / "inputs")
    output = tmp_path / "plus-full-jist-unchanged"
    command = _command(
        t1, flair, mask, maseg, output,
        "--dura-transfer-mode", "jist-unchanged",
    )
    assert main(command) == 0

    manifest = json.loads((output / "plus-full-workflow.json").read_text())
    assert manifest["status"] == "complete"
    assert manifest["terminal"] is True
    assert manifest["settings"]["dura_transfer_mode"] == "jist-unchanged"
    native_manifest = json.loads(
        Path(manifest["stages"]["native_dura_transfer"]["manifest"]["path"]).read_text()
    )
    dura_manifest = json.loads(
        Path(manifest["stages"]["dura_outputs"]["manifest"]["path"]).read_text()
    )
    native_gm = nib.load(native_manifest["outputs"]["gm"]["path"])
    dura_gm = nib.load(dura_manifest["outputs"]["gm"]["path"])
    assert native_manifest["transfer_mode"] == "jist-unchanged"
    assert native_gm.shape == dura_gm.shape
    assert native_gm.shape != nib.load(str(t1)).shape
    np.testing.assert_array_equal(native_gm.affine, nib.load(str(t1)).affine)
    assert nib.aff2axcodes(native_gm.affine) == ("R", "A", "S")

    surface_manifest = json.loads(
        Path(manifest["stages"]["plus_surface_workflow"]["manifest"]["path"]).read_text()
    )
    crop_record = surface_manifest["stages"]["transfer_crop"]
    cropped_gm_path = Path(crop_record["outputs"]["toads_gm"]["path"])
    crop_json = json.loads((cropped_gm_path.parent / "crop.json").read_text())
    from pymacruise.cropping import CropParameters, crop_image_with_reference_bounds

    crop = CropParameters(
        original_shape=tuple(crop_json["original_shape"]),
        bounds=tuple(tuple(pair) for pair in crop_json["bounds"]),
        background=float(crop_json["background"]),
        original_affine=np.asarray(crop_json["original_affine"], dtype=np.float64),
    )
    expected_crop = crop_image_with_reference_bounds(native_gm, crop)
    np.testing.assert_array_equal(nib.load(str(cropped_gm_path)).dataobj, expected_crop.dataobj)
    assert surface_manifest["status"] == "complete"
    assert all(Path(record["path"]).is_file() for record in manifest["final_outputs"].values())


def test_plus_full_run_completes_on_larger_asymmetric_folded_ras_grid(tmp_path: Path) -> None:
    shape = (36, 45, 34)
    t1, flair, mask, maseg = tuple(
        LARGER_GRID / name
        for name in ("t1.nii.gz", "flair.nii.gz", "brain-mask.nii.gz", "maseg.nii.gz")
    )
    fixture_record = json.loads((LARGER_GRID / "provenance.json").read_text())
    for path in (t1, flair, mask, maseg):
        assert path.is_file()
        record = fixture_record["images"][path.name]
        assert _sha256(path) == record["sha256"]
        image = nib.load(str(path))
        assert list(image.shape) == record["shape"] == list(shape)
        assert nib.aff2axcodes(image.affine) == tuple(record["axis_codes"]) == ("R", "A", "S")
    output = tmp_path / "plus-larger-grid"
    command = _command(t1, flair, mask, maseg, output)
    # Exercise the complete-run defaults: a 99-pass ceiling and 1e-4
    # convergence tolerance, rather than the helper's bounded two-pass setup.
    threshold_index = command.index("--loop-convergence-threshold")
    del command[threshold_index:threshold_index + 2]
    command[command.index("--n-iter-max") + 1] = "50"
    assert main(command) == 0

    manifest_path = output / "plus-full-workflow.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["status"] == "complete"
    assert manifest["terminal"] is True
    assert manifest["settings"]["max_passes"] == 99
    assert manifest["settings"]["max_distance"] == pytest.approx(1e-4)
    lesion_workflow = json.loads(
        Path(manifest["stages"]["lesion_workflow"]["manifest"]["path"]).read_text()
    )
    registered_loop = json.loads(
        Path(lesion_workflow["stages"]["registered_loop"]["manifest"]["path"]).read_text()
    )
    loop_summary = json.loads(
        Path(registered_loop["outputs"]["summary"]["path"]).read_text()
    )
    assert loop_summary["stop_reason"] == "convergence"
    assert loop_summary["completed_iteration"] == 5
    assert loop_summary["energy_ratio"] < 1e-4
    assert manifest["stage_order"][-1] == "plus_surface_workflow"
    for stage in manifest["stages"].values():
        child = Path(stage["manifest"]["path"])
        assert _sha256(child) == stage["manifest"]["sha256"]
        child_data = json.loads(child.read_text())
        for record in child_data.get("outputs", {}).values():
            path = Path(record["path"])
            assert _sha256(path) == record["sha256"]
            if path.name.endswith(".nii.gz"):
                image = nib.load(str(path))
                assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
                assert np.isfinite(np.asarray(image.dataobj)).all()

    recrop = json.loads(Path(manifest["stages"]["registered_recrop"]["manifest"]["path"]).read_text())
    recrop_image = nib.load(recrop["outputs"]["modalities"]["path"])
    assert nib.aff2axcodes(recrop_image.affine) == ("R", "A", "S")
    # LesionTOADS intentionally recrops to registered topology support.
    assert all(0 < dim <= full_dim for dim, full_dim in zip(recrop_image.shape[:3], shape, strict=True))
    final_outputs = manifest["final_outputs"]
    for name in (
        "native_gm", "native_wm_fill", "native_wm_mask", "plus_inner_level_set",
        "plus_central_level_set", "plus_outer_level_set", "plus_thickness_volume",
    ):
        image = nib.load(final_outputs[name]["path"])
        assert image.shape == shape
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert main([*command, "--resume"]) == 0
    resumed = json.loads(manifest_path.read_text())
    assert resumed["status"] == "complete"
    assert all(stage["resumed"] for stage in resumed["stages"].values())


def test_workflow_run_config_executes_and_resumes_folded_synthetic_plus(
    tmp_path: Path, monkeypatch
) -> None:
    t1, flair, mask, maseg = _inputs(tmp_path / "inputs")
    output = tmp_path / "configured-plus"
    config_path = tmp_path / "workflow.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "profile": "plus",
                "t1_path": str(t1),
                "flair_path": str(flair),
                "brain_mask_path": str(mask),
                "maseg_path": str(maseg),
                "output_dir": str(output),
                "plus_options": {"boundary_overlook_labels": [120, 121]},
                "run_options": {
                    "registration": {
                        "iterations": [2],
                        "shrink_factors": [1],
                        "smoothing_sigmas": [0],
                        "sampling_rate": 1,
                    },
                    "n4": {
                        "shrink_factor": 1,
                        "iterations": [2],
                        "spline_mesh": [1, 1, 1],
                    },
                    "plus": {
                        "crop_padding": 1,
                        "init_align_iterations": 0,
                        "pyramid_levels": 2,
                        "lesion_iterations": 1,
                        "first_thinning_limit": 1000,
                        "source_zero_variance_nan": False,
                        "loop_convergence_threshold": 0,
                        "n_iter_max": 2,
                        "loop_thinning_limit": 1000,
                    },
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(workflow, "get_bundled_atlas_descriptor", lambda _name: ATLAS)
    monkeypatch.setattr(workflow, "_ants_available", lambda: True)

    assert main(["workflow-preflight", "--config", str(config_path)]) == 0
    assert (output / "run.json").is_file()
    assert main(["workflow-run", "--config", str(config_path)]) == 0
    assert (output / "preflight.json").is_file()
    run_manifest = output / "plus-full-workflow.json"
    manifest = json.loads(run_manifest.read_text(encoding="utf-8"))
    assert manifest["status"] == "complete"
    assert manifest["terminal"] is True
    assert manifest["stage_order"][-1] == "plus_surface_workflow"
    assert manifest["settings"]["pipeline_options"]["boundary_overlook_labels"] == [120, 121]
    for record in manifest["stages"].values():
        child_manifest_path = Path(record["manifest"]["path"])
        child_manifest = json.loads(child_manifest_path.read_text(encoding="utf-8"))
        for output_record in child_manifest.get("outputs", {}).values():
            artifact = Path(output_record["path"])
            assert artifact.is_file()
            assert _sha256(artifact) == output_record["sha256"]
            if artifact.name.endswith(".nii.gz"):
                assert nib.aff2axcodes(nib.load(str(artifact)).affine) == ("R", "A", "S")

    assert main(["workflow-run", "--config", str(config_path), "--resume"]) == 0
    resumed = json.loads(run_manifest.read_text(encoding="utf-8"))
    assert resumed["status"] == "complete"
    assert resumed["terminal"] is True


def test_plus_full_run_uses_new_preserved_generation_after_code_change(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    t1, flair, mask, maseg = _inputs(tmp_path / "inputs")
    output = tmp_path / "plus-full-code-change"
    command = _command(t1, flair, mask, maseg, output)
    assert main(command) == 0

    manifest_path = output / "plus-full-workflow.json"
    before = json.loads(manifest_path.read_text())
    old_root = Path(before["stage_root"])
    old_manifests = {
        name: Path(record["manifest"]["path"])
        for name, record in before["stages"].items()
    }
    assert before["implementation_signature"]
    assert old_root == output
    assert all(path.is_file() for path in old_manifests.values())

    changed_signature = "a" * 64
    monkeypatch.setattr(
        plus_full_workflow_module, "_implementation_signature", lambda: changed_signature
    )
    assert main([*command, "--resume"]) == 0

    updated = json.loads(manifest_path.read_text())
    new_root = output / "code-revisions" / changed_signature[:16]
    assert updated["status"] == "complete"
    assert updated["previous_implementation_signature"] == before["implementation_signature"]
    assert updated["implementation_signature"] == changed_signature
    assert Path(updated["stage_root"]) == new_root
    assert old_root.is_dir() and all(path.is_file() for path in old_manifests.values())
    assert all(not record["resumed"] for record in updated["stages"].values())
    assert all(
        Path(record["manifest"]["path"]).is_relative_to(new_root)
        for record in updated["stages"].values()
    )

    assert main([*command, "--resume"]) == 0
    resumed = json.loads(manifest_path.read_text())
    assert all(record["resumed"] for record in resumed["stages"].values())
