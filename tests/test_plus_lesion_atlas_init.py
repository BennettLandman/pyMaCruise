from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from pymacruise.image_io import save_nifti_ras
from pymacruise.lesion_initial_atlas_alignment import initialize_lesion_atlas_alignment_from_descriptor
import pymacruise.plus_lesion_atlas_init as alignment_stage
from pymacruise.plus_lesion_crop import run_plus_lesion_crop_checkpoint


ATLAS = Path(__file__).parents[1] / "data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"
SYNTHETIC = Path(__file__).parent / "data/synthetic/plus_lesion_atlas_init"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_plus_lesion_file_crop_flows_into_13_class_initial_atlas_alignment(
    tmp_path: Path, monkeypatch
) -> None:
    modalities_path = SYNTHETIC / "synthetic_t1_flair.nii.gz"
    topology_path = SYNTHETIC / "synthetic_topology.nii.gz"
    crop = run_plus_lesion_crop_checkpoint(
        modalities_path,
        topology_path,
        tmp_path / "crop",
        padding=1,
    )

    result = alignment_stage.run_plus_lesion_initial_atlas_alignment_checkpoint(
        crop.cropped_path,
        crop.manifest_path,
        ATLAS,
        tmp_path / "initial-alignment",
        init_align_iterations=0,
        pyramid_levels=2,
    )
    assert not result.resumed
    memberships = nib.load(str(result.memberships_path))
    cropped = nib.load(str(crop.cropped_path))
    assert memberships.shape == (*cropped.shape[:3], 13)
    assert nib.aff2axcodes(memberships.affine) == ("R", "A", "S")
    np.testing.assert_allclose(memberships.affine, cropped.affine, rtol=0, atol=1e-6)
    actual = np.asarray(memberships.dataobj)
    assert np.isfinite(actual).all()

    manifest = json.loads(result.manifest_path.read_text())
    summary = json.loads(result.summary_path.read_text())
    assert manifest["status"] == "complete"
    shape_assets = manifest["inputs"]["shape_atlases"]
    assert len(shape_assets) == 13
    assert all(Path(record["path"]).is_file() and record["sha256"] for record in shape_assets.values())
    assert manifest["source_order"] == [
        "preprocess.initialCentroids",
        "preprocess.initialMemberships",
        "atlas.initShapeRegistration",
        "atlas.setTransform(preprocess.getTransform())",
        "atlas.updateShapeRegistration",
        "atlas.registerShapesPyramid",
    ]
    assert summary["modalities"] == ["T1", "FLAIR"]
    assert len(summary["class_names"]) == 13
    assert summary["membership_shape_xyz_classes"] == list(memberships.shape)

    crop_data = json.loads(crop.manifest_path.read_text())
    transform = np.asarray(crop_data["crop"]["post_crop_transform"], dtype=np.float32)
    bounds = tuple(crop_data["crop"]["bounds"])
    direct = initialize_lesion_atlas_alignment_from_descriptor(
        str(ATLAS),
        np.asarray(cropped.dataobj, dtype=np.float32).transpose(3, 0, 1, 2),
        transform,
        crop_bounds=(
            (bounds[0], bounds[1] + 1),
            (bounds[2], bounds[3] + 1),
            (bounds[4], bounds[5] + 1),
        ),
        pad_border=1,
        image_resolutions_mm=(0.9375,) * 3,
        template_resolutions_mm=(0.9375,) * 3,
        init_align_iterations=0,
        pyramid_levels=2,
    )
    np.testing.assert_array_equal(
        actual.view(np.uint32),
        direct.alignment.memberships.memberships.view(np.uint32),
    )
    np.testing.assert_array_equal(
        np.asarray(result.transform, dtype=np.float32).view(np.uint32),
        direct.alignment.registration.transform.view(np.uint32),
    )

    def unexpected(*args, **kwargs):
        raise AssertionError("valid resume should not repeat atlas initialization")

    monkeypatch.setattr(
        alignment_stage,
        "initialize_lesion_atlas_alignment_from_descriptor",
        unexpected,
    )
    resumed = alignment_stage.run_plus_lesion_initial_atlas_alignment_checkpoint(
        crop.cropped_path,
        crop.manifest_path,
        ATLAS,
        tmp_path / "initial-alignment",
        init_align_iterations=0,
        pyramid_levels=2,
        resume=True,
    )
    assert resumed.resumed

    original_hash = alignment_stage._sha256
    shape_paths = {
        path.resolve()
        for path in alignment_stage.lesion_atlas_shape_asset_paths(ATLAS).values()
    }
    monkeypatch.setattr(
        alignment_stage,
        "_sha256",
        lambda path: "f" * 64 if Path(path).resolve() in shape_paths else original_hash(path),
    )
    with pytest.raises(ValueError, match="inputs or settings differ"):
        alignment_stage.run_plus_lesion_initial_atlas_alignment_checkpoint(
            crop.cropped_path,
            crop.manifest_path,
            ATLAS,
            tmp_path / "initial-alignment",
            init_align_iterations=0,
            pyramid_levels=2,
            resume=True,
        )


def test_plus_lesion_initial_align_cli_routes_registration_settings(tmp_path: Path, monkeypatch) -> None:
    from pymacruise.cli import main

    captured = {}

    def fake_align(cropped, crop_manifest, descriptor, output, **kwargs):
        captured.update(cropped=cropped, crop_manifest=crop_manifest, descriptor=descriptor, output=output, **kwargs)
        return alignment_stage.PlusLesionAtlasInitResult(
            tmp_path / "manifest.json",
            tmp_path / "memberships.nii.gz",
            tmp_path / "summary.json",
            (0.0,) * 6,
            resumed=True,
        )

    monkeypatch.setattr(
        alignment_stage,
        "run_plus_lesion_initial_atlas_alignment_checkpoint",
        fake_align,
    )
    assert main(
        [
            "plus-lesion-initial-align",
            "--cropped-modalities", "cropped.nii.gz",
            "--crop-manifest", "crop.json",
            "--atlas-descriptor", str(ATLAS),
            "--output-dir", "out",
            "--init-align-iterations", "4",
            "--pyramid-levels", "3",
            "--resume",
        ]
    ) == 0
    assert captured["init_align_iterations"] == 4
    assert captured["pyramid_levels"] == 3
    assert captured["resume"] is True
