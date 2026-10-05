from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_atlas_priors import load_lesion_atlas_topology
from pymacruise.lesion_initial_atlas_alignment import initialize_lesion_atlas_alignment_from_descriptor
from pymacruise.lesion_preprocessing_crop import crop_lesion_toads_preprocessing_inputs
from pymacruise.lesion_registered_recrop import run_registered_lesion_recrop
import pymacruise.plus_lesion_registered_recrop as recrop_stage
from pymacruise.plus_lesion_atlas_init import run_plus_lesion_initial_atlas_alignment_checkpoint
from pymacruise.plus_lesion_crop import run_plus_lesion_crop_checkpoint


ROOT = Path(__file__).parent
SYNTHETIC = ROOT / "data/synthetic/plus_lesion_atlas_init"
ATLAS = SYNTHETIC / "mini-lesion-toads-atlas.txt"


def _run_upstream(tmp_path: Path):
    crop = run_plus_lesion_crop_checkpoint(
        SYNTHETIC / "synthetic_t1_flair.nii.gz",
        SYNTHETIC / "synthetic_topology.nii.gz",
        tmp_path / "crop",
        padding=1,
    )
    alignment = run_plus_lesion_initial_atlas_alignment_checkpoint(
        crop.cropped_path,
        crop.manifest_path,
        ATLAS,
        tmp_path / "alignment",
        init_align_iterations=0,
        pyramid_levels=2,
    )
    return crop, alignment


def test_registered_recrop_file_checkpoint_matches_array_adapter_and_resumes(
    tmp_path: Path, monkeypatch
) -> None:
    crop, alignment = _run_upstream(tmp_path)
    result = recrop_stage.run_plus_lesion_registered_recrop_checkpoint(
        crop.normalized_path,
        crop.manifest_path,
        alignment.manifest_path,
        ATLAS,
        tmp_path / "registered-recrops",
    )
    modalities = nib.load(str(result.modalities_path))
    topology = nib.load(str(result.topology_path))
    shapes = nib.load(str(result.shapes_path))
    assert len(modalities.shape) == 4 and modalities.shape[-1] == 2
    assert topology.shape == modalities.shape[:3]
    assert shapes.shape == (*modalities.shape[:3], 13)
    assert all(nib.aff2axcodes(image.affine) == ("R", "A", "S") for image in (modalities, topology, shapes))
    assert np.isfinite(np.asarray(modalities.dataobj)).all()
    assert np.isfinite(np.asarray(shapes.dataobj)).all()
    np.testing.assert_allclose(modalities.affine, topology.affine, rtol=0, atol=1e-6)
    np.testing.assert_allclose(modalities.affine, shapes.affine, rtol=0, atol=1e-6)

    normalized = load_nifti_ras(crop.normalized_path)
    initial_topology = load_nifti_ras(SYNTHETIC / "synthetic_topology.nii.gz")
    spacing = (0.9375,) * 3
    crop_arrays = crop_lesion_toads_preprocessing_inputs(
        np.asarray(normalized.dataobj, dtype=np.float32).transpose(3, 0, 1, 2),
        np.asarray(initial_topology.dataobj),
        np.zeros(6, dtype=np.float32),
        image_resolution_mm=spacing,
        template_resolution_mm=spacing,
        padding=1,
    )
    direct_alignment = initialize_lesion_atlas_alignment_from_descriptor(
        str(ATLAS),
        crop_arrays.crop.images,
        crop_arrays.crop.transform,
        crop_bounds=tuple((crop.crop_bounds[i], crop.crop_bounds[i + 1] + 1) for i in (0, 2, 4)),
        pad_border=1,
        image_resolutions_mm=spacing,
        template_resolutions_mm=spacing,
        init_align_iterations=0,
        pyramid_levels=2,
    )
    direct = run_registered_lesion_recrop(
        direct_alignment,
        crop_arrays,
        load_lesion_atlas_topology(ATLAS),
        (11, 11, 11),
        image_resolutions_mm=spacing,
        template_resolutions_mm=spacing,
    )
    np.testing.assert_array_equal(
        np.asarray(modalities.dataobj).view(np.uint32),
        direct.preprocessing_crop.crop.images.transpose(1, 2, 3, 0).view(np.uint32),
    )
    np.testing.assert_array_equal(np.asarray(topology.dataobj), direct.transformed_template)
    np.testing.assert_array_equal(
        np.asarray(shapes.dataobj).view(np.uint32), direct.transformed_shapes.view(np.uint32)
    )

    manifest = json.loads(result.manifest_path.read_text())
    summary = json.loads(result.summary_path.read_text())
    assert manifest["status"] == "complete"
    assert summary["spatial_convention"] == "RAS"
    assert summary["class_names"] == [
        "Background", "Sulcal-CSF", "Ventricles", "Cerebellum-GM", "Cerebrum-GM",
        "Caudate", "Thalamus", "Putamen", "Brainstem", "Cerebellum-WM",
        "OuterVentricular-WM", "InterVentricular-WM", "Frontal-WM",
    ]
    assert summary["stage_order"].index("preprocess.updateTransformedTemplate") > summary["stage_order"].index(
        "atlas.computeTransformedTemplate:optimizedAtlasTransform"
    )

    def unexpected(*args, **kwargs):
        raise AssertionError("valid resume must skip registration and recropping")

    monkeypatch.setattr(recrop_stage, "run_registered_lesion_recrop", unexpected)
    resumed = recrop_stage.run_plus_lesion_registered_recrop_checkpoint(
        crop.normalized_path,
        crop.manifest_path,
        alignment.manifest_path,
        ATLAS,
        tmp_path / "registered-recrops",
        resume=True,
    )
    assert resumed.resumed
    monkeypatch.undo()
    result.shapes_path.write_bytes(b"corrupt")
    repaired = recrop_stage.run_plus_lesion_registered_recrop_checkpoint(
        crop.normalized_path,
        crop.manifest_path,
        alignment.manifest_path,
        ATLAS,
        tmp_path / "registered-recrops",
        resume=True,
    )
    assert not repaired.resumed
    assert nib.load(str(repaired.shapes_path)).shape[-1] == 13


def test_registered_recrop_cli_routes_manifest_inputs(tmp_path: Path, monkeypatch) -> None:
    from pymacruise.cli import main

    captured = {}

    def fake_recrop(normalized, crop_manifest, alignment_manifest, descriptor, output, **kwargs):
        captured.update(
            normalized=normalized,
            crop_manifest=crop_manifest,
            alignment_manifest=alignment_manifest,
            descriptor=descriptor,
            output=output,
            **kwargs,
        )
        return recrop_stage.PlusLesionRecropResult(
            tmp_path / "manifest.json",
            tmp_path / "modalities.nii.gz",
            tmp_path / "topology.nii.gz",
            tmp_path / "shapes.nii.gz",
            tmp_path / "summary.json",
            True,
        )

    monkeypatch.setattr(recrop_stage, "run_plus_lesion_registered_recrop_checkpoint", fake_recrop)
    assert main(
        [
            "plus-lesion-registered-recrop",
            "--normalized-modalities", "normalized.nii.gz",
            "--crop-manifest", "crop.json",
            "--initial-alignment-manifest", "alignment.json",
            "--atlas-descriptor", str(ATLAS),
            "--output-dir", "out",
            "--resume",
        ]
    ) == 0
    assert captured["resume"] is True
    assert str(captured["alignment_manifest"]) == "alignment.json"
