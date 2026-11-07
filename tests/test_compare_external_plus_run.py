from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil

import nibabel as nib
import numpy as np
from scipy.io import savemat

from pymacruise.surfaces import extract_levelset_surface, write_legacy_vtk
from scripts.compare_external_plus_run import PRODUCTS, compare


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(path: Path) -> dict[str, str]:
    return {"path": str(path), "sha256": _sha(path)}


def test_plus_comparator_accepts_complete_synthetic_self_comparison(tmp_path: Path) -> None:
    fixture = Path(__file__).parent / "data/synthetic/surfaces/sphere.nii.gz"
    source_image = nib.load(str(fixture))
    shape, affine = source_image.shape, source_image.affine
    level = np.asarray(source_image.dataobj, dtype=np.float32)
    labels = (level <= 0).astype(np.int16) * 100
    t1 = tmp_path / "T1.nii.gz"
    flair = tmp_path / "FLAIR.nii.gz"
    mask = tmp_path / "brain-mask.nii.gz"
    maseg = tmp_path / "MaSeg.nii.gz"
    for path, values in ((t1, level), (flair, level), (mask, (level <= 0).astype(np.uint8)), (maseg, labels)):
        image = nib.Nifti1Image(values, affine)
        image.header.set_xyzt_units("mm")
        nib.save(image, path)

    python_root = tmp_path / "python"
    surface_root = python_root / "plus-surface"
    refinement_root = python_root / "refinement"
    atlasem_root = python_root / "atlasem"
    surface_root.mkdir(parents=True)
    refinement_root.mkdir(parents=True)
    atlasem_root.mkdir(parents=True)
    atlas_descriptor = atlasem_root / "ibsr8obj.txt"
    atlas_descriptor.write_text("synthetic descriptor")
    atlas_payload = atlasem_root / "atlas.nii.gz"
    nib.save(nib.Nifti1Image(np.ones(shape, dtype=np.float32), affine), atlas_payload)
    lesion_descriptor = atlasem_root / "lesion-atlas.txt"
    lesion_descriptor.write_text("synthetic LesionTOADS descriptor")
    lesion_topology = atlasem_root / "lesion-topology.nii.gz"
    nib.save(nib.Nifti1Image(labels, affine), lesion_topology)
    lesion_shape = atlasem_root / "lesion-shape.nii.gz"
    nib.save(nib.Nifti1Image(np.ones(shape, dtype=np.float32), affine), lesion_shape)
    atlasem_manifest = atlasem_root / "atlasem-dura-handoff.json"
    atlasem_manifest.write_text(json.dumps({
        "status": "complete", "inputs": {
            "atlasem_descriptor": _record(atlas_descriptor),
            "atlasem_assets": {str(atlas_payload): _sha(atlas_payload)},
        }, "outputs": {},
    }))
    refinement = refinement_root / "SEG_refine.nii.gz"
    refined_image = nib.Nifti1Image(labels, affine)
    refined_image.header.set_xyzt_units("mm")
    nib.save(refined_image, refinement)
    thickness = python_root / "thickness.nii.gz"
    thickness_image = nib.Nifti1Image(np.ones(shape, dtype=np.float32), affine)
    thickness_image.header.set_xyzt_units("mm")
    nib.save(thickness_image, thickness)
    outputs: dict[str, Path] = {
        "inner_level_set": python_root / "inner_level_set.nii.gz",
        "central_level_set": python_root / "central_level_set.nii.gz",
        "outer_level_set": python_root / "outer_level_set.nii.gz",
        "thickness_volume": thickness,
    }
    for key in ("inner_level_set", "central_level_set", "outer_level_set"):
        image = nib.Nifti1Image(level, affine)
        image.header.set_xyzt_units("mm")
        nib.save(image, outputs[key])
    mesh = extract_levelset_surface(source_image)
    for name in ("inner_surface", "central_surface", "outer_surface"):
        path = surface_root / f"{name}.vtk"
        write_legacy_vtk(mesh, path)
        outputs[name] = path
    metrics = surface_root / "thickness.json"
    metrics.write_text(json.dumps({"surface": {"mean_mm": 1., "standard_deviation_mm": 0., "nonzero_vertex_count": len(mesh.vertices)}}))
    outputs["thickness_metrics"] = metrics
    surface_manifest = surface_root / "run.json"
    surface_manifest.write_text(json.dumps({
        "status": "complete",
        "stages": {"surfaces_thickness": {"status": "complete", "outputs": {
            name: _record(path) for name, path in outputs.items()
        }}},
    }))
    refinement_manifest = refinement_root / "plus-segmentation-refinement.json"
    refinement_manifest.write_text(json.dumps({
        "status": "complete", "outputs": {"SEG_refine.nii.gz": _record(refinement)},
    }))

    stage_names = (
        "preprocessing", "atlasem_dura", "crop", "initial_alignment", "registered_recrop", "lesion_workflow",
        "dura_outputs", "native_dura_transfer", "maseg_cortical_separation",
        "plus_surface_workflow", "maseg_refinement",
    )
    stages = {}
    for name in stage_names:
        child = surface_manifest if name == "plus_surface_workflow" else (
            refinement_manifest if name == "maseg_refinement" else
            atlasem_manifest if name == "atlasem_dura" else python_root / f"{name}.json"
        )
        if not child.exists():
            content = {"status": "complete", "outputs": {}}
            if name == "initial_alignment":
                content["inputs"] = {"shape_atlases": {"class0": _record(lesion_shape)}}
            child.write_text(json.dumps(content))
        stages[name] = {"status": "complete", "manifest": _record(child)}
    inputs = {name: _record(path) for name, path in {
        "t1": t1, "flair": flair, "brain_mask": mask, "maseg": maseg,
    }.items()}
    inputs["atlasem_descriptor"] = _record(atlas_descriptor)
    inputs["atlas_descriptor"] = _record(lesion_descriptor)
    inputs["topology"] = _record(lesion_topology)
    (python_root / "plus-full-workflow.json").write_text(json.dumps({
        "status": "complete", "terminal": True, "stage_order": list(stage_names),
        "settings": {"atlasem_dura_enabled": True}, "inputs": inputs, "stages": stages,
        "final_outputs": {"SEG_refine": _record(refinement)},
    }))

    external_root = tmp_path / "external"
    external_root.mkdir()
    external = {}
    for name in PRODUCTS:
        source = refinement if name == "segmentation" else (
            metrics if name == "surface_thickness_mat" else outputs[name]
        )
        destination = external_root / ("surface-thickness.mat" if name == "surface_thickness_mat" else source.name)
        if name == "surface_thickness_mat":
            savemat(destination, {"whole_thickness": np.ones(len(mesh.vertices))})
        else:
            shutil.copyfile(source, destination)
        external[name] = destination
    product_map = tmp_path / "product-map.json"
    product_map.write_text(json.dumps({name: path.name for name, path in external.items()}))
    reference_info = tmp_path / "reference.json"
    reference_info.write_text(json.dumps({
        "status": "complete", "exit_code": 0, "runtime": "synthetic reference fixture",
        "source_revision": "synthetic", "command": "synthetic fixture generation",
        "atlasem_enabled": True, "atlasem_descriptor_sha256": _sha(atlas_descriptor),
        "atlasem_asset_sha256": {atlas_payload.name: _sha(atlas_payload)},
        "lesion_atlas_asset_sha256": {lesion_shape.name: _sha(lesion_shape)},
        "products_scope": "independently_generated_selected_plus_products",
        "preprocessing_boundary": "shared raw T1, FLAIR, mask, and MaSeg; no intermediate reuse",
        "surface_frame": "RAS+ mm", "required_segmentation_labels": [100],
        "active_product_attempts": {name: "synthetic-attempt-01" for name in PRODUCTS},
        "input_sha256": {name: _sha(path) for name, path in {
            "t1": t1, "flair": flair, "brain_mask": mask, "maseg": maseg,
            "atlasem_descriptor": atlas_descriptor,
            "atlas_descriptor": lesion_descriptor, "topology": lesion_topology,
        }.items()},
        "product_sha256": {name: _sha(path) for name, path in external.items()},
    }))
    log = tmp_path / "reference.log"
    log.write_text("synthetic fixture reference run complete\n")
    args = type("Args", (), {
        "python_run": python_root, "reference_info": reference_info, "reference_log": log,
        "t1": t1, "flair": flair, "brain_mask": mask, "maseg": maseg,
        "product_map": product_map, "external_root": external_root,
    })()
    report = compare(args)
    assert report["segmentation_and_surface_gates_passed"]
    assert report["products"]["segmentation"]["minimum_dice"] == 1.0
    assert all(report["products"][name]["hausdorff_mm"] == 0.0 for name in (
        "inner_surface", "central_surface", "outer_surface",
    ))
    assert "vertexwise_comparison" in report["products"]["surface_thickness"]


def test_plus_comparator_rejects_incomplete_active_products(tmp_path: Path) -> None:
    from scripts.compare_external_plus_run import _python_products
    import pytest

    with pytest.raises((KeyError, ValueError)):
        _python_products({"stages": {}})
