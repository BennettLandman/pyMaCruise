"""Synthetic contract tests for independent selected PLUS comparison."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import nibabel as nib
import numpy as np
import pytest
from scipy.io import savemat

from pymacruise.surfaces import SurfaceMesh, write_legacy_vtk

SCRIPT = Path(__file__).parents[1] / "scripts/compare_external_plus_run.py"
spec = importlib.util.spec_from_file_location("external_plus_comparison", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def _write_json(path: Path, value: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")
    return path


def _save_image(path: Path, values: np.ndarray, affine: np.ndarray) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    image = nib.Nifti1Image(values, affine)
    image.header.set_xyzt_units("mm")
    nib.save(image, path)
    return path


def case(tmp_path: Path, *, atlasem: bool = True, surface_offset_mm: float = 0.0):
    affine = np.eye(4)
    affine[:3, 3] = (-2, -2, -2)
    t1_values = np.arange(64, dtype=np.float32).reshape(4, 4, 4)
    flair_values = t1_values[::-1].copy()
    brain_mask = np.ones((4, 4, 4), dtype=np.uint8)
    labels = np.zeros((4, 4, 4), dtype=np.int16)
    labels[1:3, 1:3, 1:3] = 100
    inputs = {}
    for name, values in (("t1", t1_values), ("flair", flair_values),
                         ("brain_mask", brain_mask), ("maseg", labels)):
        inputs[name] = _save_image(tmp_path / "inputs" / f"{name}.nii.gz", values, affine)

    python_root = tmp_path / "python"
    external_root = tmp_path / "external"
    products_root = external_root / "outputs"
    products_root.mkdir(parents=True)
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    faces = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])

    def add_external(name: str, suffix: str, writer) -> Path:
        path = products_root / f"{name}{suffix}"
        path.parent.mkdir(parents=True, exist_ok=True)
        writer(path)
        return path

    external_products: dict[str, Path] = {}
    ext_values = {
        "segmentation": labels,
        "inner_level_set": np.where(labels > 0, -1.0, 1.0).astype(np.float32),
        "central_level_set": np.where(labels > 0, -0.5, 0.5).astype(np.float32),
        "outer_level_set": np.where(labels > 0, -2.0, 2.0).astype(np.float32),
        "thickness_volume": np.ones((4, 4, 4), dtype=np.float32),
    }
    for name, values in ext_values.items():
        external_products[name] = add_external(
            name, ".nii.gz", lambda path, values=values: _save_image(path, values, affine)
        )
    for name in module.MESHES:
        mesh_path = add_external(
            name, ".vtk", lambda path: write_legacy_vtk(SurfaceMesh(vertices, faces), path)
        )
        external_products[name] = mesh_path
    mat_path = add_external(
        "surface_thickness", ".mat",
        lambda path: savemat(path, {"whole_thickness": np.ones(4)}),
    )
    external_products["surface_thickness_mat"] = mat_path

    # Create child manifests and their products in the same schema as the
    # production full-workflow runner. Old/inactive attempts are omitted.
    stages: dict[str, dict] = {}
    stage_order = [
        "preprocessing", "crop", "initial_alignment", "registered_recrop",
        "lesion_workflow", "dura_outputs", "native_dura_transfer",
        "maseg_cortical_separation", "plus_surface_workflow", "maseg_refinement",
    ]
    if atlasem:
        stage_order.insert(1, "atlasem_dura")
    atlas_descriptor = _write_json(tmp_path / "atlas/descriptor.json", {"fixture": True})
    atlas_payload = tmp_path / "atlas/payload.nii.gz"
    _save_image(atlas_payload, t1_values, affine)
    atlas_assets = {str(atlas_payload): module.sha256(atlas_payload)}
    lesion_atlas = tmp_path / "atlas/lesion-shape-atlas.nii.gz"
    _save_image(lesion_atlas, labels.astype(np.float32), affine)
    lesion_assets = {lesion_atlas.name: module.sha256(lesion_atlas)}
    for stage in stage_order:
        stage_root = python_root / stage
        child_outputs = {}
        if stage == "plus_surface_workflow":
            output_specs = {
                "inner_level_set": ("inner_level_set.nii.gz", ext_values["inner_level_set"]),
                "central_level_set": ("central_level_set.nii.gz", ext_values["central_level_set"]),
                "outer_level_set": ("outer_level_set.nii.gz", ext_values["outer_level_set"]),
                "thickness_volume": ("thickness_volume.nii.gz", ext_values["thickness_volume"]),
            }
            for name, (filename, values) in output_specs.items():
                path = _save_image(stage_root / filename, values, affine)
                child_outputs[name] = {"path": str(path), "sha256": module.sha256(path)}
            meshes = {}
            for name in module.MESHES:
                shifted = vertices.copy()
                shifted[:, 0] += surface_offset_mm
                path = stage_root / f"{name}.vtk"
                path.parent.mkdir(parents=True, exist_ok=True)
                write_legacy_vtk(SurfaceMesh(shifted, faces), path)
                child_outputs[name] = {"path": str(path), "sha256": module.sha256(path)}
            metrics_path = _write_json(stage_root / "thickness_metrics.json", {
                "surface": {"mean_mm": 1.0, "standard_deviation_mm": 0.0,
                            "nonzero_vertex_count": 4}
            })
            child_outputs["thickness_metrics"] = {
                "path": str(metrics_path), "sha256": module.sha256(metrics_path)
            }
            surface_child = _write_json(stage_root / "run.json", {
                "status": "complete",
                "stages": {"surfaces_thickness": {
                    "status": "complete", "outputs": child_outputs
                }},
            })
        elif stage == "maseg_refinement":
            path = _save_image(stage_root / "SEG_refine.nii.gz", labels, affine)
            child_outputs["SEG_refine.nii.gz"] = {
                "path": str(path), "sha256": module.sha256(path)
            }
            surface_child = _write_json(stage_root / "run.json", {
                "status": "complete", "outputs": child_outputs,
            })
        elif stage == "atlasem_dura":
            surface_child = _write_json(stage_root / "run.json", {
                "status": "complete",
                "inputs": {
                    "atlasem_descriptor": {"sha256": module.sha256(atlas_descriptor)},
                    "atlasem_assets": atlas_assets,
                },
            })
        elif stage == "initial_alignment":
            surface_child = _write_json(stage_root / "run.json", {
                "status": "complete",
                "inputs": {"shape_atlases": {
                    lesion_atlas.name: {"path": str(lesion_atlas), "sha256": module.sha256(lesion_atlas)}
                }},
            })
        else:
            path = stage_root / f"{stage}.json"
            surface_child = _write_json(path, {"status": "complete", "stage": stage})
        stages[stage] = {
            "status": "complete",
            "manifest": {"path": str(surface_child), "sha256": module.sha256(surface_child)},
        }

    input_hashes = {name: module.sha256(path) for name, path in inputs.items()}
    python_manifest = {
        "status": "complete", "terminal": True,
        "stage_order": stage_order, "stages": stages,
        "settings": {"atlasem_dura_enabled": atlasem},
        "inputs": {name: {"path": str(path), "sha256": input_hashes[name]}
                   for name, path in inputs.items()},
    }
    run_path = _write_json(python_root / "plus-full-workflow.json", python_manifest)

    info_path = _write_json(tmp_path / "reference-info.json", {
        "status": "complete", "exit_code": 0, "runtime": "synthetic fixture",
        "source_revision": "test-only", "command": "not-an-external-run",
        "input_sha256": input_hashes,
        "atlasem_enabled": atlasem,
        "atlasem_descriptor_sha256": module.sha256(atlas_descriptor),
        "atlasem_asset_sha256": {Path(name).name: digest for name, digest in atlas_assets.items()},
        "lesion_atlas_asset_sha256": lesion_assets,
        "preprocessing_boundary": "same supplied T1/FLAIR/brain-mask/MaSeg files",
        "surface_frame": "RAS+ mm",
        "required_segmentation_labels": [100],
        "active_product_attempts": {name: "attempt-01" for name in module.PRODUCTS},
        "products_scope": "independently_generated_selected_plus_products",
        "product_sha256": {name: module.sha256(path) for name, path in external_products.items()},
    })
    log_path = tmp_path / "reference.log"
    log_path.write_text("Synthetic comparator test; no external execution.\n")
    product_map = _write_json(tmp_path / "product-map.json", {
        name: str(path.relative_to(products_root)) for name, path in external_products.items()
    })
    for name, filename in (("t1", "T1.nii.gz"), ("flair", "FLAIR.nii.gz"),
                           ("brain_mask", "BrainMask.nii.gz"), ("maseg", "MaSeg.nii.gz")):
        # Keep a stable descriptive alias for direct CLI argument paths.
        alias = tmp_path / filename
        alias.write_bytes(inputs[name].read_bytes())
        inputs[name] = alias
    # The hashes remain identical after writing aliases.
    args = SimpleNamespace(
        external_root=products_root, python_run=python_root,
        t1=inputs["t1"], flair=inputs["flair"], brain_mask=inputs["brain_mask"],
        maseg=inputs["maseg"], reference_info=info_path, reference_log=log_path,
        product_map=product_map, report=tmp_path / "comparison.json",
    )
    return args


def test_identical_terminal_products_pass_segmentation_and_surface_gates(tmp_path: Path):
    result = module.compare(case(tmp_path))
    assert result["segmentation_and_surface_gates_passed"]
    assert result["products"]["segmentation"]["minimum_dice"] == 1.0
    assert result["overall_port_acceptance"] == "unverified"


def test_atlasem_enabled_run_requires_and_accepts_atlas_stage(tmp_path: Path):
    args = case(tmp_path, atlasem=True)
    result = module.compare(args)
    assert result["atlasem_enabled"]
    args.python_run.joinpath("plus-full-workflow.json").write_text(json.dumps({
        **json.loads(args.python_run.joinpath("plus-full-workflow.json").read_text()),
        "stage_order": [stage for stage in json.loads(
            args.python_run.joinpath("plus-full-workflow.json").read_text()
        )["stage_order"] if stage != "atlasem_dura"],
    }))
    with pytest.raises(ValueError, match="stage order differs"):
        module.compare(args)


def test_nonterminal_or_wrong_flair_input_is_rejected(tmp_path: Path):
    args = case(tmp_path)
    path = args.python_run / "plus-full-workflow.json"
    manifest = json.loads(path.read_text())
    manifest["terminal"] = False
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="complete and terminal"):
        module.compare(args)


def test_surface_disagreement_fails_gate_but_keeps_metrics(tmp_path: Path):
    args = case(tmp_path, surface_offset_mm=0.05)
    result = module.compare(args)
    assert not result["segmentation_and_surface_gates_passed"]
    assert result["products"]["inner_surface"]["mean_absolute_distance_mm"] > 0.01


def test_input_identity_mismatch_and_corrupt_products_are_rejected(tmp_path: Path):
    args = case(tmp_path)
    wrong_flair = tmp_path / "wrong-flair.nii.gz"
    _save_image(wrong_flair, np.zeros((4, 4, 4), dtype=np.float32), np.diag([1, 1, 1, 1]))
    args.flair = wrong_flair
    with pytest.raises(ValueError, match="input hashes do not match"):
        module.compare(args)

    args = case(tmp_path / "corrupt")
    python_path = args.python_run / "maseg_refinement/SEG_refine.nii.gz"
    python_path.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="hash mismatch"):
        module.compare(args)


def test_cli_preserves_failure_report_for_missing_product(tmp_path: Path):
    args = case(tmp_path)
    args.product_map.write_text("{}")
    argv = []
    for key, value in vars(args).items():
        argv.extend(["--" + key.replace("_", "-"), str(value)])
    assert module.main(argv) == 2
    assert json.loads(args.report.read_text())["overall_port_acceptance"] == "unverified"


def test_exploratory_atlas_policy_reports_differences_without_waiving_inputs(tmp_path: Path):
    args = case(tmp_path)
    info = json.loads(args.reference_info.read_text())
    info['atlasem_descriptor_sha256'] = 'native-java-format-differs'
    args.reference_info.write_text(json.dumps(info))
    with pytest.raises(ValueError, match='AtlasEM descriptor or payload hashes differ'):
        module.compare(args)
    args.atlas_policy = 'report-differences'
    report = module.compare(args)
    assert report['segmentation_and_surface_gates_passed']
    assert report['atlas_identity_differences']
    assert report['overall_port_acceptance'] == 'unverified'
    info['input_sha256']['flair'] = 'wrong-scan'
    args.reference_info.write_text(json.dumps(info))
    with pytest.raises(ValueError, match='input hashes do not match'):
        module.compare(args)


def test_explicit_voxel_mesh_coordinates_are_mapped_with_t1_affine(tmp_path: Path):
    from pymacruise.surfaces import SurfaceMesh, read_legacy_vtk, write_legacy_vtk
    args = case(tmp_path)
    products = json.loads(args.product_map.read_text())
    affine = nib.load(args.t1).affine
    for name in module.MESHES:
        path = args.external_root / products[name]
        mesh = read_legacy_vtk(path, input_frame='RAS')
        vertices = nib.affines.apply_affine(np.linalg.inv(affine), mesh.vertices)
        write_legacy_vtk(SurfaceMesh(vertices, mesh.faces, mesh.coordinate_frame), path)
    info = json.loads(args.reference_info.read_text())
    info['surface_frame'] = 'T1 voxel indices'
    info['product_sha256'] = {name: module.sha256(args.external_root / path) for name, path in products.items()}
    args.reference_info.write_text(json.dumps(info))
    args.source_mesh_frame = 'voxel'
    result = module.compare(args)
    assert result['segmentation_and_surface_gates_passed']


def test_python_only_audit_accepts_terminal_tree_and_rejects_changed_final_product(tmp_path, monkeypatch):
    import importlib
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / 'scripts'))
    audit = importlib.import_module('audit_python_plus_run').audit
    args = case(tmp_path)
    result = audit(args.python_run)
    assert result['status'] == 'passed'
    assert result['stage_count'] == 11
    assert result['verified_file_count'] > 10
    assert len(result['meshes']) == 3
    product = module._python_products(json.loads((args.python_run / 'plus-full-workflow.json').read_text()))['segmentation']
    product.write_bytes(b'corrupt')
    with pytest.raises(ValueError, match='hash (verification failed|mismatch)'):
        audit(args.python_run)


def test_record_original_run_uses_original_outputs_and_preserves_atlas_differences(tmp_path, monkeypatch):
    import importlib
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / 'scripts'))
    recorder = importlib.import_module('record_plus_reference')
    args = case(tmp_path)
    native = tmp_path / 'native-atlas.txt'
    raw = tmp_path / 'atlas.raw'
    raw.write_bytes(b'original-native-atlas-payload')
    native.write_text('File:\tatlas.raw\n')
    kit = tmp_path / 'kit.json'
    kit.write_text(json.dumps({
        'inputs': {name: {'path': str(getattr(args, name)), 'sha256': module.sha256(getattr(args, name))}
                   for name in ('t1', 'flair', 'brain_mask', 'maseg')},
        'source_atlasem_descriptor': str(native),
        'source_atlasem_descriptor_sha256': module.sha256(native),
    }))
    source_map = json.loads(args.product_map.read_text())
    outputs = {('refined_segmentation' if name == 'segmentation' else name): str(args.external_root / path)
               for name, path in source_map.items()}
    outputs.update(input_dir=str(args.t1.parent), atlasem_descriptor=str(native))
    (args.external_root / 'driver_outputs.json').write_text(json.dumps(outputs))
    record_args = SimpleNamespace(kit=kit, source_root=args.external_root, lesion_descriptor=native,
        log=args.reference_log, exit_code=0, runtime='synthetic original runtime',
        source_revision='synthetic source', command='synthetic independent run',
        mesh_frame='ras-mm', attempt='fixture-01')
    info_path = recorder.record(record_args)
    args.reference_info = info_path
    args.product_map = args.external_root / 'product-map.json'
    args.atlas_policy = 'report-differences'
    report = module.compare(args)
    assert report['segmentation_and_surface_gates_passed']
    assert report['atlas_identity_differences']
    # A changed starting scan cannot be papered over by a new run record.
    args.flair.write_bytes(b'changed')
    with pytest.raises(ValueError, match='input directory/hash differs'):
        recorder.record(record_args)
