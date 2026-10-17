from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.cli import main


def _write_synthetic_descriptor(directory: Path) -> Path:
    shape = (12, 12, 12)
    x, y, z = np.indices(shape, dtype=np.float32)
    radius = np.sqrt((x - 5.5) ** 2 + (y - 5.5) ** 2 + (z - 5.5) ** 2)
    topology = np.where(radius < 2.4, 3, np.where(radius < 4.6, 2, 1)).astype(np.uint8)
    background = np.ones(shape, dtype=np.float32)
    gm = np.exp(-((radius - 3.5) ** 2) / 2.5).astype(np.float32)
    wm = np.exp(-(radius**2) / 18.0).astype(np.float32)
    affine = np.eye(4)
    for name, data in (("topology", topology), ("background", background), ("gm", gm), ("wm", wm)):
        nib.save(nib.Nifti1Image(data, affine), directory / f"{name}.nii.gz")
    descriptor = directory / "atlas.txt"
    descriptor.write_text(
        "\n".join(
            (
                "Structure Atlas File (edit at your own risks)", "", "Structures:\t3",
                "Background\t1\tmask", "Cerebrum-GM\t2\tclosed", "Cerebrum-WM\t3\tclosed",
                "", "Topology Atlas:", "File:\ttopology.nii.gz", "Dimensions:\t12\t12\t12",
                "Resolutions:\t1\t1\t1", "", "Shape Atlas:\t4", "Dimensions:\t12\t12\t12",
                "Resolutions:\t1\t1\t1", "Orientations:\t1\t1\t6\t3",
                "Structure:\tBackground", "File:\tbackground.nii.gz", "Structure:\tCerebrum-GM",
                "File:\tgm.nii.gz", "Structure:\tCerebrum-WM", "File:\twm.nii.gz",
                "Registered Shapes:", "0\t1\t1", "0\t1\t1", "", "Intensity Atlas:\t10",
                "T1_MPRAGE\t0.2\t0.5\t0.8",
            )
        ) + "\n",
        encoding="utf-8",
    )
    return descriptor


def _write_eight_object_descriptor(directory: Path) -> Path:
    """Build an 8-class atlas-shaped case without relying on external assets."""
    shape = (12, 12, 12)
    x, y, z = np.indices(shape, dtype=np.float32)
    centers = ((5.5, 5.5, 5.5), (3.5, 5.5, 5.5), (7.5, 5.5, 5.5),
               (5.5, 3.5, 5.5), (5.5, 7.5, 5.5), (5.5, 5.5, 3.5),
               (5.5, 5.5, 7.5))
    names = (
        "Background", "SulcalCSF", "CorticalGM", "Cerebrum-WM",
        "Cerebellum-GM", "Cerebellum-WM", "Brainstem", "Ventricle",
    )
    topology = np.ones(shape, dtype=np.uint8)
    paths = []
    background = np.ones(shape, dtype=np.float32)
    nib.save(nib.Nifti1Image(background, np.eye(4)), directory / "obj0.nii.gz")
    paths.append("obj0.nii.gz")
    for class_index, (cx, cy, cz) in enumerate(centers, start=1):
        prior = np.exp(-((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2) / 20.0)
        nib.save(nib.Nifti1Image(prior.astype(np.float32), np.eye(4)), directory / f"obj{class_index}.nii.gz")
        paths.append(f"obj{class_index}.nii.gz")
        topology[prior > 0.55] = class_index + 1
    nib.save(nib.Nifti1Image(topology, np.eye(4)), directory / "topology8.nii.gz")
    descriptor = directory / "atlas8.txt"
    rows = [
        "Structure Atlas File (edit at your own risks)", "", "Structures:\t8",
        *[f"{name}\t{index + 1}\t{'mask' if index == 0 else 'closed'}" for index, name in enumerate(names)],
        "", "Topology Atlas:", "File:\ttopology8.nii.gz", "Dimensions:\t12\t12\t12",
        "Resolutions:\t1\t1\t1", "", "Shape Atlas:\t4", "Dimensions:\t12\t12\t12",
        "Resolutions:\t1\t1\t1", "Orientations:\t1\t1\t6\t3",
    ]
    for name, payload in zip(names, paths, strict=True):
        rows.extend((f"Structure:\t{name}", f"File:\t{payload}"))
    rows.extend(("Registered Shapes:", "0\t1\t1\t1\t1\t1\t1\t1", "0\t1\t1\t1\t1\t1\t1\t1",
                 "", "Intensity Atlas:\t10", "T1_MPRAGE\t0.1\t0.2\t0.3\t0.4\t0.5\t0.6\t0.7\t0.8"))
    descriptor.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return descriptor


def test_atlasem_run_cli_executes_synthetic_descriptor_workflow(tmp_path: Path) -> None:
    shape = (12, 12, 12)
    x, y, z = np.indices(shape, dtype=np.float32)
    radius = np.sqrt((x - 5.5) ** 2 + (y - 5.5) ** 2 + (z - 5.5) ** 2)
    t1_data = np.where(radius < 4.6, 0.35 + 0.025 * x + 0.015 * y, 0.0).astype(np.float32)
    t1_path = tmp_path / "t1.nii.gz"
    mask_path = tmp_path / "brain_mask.nii.gz"
    nib.save(nib.Nifti1Image(t1_data, np.eye(4)), t1_path)
    nib.save(nib.Nifti1Image((radius < 4.6).astype(np.uint8), np.eye(4)), mask_path)
    descriptor = _write_synthetic_descriptor(tmp_path)
    output = tmp_path / "out"
    provenance = tmp_path / "run.json"
    command = [
        "atlasem-run", "--t1", str(t1_path), "--brain-mask", str(mask_path),
        "--descriptor", str(descriptor), "--output-dir", str(output),
        "--crop-border", "2", "--levels", "2",
        "--initial-registration-iterations", "1", "--post-recrop-registration-iterations", "1",
        "--max-em-iterations", "1", "--max-dura-iterations", "1",
        "--provenance", str(provenance),
    ]

    status = main(command)

    assert status == 0
    expected = {
        "atlasem_memberships.nii.gz", "atlasem_segmentation.nii.gz",
        "atlasem_internal_mask.nii.gz", "atlasem_dura_inputs.nii.gz",
        "prelesion_dura_inputs.nii.gz", "prelesion_dura_segmentation.nii.gz",
        "dura_stripped_t1.nii.gz", "dura_corrected_white_matter.nii.gz",
        "dura_corrected_gray_matter.nii.gz", "dura_brain_mask_level_set.nii.gz",
        "dura_brain_mask.nii.gz", "dura_wm_topology.nii.gz", "dura_gm_topology.nii.gz",
    }
    assert {path.name for path in output.glob("*.nii.gz")} == expected
    assert all(nib.aff2axcodes(nib.load(str(output / name)).affine) == ("R", "A", "S") for name in expected)
    assert all(np.isfinite(np.asarray(nib.load(str(output / name)).dataobj)).all() for name in expected)
    manifest = json.loads(provenance.read_text(encoding="utf-8"))
    assert manifest["command"] == "atlasem-run"
    assert manifest["parameters"]["em_iterations_completed"] == 1
    assert manifest["parameters"]["dura_outer_passes_completed"] == 1
    stage_manifest_path = output / "run.json"
    first_run = json.loads(stage_manifest_path.read_text(encoding="utf-8"))
    assert first_run["status"] == "complete"
    assert [event["stage"] for event in first_run["history"]] == ["atlasem", "prelesion_dura"]

    assert main([*command, "--resume"]) == 0
    resumed = json.loads(stage_manifest_path.read_text(encoding="utf-8"))
    assert resumed["history"] == first_run["history"]

    damaged = output / "dura_brain_mask.nii.gz"
    damaged.write_bytes(b"corrupt stage output")
    assert main([*command, "--resume"]) == 0
    repaired = json.loads(stage_manifest_path.read_text(encoding="utf-8"))
    assert [event["stage"] for event in repaired["history"]] == [
        "atlasem", "prelesion_dura", "prelesion_dura"
    ]

    mask_data = nib.load(str(mask_path)).get_fdata().astype(np.uint8)
    mask_data[0, 0, 0] = 1
    nib.save(nib.Nifti1Image(mask_data, np.eye(4)), mask_path)
    with np.testing.assert_raises_regex(ValueError, "inputs, settings, atlas assets, or implementation"):
        main([*command, "--resume"])


def test_atlasem_run_cli_accepts_eight_object_synthetic_atlas(tmp_path: Path) -> None:
    shape = (12, 12, 12)
    x, y, z = np.indices(shape, dtype=np.float32)
    radius = np.sqrt((x - 5.5) ** 2 + (y - 5.5) ** 2 + (z - 5.5) ** 2)
    affine = np.eye(4)
    t1 = np.where(radius < 5.0, 0.25 + 0.01 * x + 0.005 * y, 0.0).astype(np.float32)
    nib.save(nib.Nifti1Image(t1, affine), tmp_path / "t1.nii.gz")
    nib.save(nib.Nifti1Image((radius < 5.0).astype(np.uint8), affine), tmp_path / "mask.nii.gz")
    descriptor = _write_eight_object_descriptor(tmp_path)
    output = tmp_path / "out8"
    status = main([
        "atlasem-run", "--t1", str(tmp_path / "t1.nii.gz"),
        "--brain-mask", str(tmp_path / "mask.nii.gz"),
        "--descriptor", str(descriptor), "--output-dir", str(output),
        "--crop-border", "2", "--levels", "2",
        "--initial-registration-iterations", "1", "--post-recrop-registration-iterations", "1",
        "--max-em-iterations", "2", "--max-dura-iterations", "1",
    ])

    assert status == 0
    manifest = json.loads((output / "run.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "complete"
    assert manifest["stages"]["atlasem"]["em_iterations_completed"] == 2
    memberships = nib.load(str(output / "atlasem_memberships.nii.gz"))
    assert memberships.shape == (*shape, 8)
    values = np.asarray(memberships.dataobj)
    assert np.isfinite(values).all()
    assert nib.aff2axcodes(memberships.affine) == ("R", "A", "S")
    assert {path.name for path in output.glob("*.nii.gz")} == {
        "atlasem_memberships.nii.gz", "atlasem_segmentation.nii.gz",
        "atlasem_internal_mask.nii.gz", "atlasem_dura_inputs.nii.gz",
        "prelesion_dura_inputs.nii.gz", "prelesion_dura_segmentation.nii.gz",
        "dura_stripped_t1.nii.gz", "dura_corrected_white_matter.nii.gz",
        "dura_corrected_gray_matter.nii.gz", "dura_brain_mask_level_set.nii.gz",
        "dura_brain_mask.nii.gz", "dura_wm_topology.nii.gz", "dura_gm_topology.nii.gz",
    }
