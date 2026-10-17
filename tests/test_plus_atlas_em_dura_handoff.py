from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.plus_atlas_em_dura_handoff import (
    run_plus_atlas_em_dura_handoff_checkpoint,
)


def _synthetic_inputs(root: Path) -> tuple[Path, Path, Path, Path]:
    root.mkdir(parents=True, exist_ok=True)
    shape = (12, 12, 12)
    x, y, z = np.indices(shape, dtype=np.float32)
    radius = np.sqrt((x - 5.5) ** 2 + (y - 5.5) ** 2 + (z - 5.5) ** 2)
    brain = radius < 4.6
    t1_values = np.where(brain, 0.35 + 0.025 * x + 0.015 * y, 0.0).astype(np.float32)
    flair_values = np.where(brain, 0.45 + 0.01 * z, 0.0).astype(np.float32)
    affine = np.diag((1.0, 1.0, 1.0, 1.0))
    for name, values in (
        ("t1_brain.nii.gz", t1_values),
        ("flair_brain.nii.gz", flair_values),
        ("brain_mask.nii.gz", brain.astype(np.uint8)),
    ):
        nib.save(nib.Nifti1Image(values, affine), root / name)

    gm = np.exp(-((radius - 3.5) ** 2) / 2.5).astype(np.float32)
    wm = np.exp(-(radius**2) / 18.0).astype(np.float32)
    background = np.ones(shape, dtype=np.float32)
    nib.save(
        nib.Nifti1Image(np.where(radius < 2.4, 3, np.where(radius < 4.6, 2, 1)).astype(np.uint8), affine),
        root / "topology.nii.gz",
    )
    nib.save(nib.Nifti1Image(background, affine), root / "background.nii.gz")
    nib.save(nib.Nifti1Image(gm, affine), root / "gm.nii.gz")
    nib.save(nib.Nifti1Image(wm, affine), root / "wm.nii.gz")
    descriptor = root / "atlasem.txt"
    descriptor.write_text(
        "\n".join((
            "Structure Atlas File (edit at your own risks)", "", "Structures:\t3",
            "Background\t1\tmask", "Cerebrum-GM\t2\tclosed", "Cerebrum-WM\t3\tclosed",
            "", "Topology Atlas:", "File:\ttopology.nii.gz", "Dimensions:\t12\t12\t12",
            "Resolutions:\t1\t1\t1", "", "Shape Atlas:\t4", "Dimensions:\t12\t12\t12",
            "Resolutions:\t1\t1\t1", "Orientations:\t1\t1\t6\t3",
            "Structure:\tBackground", "File:\tbackground.nii.gz", "Structure:\tCerebrum-GM",
            "File:\tgm.nii.gz", "Structure:\tCerebrum-WM", "File:\twm.nii.gz",
            "Registered Shapes:", "0\t1\t1", "0\t1\t1", "",
            "Intensity Atlas:\t10", "T1_MPRAGE\t0.2\t0.5\t0.8",
        )) + "\n",
        encoding="utf-8",
    )
    return (root / "t1_brain.nii.gz", root / "flair_brain.nii.gz",
            root / "brain_mask.nii.gz", descriptor)


def test_synthetic_prelesion_atlasem_dura_handoff_and_resume(tmp_path: Path) -> None:
    inputs = _synthetic_inputs(tmp_path / "inputs")
    t1, flair, mask, descriptor = inputs
    out = tmp_path / "atlasem-dura-stage"
    result = run_plus_atlas_em_dura_handoff_checkpoint(
        t1, flair, mask, descriptor, out,
        levels=2,
        initial_registration_iterations=1,
        post_recrop_registration_iterations=1,
        max_em_iterations=1,
        max_dura_iterations=1,
    )
    manifest = json.loads(result.manifest_path.read_text())
    assert manifest["status"] == "complete"
    assert manifest["source_order"] == [
        "AtlasEM on prepared T1", "selected AtlasEM-output Dura removal",
        "Dura-mask registered FLAIR", "form Dura-conditioned T1/FLAIR handoff",
    ]
    assert set(manifest["outputs"]) >= {
        "memberships", "segmentation", "dura_brain_mask", "dura_stripped_t1",
        "dura_masked_flair", "lesion_modalities",
    }
    for path in result.output_paths.values():
        assert path.name.endswith(".nii.gz")
        image = nib.load(str(path))
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert image.shape[:3] == (12, 12, 12)
        assert np.isfinite(np.asarray(image.dataobj)).all()

    modalities = np.asarray(nib.load(str(result.modalities_path)).dataobj)
    dura_t1 = np.asarray(nib.load(str(result.dura_stripped_t1_path)).dataobj)
    dura_flair = np.asarray(nib.load(str(result.dura_masked_flair_path)).dataobj)
    dura_mask = np.asarray(nib.load(str(result.output_paths["dura_brain_mask"])).dataobj) > 0
    np.testing.assert_array_equal(modalities[..., 0], dura_t1)
    np.testing.assert_array_equal(modalities[..., 1], dura_flair)
    assert np.all(dura_flair[~dura_mask] == 0)
    assert np.all(dura_flair >= 0)

    resumed = run_plus_atlas_em_dura_handoff_checkpoint(
        t1, flair, mask, descriptor, out,
        levels=2,
        initial_registration_iterations=1,
        post_recrop_registration_iterations=1,
        max_em_iterations=1,
        max_dura_iterations=1,
        resume=True,
    )
    assert resumed.resumed is True
    assert resumed.modalities_path == result.modalities_path

