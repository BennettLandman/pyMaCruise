from __future__ import annotations

import json

import nibabel as nib
import numpy as np

from pymacruise.cli import main


def _write_descriptor(path):
    path.write_text(
        "\n".join(
            (
                "Structure Atlas File (edit at your own risks)", "", "Structures:\t2",
                "Background\t1\tmask", "Cortex\t2\tclosed", "", "Topology Atlas:",
                "File:\ttopology.nii.gz", "Dimensions:\t5\t6\t7", "Resolutions:\t1\t1\t1",
                "", "Shape Atlas:\t4", "Dimensions:\t5\t6\t7", "Resolutions:\t1\t1\t1",
                "Orientations:\t1\t1\t6\t3", "Structure:\tBackground", "File:\tbackground.nii.gz",
                "Structure:\tCortex", "File:\tcortex.nii.gz",
            )
        ) + "\n",
        encoding="utf-8",
    )


def test_atlasem_inventory_reports_descriptor_resources_without_loading_voxels(tmp_path, capsys):
    descriptor = tmp_path / "atlas.txt"
    _write_descriptor(descriptor)
    image = nib.Nifti1Image(np.ones((5, 6, 7), dtype=np.float32), np.eye(4))
    for name in ("topology.nii.gz", "background.nii.gz", "cortex.nii.gz"):
        nib.save(image, tmp_path / name)

    exit_code = main(["atlasem-atlas-inventory", "--descriptor", str(descriptor)])
    inventory = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert inventory["resources_present_as_nifti_gz"] is True
    assert inventory["voxel_data_loaded"] is False
    assert [item["name"] for item in inventory["class_order"]] == ["Background", "Cortex"]
    assert all(item["exists"] and item["compressed_nifti"] for item in inventory["assets"])


def test_atlasem_inventory_exits_nonzero_for_missing_payload(tmp_path, capsys):
    descriptor = tmp_path / "atlas.txt"
    _write_descriptor(descriptor)
    nib.save(
        nib.Nifti1Image(np.ones((5, 6, 7), dtype=np.uint8), np.eye(4)),
        tmp_path / "topology.nii.gz",
    )

    exit_code = main(["atlasem-atlas-inventory", "--descriptor", str(descriptor)])
    inventory = json.loads(capsys.readouterr().out)

    assert exit_code == 1
    assert inventory["resources_present_as_nifti_gz"] is False
    assert [item["role"] for item in inventory["assets"] if not item["exists"]] == [
        "shape:Background", "shape:Cortex"
    ]
