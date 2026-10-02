from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.cli import main
from pymacruise.toads_cerebellum import (
    CEREBELLUM_AND_BRAINSTEM_LABELS,
    remove_brainstem_and_cerebellum,
)


FIXTURE = Path(__file__).parent / "data" / "synthetic" / "toads_cerebellum"


def test_remove_brainstem_and_cerebellum_transfers_labels_and_preserves_map() -> None:
    segmentation = nib.load(str(FIXTURE / "maseg.nii.gz"))
    expected_removed = len(CEREBELLUM_AND_BRAINSTEM_LABELS)
    for filename in ("gm.nii.gz", "wmfill.nii.gz", "wmmask.nii.gz"):
        source = nib.load(str(FIXTURE / filename))
        result = remove_brainstem_and_cerebellum(source, segmentation)
        source_values = np.asarray(source.dataobj)
        output_values = np.asarray(result.image.dataobj)
        labels = np.asarray(result.transferred_segmentation.dataobj)
        remove = np.isin(labels, tuple(CEREBELLUM_AND_BRAINSTEM_LABELS))
        assert result.removed_voxels == expected_removed
        assert result.image.shape == source.shape
        assert nib.aff2axcodes(result.image.affine) == ("R", "A", "S")
        assert result.image.get_data_dtype() == source.get_data_dtype()
        np.testing.assert_array_equal(output_values[remove], 0)
        np.testing.assert_array_equal(output_values[~remove], source_values[~remove])


def test_remove_cerebellum_cli_writes_three_native_ras_maps_and_provenance(tmp_path) -> None:
    output_dir = tmp_path / "cleaned"
    provenance = tmp_path / "run.json"
    result = main([
        "remove-cerebellum",
        "--segmentation", str(FIXTURE / "maseg.nii.gz"),
        "--gm", str(FIXTURE / "gm.nii.gz"),
        "--wm-fill", str(FIXTURE / "wmfill.nii.gz"),
        "--wm-mask", str(FIXTURE / "wmmask.nii.gz"),
        "--output-dir", str(output_dir),
        "--provenance", str(provenance),
    ])
    assert result == 0
    assert provenance.is_file()
    for filename in (
        "gm_cerebellum_removed.nii.gz",
        "wmfill_cerebellum_removed.nii.gz",
        "wmmask_cerebellum_removed.nii.gz",
    ):
        image = nib.load(str(output_dir / filename))
        assert image.shape == (9, 8, 7)
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert image.get_data_dtype() == np.dtype(np.float32)
        assert np.count_nonzero(np.asarray(image.dataobj) == 0) == 8
