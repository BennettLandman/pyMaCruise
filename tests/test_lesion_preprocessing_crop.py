import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_preprocessing_crop import crop_lesion_toads_preprocessing_inputs


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_preprocessing_crop"


def test_lesion_preprocessing_crop_matches_selected_java_excerpt():
    metadata = json.loads((ROOT / "inputs.json").read_text())
    image = nib.load(str(ROOT / metadata["image_nifti"]))
    topology_image = nib.load(str(ROOT / metadata["topology_nifti"]))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(topology_image.affine) == ("R", "A", "S")
    image_values = np.asarray(image.dataobj, dtype=np.float32)
    normalized = np.transpose(image_values, (3, 0, 1, 2))
    topology = np.asarray(topology_image.dataobj)
    transform = np.asarray(metadata["initial_transform"], dtype=np.float32)
    before = (normalized.copy(), topology.copy(), transform.copy())
    reference = json.loads((ROOT / "java_reference.json").read_text())

    result = crop_lesion_toads_preprocessing_inputs(
        normalized,
        topology,
        transform,
        image_resolution_mm=tuple(metadata["image_resolution_mm"]),
        template_resolution_mm=tuple(metadata["template_resolution_mm"]),
        orientation_signs=tuple(metadata["orientation_signs"]),
        padding=metadata["padding"],
    )

    assert result.crop.bounds == tuple(reference["bounds_inclusive"])
    assert result.padding == reference["padding"] == 10
    assert result.background_ratio.view(np.int32) == reference["background_ratio_bits"]
    assert result.crop.images.shape[1:] == tuple(reference["cropped_dims_xyz"])
    assert list(result.crop.transform.view(np.int32)) == reference["cropped_transform_bits"]
    expected = np.asarray(reference["cropped_image_bits"], dtype=np.int32).view(np.float32)
    expected = expected.reshape(result.crop.images.shape)
    np.testing.assert_array_equal(result.crop.images.view(np.int32), expected.view(np.int32))
    assert result.crop.bounds[0] < 3 and result.crop.bounds[1] > 7
    assert np.array_equal(normalized, before[0])
    assert np.array_equal(topology, before[1])
    assert np.array_equal(transform, before[2])
