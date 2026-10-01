import json
import hashlib
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_preprocessing_crop import crop_lesion_toads_preprocessing_inputs


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_preprocessing_crop_rigid"


def test_lesion_preprocessing_crop_maps_topology_under_rigid_transform():
    metadata = json.loads((ROOT / "inputs.json").read_text())
    image = nib.load(str(ROOT / metadata["image_nifti"]))
    topology_image = nib.load(str(ROOT / metadata["topology_nifti"]))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert nib.aff2axcodes(topology_image.affine) == ("R", "A", "S")
    normalized = np.transpose(np.asarray(image.dataobj, dtype=np.float32), (3, 0, 1, 2))
    topology = np.asarray(topology_image.dataobj)
    transform = np.asarray(metadata["initial_transform"], dtype=np.float32)
    assert np.any(transform[:3] != 0) and np.any(transform[3:] != 0)
    before = (normalized.copy(), topology.copy(), transform.copy())
    reference = json.loads((ROOT / "java_reference.json").read_text())
    harness = Path(__file__).parent / "java_reference" / "LesionPreprocessingCropRigidReferenceHarness.java"
    assert hashlib.sha256(harness.read_bytes()).hexdigest() == "43c1846cb48d2dec52dad9de1c85c805b7666bfeba5213335b950aaaaebda4e0"
    assert reference["source_sha256"] == "72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08"

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
    assert result.crop.images.shape[1:] == tuple(reference["cropped_dims_xyz"])
    assert list(result.crop.transform.view(np.int32)) == reference["cropped_transform_bits"]
    expected = np.asarray(reference["cropped_image_bits"], dtype=np.int32).view(np.float32)
    expected = expected.reshape(result.crop.images.shape)
    np.testing.assert_array_equal(result.crop.images.view(np.int32), expected.view(np.int32))
    assert np.any(result.crop.transform[:3] != 0)  # rotation is handed through unchanged
    assert result.crop.bounds[1] >= normalized.shape[1]
    assert result.crop.bounds[3] >= normalized.shape[2]
    assert result.crop.bounds[5] >= normalized.shape[3]
    assert np.all(result.crop.images[:, -1] == 0)  # mapped topology extends beyond all high image edges
    assert np.array_equal(normalized, before[0])
    assert np.array_equal(topology, before[1])
    assert np.array_equal(transform, before[2])
